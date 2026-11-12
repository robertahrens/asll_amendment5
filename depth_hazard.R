###############################################################################
# American Samoa longline: change in interaction risk if Amendment 5 gear
# requirements (70 m blank mainline at floats) were removed
#
# Approach
#   1. Rebuild hook EFFORT by depth from set-level gear (every hook position
#      in every basket), not from catch records.
#   2. Count catches by set x hook position for each species.
#   3. Fit a depth-specific hazard (catch per hook as a smooth function of
#      hook depth) shared across regulatory periods, with trip and year random
#      effects and a period fixed effect absorbing non-depth differences.
#   4. Apply the fitted depth hazard to the current (post-2011) fleet's sets
#      under counterfactual gear scenarios and report the risk ratio relative
#      to the gear as fished, with uncertainty from the depth smooth.
#
# Key assumption: depth-specific catchability estimated from the observed
# record (2006 onward) holds in the future.
###############################################################################

library(dplyr)
library(tidyr)
library(purrr)
library(mgcv)
library(ggplot2)

# ---------------------------------------------------------------- settings --
working_dir <- "set your working dir"

data_file  <- file.path(working_dir,"as_data.Rdata")
out_dir    <- file.path(working_dir,"dereg_risk")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

phi_deg     <- 68
# Set "shoaled" to the multiplier used to create PRED_DEPTH_68_SHOAL
depth_mults <- c(catenary = 1, shoaled = 0.70)
blank_m     <- 70
species     <- c("ALBACORE", "WAHOO", "OCEANIC WHITETIP SHARK",
                 "MAHIMAHI", "GREEN SEA TURTLE")
rare_species  <- "GREEN SEA TURTLE"
current_years <- NULL      # e.g. 2019:2025 to define the "current fleet"
n_draws       <- 2000
set.seed(2026)

# ------------------------------------------------- catenary depth function --
yoshihara_vectorized_v2 <- function(hook_number, total_basket_length,
                                    catenary_angle_deg, floatline_length,
                                    branchline_length, n_actual_hooks,
                                    blank_length, depth_mult = 1) {
  phi   <- catenary_angle_deg * pi / 180
  cotp  <- 1 / tan(phi)
  L     <- total_basket_length / 2
  S     <- (total_basket_length - 2 * blank_length) / (n_actual_hooks + 1)
  n_bl  <- blank_length / S
  denom <- n_actual_hooks + 2 * n_bl + 1
  j     <- hook_number + n_bl
  sag      <- L * (sqrt(1 + cotp^2) - sqrt((1 - 2 * j / denom)^2 + cotp^2))
  d0       <- floatline_length + branchline_length
  sag_apex <- L * tan(phi / 2)
  s        <- pmax((depth_mult * (d0 + sag_apex) - d0) / sag_apex, 0)
  depths <- d0 + s * sag
  depths[is.na(hook_number)] <- NA
  round(depths, 2)
}

# -------------------------------------------------------------------- data --
load(data_file)
# data file is the individual record of an interaction with the associated gear configuration derived from data extract and data_estim_hk_depth. 
# 1. Data checks ------------------------------------------------------------
cat("\n== Hook number checks ==\n")
print(as_data %>% summarise(
  n_catch            = n(),
  hk_na              = sum(is.na(HK_NUM)),
  hk_noninteger      = sum(HK_NUM != round(HK_NUM), na.rm = TRUE),
  hk_gt_hks_per_flt  = sum(HK_NUM > HKS_PER_FLT, na.rm = TRUE)))

gear_conflict <- as_data %>%
  distinct(TRIP_NUM, SET_NUM, CALC_BSK_LEN, FLTLN_LEN, BRNCHLN_LEN,
           HKS_PER_FLT, NUM_FLTS, NUM_HKS_SET, REG_TURTLE) %>%
  count(TRIP_NUM, SET_NUM) %>% filter(n > 1)
if (nrow(gear_conflict) > 0)
  warning(nrow(gear_conflict), " sets have more than one gear record.")

# 2. Set-level gear table ---------------------------------------------------
sets_all <- as_data %>%
  distinct(TRIP_NUM, SET_NUM, .keep_all = TRUE) %>%
  transmute(TRIP_NUM, SET_NUM, PERMIT_NUM, REG_TURTLE,
            year       = as.integer(format(HAUL_BEGIN_DTM, "%Y")),
            mainline_m = MAINLN_LEN_RPTD_M,
            bsk_len_flts = CALC_BSK_LEN,          # mainline / (floats - 1)
            fl         = FLTLN_LEN,
            bl         = BRNCHLN_LEN,
            hpf        = round(HKS_PER_FLT),
            n_flts     = NUM_FLTS,
            hks_set    = NUM_HKS_SET,
            n_bsk_flts = n_flts - 1,
            n_bsk_hks  = hks_set / hpf,
            # baskets from reported hooks where available, else floats - 1;
            # basket length recomputed from the same basket count
            n_bsk      = coalesce(n_bsk_hks, n_bsk_flts),
            bsk_len    = mainline_m / n_bsk,
            apex_catenary = yoshihara_vectorized_v2((hpf + 1) / 2, bsk_len,
                              phi_deg, fl, bl, hpf, REG_TURTLE, 1))

cat("\n== Reported hooks vs (floats - 1) x hooks per float ==\n")
print(summary(with(sets_all, n_bsk_hks / n_bsk_flts)))

cat("\n== Gear distributions (check the tails) ==\n")
gear_vars <- c("mainline_m", "n_flts", "hpf", "hks_set", "bsk_len_flts",
               "bsk_len", "fl", "bl", "apex_catenary")
print(round(sapply(sets_all[gear_vars], quantile,
                   probs = c(0, .01, .05, .5, .95, .99, 1), na.rm = TRUE), 1))

# Plausibility limits: review the tables above and the flagged sets before
# accepting these values
max_apex_catenary <- 500     # m, unshoaled catenary depth of the apex hook
bsk_len_range     <- c(200, 2000)   # m
basket_mismatch   <- 0.15    # max |hooks-based / floats-based basket count - 1|
# Basket length uses the hooks-based basket count, so a float-count mismatch
# does not by itself distort depth. Set TRUE to drop those sets anyway.
drop_basket_mismatch <- FALSE

sets_all <- sets_all %>%
  mutate(flag = case_when(
    is.na(bsk_len) | is.na(fl) | is.na(bl) | is.na(hpf) | is.na(year) ~ "missing gear",
    hpf <= 0                                                           ~ "hpf <= 0",
    bsk_len < bsk_len_range[1] | bsk_len > bsk_len_range[2]            ~ "basket length",
    apex_catenary > max_apex_catenary                                  ~ "apex too deep",
    !is.na(n_bsk_hks) & abs(n_bsk_hks / n_bsk_flts - 1) > basket_mismatch
                                                                       ~ "hooks vs floats",
    TRUE ~ "ok"))

cat("\n== Set flags ==\n")
print(count(sets_all, REG_TURTLE, flag))
flagged <- filter(sets_all, flag != "ok") %>% arrange(desc(apex_catenary))
write.csv(flagged, file.path(out_dir, "flagged_sets.csv"), row.names = FALSE)
cat("\nDeepest flagged sets (full list in flagged_sets.csv):\n")
print(flagged %>%
        select(TRIP_NUM, SET_NUM, flag, mainline_m, n_flts, hpf, hks_set,
               bsk_len_flts, bsk_len, apex_catenary) %>%
        head(15))

keep_flags <- if (drop_basket_mismatch) "ok" else c("ok", "hooks vs floats")
sets <- filter(sets_all, flag %in% keep_flags)
cat("\nSets retained:", nrow(sets), "of", nrow(sets_all), "\n")
print(count(sets, REG_TURTLE, name = "n_sets"))

# 3. Catch by set x hook position -------------------------------------------
catch <- as_data %>%
  filter(ENGLISH_NAME %in% species, !is.na(HK_NUM)) %>%
  semi_join(sets, by = c("TRIP_NUM", "SET_NUM")) %>%
  mutate(hk = as.integer(round(HK_NUM))) %>%
  left_join(select(sets, TRIP_NUM, SET_NUM, hpf), by = c("TRIP_NUM", "SET_NUM"))

cat("\n== Catches dropped for hook position outside 1..HKS_PER_FLT ==\n")
print(catch %>% group_by(ENGLISH_NAME) %>%
        summarise(n = n(), dropped = sum(hk < 1 | hk > hpf)))

catch <- catch %>%
  filter(hk >= 1, hk <= hpf) %>%
  count(TRIP_NUM, SET_NUM, hk, ENGLISH_NAME, name = "y")

# 4. Helpers ----------------------------------------------------------------
# Expand sets to one row per hook position; n_hooks = baskets per set.
# `s` must contain a `blank` column (m of blank mainline at each float).
make_cells <- function(s, mult) {
  s %>%
    mutate(hk = map(hpf, seq_len)) %>%
    unnest(hk) %>%
    mutate(n_hooks = n_bsk,
           depth   = yoshihara_vectorized_v2(hk, bsk_len, phi_deg, fl, bl,
                                             hpf, blank, mult))
}

fit_species <- function(sp, cells) {
  d <- cells %>%
    left_join(filter(catch, ENGLISH_NAME == sp) %>% select(-ENGLISH_NAME),
              by = c("TRIP_NUM", "SET_NUM", "hk")) %>%
    mutate(y = coalesce(y, 0L))
  rare <- sp %in% rare_species
  f <- as.formula(paste0(
    "y ~ s(depth, k = ", if (rare) 4 else 8, ") + regime",
    " + s(year_f, bs = 're') + s(trip_f, bs = 're') + offset(log(n_hooks))"))
  fit <- bam(f, family = if (rare) poisson() else nb(),
             data = d, method = "fREML", discrete = TRUE)
  list(fit = fit, data = d)
}

# Depth-only relative hazard on a 1 m grid, point estimate and draws.
# Intercept, period, trip and year terms are common to all scenarios and
# cancel in the ratios, so only the s(depth) coefficients are needed.
depth_hazard <- function(fit, depth_grid, template) {
  nd <- template[rep(1, length(depth_grid)), ]
  nd$depth <- depth_grid
  X    <- predict(fit, newdata = nd, type = "lpmatrix")
  cols <- grep("^s\\(depth\\)", colnames(X))
  Xd   <- X[, cols, drop = FALSE]
  b    <- coef(fit)[cols]
  B    <- rmvn(n_draws, b, vcov(fit)[cols, cols])
  list(point = exp(drop(Xd %*% b)),
       draws = exp(Xd %*% t(B)),           # depth x draws
       se    = sqrt(rowSums((Xd %*% vcov(fit)[cols, cols]) * Xd)),
       eta   = drop(Xd %*% b))
}

# Hooks aggregated to the 1 m grid for a scenario
exposure <- function(cells, depth_grid) {
  e <- cells %>% mutate(d = round(depth)) %>%
    group_by(d) %>% summarise(h = sum(n_hooks), .groups = "drop")
  out <- numeric(length(depth_grid))
  out[match(e$d, depth_grid)] <- e$h
  out
}

# 5. Scenario gear ----------------------------------------------------------
pre_fl_median <- median(sets$fl[sets$REG_TURTLE == 0])
cur <- sets %>% filter(REG_TURTLE == 70)
if (!is.null(current_years)) cur <- filter(cur, year %in% current_years)
pre <- sets %>% filter(REG_TURTLE == 0)

add_hooks <- function(s) {      # hooks that fit in the 70 m at current spacing
  S <- (s$bsk_len - 2 * blank_m) / (s$hpf + 1)
  2 * floor(blank_m / S)
}

scenarios <- list(
  base_regulated      = mutate(cur, blank = blank_m),
  S1_redistribute     = mutate(cur, blank = 0),
  S2_fill_blank       = mutate(cur, blank = 0, n_add = add_hooks(cur),
                               hks_set = hks_set + n_add * n_bsk,
                               hpf = hpf + n_add) %>% select(-n_add),
  S3_fill_short_float = mutate(cur, blank = 0, n_add = add_hooks(cur),
                               hks_set = hks_set + n_add * n_bsk,
                               hpf = hpf + n_add,
                               fl = pre_fl_median) %>% select(-n_add),
  S4_pre2011_gear     = mutate(pre, blank = 0)
)
cat("\nPre-2011 median floatline length (m):", pre_fl_median, "\n")

# 6. Fit and project --------------------------------------------------------
results  <- list()
rates    <- list()
decomp   <- list()
fit_summ <- list()
hazards  <- list()
diag_tab <- list()
exposures <- list()

for (dm in names(depth_mults)) {
  mult <- depth_mults[[dm]]
  cat("\n================ depth model:", dm, "(multiplier", mult, ") ====\n")

  fit_cells <- sets %>%
    mutate(blank  = REG_TURTLE,
           regime = factor(REG_TURTLE, levels = c(0, 70)),
           year_f = factor(year),
           trip_f = factor(TRIP_NUM)) %>%
    make_cells(mult)

  scen_cells <- map(scenarios, make_cells, mult = mult)

  obs_rng  <- range(fit_cells$depth, na.rm = TRUE)
  all_d    <- c(fit_cells$depth, unlist(map(scen_cells, "depth")))
  depth_grid <- seq(floor(min(all_d, na.rm = TRUE)),
                    ceiling(max(all_d, na.rm = TRUE)), by = 1)

  # extrapolation check
  ex <- map_dfr(scen_cells, ~ tibble(
    pct_hooks_outside_fit_range =
      100 * sum(.x$n_hooks[.x$depth < obs_rng[1] | .x$depth > obs_rng[2]]) /
      sum(.x$n_hooks)), .id = "scenario")
  cat("\nFitted depth range:", round(obs_rng, 1), "m\n")
  print(ex)

  expo <- map(scen_cells, exposure, depth_grid = depth_grid)
  exposures[[dm]] <- map_dfr(expo, ~ tibble(depth = depth_grid, hooks = .x),
                             .id = "scenario") %>% mutate(depth_model = dm)
  H <- map_dbl(expo, sum)

  for (sp in species) {
    cat("  fitting", sp, "\n")
    fs  <- fit_species(sp, fit_cells)
    fit <- fs$fit
    fit_summ[[length(fit_summ) + 1]] <- tibble(
      depth_model = dm, species = sp, n_cells = nrow(fs$data),
      intercept    = coef(fit)[["(Intercept)"]],
      intercept_se = sqrt(vcov(fit)["(Intercept)", "(Intercept)"]),
      n_interactions = sum(fs$data$y),
      regime70     = coef(fit)[["regime70"]],
      regime70_se  = sqrt(vcov(fit)["regime70", "regime70"]),
      theta = if (inherits(fit$family, "extended.family"))
                fit$family$getTheta(TRUE) else NA_real_,
      dev_expl  = summary(fit)$dev.expl,
      scale_est = fit$sig2)
    hz  <- depth_hazard(fit, depth_grid, fit_cells[1, ])

    # Expected interactions (relative units) per scenario, point and draws
    E_pt <- map_dbl(expo, ~ sum(.x * hz$point))
    E_dr <- map(expo, ~ colSums(.x * hz$draws))

    bad <- c(!is.finite(E_pt), map_lgl(E_dr, ~ any(!is.finite(.x))))
    if (any(bad)) {
      used <- depth_grid[Reduce(`|`, map(expo, ~ .x > 0))]
      stop(sp, " (", dm, "): non-finite expected interactions. ",
           "Depth curve range over hook depths in use: eta ",
           paste(round(range(hz$eta[depth_grid %in% used]), 1), collapse = " to "),
           ". Check gear outliers (flagged_sets.csv) and the s(depth) fit.")
    }

    base <- "base_regulated"

    # Observed baseline: interactions per set and per 1000 hooks in the
    # regulated-period sets, with a trip-level bootstrap for uncertainty
    trip_tab <- fs$data %>%
      semi_join(cur, by = c("TRIP_NUM", "SET_NUM")) %>%
      group_by(TRIP_NUM) %>%
      summarise(y = sum(y), sets = n_distinct(SET_NUM), hooks = sum(n_hooks),
                .groups = "drop")
    rate_set_obs  <- sum(trip_tab$y) / sum(trip_tab$sets)
    rate_hook_obs <- 1000 * sum(trip_tab$y) / sum(trip_tab$hooks)
    boot_idx <- replicate(n_draws, sample.int(nrow(trip_tab), replace = TRUE),
                          simplify = FALSE)
    rate_set_b  <- map_dbl(boot_idx, ~ sum(trip_tab$y[.x]) / sum(trip_tab$sets[.x]))
    rate_hook_b <- map_dbl(boot_idx, ~ 1000 * sum(trip_tab$y[.x]) / sum(trip_tab$hooks[.x]))

    add_rate <- function(sc, set_pt, set_draws, hook_pt, hook_draws) {
      tibble(depth_model = dm, species = sp, scenario = sc,
             n_interactions_base = sum(trip_tab$y),
             n_sets_base = sum(trip_tab$sets), n_trips_base = nrow(trip_tab),
             per_set    = set_pt,
             per_set_lo = unname(quantile(set_draws, 0.025)),
             per_set_hi = unname(quantile(set_draws, 0.975)),
             per_1000_hooks    = hook_pt,
             per_1000_hooks_lo = unname(quantile(hook_draws, 0.025)),
             per_1000_hooks_hi = unname(quantile(hook_draws, 0.975)))
    }
    rates[[length(rates) + 1]] <- add_rate(base, rate_set_obs, rate_set_b,
                                           rate_hook_obs, rate_hook_b)

    # Retrospective decomposition of the observed pre -> post change per set:
    #   gear (depth) component = regulated gear / pre-2011 gear, same hazard
    #   period component       = exp(regime coefficient), change in rate at
    #                            a given depth not explained by hook depth
    pre_tab <- fs$data %>%
      semi_join(pre, by = c("TRIP_NUM", "SET_NUM")) %>%
      group_by(TRIP_NUM) %>%
      summarise(y = sum(y), sets = n_distinct(SET_NUM), .groups = "drop")
    obs_pre  <- sum(pre_tab$y) / sum(pre_tab$sets)
    boot_pre <- map_dbl(seq_len(n_draws), function(i) {
      k <- sample.int(nrow(pre_tab), replace = TRUE)
      sum(pre_tab$y[k]) / sum(pre_tab$sets[k])
    })
    b_reg  <- coef(fit)[["regime70"]]
    se_reg <- sqrt(vcov(fit)["regime70", "regime70"])
    gear_pt <- (E_pt[[base]] / nrow(scenarios[[base]])) /
               (E_pt[["S4_pre2011_gear"]] / nrow(scenarios[["S4_pre2011_gear"]]))
    gear_dr <- (E_dr[[base]] / nrow(scenarios[[base]])) /
               (E_dr[["S4_pre2011_gear"]] / nrow(scenarios[["S4_pre2011_gear"]]))
    per_dr  <- exp(rnorm(n_draws, b_reg, se_reg))
    q <- function(x) unname(quantile(x, c(0.025, 0.975), na.rm = TRUE))
    decomp[[length(decomp) + 1]] <- tibble(
      depth_model = dm, species = sp,
      obs_pre_per_set  = obs_pre,
      obs_post_per_set = rate_set_obs,
      obs_ratio        = rate_set_obs / obs_pre,
      obs_ratio_lo     = q(rate_set_b / boot_pre)[1],
      obs_ratio_hi     = q(rate_set_b / boot_pre)[2],
      gear_component    = gear_pt,
      gear_lo           = q(gear_dr)[1],
      gear_hi           = q(gear_dr)[2],
      period_component  = exp(b_reg),
      period_lo         = exp(b_reg - 1.96 * se_reg),
      period_hi         = exp(b_reg + 1.96 * se_reg),
      model_ratio       = gear_pt * exp(b_reg),
      model_lo          = q(gear_dr * per_dr)[1],
      model_hi          = q(gear_dr * per_dr)[2])

    for (sc in setdiff(names(scenarios), base)) {
      # draws (named differently from the tibble columns, which tibble()
      # evaluates sequentially)
      set_dr  <- (E_dr[[sc]] / nrow(scenarios[[sc]])) /
                 (E_dr[[base]] / nrow(scenarios[[base]]))
      hook_dr <- (E_dr[[sc]] / H[[sc]]) / (E_dr[[base]] / H[[base]])
      results[[length(results) + 1]] <- tibble(
        depth_model = dm, species = sp, scenario = sc,
        per_set_ratio   = (E_pt[[sc]] / nrow(scenarios[[sc]])) /
                          (E_pt[[base]] / nrow(scenarios[[base]])),
        per_set_lo      = unname(quantile(set_dr, 0.025)),
        per_set_hi      = unname(quantile(set_dr, 0.975)),
        per_hook_ratio  = (E_pt[[sc]] / H[[sc]]) / (E_pt[[base]] / H[[base]]),
        per_hook_lo     = unname(quantile(hook_dr, 0.025)),
        per_hook_hi     = unname(quantile(hook_dr, 0.975)),
        p_increase      = mean(set_dr > 1))

      # Scenario rates = observed baseline x depth-driven ratio; baseline
      # bootstrap and depth-curve draws are combined draw by draw
      rates[[length(rates) + 1]] <- add_rate(
        sc,
        rate_set_obs  * (E_pt[[sc]] / nrow(scenarios[[sc]])) /
                        (E_pt[[base]] / nrow(scenarios[[base]])),
        rate_set_b * set_dr,
        rate_hook_obs * (E_pt[[sc]] / H[[sc]]) / (E_pt[[base]] / H[[base]]),
        rate_hook_b * hook_dr)
    }

    # relative hazard curve, scaled to 1 at its maximum within fitted range
    in_rng <- depth_grid >= obs_rng[1] & depth_grid <= obs_rng[2]
    ref    <- max(hz$eta[in_rng])
    hazards[[length(hazards) + 1]] <- tibble(
      depth_model = dm, species = sp, depth = depth_grid,
      rel   = exp(hz$eta - ref),
      lo    = exp(hz$eta - 1.96 * hz$se - ref),
      hi    = exp(hz$eta + 1.96 * hz$se - ref),
      in_fit_range = in_rng)

    # observed vs fitted catch by 20 m depth bin and period
    diag_tab[[length(diag_tab) + 1]] <- fs$data %>%
      mutate(fitted = fitted(fit), bin = floor(depth / 20) * 20 + 10) %>%
      group_by(regime, bin) %>%
      summarise(observed = sum(y), fitted = sum(fitted), .groups = "drop") %>%
      mutate(depth_model = dm, species = sp)

    print(summary(fit)$s.table)
  }
}

results  <- bind_rows(results)
hazards  <- bind_rows(hazards)
diag_tab <- bind_rows(diag_tab)
exposures <- bind_rows(exposures)

# 7. Output -----------------------------------------------------------------
write.csv(results, file.path(out_dir, "dereg_risk_ratios.csv"), row.names = FALSE)
write.csv(bind_rows(fit_summ), file.path(out_dir, "model_summary.csv"), row.names = FALSE)
write.csv(diag_tab, file.path(out_dir, "obs_vs_fitted_by_depth.csv"), row.names = FALSE)

decomp <- bind_rows(decomp)
write.csv(decomp, file.path(out_dir, "amend5_decomposition.csv"), row.names = FALSE)
cat("\n== Observed pre -> post change per set and its decomposition ==\n")
print(decomp %>%
        transmute(depth_model, species,
                  observed = sprintf("%.2f (%.2f-%.2f)", obs_ratio, obs_ratio_lo, obs_ratio_hi),
                  gear     = sprintf("%.2f (%.2f-%.2f)", gear_component, gear_lo, gear_hi),
                  period   = sprintf("%.2f (%.2f-%.2f)", period_component, period_lo, period_hi),
                  model    = sprintf("%.2f (%.2f-%.2f)", model_ratio, model_lo, model_hi)),
      n = Inf)

rates <- bind_rows(rates)
write.csv(rates, file.path(out_dir, "dereg_interaction_rates.csv"), row.names = FALSE)

cat("\n== Interactions per set and per 1000 hooks (observed regulated baseline x ratio) ==\n")
print(rates %>%
        transmute(depth_model, species, scenario,
                  per_set = sprintf("%.3f (%.3f-%.3f)", per_set, per_set_lo, per_set_hi),
                  per_1000_hooks = sprintf("%.4f (%.4f-%.4f)", per_1000_hooks,
                                           per_1000_hooks_lo, per_1000_hooks_hi)),
      n = Inf)

cat("\n== Risk ratios relative to regulated gear (per set) ==\n")
print(results %>%
        transmute(depth_model, species, scenario,
                  per_set = sprintf("%.2f (%.2f-%.2f)",
                                    per_set_ratio, per_set_lo, per_set_hi),
                  per_hook = sprintf("%.2f (%.2f-%.2f)",
                                     per_hook_ratio, per_hook_lo, per_hook_hi),
                  p_increase = round(p_increase, 3)),
      n = Inf)

p_hz <- ggplot(filter(hazards, in_fit_range),
               aes(depth, rel, colour = depth_model, fill = depth_model)) +
  geom_ribbon(aes(ymin = lo, ymax = pmin(hi, 3)), alpha = 0.2, colour = NA) +
  geom_line() +
  facet_wrap(~ species, scales = "free_y") +
  labs(x = "Hook depth (m)", y = "Relative catch per hook",
       colour = "Depth model", fill = "Depth model") +
  theme_bw()
ggsave(file.path(out_dir, "depth_hazard_curves.png"), p_hz, width = 10, height = 6)

p_ex <- exposures %>%
  filter(scenario %in% c("base_regulated", "S1_redistribute", "S4_pre2011_gear")) %>%
  group_by(depth_model, scenario) %>%
  mutate(prop = hooks / sum(hooks)) %>%
  ggplot(aes(depth, prop, colour = scenario)) +
  geom_line() +
  facet_wrap(~ depth_model) +
  labs(x = "Hook depth (m)", y = "Proportion of hooks", colour = "Scenario") +
  theme_bw()
ggsave(file.path(out_dir, "hook_depth_exposure.png"), p_ex, width = 10, height = 4)

p_dg <- ggplot(diag_tab, aes(bin, observed, colour = regime)) +
  geom_point() +
  geom_line(aes(y = fitted)) +
  facet_grid(species ~ depth_model, scales = "free_y") +
  labs(x = "Hook depth bin (m)", y = "Catch (points observed, lines fitted)",
       colour = "Blank (m)") +
  theme_bw()
ggsave(file.path(out_dir, "obs_vs_fitted_by_depth.png"), p_dg, width = 9, height = 11)

###############################################################################
# Simulated observer data for the American Samoa longline depth hazard
# analysis
#
# Produces `as_data`, one row per animal caught, with the same columns the
# analysis uses. All values are SIMULATED: trip, set and permit identifiers
# are fictitious and no record corresponds to a real vessel, trip or catch.
# The data exist so the analysis code can be run and checked end to end.
#
# Structure
#   - 109 pre-regulation trips (2006 to Sep 2011) and 154 regulated trips
#     (Sep 2011 to 2025), 20-38 sets per trip
#   - Set gear drawn from ranges typical of deep-set tuna longlining;
#     regulated sets use a 70 m blank mainline at each float and float lines
#     of at least 30 m
#   - Catches generated hook by hook from known depth hazards (below), with
#     trip-level random variation, so the analysis's estimates can be compared
#     with the truth
#   - A small number of sets with corrupted float counts or filled mainline
#     values to exercise the gear screening step
#
# Output: sim_as_data.Rdata (object `as_data`) and sim_as_data.csv
###############################################################################

set.seed(20261110)
working_dir <- "set your working dir"

# Yoshihara catenary depth with blank mainline (same as the analysis)
hook_depth <- function(hook, bsk_len, phi_deg, fl, bl, n_hooks, blank) {
  phi  <- phi_deg * pi / 180
  cotp <- 1 / tan(phi)
  L    <- bsk_len / 2
  S    <- (bsk_len - 2 * blank) / (n_hooks + 1)
  n_bl <- blank / S
  j    <- hook + n_bl
  fl + bl + L * (sqrt(1 + cotp^2) -
                 sqrt((1 - 2 * j / (n_hooks + 2 * n_bl + 1))^2 + cotp^2))
}

# True depth hazards: expected catch per hook at depth d (m)
true_hazard <- list(
  "ALBACORE"               = function(d) 0.012  * exp(-0.5 * ((d - 230) / 70)^2),
  "WAHOO"                  = function(d) 0.004  * exp(-0.5 * ((d - 110) / 60)^2),
  "OCEANIC WHITETIP SHARK" = function(d) 0.0012 * exp(-0.5 * ((d - 90)  / 50)^2),
  "MAHIMAHI"               = function(d) 0.003  * exp(-0.5 * ((d - 60)  / 40)^2),
  "GREEN SEA TURTLE"       = function(d) 0.00006 * exp(-d / 50)
)
trip_sd  <- c("ALBACORE" = 0.3, "WAHOO" = 0.3, "OCEANIC WHITETIP SHARK" = 0.3,
              "MAHIMAHI" = 0.5, "GREEN SEA TURTLE" = 0)
reg_date <- as.POSIXct("2011-09-23", tz = "UTC")
nm_to_m  <- 1852
vessels  <- sprintf("SIM%02d", 1:30)

rows <- list()
trip_id <- 0
for (period in c("pre", "regulated")) {
  n_trips <- if (period == "pre") 109 else 154
  for (t in seq_len(n_trips)) {
    trip_id <- trip_id + 1
    permit  <- sample(vessels, 1)
    start   <- if (period == "pre")
      as.POSIXct(runif(1, as.numeric(as.POSIXct("2006-01-15", tz = "UTC")),
                          as.numeric(reg_date) - 40 * 86400), origin = "1970-01-01", tz = "UTC")
    else
      as.POSIXct(runif(1, as.numeric(reg_date) + 86400,
                          as.numeric(as.POSIXct("2025-11-30", tz = "UTC"))),
                 origin = "1970-01-01", tz = "UTC")
    trip_eff <- sapply(trip_sd, function(s) rnorm(1, 0, s))
    n_sets   <- sample(20:38, 1)

    for (s in seq_len(n_sets)) {
      haul   <- start + (s - 1) * 86400 + 16 * 3600
      blank  <- if (haul >= reg_date) 70 else 0
      hpf    <- sample(25:35, 1)
      n_flts <- sample(80:120, 1)
      ml_nm  <- round(runif(1, 30, 50), 1)
      fl     <- round(if (blank == 70) runif(1, 30, 40) else runif(1, 15, 30), 1)
      bl     <- round(runif(1, 9, 13), 1)
      n_bsk  <- n_flts - 1
      hks    <- as.integer(n_bsk * hpf)
      bsk    <- ml_nm * nm_to_m / n_bsk
      d      <- hook_depth(1:hpf, bsk, 68, fl, bl, hpf, blank)

      for (sp in names(true_hazard)) {
        # catch per basket at each hook position, summed over baskets
        lam <- n_bsk * true_hazard[[sp]](d) * exp(trip_eff[[sp]])
        y   <- rpois(hpf, lam)
        if (sum(y) == 0) next
        hk  <- rep(seq_len(hpf), y)
        rows[[length(rows) + 1]] <- data.frame(
          TRIP_NUM          = sprintf("SIM%04d", trip_id),
          SET_NUM           = s,
          PERMIT_NUM        = permit,
          HAUL_BEGIN_DTM    = haul,
          MAINLN_LEN_RPTD   = ml_nm,
          FLTLN_LEN         = fl,
          BRNCHLN_LEN       = bl,
          NUM_FLTS          = n_flts,
          HKS_PER_FLT       = hpf,
          NUM_HKS_SET       = hks,
          FLT_NUM           = sample.int(n_bsk, length(hk), replace = TRUE),
          HK_NUM            = hk,
          ENGLISH_NAME      = sp,
          REG_TURTLE        = as.integer(blank),
          MAINLN_LEN_RPTD_M = ml_nm * nm_to_m,
          stringsAsFactors  = FALSE)
      }
    }
  }
}
as_data <- do.call(rbind, rows)

# Corrupt a few sets to exercise the gear screening:
#  (a) float count recorded far too low -> flagged "hooks vs floats" only,
#      since basket length is computed from reported hooks
#  (b) filled mainline length with very few reported hooks -> huge baskets,
#      flagged "basket length" or "apex too deep" and excluded
set_keys <- unique(as_data[c("TRIP_NUM", "SET_NUM")])
pick <- set_keys[sample(nrow(set_keys), 40), ]
key  <- paste(as_data$TRIP_NUM, as_data$SET_NUM)
a <- paste(pick$TRIP_NUM, pick$SET_NUM)[1:20]
b <- paste(pick$TRIP_NUM, pick$SET_NUM)[21:40]
as_data$NUM_FLTS[key %in% a] <- 20
fill_m <- 73149.19
as_data$MAINLN_LEN_RPTD[key %in% b]   <- NA
as_data$MAINLN_LEN_RPTD_M[key %in% b] <- fill_m
as_data$NUM_HKS_SET[key %in% b] <- as.integer(round(
  as_data$NUM_HKS_SET[key %in% b] * runif(sum(key %in% b), 0.1, 0.25)))

rownames(as_data) <- NULL
save(as_data, file = "sim_as_data.Rdata")
write.csv(as_data, "sim_as_data.csv", row.names = FALSE)

cat("Simulated catch records:", nrow(as_data), "\n")
cat("Sets:", nrow(unique(as_data[c("TRIP_NUM", "SET_NUM")])), "\n")
print(table(as_data$ENGLISH_NAME, as_data$REG_TURTLE))

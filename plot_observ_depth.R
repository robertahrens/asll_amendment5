library(dplyr)
library(ggplot2)


yoshihara_vectorized <- function(hook_number,
                                 total_basket_length, 
                                 catenary_angle_deg, 
                                 floatline_length, 
                                 branchline_length, 
                                 n_actual_hooks, 
                                 blank_length) {
  
  # Convert catenary angle from degrees to radians
  phi_rad <- catenary_angle_deg * pi / 180
  cot_phi <- 1 / tan(phi_rad)
  L <- total_basket_length / 2
  
  # Structural calculations for the basket
  active_length <- total_basket_length - (2 * blank_length)
  hook_spacing <- active_length / (n_actual_hooks + 1)
  n_blank_hooks <- blank_length / hook_spacing
  n_virtual_total <- n_actual_hooks + (2 * n_blank_hooks)
  denominator <- n_virtual_total + 1
  
  # Map actual hook numbers to virtual positions (vectorized math)
  j <- hook_number + n_blank_hooks
  
  # CORRECTED Yoshihara geometric calculation (sin_phi removed from term_1)
  term_1 <- sqrt(1 + cot_phi^2)
  term_2 <- sqrt((1 - (2 * j) / denominator)^2 + cot_phi^2)
  
  # Calculate depths across the vectors
  depths <- floatline_length + branchline_length + L * (term_1 - term_2)
  
  # Explicitly force indices that were originally NA in hook_number back to NA
  depths[is.na(hook_number)] <- NA
  
  return(round(depths, 2))
}

working_dir <- "set your working dir"
setwd(working_dir)

load("as_data.R")

as_data <- as_data %>%
  mutate(
    PRED_DEPTH = yoshihara_vectorized(
      hook_number = HK_NUM,
      total_basket_length = CALC_BSK_LEN, # Dynamic input
      catenary_angle_deg = 68,                          # Assumed fleet average
      floatline_length = FLTLN_LEN,
      branchline_length = BRNCHLN_LEN,
      n_actual_hooks = HKS_PER_FLT,
      blank_length = REG_TURTLE
    )
  )

plot_dat <- as_data %>%
  filter(
    !is.na(ENGLISH_NAME),
    !is.na(PRED_DEPTH),
    REG_TURTLE %in% c(0, 70)
  ) %>%
  mutate(
    REG_TURTLE = factor(
      REG_TURTLE,
      levels = c(0, 70),
      labels = c("0", "70")
    )
  )

set.seed(123)

jitter_dat <- plot_dat %>%
  group_by(REG_TURTLE, ENGLISH_NAME) %>%
  slice_sample(prop = 1) %>%
  slice_head(n = 500) %>%
  ungroup()

p1 <- ggplot(
  plot_dat,
  aes(
    x = ENGLISH_NAME,
    y = PRED_DEPTH,
    fill = REG_TURTLE
  )
) +
  geom_violin(
    position = position_dodge(width = 0.9),
    trim = FALSE,
    alpha = 0.5
  ) +
  geom_boxplot(
    aes(group = interaction(ENGLISH_NAME, REG_TURTLE)),
    position = position_dodge(width = 0.9),
    width = 0.12,
    outlier.shape = NA
  ) +
  geom_point(
    data = jitter_dat,
    aes(group = REG_TURTLE),
    position = position_jitterdodge(
      jitter.width = 0.08,
      dodge.width = 0.9
    ),
    alpha = 0.25,
    size = 1
  ) +
  scale_y_reverse(
    limits = c(500, 0)
  ) +
  scale_fill_discrete(
    name = "Turtle Reg.",
    labels = c(
      "0" = "None",
      "70" = "100m"
    )
  ) +
  labs(
    x = "Species",
    y = "Predicted interaction depth (m)"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    )
  )


ggplot2::ggsave(
  filename = "interaction_depth.png",
  plot = p1,
  width = 10,
  height = 7,
  units = "in",
  dpi = 300,
  bg = "white"
)

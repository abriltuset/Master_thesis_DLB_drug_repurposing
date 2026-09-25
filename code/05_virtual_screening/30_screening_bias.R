# ============================================================
# 30_screening_bias.R
# Minimal bias assessment of AutoDock Vina screening
# ============================================================

library(dplyr)

# -------------------------
# 1. Load data
# -------------------------

df <- read.csv(
  "vina_screening_results.csv",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# -------------------------
# 2. Spearman correlations
# -------------------------

cor_mw <- cor.test(
  df$vina_affinity_kcal_mol,
  df$mw_freebase,
  method = "spearman",
  exact = FALSE
)

cor_rot <- cor.test(
  df$vina_affinity_kcal_mol,
  df$rotatable_bonds,
  method = "spearman",
  exact = FALSE
)

cat("=== SPEARMAN CORRELATIONS ===\n\n")

cat("Vina score vs molecular weight\n")
cat("rho =", round(unname(cor_mw$estimate), 3), "\n")
cat("p =", format.pval(cor_mw$p.value, digits = 3), "\n\n")

cat("Vina score vs rotatable bonds\n")
cat("rho =", round(unname(cor_rot$estimate), 3), "\n")
cat("p =", format.pval(cor_rot$p.value, digits = 3), "\n")

# -------------------------
# 3. Additional MW analyses and save correlation summary
# -------------------------

# Compounds with negative Vina scores
df_negative <- df %>%
  filter(vina_affinity_kcal_mol < 0)

rho_negative <- cor(
  df_negative$vina_affinity_kcal_mol,
  df_negative$mw_freebase,
  method = "spearman",
  use = "complete.obs"
)

# Top 10% of the Vina ranking
n_top10 <- ceiling(nrow(df) * 0.10)

df_top10 <- df %>%
  arrange(vina_affinity_kcal_mol) %>%
  slice_head(n = n_top10)

rho_top10 <- cor(
  df_top10$vina_affinity_kcal_mol,
  df_top10$mw_freebase,
  method = "spearman",
  use = "complete.obs"
)

bias_summary <- data.frame(
  comparison = c(
    "Vina score vs molecular weight",
    "Vina score vs rotatable bonds",
    "Vina score vs molecular weight",
    "Vina score vs molecular weight"
  ),
  subset = c(
    "All compounds",
    "All compounds",
    "Negative Vina scores",
    "Top 10%"
  ),
  spearman_rho = c(
    unname(cor_mw$estimate),
    unname(cor_rot$estimate),
    rho_negative,
    rho_top10
  ),
  p_value = c(
    cor_mw$p.value,
    cor_rot$p.value,
    NA,
    NA
  )
)

cat("\nMW vs Vina - negative scores:",
    round(rho_negative, 3), "\n")

cat("MW vs Vina - top 10%:",
    round(rho_top10, 3), "\n")

write.csv(
  bias_summary,
  "screening_bias_summary.csv",
  row.names = FALSE
)

# -------------------------
# 4. Vina score distribution
# -------------------------

# Reference score: mean redocking score of M-MPEP
mmpep_score <- -10.906

# TFM palette
col_blue_main   <- "#355C7D"
col_blue_light  <- "#DCE6EF"
col_grey_light  <- "#F5F7F9"
col_grey_text   <- "#6E7B87"
col_text_main   <- "#1F2A33"
col_orange      <- "#C47A3A"
col_orange_light<- "#EED8C6"

png(
  "vina_score_distribution.png",
  width = 2400,
  height = 1500,
  res = 300,
  bg = "white"
)

# More left margin so Y label is not cut
par(
  mar = c(5.2, 6.2, 1.2, 1.2),
  mgp = c(3.0, 0.8, 0),
  family = "sans",
  col.axis = col_text_main,
  col.lab = col_text_main,
  col.main = col_text_main,
  fg = col_text_main,
  xaxs = "i",
  yaxs = "i"
)

# Histogram
h <- hist(
  df$vina_affinity_kcal_mol,
  breaks = 60,
  col = "#BFD0E0",
  border = "white",
  main = "",
  xlab = "",
  ylab = "",
  axes = FALSE
)

# Axes
axis(
  1,
  col = col_grey_text,
  col.axis = col_text_main,
  cex.axis = 0.95,
  lwd = 0.8
)

axis(
  2,
  col = col_grey_text,
  col.axis = col_text_main,
  cex.axis = 0.95,
  lwd = 0.8,
  las = 1
)

# Axis labels
mtext(
  "AutoDock Vina score (kcal/mol)",
  side = 1,
  line = 3.3,
  cex = 1.05,
  col = col_text_main
)

mtext(
  "Number of compounds",
  side = 2,
  line = 4.2,
  cex = 1.05,
  col = col_text_main
)

# Clean scientific-style box
box(bty = "l", col = col_grey_text, lwd = 0.9)

# Reference line for M-MPEP

abline(
  v = mmpep_score,
  col = col_orange,
  lty = 2,
  lwd = 2
)

# Right-side legend style annotation
usr <- par("usr")

x_leg_line_start <- usr[2] - 48
x_leg_line_end   <- usr[2] - 42
y_leg            <- usr[4] * 0.92

segments(
  x0 = x_leg_line_start,
  y0 = y_leg,
  x1 = x_leg_line_end,
  y1 = y_leg,
  col = col_orange,
  lty = 2,
  lwd = 2
)

text(
  x = usr[2] - 40.5,
  y = y_leg,
  labels = "M-MPEP = -10.91 kcal/mol",
  adj = c(0, 0.5),
  cex = 0.95,
  col = col_text_main
)

dev.off()

cat("\nGenerated files:\n")
cat("- screening_bias_summary.csv\n")
cat("- vina_score_distribution.png\n")


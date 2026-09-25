# ============================================================
# 29_screening_qc.R
# Quality control of the complete AutoDock Vina screening
# ============================================================

library(dplyr)

# -------------------------
# 1. Load screening results
# -------------------------

input_file <- "vina_screening_results.csv"

df <- read.csv(
  input_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

cat("=== SCREENING QC ===\n\n")

# -------------------------
# 2. Basic completeness
# -------------------------

cat("Total rows:", nrow(df), "\n")

cat("\nStatus:\n")
print(table(df$status, useNA = "ifany"))

cat("\nMissing Vina scores:",
    sum(is.na(df$vina_affinity_kcal_mol)), "\n")

cat("Missing ChEMBL IDs:",
    sum(is.na(df$parent_chembl_id) | df$parent_chembl_id == ""), "\n")

cat("Duplicated ChEMBL IDs:",
    sum(duplicated(df$parent_chembl_id)), "\n")

# -------------------------
# 3. Check ranking consistency
# -------------------------

df_sorted <- df %>%
  arrange(vina_affinity_kcal_mol)

expected_rank <- seq_len(nrow(df_sorted))

rank_mismatches <- sum(df_sorted$vina_rank != expected_rank)

cat("\nRank inconsistencies:", rank_mismatches, "\n")

# -------------------------
# 4. Score distribution
# -------------------------

scores <- df$vina_affinity_kcal_mol

cat("\n=== VINA SCORE DISTRIBUTION ===\n")

print(summary(scores))

cat("\nNegative scores:", sum(scores < 0), "\n")
cat("Zero scores:", sum(scores == 0), "\n")
cat("Positive scores:", sum(scores > 0), "\n")

cat(
  "Positive scores (%):",
  round(100 * mean(scores > 0), 1),
  "%\n"
)

# -------------------------
# 5. Ranking percentiles
# -------------------------

cat("\n=== SCORE QUANTILES ===\n")

print(
  quantile(
    scores,
    probs = c(0, 0.01, 0.025, 0.05, 0.10, 0.25, 0.50, 0.75, 1),
    na.rm = TRUE
  )
)

# -------------------------
# 6. Top 20 screening hits
# -------------------------

top20 <- df_sorted %>%
  slice_head(n = 20) %>%
  select(
    vina_rank,
    parent_chembl_id,
    name,
    vina_affinity_kcal_mol,
    bbb_martins_probability,
    withdrawn_flag,
    mw_freebase,
    alogp,
    psa,
    heavy_atoms,
    rotatable_bonds
  )

cat("\n=== TOP 20 ===\n")
print(top20, row.names = FALSE)

write.csv(
  top20,
  "top20_vina_hits.csv",
  row.names = FALSE
)

# -------------------------
# 7. Save QC summary
# -------------------------

qc_summary <- data.frame(
  metric = c(
    "Total compounds",
    "Completed",
    "Missing scores",
    "Duplicated ChEMBL IDs",
    "Negative scores",
    "Positive scores",
    "Minimum Vina score",
    "Median Vina score",
    "Maximum Vina score"
  ),
  value = c(
    nrow(df),
    sum(df$status == "completed"),
    sum(is.na(scores)),
    sum(duplicated(df$parent_chembl_id)),
    sum(scores < 0),
    sum(scores > 0),
    min(scores, na.rm = TRUE),
    median(scores, na.rm = TRUE),
    max(scores, na.rm = TRUE)
  )
)

write.csv(
  qc_summary,
  "screening_qc_summary.csv",
  row.names = FALSE
)

cat("\nFiles generated:\n")
cat("- top20_vina_hits.csv\n")
cat("- screening_qc_summary.csv\n")

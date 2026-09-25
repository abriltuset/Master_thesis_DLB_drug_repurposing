# ============================================================
# 31_select_high_scoring_hits.R
# Selection of high-scoring compounds for pharmacological review
# ============================================================

library(dplyr)

# -------------------------
# 1. Load complete screening
# -------------------------

df <- read.csv(
  "vina_screening_results.csv",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# -------------------------
# 2. Define structural cutoff
# -------------------------

vina_cutoff <- -8.5

# -------------------------
# 3. Select high-scoring compounds
# -------------------------

high_scoring_hits <- df %>%
  filter(vina_affinity_kcal_mol <= vina_cutoff) %>%
  arrange(vina_affinity_kcal_mol)

# -------------------------
# 4. Keep relevant existing variables
# -------------------------

high_scoring_table <- high_scoring_hits %>%
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
    rotatable_bonds,
    first_approval
  )

# -------------------------
# 5. Export structural shortlist
# -------------------------

write.csv(
  high_scoring_table,
  "vina_high_scoring_hits.csv",
  row.names = FALSE,
  na = ""
)

# -------------------------
# 6. Summary
# -------------------------

cat("=== HIGH-SCORING VINA HITS ===\n\n")

cat("Vina cutoff:", vina_cutoff, "kcal/mol\n")
cat("Total screened compounds:", nrow(df), "\n")
cat("Compounds passing cutoff:", nrow(high_scoring_table), "\n")

cat(
  "Percentage of screened library:",
  round(100 * nrow(high_scoring_table) / nrow(df), 2),
  "%\n"
)

cat(
  "Score range:",
  min(high_scoring_table$vina_affinity_kcal_mol),
  "to",
  max(high_scoring_table$vina_affinity_kcal_mol),
  "kcal/mol\n"
)

cat(
  "Withdrawn flag = 1:",
  sum(high_scoring_table$withdrawn_flag == 1, na.rm = TRUE),
  "\n"
)

cat("\nGenerated file:\n")
cat("- vina_high_scoring_hits.csv\n")

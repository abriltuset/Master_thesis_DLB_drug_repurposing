###############################################################################
# First iteration of the therapeutic target selection pipeline for DLB
#
# Reproducible implementation of the first target-prioritization iteration.
#
# Conceptual criteria preserved:
#   1) Control_vs_DLB
#   2) Excitatory neurons
#   3) |coef| > 2
#   4) Independent Braak evidence REQUIRED
#   5) Directional concordance REQUIRED
#   6) Protein-coding genes
#   7) Functional class compatible with small-molecule modulation
#   8) Additional biological/pharmacological/structural evidence
#
# Important difference vs SECOND iteration:
#   - Braak evidence is a HARD FILTER here.
#   - TDL is ANNOTATED but NOT used as a Tchem/Tclin hard filter.
#   - Open Targets LBD score is supportive evidence, NOT the ranking rule.
###############################################################################

# == 0. Setup ================================================================
library(readxl)
library(httr)
library(jsonlite)
library(dplyr)
library(stringr)

set.seed(42)
options(timeout = 300)

# == 1. Load transcriptomic datasets ==========================================

# Main DLB dataset: GSE178146
# Supplementary Table 5 contains the cell-type-specific DEGs reported by
# Feleke et al. after applying their original significance criteria
# (FDR < 0.05 and absolute fold change > 1.5).
# The filters below are additional prioritization criteria applied in this study.

cat("1. Loading transcriptomic datasets...\n")

# Main DLB dataset: GSE178146
results <- read_excel(
  "data/raw/main_snRNA_table.xlsx",
  sheet = 2
)

cat(sprintf(
  "   GSE178146: %d rows, %d columns\n",
  nrow(results),
  ncol(results)
))

# Independent Braak dataset: GSE216281
# Table contains significant DEGs for Braak 0 vs 5.
# skip = 1 because the first row contains the supplementary table title.
braak_degs <- read_excel(
  "data/raw/braak_table.xlsx",
  sheet = 1,
  skip = 1
)

cat(sprintf(
  "   GSE216281 Braak 0 vs 5: %d rows, %d columns\n",
  nrow(braak_degs),
  ncol(braak_degs)
))

# Human Protein Atlas
hpa <- read.csv(
  "data/raw/proteinatlas.tsv",
  sep = "\t",
  stringsAsFactors = FALSE,
  quote = "\""
)
cat(sprintf("   HPA proteinatlas: %d rows\n", nrow(hpa)))

# TCRD / Pharos
tcrd <- read.csv(
  "data/raw/PharosTCRD_UniProt_Mapping.tsv",
  sep = "\t",
  stringsAsFactors = FALSE
)
cat(sprintf("   TCRD/Pharos mapping: %d entries\n", nrow(tcrd)))

# Open Targets associations
ot_assoc <- read.csv(
  "data/derived/ot_associations_lbd_pd_direct.csv",
  stringsAsFactors = FALSE
)
cat(sprintf("   OT associations LBD+PD: %d rows\n", nrow(ot_assoc)))

# Open Targets targets
ot_targets <- read.csv(
  "data/derived/ot_targets_all.csv",
  stringsAsFactors = FALSE
)
cat(sprintf("   OT targets: %d rows\n", nrow(ot_targets)))


# == 2. Transcriptomic filters ================================================

cat("\n2. Applying transcriptomic filters...\n")

# Step 1: Control vs DLB
step1 <- results %>%
  filter(comparison == "Control_vs_DLB")

cat(sprintf(
  "   Control_vs_DLB: %d records\n",
  nrow(step1)
))

# Step 2: Excitatory neurons
step2 <- step1 %>%
  filter(cell_type == "Excitatory neuron")

cat(sprintf(
  "   Excitatory neuron: %d records\n",
  nrow(step2)
))

# Step 3: absolute coefficient > 2
step3 <- step2 %>%
  filter(abs(coef) > 2)

cat(sprintf(
  "   |coef| > 2: %d records (%d unique genes)\n",
  nrow(step3),
  n_distinct(step3$hgnc_symbol)
))

candidates <- step3 %>%
  select(
    hgnc_symbol,
    coef,
    fdr
  ) %>%
  distinct()

duplicate_candidates <- candidates %>%
  count(hgnc_symbol) %>%
  filter(n > 1)

if (nrow(duplicate_candidates) > 0) {
  stop(
    "Duplicated genes detected after transcriptomic filtering. ",
    "Inspect GSE178146 before continuing."
  )
}

cat(sprintf(
  "   Unique candidate genes: %d\n",
  nrow(candidates)
))


# == 3. Braak independent transcriptomic validation ===========================

cat("\n3. Applying Braak validation from GSE216281...\n")

braak_lookup <- braak_degs %>%
  select(
    hgnc_symbol,
    log2FoldChange,
    padj
  ) %>%
  mutate(
    log2FoldChange = as.numeric(log2FoldChange),
    padj = as.numeric(padj)
  ) %>%
  filter(
    !is.na(hgnc_symbol),
    hgnc_symbol != ""
  )

braak_duplicates <- braak_lookup %>%
  count(hgnc_symbol) %>%
  filter(n > 1)

if (nrow(braak_duplicates) > 0) {
  stop(
    "Duplicated genes detected in the Braak dataset. ",
    "Inspect GSE216281 before joining."
  )
}

candidates <- candidates %>%
  left_join(
    braak_lookup,
    by = "hgnc_symbol",
    relationship = "many-to-one"
  ) %>%
  mutate(
    Braak_DEG = if_else(
      !is.na(log2FoldChange),
      "YES",
      "NO"
    ),
    Braak_Direction = case_when(
      is.na(log2FoldChange) ~ NA_character_,
      (coef > 0 & log2FoldChange > 0) |
        (coef < 0 & log2FoldChange < 0) ~ "consistent",
      TRUE ~ "opposite"
    )
  ) %>%
  rename(
    Braak_log2FC = log2FoldChange,
    Braak_FDR = padj
  )

n_braak_yes <- sum(
  candidates$Braak_DEG == "YES",
  na.rm = TRUE
)

n_braak_consistent <- sum(
  candidates$Braak_Direction == "consistent",
  na.rm = TRUE
)

cat(sprintf(
  "   Braak significant DEGs: %d genes\n",
  n_braak_yes
))

cat(sprintf(
  "   Braak direction consistent: %d genes\n",
  n_braak_consistent
))

# FIRST-ITERATION DIFFERENCE:
# Braak evidence and consistent direction are REQUIRED.
candidates_braak <- candidates %>%
  filter(
    Braak_DEG == "YES",
    Braak_Direction == "consistent"
  )

cat(sprintf(
  "   After Braak hard filter: %d genes\n",
  nrow(candidates_braak)
))


# == 4. Protein-coding filter =================================================

cat("\n4. Filtering protein-coding genes (Open Targets biotype)...\n")

ot_protein_coding <- ot_targets %>%
  filter(biotype == "protein_coding") %>%
  pull(approvedSymbol) %>%
  unique()

candidates_braak <- candidates_braak %>%
  mutate(
    protein_coding = ifelse(
      hgnc_symbol %in% ot_protein_coding,
      "YES",
      "NO"
    )
  )

n_pc <- sum(
  candidates_braak$protein_coding == "YES"
)

n_not_pc <- sum(
  candidates_braak$protein_coding == "NO"
)

cat(sprintf(
  "   Protein-coding (OT biotype): %d genes\n",
  n_pc
))

cat(sprintf(
  "   Non-protein-coding (excluded): %d genes\n",
  n_not_pc
))

candidates_pc <- candidates_braak %>%
  filter(protein_coding == "YES")


# == 5. Curated functional classification ====================================

cat("\n5. Applying curated functional classification from first iteration...\n")

# HPA is still used to obtain UniProt and Ensembl identifiers.
# Functional classification itself corresponds to the curated review
# performed during the exploratory first iteration and is therefore
# NOT derived automatically from HPA Protein.class.

hpa_select <- hpa %>%
  select(
    Gene,
    Protein.class,
    Uniprot,
    Ensembl
  ) %>%
  filter(!is.na(Gene)) %>%
  group_by(Gene) %>%
  summarise(
    Protein.class = first(Protein.class),
    Uniprot = first(Uniprot),
    Ensembl = first(Ensembl),
    .groups = "drop"
  )

candidates_pc <- candidates_pc %>%
  left_join(
    hpa_select,
    by = c("hgnc_symbol" = "Gene"),
    relationship = "many-to-one"
  )


# ---------------------------------------------------------------------------
# Curated functional classification used in the first iteration
# ---------------------------------------------------------------------------

functional_review <- data.frame(
  hgnc_symbol = c(
    "ANK1",
    "CCDC24",
    "HP1BP3",
    "APLP1",
    "C3orf18",
    "ARHGAP9",
    "PNKD",
    "ATP4A",
    "JRK",
    "ABCB9",
    "ESYT1",
    "TSPAN9",
    "GPI",
    "PSMA7",
    "PLEKHJ1",
    "FLYWCH2",
    "CDC42BPA"
  ),
  
  function_type = c(
    "structural",
    "structural",
    "transcription regulation",
    "signaling",
    "poorly characterized",
    "signaling",
    "enzyme",
    "transporter",
    "transcription regulation",
    "transporter",
    "signaling",
    "membrane organization",
    "enzyme",
    "protein degradation",
    "signaling",
    "transcription regulation",
    "enzyme"
  ),
  
  stringsAsFactors = FALSE
)


# Check that every protein-coding candidate reaching this stage
# has a curated functional annotation.
missing_function_review <- setdiff(
  candidates_pc$hgnc_symbol,
  functional_review$hgnc_symbol
)

if (length(missing_function_review) > 0) {
  stop(
    paste0(
      "Missing curated functional classification for: ",
      paste(missing_function_review, collapse = ", ")
    )
  )
}


# Add curated functional annotation
candidates_pc <- candidates_pc %>%
  left_join(
    functional_review,
    by = "hgnc_symbol",
    relationship = "many-to-one"
  )


cat("\n   Curated function type distribution:\n")

print(
  table(
    candidates_pc$function_type,
    useNA = "ifany"
  )
)

# Retain proteins classified as signaling-related, transporters or enzymes.
candidates_func <- candidates_pc %>%
  filter(
    function_type %in% c(
      "signaling",
      "transporter",
      "enzyme"
    )
  )


cat(sprintf(
  "\n   After curated function type filter: %d genes\n",
  nrow(candidates_func)
))

cat("\n   Retained candidates:\n")

print(
  candidates_func %>%
    select(
      hgnc_symbol,
      coef,
      Braak_log2FC,
      function_type
    ) %>%
    arrange(desc(abs(coef)))
)


# == 6. Druggability annotation ===============================================

cat("\n6. Annotating druggability (TCRD/Pharos TDL)...\n")

tdl_rank <- c(
  Tdark = 1,
  Tbio = 2,
  Tchem = 3,
  Tclin = 4
)

tcrd_tdl <- tcrd %>%
  select(
    Uniprot = UniProt_accession,
    TDL
  ) %>%
  filter(
    !is.na(Uniprot),
    Uniprot != ""
  ) %>%
  mutate(
    TDL_rank = unname(tdl_rank[TDL])
  ) %>%
  group_by(Uniprot) %>%
  arrange(desc(TDL_rank), .by_group = TRUE) %>%
  slice(1) %>%
  ungroup() %>%
  select(-TDL_rank)

candidates_func <- candidates_func %>%
  left_join(
    tcrd_tdl,
    by = "Uniprot",
    relationship = "many-to-one"
  )

cat("\n   TDL distribution (annotation only; NOT a filter):\n")

print(
  table(
    candidates_func$TDL,
    useNA = "ifany"
  )
)

# FIRST-ITERATION DIFFERENCE:
# No Tchem/Tclin hard filter is applied here.


# == 7. Intermediate candidate table ==========================================

cat("\n7. Generating first-iteration filtered candidates table...\n")

table_intermediate <- candidates_func %>%
  arrange(desc(abs(coef))) %>%
  select(
    Gene = hgnc_symbol,
    Coef = coef,
    FDR = fdr,
    `Braak log2FC` = Braak_log2FC,
    `Braak FDR` = Braak_FDR,
    `Braak Direction` = Braak_Direction,
    `Function Type` = function_type,
    Druggability = TDL
  ) %>%
  mutate(
    Coef = round(Coef, 2),
    FDR = signif(FDR, 2),
    `Braak log2FC` = round(`Braak log2FC`, 3),
    `Braak FDR` = signif(`Braak FDR`, 2)
  )

write.csv(
  table_intermediate,
  "results/first_iteration_candidates_filtered.csv",
  row.names = FALSE
)

cat(sprintf(
  "   Saved: first_iteration_candidates_filtered.csv (%d genes)\n",
  nrow(table_intermediate)
))


# == 8. Additional evidence ===================================================

cat("\n8. Gathering additional evidence for first-iteration candidates...\n")

final_candidates <- candidates_func


# ---------------------------------------------------------------------------
# 8a. Open Targets association scores
# ---------------------------------------------------------------------------

cat("\n   8a. Open Targets association scores...\n")

ot_lbd <- ot_assoc %>%
  filter(diseaseId == "MONDO_0007488") %>%
  group_by(targetId) %>%
  summarise(
    OT_LBD_score = max(
      associationScore,
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  rename(
    Ensembl = targetId
  )

ot_pd <- ot_assoc %>%
  filter(diseaseId == "MONDO_0005180") %>%
  group_by(targetId) %>%
  summarise(
    OT_PD_score = max(
      associationScore,
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  rename(
    Ensembl = targetId
  )

final_candidates <- final_candidates %>%
  left_join(
    ot_lbd,
    by = "Ensembl",
    relationship = "many-to-one"
  ) %>%
  left_join(
    ot_pd,
    by = "Ensembl",
    relationship = "many-to-one"
  )

cat(sprintf(
  "   OT LBD found: %d / %d\n",
  sum(!is.na(final_candidates$OT_LBD_score)),
  nrow(final_candidates)
))

cat(sprintf(
  "   OT PD found: %d / %d\n",
  sum(!is.na(final_candidates$OT_PD_score)),
  nrow(final_candidates)
))


# ---------------------------------------------------------------------------
# 8b. Tissue specificity
# ---------------------------------------------------------------------------

cat("\n   8b. HPA tissue specificity...\n")

format_tissue_specificity <- function(specificity, specific_ntpm) {
  
  if (is.na(specificity) || specificity == "") {
    return(NA_character_)
  }
  
  if (
    is.na(specific_ntpm) ||
    specific_ntpm == "" ||
    !specificity %in% c(
      "Tissue enriched",
      "Tissue enhanced",
      "Group enriched"
    )
  ) {
    return(specificity)
  }
  
  tissues <- str_split(
    specific_ntpm,
    ";"
  )[[1]] %>%
    str_trim() %>%
    str_remove(":.*$") %>%
    str_to_title()
  
  paste0(
    specificity,
    " (",
    paste(tissues, collapse = ", "),
    ")"
  )
}

hpa_tissue <- hpa %>%
  select(
    Gene,
    RNA.tissue.specificity,
    RNA.tissue.specific.nTPM
  ) %>%
  filter(
    !is.na(Gene),
    Gene != ""
  ) %>%
  group_by(Gene) %>%
  summarise(
    RNA.tissue.specificity = first(RNA.tissue.specificity),
    RNA.tissue.specific.nTPM = first(RNA.tissue.specific.nTPM),
    .groups = "drop"
  ) %>%
  rowwise() %>%
  mutate(
    Tissue_Specificity = format_tissue_specificity(
      RNA.tissue.specificity,
      RNA.tissue.specific.nTPM
    )
  ) %>%
  ungroup() %>%
  select(
    Gene,
    Tissue_Specificity
  )

final_candidates <- final_candidates %>%
  left_join(
    hpa_tissue,
    by = c("hgnc_symbol" = "Gene"),
    relationship = "many-to-one"
  )

cat("   Tissue specificity distribution:\n")

print(
  table(
    final_candidates$Tissue_Specificity,
    useNA = "ifany"
  )
)


# ---------------------------------------------------------------------------
# 8c. PubMed literature evidence
# ---------------------------------------------------------------------------

cat("\n   8c. PubMed literature evidence...\n")


# ---------------------------------------------------------------------------
# Add official gene names
# ---------------------------------------------------------------------------

ot_gene_names <- ot_targets %>%
  select(
    hgnc_symbol = approvedSymbol,
    approved_name = approvedName
  ) %>%
  filter(
    !is.na(hgnc_symbol),
    hgnc_symbol != ""
  ) %>%
  group_by(hgnc_symbol) %>%
  summarise(
    approved_name = first(
      approved_name[
        !is.na(approved_name) &
          approved_name != ""
      ],
      default = NA_character_
    ),
    .groups = "drop"
  )


final_candidates <- final_candidates %>%
  left_join(
    ot_gene_names,
    by = "hgnc_symbol",
    relationship = "many-to-one"
  )


# ---------------------------------------------------------------------------
# PubMed query
# ---------------------------------------------------------------------------

query_pubmed <- function(gene_sym, gene_name = NA_character_) {
  
  # GPI is a highly ambiguous symbol (e.g. glycosylphosphatidylinositol).
  # Therefore, only its official full gene name is used.
  if (
    gene_sym == "GPI" &&
    !is.na(gene_name) &&
    gene_name != ""
  ) {
    
    gene_terms <- sprintf(
      '"%s"[Title/Abstract]',
      gene_name
    )
    
  } else if (
    !is.na(gene_name) &&
    gene_name != ""
  ) {
    
    # Default strategy:
    # official HGNC symbol OR official full gene name
    gene_terms <- sprintf(
      '("%s"[Title/Abstract] OR "%s"[Title/Abstract])',
      gene_sym,
      gene_name
    )
    
  } else {
    
    # Fallback if no official full name is available
    gene_terms <- sprintf(
      '"%s"[Title/Abstract]',
      gene_sym
    )
  }
  
  
  # Disease/context block:
  # DLB + Parkinson disease + synucleinopathies + alpha-synuclein.
  # MeSH terms are combined with free-text terms to improve sensitivity.
  disease_terms <- paste0(
    "(",
    
    '"Lewy Body Disease"[MeSH Terms]',
    ' OR "dementia with Lewy bodies"[Title/Abstract]',
    ' OR "Lewy body dementia"[Title/Abstract]',
    
    ' OR "Parkinson Disease"[MeSH Terms]',
    ' OR "Parkinson disease"[Title/Abstract]',
    ' OR "Parkinson\'s disease"[Title/Abstract]',
    
    ' OR "Synucleinopathies"[MeSH Terms]',
    ' OR synucleinopath*[Title/Abstract]',
    
    ' OR "alpha-Synuclein"[MeSH Terms]',
    ' OR "alpha-synuclein"[Title/Abstract]',
    
    ")"
  )
  
  
  search_term <- paste(
    gene_terms,
    "AND",
    disease_terms
  )
  
  
  # PubMed E-utilities
  url <- paste0(
    "https://eutils.ncbi.nlm.nih.gov/",
    "entrez/eutils/esearch.fcgi"
  )
  
  
  params <- list(
    db = "pubmed",
    term = search_term,
    retmode = "json",
    retmax = 0,
    tool = "DLB_TFM_pipeline"
  )
  
  
  tryCatch({
    
    r <- GET(
      url,
      query = params,
      timeout(30)
    )
    
    
    if (status_code(r) == 200) {
      
      data <- content(
        r,
        "parsed"
      )
      
      return(
        list(
          hits = as.integer(
            data$esearchresult$count
          ),
          query = search_term
        )
      )
    }
    
    
    list(
      hits = NA_integer_,
      query = search_term
    )
    
    
  }, error = function(e) {
    
    list(
      hits = NA_integer_,
      query = search_term
    )
  })
}


# ---------------------------------------------------------------------------
# Run PubMed queries
# ---------------------------------------------------------------------------

pubmed_results <- lapply(
  seq_len(nrow(final_candidates)),
  function(i) {
    
    gene <- final_candidates$hgnc_symbol[i]
    gene_name <- final_candidates$approved_name[i]
    
    
    cat(sprintf(
      "      PubMed: %s\n",
      gene
    ))
    
    
    res <- query_pubmed(
      gene,
      gene_name
    )
    
    
    Sys.sleep(0.4)
    
    res
  }
)


# Number of retrieved publications
final_candidates$PubMed_Hits <- sapply(
  pubmed_results,
  `[[`,
  "hits"
)


# Keep exact query for traceability
final_candidates$PubMed_Query <- sapply(
  pubmed_results,
  `[[`,
  "query"
)


cat(sprintf(
  "   PubMed hits gathered for %d genes\n",
  nrow(final_candidates)
))


# ---------------------------------------------------------------------------
# 8d. Structural evidence: experimental / computational / none
# ---------------------------------------------------------------------------

cat("\n   8d. Structural evidence...\n")

query_pdb_experimental <- function(uniprot_acc) {
  
  if (
    is.na(uniprot_acc) ||
    uniprot_acc == ""
  ) {
    return(0L)
  }
  
  query_list <- list(
    type = "terminal",
    service = "text",
    parameters = list(
      attribute = paste0(
        "rcsb_polymer_entity_container_identifiers.",
        "reference_sequence_identifiers.database_accession"
      ),
      operator = "exact_match",
      value = uniprot_acc
    )
  )
  
  body <- list(
    query = query_list,
    return_type = "entry",
    request_options = list(
      return_all_hits = TRUE
    )
  )
  
  tryCatch({
    
    r <- POST(
      "https://search.rcsb.org/rcsbsearch/v2/query",
      body = body,
      encode = "json",
      timeout(30)
    )
    
    if (status_code(r) == 200) {
      
      data <- content(
        r,
        "parsed"
      )
      
      return(
        as.integer(data$total_count)
      )
    }
    
    if (status_code(r) == 204) {
      return(0L)
    }
    
    NA_integer_
    
  }, error = function(e) {
    
    NA_integer_
  })
}

query_alphafold <- function(uniprot_acc) {
  
  if (
    is.na(uniprot_acc) ||
    uniprot_acc == ""
  ) {
    return(FALSE)
  }
  
  url <- paste0(
    "https://alphafold.ebi.ac.uk/api/prediction/",
    uniprot_acc
  )
  
  tryCatch({
    
    r <- GET(
      url,
      timeout(30)
    )
    
    if (status_code(r) == 200) {
      
      data <- content(
        r,
        "parsed"
      )
      
      return(
        length(data) > 0
      )
    }
    
    FALSE
    
  }, error = function(e) {
    
    FALSE
  })
}

classify_structure <- function(uniprot_acc) {
  
  if (
    is.na(uniprot_acc) ||
    uniprot_acc == ""
  ) {
    return("none")
  }
  
  n_pdb <- query_pdb_experimental(
    uniprot_acc
  )
  
  Sys.sleep(0.4)
  
  if (
    !is.na(n_pdb) &&
    n_pdb > 0
  ) {
    return("experimental")
  }
  
  has_af <- query_alphafold(
    uniprot_acc
  )
  
  Sys.sleep(0.4)
  
  if (isTRUE(has_af)) {
    return("computational")
  }
  
  return("none")
}

final_candidates$structure_evidence <- sapply(
  final_candidates$Uniprot,
  function(up) {
    
    cat(sprintf(
      "      Structure: %s\n",
      up
    ))
    
    classify_structure(up)
  }
)

cat("\n   Structural classification distribution:\n")

print(
  table(
    final_candidates$structure_evidence,
    useNA = "ifany"
  )
)


# == 9. Generate final first-iteration table ==================================

cat("\n9. Generating first-iteration annotated candidates table...\n")

# IMPORTANT:
# No single quantitative ranking is imposed in the first iteration.
# Sorting by |coef| is only for presentation.
final_table <- final_candidates %>%
  arrange(desc(abs(coef))) %>%
  select(
    Gene = hgnc_symbol,
    Coef = coef,
    FDR = fdr,
    `Braak log2FC` = Braak_log2FC,
    `Braak FDR` = Braak_FDR,
    `Braak Direction` = Braak_Direction,
    `Function Type` = function_type,
    Druggability = TDL,
    `OT LBD Score` = OT_LBD_score,
    `OT PD Score` = OT_PD_score,
    `PubMed Hits` = PubMed_Hits,
    `Tissue Specificity` = Tissue_Specificity,
    `Structure evidence` = structure_evidence
  ) %>%
  mutate(
    Coef = round(Coef, 2),
    FDR = signif(FDR, 2),
    `Braak log2FC` = round(
      `Braak log2FC`,
      3
    ),
    `Braak FDR` = signif(
      `Braak FDR`,
      2
    ),
    `OT LBD Score` = round(
      `OT LBD Score`,
      4
    ),
    `OT PD Score` = round(
      `OT PD Score`,
      4
    )
  )

write.csv(
  final_table,
  "results/first_iteration_candidates_annotated.csv",
  row.names = FALSE
)

cat(sprintf(
  "   Saved: first_iteration_candidates_annotated.csv (%d genes)\n",
  nrow(final_table)
))

cat("\n   Final first-iteration candidates (ordered by |coef| for display only):\n")

for (i in seq_len(nrow(final_table))) {
  
  cat(sprintf(
    "   %2d. %-12s | coef=%6.2f | Braak=%7.3f | func=%s | TDL=%s | OT_LBD=%s\n",
    i,
    final_table$Gene[i],
    final_table$Coef[i],
    final_table$`Braak log2FC`[i],
    final_table$`Function Type`[i],
    ifelse(
      is.na(final_table$Druggability[i]),
      "NA",
      final_table$Druggability[i]
    ),
    ifelse(
      is.na(final_table$`OT LBD Score`[i]),
      "NA",
      sprintf("%.4f", final_table$`OT LBD Score`[i])
    )
  ))
}


# == 10. Pipeline summary =====================================================

cat("\n=== PIPELINE SUMMARY (FIRST ITERATION) ===\n")

summary_table <- data.frame(
  Stage = c(
    "Control_vs_DLB",
    "Excitatory neuron",
    "|coef| > 2",
    "Braak DEG",
    "Braak direction consistent",
    "Braak hard filter",
    "Protein-coding (OT biotype)",
    "Function type filter",
    "Final candidates"
  ),
  N = c(
    nrow(step1),
    nrow(step2),
    nrow(candidates),
    n_braak_yes,
    n_braak_consistent,
    nrow(candidates_braak),
    nrow(candidates_pc),
    nrow(candidates_func),
    nrow(final_table)
  )
)

print(
  summary_table,
  row.names = FALSE
)

write.csv(
  summary_table,
  "results/first_iteration_pipeline_summary.csv",
  row.names = FALSE
)

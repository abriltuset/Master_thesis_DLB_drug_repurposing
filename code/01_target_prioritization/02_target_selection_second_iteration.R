###############################################################################
# Second iteration of the therapeutic target selection pipeline for DLB
#
# Key differences from the first iteration:
#   - Braak evidence is retained as supporting evidence, not as a hard filter.
#   - Functional classes are derived reproducibly from HPA Protein.class.
#   - TDL Tchem/Tclin is used as a mandatory druggability filter.
#   - Candidates are ranked by Open Targets LBD association score,
#     using |coef| only as a secondary criterion.
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
# The input corresponds to the cell-type-specific DEG table reported by
# Feleke et al. The |coef| > 2 threshold below is an additional
# prioritization criterion applied in this study.

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
# skip = 1 because the first row contains the supplementary table title
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

# HPA
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

# == 2. Transcriptomic filters ===============================================

cat("\n2. Applying transcriptomic filters...\n")

# Step 1: Control vs DLB
step1 <- results %>%
  filter(comparison == "Control_vs_DLB")

cat(sprintf(
  "   Control_vs_DLB: %d records\n",
  nrow(step1)
))


# Step 2: Excitatory neuron
step2 <- step1 %>%
  filter(cell_type == "Excitatory neuron")

cat(sprintf(
  "   Excitatory neuron: %d records\n",
  nrow(step2)
))


# Step 3: additional prioritization threshold |log2FC| > 2
step3 <- step2 %>%
  filter(abs(coef) > 2)

cat(sprintf(
  "   |coef| > 2: %d records (%d unique genes)\n",
  nrow(step3),
  n_distinct(step3$hgnc_symbol)
))


# Candidate table from GSE178146
candidates <- step3 %>%
  select(
    hgnc_symbol,
    coef,
    fdr
  ) %>%
  distinct()


# Verify that each gene appears only once
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

# == 3. Braak external transcriptomic evidence ===============================

cat("\n3. Adding Braak evidence from GSE216281...\n")


# Keep the variables needed from the independent dataset
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


# Check that there is only one row per gene
braak_duplicates <- braak_lookup %>%
  count(hgnc_symbol) %>%
  filter(n > 1)

if (nrow(braak_duplicates) > 0) {
  stop(
    "Duplicated genes detected in the Braak dataset. ",
    "Inspect GSE216281 before joining."
  )
}

# Join Braak information to GSE178146 candidates
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

# == 4. Protein-coding filter =================================================
cat("\n4. Filtering protein-coding genes (Open Targets biotype)...\n")

ot_protein_coding <- ot_targets %>%
  filter(biotype == "protein_coding") %>%
  pull(approvedSymbol) %>%
  unique()

candidates <- candidates %>%
  mutate(
    protein_coding = ifelse(
      hgnc_symbol %in% ot_protein_coding,
      "YES",
      "NO"
    )
  )

n_pc <- sum(candidates$protein_coding == "YES")
n_not_pc <- sum(candidates$protein_coding == "NO")

cat(sprintf(
  "   Protein-coding (OT biotype): %d genes\n",
  n_pc
))

cat(sprintf(
  "   Non-protein-coding (excluded): %d genes\n",
  n_not_pc
))

candidates_pc <- candidates %>%
  filter(protein_coding == "YES")


# == 5. Function type classification =========================================
cat("\n5. Classifying function type from HPA Protein class...\n")

# Guarantee one HPA row per gene before joining
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

# Multi-label classification
classify_function_multilabel <- function(pc_str) {
  
  if (
    is.na(pc_str) ||
    pc_str == "" ||
    pc_str == "None"
  ) {
    return("")
  }
  
  s <- tolower(pc_str)
  cats <- c()
  
  if (str_detect(s, "g-protein coupled receptors")) {
    cats <- c(cats, "signaling")
  }
  
  if (str_detect(s, "nuclear receptors")) {
    cats <- c(cats, "signaling")
  }
  
  if (str_detect(s, "cd markers")) {
    cats <- c(cats, "signaling")
  }
  
  if (str_detect(s, "transporters")) {
    cats <- c(cats, "transporter")
  }
  
  if (str_detect(s, "voltage-gated ion channels")) {
    cats <- c(cats, "transporter")
  }
  
  if (str_detect(s, "enzymes")) {
    cats <- c(cats, "enzyme")
  }
  
  if (length(cats) == 0) {
    return("other")
  }
  
  paste(
    sort(unique(cats)),
    collapse = ", "
  )
}

candidates_pc <- candidates_pc %>%
  mutate(
    function_type = sapply(
      Protein.class,
      classify_function_multilabel
    )
  )

cat("\n   Function type distribution (all protein-coding):\n")

print(
  table(
    candidates_pc$function_type,
    useNA = "ifany"
  )
)

candidates_func <- candidates_pc %>%
  filter(
    function_type != "other",
    function_type != ""
  )

cat(sprintf(
  "\n   After function type filter: %d genes\n",
  nrow(candidates_func)
))


# == 6. Druggability filter ===================================================
cat("\n6. Filtering by druggability (Tchem/Tclin)...\n")

# TDL hierarchy used only to guarantee one value per UniProt accession
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

candidates_pc <- candidates_pc %>%
  left_join(
    tcrd_tdl,
    by = "Uniprot",
    relationship = "many-to-one"
  )

candidates_func <- candidates_func %>%
  left_join(
    tcrd_tdl,
    by = "Uniprot",
    relationship = "many-to-one"
  )

cat("\n   TDL distribution (after function type filter):\n")

print(
  table(
    candidates_func$TDL,
    useNA = "ifany"
  )
)

# Strict Tchem/Tclin filter
candidates_tchem <- candidates_func %>%
  filter(
    TDL %in% c("Tchem", "Tclin")
  )

cat(sprintf(
  "\n   Tchem/Tclin candidates: %d genes\n",
  nrow(candidates_tchem)
))


# == 7. Intermediate table ====================================================
cat("\n7. Generating second_iteration_candidates_protein_coding table...\n")

table_all_pc <- candidates_pc %>%
  arrange(desc(abs(coef))) %>%
  select(
    Gene = hgnc_symbol,
    Coef = coef,
    FDR = fdr,
    `Braak DEG` = Braak_DEG,
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
    `Braak FDR` = signif(`Braak FDR`, 2),
  )

write.csv(
  table_all_pc,
  "results/second_iteration_candidates_protein_coding.csv",
  row.names = FALSE
)

cat(sprintf(
  "   Saved: second_iteration_candidates_protein_coding.csv (%d genes)\n",
  nrow(table_all_pc)
))


# == 8. Additional evidence ===================================================
cat("\n8. Gathering additional evidence for Tchem/Tclin candidates...\n")

final_candidates <- candidates_tchem


# ---------------------------------------------------------------------------
# 8a. Open Targets association scores
# ---------------------------------------------------------------------------
cat("\n   8a. Open Targets association scores...\n")

# Aggregate to one row per target to prevent row multiplication
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

# HPA reports the specificity category and the associated tissue(s)
# separately. For example, GRM5:
#   RNA tissue specificity     = "Tissue enriched"
#   RNA tissue specific nTPM   = "brain: 13.9"
# is reported here as "Tissue enriched (Brain)".

format_tissue_specificity <- function(specificity, specific_ntpm) {
  
  if (is.na(specificity) || specificity == "") {
    return(NA_character_)
  }
  
  # Categories without a specific enriched tissue
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
  
  # Extract tissue names from strings such as:
  # "brain: 13.9"
  # or "brain: 20.1; retina: 14.3"
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


# Query experimental PDB structures associated with a UniProt accession
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


# Check whether an AlphaFold model exists
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
  
  # Experimental PDB structure
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
  
  # Computational AlphaFold model
  has_af <- query_alphafold(
    uniprot_acc
  )
  
  Sys.sleep(0.4)
  
  if (isTRUE(has_af)) {
    return("computational")
  }
  
  # Neither experimental nor AlphaFold structure
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


# == 9. Generate final table ==================================================
cat("\n9. Generating second_iteration_candidates_tchem_tclin table...\n")

# IMPORTANT:
# Sorting uses the original unrounded values.
# Rounding is applied only for presentation/output.
final_table <- final_candidates %>%
  arrange(desc(abs(coef))) %>%
  select(
    Gene = hgnc_symbol,
    Coef = coef,
    FDR = fdr,
    `Braak DEG` = Braak_DEG,
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
  "results/second_iteration_candidates_tchem_tclin.csv",
  row.names = FALSE
)

cat(sprintf(
  "   Saved: second_iteration_candidates_tchem_tclin.csv (%d genes)\n",
  nrow(final_table)
))


# == 10. Shortlist ============================================================
cat("\n10. Creating shortlist (top 10 by OT_LBD score)...\n")

# Ranking uses full-precision values.
shortlist <- final_candidates %>%
  mutate(
    OT_LBD_score = coalesce(
      OT_LBD_score,
      0
    )
  ) %>%
  arrange(
    desc(OT_LBD_score),
    desc(abs(coef))
  ) %>%
  slice_head(n = 10) %>%
  select(
    Gene = hgnc_symbol,
    Coef = coef,
    FDR = fdr,
    `Braak DEG` = Braak_DEG,
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
  shortlist,
  "results/second_iteration_shortlist_top10.csv",
  row.names = FALSE
)

cat(sprintf(
  "   Saved: second_iteration_shortlist_top10.csv (%d genes)\n",
  nrow(shortlist)
))

cat("\n   Shortlist (top 10 by OT_LBD score):\n")

for (i in seq_len(nrow(shortlist))) {
  
  cat(sprintf(
    "   %2d. %-12s | OT_LBD=%.4f | coef=%.2f | func=%s | TDL=%s\n",
    i,
    shortlist$Gene[i],
    shortlist$`OT LBD Score`[i],
    shortlist$Coef[i],
    shortlist$`Function Type`[i],
    shortlist$Druggability[i]
  ))
}


# == 11. Pipeline summary =====================================================

cat("\n=== PIPELINE SUMMARY (SECOND ITERATION) ===\n")

summary_table <- data.frame(
  Stage = c(
    "Control_vs_DLB",
    "Excitatory neuron",
    "|coef| > 2",
    "Braak DEG",
    "Protein-coding (OT biotype)",
    "Function type filter",
    "Druggability (Tchem/Tclin)",
    "Shortlist (top 10 OT_LBD)"
  ),
  N = c(
    nrow(step1),
    nrow(step2),
    nrow(candidates),
    n_braak_yes,
    nrow(candidates_pc),
    nrow(candidates_func),
    nrow(candidates_tchem),
    nrow(shortlist)
  )
)

print(
  summary_table,
  row.names = FALSE
)

write.csv(
  summary_table,
  "results/second_iteration_pipeline_summary.csv",
  row.names = FALSE
)


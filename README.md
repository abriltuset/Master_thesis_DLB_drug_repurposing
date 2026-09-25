# From transcriptomics to drug repurposing in Dementia with Lewy Bodies

Computational pipeline for therapeutic target prioritization and drug repurposing in **Dementia with Lewy Bodies (DLB)**.

This repository contains the code, selected intermediate outputs, structural files and final results generated for a Master's Thesis in Bioinformatics.

## Project at a glance

This project aimed to develop a computational workflow for **therapeutic target prioritization and approved-drug repurposing in Dementia with Lewy Bodies (DLB)**. The approach integrates **human brain transcriptomics, target prioritization, structural bioinformatics, cheminformatics and molecular docking** to identify druggable targets and generate experimentally testable repurposing hypotheses.

## Key results

- Prioritized **GRM5 / mGluR5 from human DLB transcriptomic data**, followed by functional, pharmacological and structural assessment.

- Reduced an initial set of **3,475 approved small molecules** to a final library of **1,477 docking-ready compounds with predicted BBB permeability** (`P(BBB) ≥ 0.50`), after chemical and technical filtering.

- Validated the docking protocol by redocking **M-MPEP**, with a mean heavy-atom RMSD of approximately **0.27 Å** across four seeds.

- Screened all **1,477 compounds** against the experimentally characterized mGluR5 NAM pocket.

- Selected **37 structural hits** and reduced them to **8 final candidates** after pharmacological and structural review.

## Tech stack

**Languages:** Python, R, Bash\
**Transcriptomics & data integration:** GEO, Open Targets, Human Protein Atlas\
**Cheminformatics:** RDKit, ChEMBL API, ADMET-AI, Molscrub\
**Structural bioinformatics:** AlphaFold, P2Rank, fpocket, UCSF Chimera\
**Molecular docking:** Meeko, AutoDock Vina\
**Reproducibility:** Conda environments, scripted multi-stage workflow

## Workflow

### 1. Target prioritization

Human transcriptomic data were obtained from public GEO datasets:

- **GSE178146**: single-nucleus RNA-seq from human DLB brain tissue

- **GSE216281**: independent bulk transcriptomic dataset used as supporting evidence

Genes were progressively filtered according to differential expression, functional class, druggability and disease-related evidence.

A first prioritization iteration identified **APLP1**, but structural assessment did not support a sufficiently defensible small-molecule binding site. After refinement of the prioritization strategy, **GRM5**, encoding the metabotropic glutamate receptor 5 (**mGluR5**), was selected as the final target.

### 2. Structural assessment

The first candidate, APLP1, was evaluated using AlphaFold structural confidence together with pocket predictions from P2Rank and fpocket.

For GRM5, the human mGluR5 crystal structure **PDB 6FFI** was selected. The receptor was prepared and the experimentally characterized negative allosteric modulator binding site was retained for docking.

### 3. Docking validation

The co-crystallized ligand **M-MPEP** was redocked into the mGluR5 NAM pocket.

The protocol was evaluated using heavy-atom RMSD relative to the crystallographic pose and repeated using four random seeds.

Mean RMSD across replicates was approximately **0.27 Å**, supporting recovery of the experimental binding geometry.

### 4. Approved-drug library

Approved small molecules were retrieved from **ChEMBL 37** and normalized to parent compounds.

The library was subsequently filtered using:

- structural availability
- predicted blood-brain barrier permeability
- chemical compatibility with the docking workflow
- successful 3D preparation
- successful PDBQT generation

The final screening library contained **1,477 compounds**.

Blood-brain barrier predictions were generated using **ADMET-AI**, applying:

`P(BBB) >= 0.50`

as the operational prefiltering threshold.

### 5. Virtual screening

Docking was performed with **AutoDock Vina 1.2.7** against the mGluR5 NAM binding site.

Main screening parameters:

- Receptor: `6FFI_dry.pdbqt`
- Center: `(-23.720, -4.974, 41.888)`
- Box size: `(13.367, 18.945, 18.362) Å`
- Exhaustiveness: `32`
- Number of modes: `20`
- Energy range: `3 kcal/mol`
- Grid spacing: `0.375 Å`
- Random seed: `2026`

All 1,477 compounds were successfully screened.

Compounds with:

`Vina score <= -8.5 kcal/mol`

were retained as structural hits, resulting in **37 compounds**.

### 6. Pharmacological prioritization

The 37 docking hits underwent additional pharmacological review including:

- withdrawal status
- systemic exposure of the docked molecular entity
- compatibility with potential chronic use in DLB
- evidence of CNS exposure

Eight compounds remained after this review.

Their docking poses were subsequently assessed for occupation of the known mGluR5 NAM pocket and for abnormal short contacts.

Final prioritization was based primarily on:

- evidence of brain exposure
- suitability for potential chronic use in DLB

Among the analysed compounds, **fluoxetine showed the highest relative priority within the final candidate set**.

This prioritization should be interpreted as hypothesis generation rather than evidence of direct mGluR5 binding or therapeutic efficacy.

## Repository structure

``` text
.
├── code/
│   ├── 01_target_prioritization/
│   ├── 02_aplp1_structural_feasibility/
│   ├── 03_grm5_structure_and_redocking/
│   ├── 04_drug_library_preparation/
│   └── 05_virtual_screening/
│
├── structures/
│   ├── 6FFI_dry.pdbqt
│   ├── 6ffi_J_D8B.sdf
│   └── 6FFI_dry.box.txt
│
├── results/
│   ├── target_prioritization/
│   ├── docking_validation/
│   ├── virtual_screening/
│   └── pharmacological_review/
│
├── environment/
│   ├── drug_repurposing.yml
│   └── admet_ai.yml
│
└── README.md
```

## Environments

Two Conda environments were used during the workflow.

### Main drug-repurposing environment

``` bash
conda env create -f environment/drug_repurposing.yml
conda activate drug_repurposing
```

This environment contains the main cheminformatics and docking tools, including RDKit, Meeko, Molscrub and AutoDock Vina.

### ADMET-AI environment

``` bash
conda env create -f environment/admet_ai.yml
conda activate admet_ai
```

This environment was used for blood-brain barrier prediction with ADMET-AI.

## Running the workflow

Scripts are numbered according to the order of the computational workflow and are grouped by analysis stage.

The repository reflects the workflow used for the thesis. Some scripts use relative input/output paths corresponding to the original working environment, so local path adjustment may be required when reproducing individual steps.

Large raw datasets, the complete PDBQT ligand library, docking poses and individual Vina log files are not included in the repository.

Selected intermediate and final outputs are provided under `results/` to document the analyses reported in the thesis.

## Main outputs

Key files include:

- `first_iteration_candidates_annotated.csv`
- `second_iteration_shortlist_top10.csv`
- `MMPEP_redocking_dry_rmsd.csv`
- `MMPEP_redocking_replicates_summary.csv`
- `vina_screening_results.csv`
- `screening_qc_summary.csv`
- `screening_bias_summary.csv`
- `vina_high_scoring_hits.csv`
- `pose_review_summary_8hits_dry.csv`
- `vina_pharmacological_review.xlsx`

## Interpretation and limitations

The docking scores reported in this project are used for **relative structural prioritization** and should not be interpreted as experimentally measured binding affinities.

Likewise, predicted BBB permeability does not demonstrate in vivo brain exposure.

The final drug candidates therefore represent computationally prioritized hypotheses requiring experimental validation of:

1.  direct interaction with mGluR5
2.  functional modulation of the receptor
3.  CNS exposure at relevant concentrations
4.  potential therapeutic relevance in DLB

## Thesis

**De la transcriptómica a la priorización farmacológica: identificación de dianas y reposicionamiento de fármacos en demencia con cuerpos de Lewy**

Master's Thesis in Bioinformatics

Abril Tuset Andreu\
2025-2026

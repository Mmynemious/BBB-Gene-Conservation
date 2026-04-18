# Project Handover: BBB Gene Conservation Analysis
*This file was written to bring Claude Code fully up to speed. Read this before doing anything.*

---

## Who You're Working With
- **Name:** Yara Shamshoum
- **Level:** Beginner in R, Python, and bioinformatics. Has a strong conceptual understanding of the biology but is new to coding and data analysis pipelines.
- **Supervision:** Dr. Clelland (Clelland Lab). Andrew Yang (Gladstone) may provide additional data later.
- **GitHub:** Mmynemious
- **Guiding principle:** Always explain *why* you're doing each step conceptually before writing any code. Yara wants to understand what each step means in the grand scheme of the project, not just copy-paste commands.

---

## The Research Question
**How similar are the genomes of human, non-human primate (rhesus macaque), and mouse across blood-brain barrier (BBB) genes?**

Specifically: compare BBB gene coding regions, non-coding regions, and regulatory regions across three species:
- *Homo sapiens* (human)
- *Macaca mulatta* (rhesus macaque)
- *Mus musculus* (mouse)

**Why this matters:** Mouse models are widely used for BBB research, but we don't know how well their BBB genes actually mirror human ones. If the BBB genes are evolutionarily conserved, mouse models are valid. If not, that has major implications for drug development and disease research.

**Control group needed:** A randomly selected set of genes (matched for length, chromosome distribution) to compare against BBB genes — to establish whether any conservation we find is unusual or just baseline genomic conservation. *Dr. Clelland must approve the control group strategy before it's finalised.*

---

## What Has Been Done Already

### ✅ Data collection complete
All six raw datasets are downloaded, renamed, and organised into the repo. Git is initialised with an initial commit. The repo still needs to be pushed to GitHub (see below).

### ✅ Repo structure established
```
BBB/
├── README.md               ← project overview
├── CONTEXT.md              ← this file
├── .gitignore
├── raw_data/
│   ├── daneman_2010/       ← 7 supplementary tables (.xls)
│   ├── munji_2019/         ← 3 supplementary files (.xlsx)
│   ├── winkler_2022/       ← 10 supplementary tables (.xlsx)
│   ├── yang_2022/          ← 8 supplementary tables (.xlsx)
│   ├── walchli_2024/       ← 1 Seurat object (.rds.gz)
│   └── giger_2010/         ← 1 GEO series matrix (.txt)
├── processed/              ← EMPTY — outputs go here
├── scripts/                ← EMPTY — analysis scripts go here
└── docs/                   ← research papers (PDF) + Notion export
```

### ✅ GitHub repo needs pushing
The local git repo is ready. Yara still needs to:
1. Create an empty repo called `BBB-Gene-Conservation` on github.com (no README, no .gitignore)
2. Run from the BBB folder:
```bash
git remote add origin https://github.com/Mmynemious/BBB-Gene-Conservation.git
git push -u origin main
```

---

## Dataset Inventory — What Each File Contains

### Daneman et al. 2010 — Mouse BBB Gene Reference
**Folder:** `raw_data/daneman_2010/`
**Species:** Mouse
**Method:** Microarray — FACS-sorted endothelial cells from brain, liver, lung

| File | Rows | Contents |
|------|------|----------|
| S1_FullExpressionMatrix_AllTissues.xls | 45,102 | All probe sets across all tissues — use only if re-deriving gene lists |
| S2_BrainVascular_vs_Parenchyma_Enriched.xls | 2,155 | Genes enriched in brain vascular vs brain parenchyma |
| **S3_CoreBBBGenes_BrainEC_Enriched.xls** | **214** | **THE KEY FILE: core BBB genes enriched in brain ECs vs BOTH liver and lung** |
| S4_BrainEC_vs_LiverEC_Enriched.xls | 1,009 | Brain EC vs liver EC only |
| S5_BrainEC_vs_LungEC_Enriched.xls | 949 | Brain EC vs lung EC only |
| S6_PostnatalBrainEC_Enriched.xls | 542 | Postnatal-specific BBB genes |
| S7_AdultBrainEC_Enriched.xls | 671 | Adult-specific BBB genes |

**Columns in S3:** Probe set, Gene Title, Gene Symbol, Brain GFP- Avg, Brain GFP+ Avg
**Note:** Gene symbols are in mouse format (e.g. *Slc2a1*, not *SLC2A1*)

---

### Munji et al. 2019 — Updated Mouse BBB Transcriptome
**Folder:** `raw_data/munji_2019/`
**Species:** Mouse
**Method:** RNA-seq — brain endothelial cells in health and multiple disease models

| File | Rows/Sheets | Contents |
|------|-------------|----------|
| S3_FullExpressionData_HealthAndDisease.xlsx | 37,992 / 4 sheets | Full expression matrix: Replicates, Averages, Health stats, Disease stats, Coexpression |
| **S5_BBBEnrichedGenes.xlsx** | **519 / 'BBB-enriched'** | **THE KEY FILE: updated mouse BBB-enriched gene list** |
| S7_PeripheralEndothelialGenes.xlsx | 1,400 / 'Peripheral Specific' | Genes that go UP when BBB breaks down — useful as negative contrast |

**Columns in S5:** Gene, Ensembl ID, Gene Description, Brain Vascular1, Brain Vascular2, Brain Endo1...
**Note:** Gene symbols are in mouse format

---

### Winkler et al. 2022 — Human Brain Vasculature Single-Cell Atlas
**Folder:** `raw_data/winkler_2022/science.abi7377_tables_s1_to_s10/`
**Species:** Human
**Method:** Single-cell RNA-seq — 181,388 cells from adult cerebrovasculature

| File | Contents |
|------|----------|
| **S10_VascularCellState_GeneSets_UCell.xlsx** | **THE KEY FILE: curated gene sets per vascular cell state for UCell scoring** |
| S1_CellTypeAnnotations.xlsx | Cell type labels and cluster annotations |
| S2–S6: Various cell type marker genes (endothelial, pericyte, SMC, fibroblast, macrophage) |
| S7–S9: AVM (arteriovenous malformation) differential genes |

**Structure of S10:** Cell State column (Artery, Capillary, Venule, Vein, Endothelium, Pericyte, SMC, Fibroblast...) + Genes column (comma-separated gene symbols in human format)

---

### Yang et al. 2022 — Human Brain Vascular Atlas (with Alzheimer's dimension)
**Folder:** `raw_data/yang_2022/`
**Species:** Human
**Method:** VINE-seq — 143,793 single-nucleus transcriptomes, hippocampus + cortex, 17 subjects (healthy + AD)

| File | Contents |
|------|----------|
| S1_PatientMetadata.xlsx | Patient demographics (age, sex, brain region, PMI, AD status) |
| **S2_DifferentialExpression_CellTypes.xlsx** | DEGs per cell type (Gene, p_val, avg_logFC, pct.1, pct.2) |
| S3_CortexVsHippocampus_DEGs.xlsx | Regional comparison |
| **S4_AverageExpression_PerCellType.xlsx** | **Average expression per cell type from Seurat AverageExpression — ready to use** |
| S5_MuralCell_SubtypeMarkers.xlsx | Pericyte/SMC subtypes |
| S6_EndothelialSubtype_Markers.xlsx | Endothelial subtypes |
| S7_AlzheimersRisk_Genes.xlsx | AD-associated genes in vasculature |
| S8_ProportionalExpression.xlsx | Cell type proportions |

---

### Wälchli et al. 2024 — Comprehensive Human Brain Vasculature Single-Cell Atlas
**Folder:** `raw_data/walchli_2024/`
**Species:** Human
**Method:** Single-cell RNA-seq — 606,380 cells across fetal, adult, and 5 CNS pathologies
**File:** `Walchli2024_AdultBrain_TemporalLobe_SortedEndothelial_SeuratObject.rds.gz`

This is a compressed Seurat object (.rds.gz). To use it:
```r
library(Seurat)
dataset <- readRDS(gzcon(file("path/to/file.rds.gz", "rb")))
avg <- AverageExpression(dataset)$RNA
write.csv(avg, "processed/Walchli2024_AverageExpression.csv")
```
**Note:** Large file. Requires R ≥ 4.2 and Seurat ≥ 4.0. This produces the average expression per cluster/cell type that feeds into downstream analysis.

---

### Giger et al. 2010 — Primate Neuronal vs Endothelial Transcriptome Evolution
**Folder:** `raw_data/giger_2010/`
**Species:** Human, Chimpanzee, Rhesus Macaque
**Method:** Affymetrix microarray (GPL570) — laser-capture microdissected neurons and endothelial cells
**File:** `Giger2010_GSE12293_Human_NeuronEndothelial_ExpressionMatrix.txt`

This is a GEO series matrix file. Format:
- Lines starting with `!` = metadata headers (skip these)
- Lines starting with `"ID_REF"` = start of data matrix
- Columns = samples (brain_endothelial_cells_rep1, rep2, rep3...)
- Rows = Affymetrix probe IDs (need to map to gene symbols using GPL570 annotation)

**Key insight from this paper:** Endothelial genes are significantly more conserved across primates than neuronal genes. This supports the hypothesis that BBB genes will be conserved.

---

## The 5 Next Steps (in order)

### Step 1: Extract and standardise the BBB gene list ← START HERE
**What:** Open Daneman S3 and Munji S5. Extract the gene symbol columns. Convert mouse gene symbols to human HGNC format (e.g. *Slc2a1* → *SLC2A1*). Remove duplicates and merge into one master BBB gene list.
**Output:** `processed/master_BBB_genelist.csv` — columns: GeneSymbol_Mouse, GeneSymbol_Human, Source (Daneman/Munji/Both)
**Tools:** R with `readxl`, `dplyr`. For mouse→human conversion: `biomaRt` package.
**Why:** This is the foundation. Every downstream step depends on having a clean, standardised list of BBB genes.

### Step 2: Process the Wälchli RDS file
**What:** Load the Seurat object in R, run `AverageExpression()` to get average expression per cluster, save as CSV.
**Output:** `processed/Walchli2024_AverageExpression.csv`
**Tools:** R with `Seurat`
**Why:** The RDS file is not human-readable. This converts it into a usable expression matrix.

### Step 3: Map BBB genes to Ensembl IDs across all three species
**What:** Take the human BBB gene symbols from Step 1 and use BioMart to find their mouse (*Mus musculus*) and macaque (*Macaca mulatta*) orthologues. Record Ensembl gene IDs for all three species.
**Output:** `processed/BBB_genes_crossspecies_ensembl.csv` — columns: HumanSymbol, HumanEnsembl, MouseEnsembl, MacaqueEnsembl
**Tools:** R with `biomaRt`
**Why:** Ensembl IDs are the universal identifier that lets us pull genome coordinates for each gene in each species.

### Step 4: Download Ensembl GTF annotation files
**What:** Download the gene annotation files (GTF format) for human (GRCh38), macaque (Mmul_10), and mouse (GRCm39) from Ensembl. These list every gene's chromosome, start, end, exons, introns, and UTRs.
**Output:** Saved in a new `genomes/` folder (not committed to git — too large)
**Tools:** `wget` or R's `AnnotationHub`
**Why:** GTF files let us extract exact genomic coordinates for each BBB gene in each species, which is needed for the sequence comparison in Step 5.

### Step 5: Define and get approval for the control gene set
**What:** Randomly sample ~3x as many genes as are in the BBB gene list, matched for gene length and chromosome distribution. Prepare a brief proposal document for Dr. Clelland.
**Output:** `processed/control_genelist_proposal.csv` + short summary for Dr. Clelland
**Tools:** R with `dplyr`
**Why:** Without a control, we can't say whether BBB gene conservation is *unusual*. Dr. Clelland must approve this before the comparison runs.

---

## Eval Function — How to Know If the Work Is Correct

Run this R script after completing Step 1 to validate the master BBB gene list:

```r
# eval_step1.R — Validate the master BBB gene list
library(readr)
library(dplyr)

eval_master_BBB_genelist <- function(filepath = "processed/master_BBB_genelist.csv") {

  cat("=== EVAL: Master BBB Gene List ===\n\n")
  df <- read_csv(filepath, show_col_types = FALSE)

  # 1. Check required columns exist
  required_cols <- c("GeneSymbol_Mouse", "GeneSymbol_Human", "Source")
  missing_cols <- setdiff(required_cols, colnames(df))
  if (length(missing_cols) > 0) {
    cat("FAIL: Missing columns:", paste(missing_cols, collapse = ", "), "\n")
  } else {
    cat("PASS: All required columns present\n")
  }

  # 2. Check gene count is in expected range
  n <- nrow(df)
  cat(sprintf("INFO: Total genes in master list: %d\n", n))
  if (n < 300) {
    cat("WARN: Fewer than 300 genes — Daneman S3 has 213 and Munji S5 has 518, expect ~400-600 after deduplication\n")
  } else if (n > 800) {
    cat("WARN: More than 800 genes — possible duplication issue\n")
  } else {
    cat("PASS: Gene count looks reasonable\n")
  }

  # 3. Check for duplicates
  n_dupes <- sum(duplicated(df$GeneSymbol_Human))
  if (n_dupes > 0) {
    cat(sprintf("FAIL: %d duplicate human gene symbols found — deduplication needed\n", n_dupes))
  } else {
    cat("PASS: No duplicate human gene symbols\n")
  }

  # 4. Check human gene symbols are uppercase (HGNC format)
  n_lowercase <- sum(df$GeneSymbol_Human != toupper(df$GeneSymbol_Human), na.rm = TRUE)
  if (n_lowercase > 0) {
    cat(sprintf("WARN: %d gene symbols are not uppercase — may not be in HGNC format\n", n_lowercase))
  } else {
    cat("PASS: All human gene symbols are uppercase (HGNC format)\n")
  }

  # 5. Check Source column has expected values
  valid_sources <- c("Daneman", "Munji", "Both")
  invalid_sources <- df %>% filter(!Source %in% valid_sources) %>% nrow()
  if (invalid_sources > 0) {
    cat(sprintf("WARN: %d rows have unexpected Source values\n", invalid_sources))
  } else {
    cat("PASS: Source column values are valid\n")
  }

  # 6. Check that known landmark BBB genes are present
  landmark_genes <- c("SLC2A1", "CLDN5", "ABCB1", "TJP1", "PECAM1", "VWF")
  missing_landmarks <- setdiff(landmark_genes, df$GeneSymbol_Human)
  if (length(missing_landmarks) > 0) {
    cat(sprintf("WARN: These known BBB landmark genes are missing: %s\n",
                paste(missing_landmarks, collapse = ", ")))
    cat("      (This may be okay — not all landmarks appear in every study)\n")
  } else {
    cat("PASS: All landmark BBB genes present (SLC2A1, CLDN5, ABCB1, TJP1, PECAM1, VWF)\n")
  }

  # 7. Summary table by source
  cat("\n--- Gene counts by source ---\n")
  print(df %>% count(Source))

  cat("\n=== EVAL COMPLETE ===\n")
}

# Run it
eval_master_BBB_genelist()
```

Save this as `scripts/eval_step1.R` and run it after Step 1 is done. All PASS = good to proceed to Step 2.

---

## Key Decisions Already Made
- **Daneman S3** is the primary mouse BBB gene list (most stringent: enriched vs BOTH liver and lung)
- **Munji S5** is the secondary/updated mouse BBB gene list to merge with Daneman
- **Winkler S10** and **Yang S4** are the primary human reference datasets
- **Wälchli RDS** needs to be processed first before it's usable
- **Giger** provides the evolutionary conservation baseline — not a gene list source
- **Large .rds.gz files** are excluded from git via .gitignore — they stay local only
- **Control group** requires Dr. Clelland approval before analysis proceeds

---

## Still Pending / Blockers
- [x] Push repo to GitHub — done at github.com/Mmynemious/BBB-Gene-Conservation
- [x] Step 1 complete — 545-gene master BBB gene list created (Daneman S3 + Munji S5)
- [ ] Dr. Clelland approval of control gene set strategy (Step 5)
- [ ] Andrew Yang contact (potential additional dataset)

---

## Future Enhancement: Macaque Input Dataset

**Decision:** Core analysis proceeds without a macaque-specific input dataset. BioMart ortholog mapping (Step 3) handles the mouse→macaque gene mapping computationally. A macaque expression dataset can be added in a later phase.

**Best candidate when ready:**
- **Paper:** "A single-cell multi-omic atlas spanning the adult rhesus macaque brain" — Science Advances, 2023
- **Why it's good:** 1 million nuclei from 28 brain regions, explicitly identifies endothelial cells and vascular subtypes, directly compared to Wälchli and Winkler human atlases
- **Data:** Available on CellxGene — https://cellxgene.cziscience.com/collections/8c4bcf0d-b4df-45c7-888c-74fb0013e9e7
- **What it would add:** Any BBB gene confirmed in mouse (Daneman/Munji) + human (Wälchli/Winkler) + macaque (this atlas) = very high-confidence conserved BBB gene
- **Action required:** Flag to Dr. Clelland before adding — it expands project scope

**Claude Code prompt to add this dataset later** (save and use when ready):

> "I want to add a macaque brain endothelial dataset to the BBB gene conservation project as an optional input layer. The dataset is the rhesus macaque brain single-cell atlas from Science Advances 2023, available on CellxGene at https://cellxgene.cziscience.com/collections/8c4bcf0d-b4df-45c7-888c-74fb0013e9e7. Please: (1) download the endothelial cell subset from CellxGene, (2) extract average expression per cluster using Seurat, (3) filter to genes that overlap with our master BBB gene list at processed/master_BBB_genelist.csv, and (4) add a new column to the master list called Macaque_Expressed (TRUE/FALSE) indicating whether each BBB gene is expressed in macaque brain endothelium above a reasonable threshold. Save the updated list as processed/master_BBB_genelist_with_macaque.csv. Explain each step conceptually as you go — I am a beginner in R and bioinformatics."

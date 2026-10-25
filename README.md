# BBB Gene Conservation Analysis
**Clelland Lab Collaboration — Yara**

## Project Goal
Do blood-brain barrier (BBB) transporters that differ in **expression** between human, rhesus macaque and mouse also differ in **coding sequence**? The project focuses on the solute carrier (SLC) and ATP binding cassette (ABC) transporter families, which decide whether a drug crosses the barrier.

Project website with the literature review: <https://bbb-gene-conservation.vercel.app>

The first pipeline, its results and its audit are archived in [`archive/first_pipeline/`](archive/first_pipeline/README.md).

## Current plan
1. Replicate Miao et al. 2026 (capillary endothelial expression, mouse vs human) as a calibration.
2. Re-run sequence conservation on the transporter set, reusing the clean steps of the first pipeline and applying the audit fixes.
3. Link expression and coding-sequence differences gene by gene.
4. Write the introduction.

---

## Repository Structure

```
BBB-Gene-Conservation/
├── raw_data/           # Original unmodified datasets from each paper
│   ├── daneman_2010/
│   ├── munji_2019/
│   ├── winkler_2022/
│   ├── yang_2022/
│   ├── walchli_2024/
│   └── hpa_liver/      # Human Protein Atlas, used for the first pipeline's control set
├── processed/          # Outputs of the first pipeline
├── scripts/            # R scripts, step1 to step5h, plus eval_step* audit scripts
├── archive/            # First pipeline: audit, results, notes
├── site/               # Quarto sources for the project website
└── docs/               # The rendered website (what Vercel serves) — no paper PDFs here
```

Paper PDFs are not stored in this repository; the website's Resources page links each paper by DOI.

---

## Datasets

### Daneman et al. (2010) — Mouse BBB Gene Reference
- **Species:** Mouse (*Mus musculus*)
- **Method:** Microarray — FACS-sorted brain, liver, lung endothelial cells
- **Key file:** `Daneman2010_S6` — the BBB-enriched gene list (`Slc2a1`, `Slco1a4`, `Slc7a5`, `Ocln`, `Abcb1a`)
- **⚠️ Correction:** the file named `Daneman2010_S3_CoreBBBGenes_BrainEC_Enriched.xls` is **not** core BBB genes — it is Daneman's pericyte-specific list. S4 and S5 are developmentally up-/down-regulated CNS vascular genes, not tissue contrasts. The pipeline used S3/S4/S5 in error; see the [audit](https://bbb-gene-conservation.vercel.app/audit.html).
- **Why it matters:** The foundational mouse BBB gene list; >1000 citations

### Munji et al. (2019) — Updated Mouse BBB Transcriptome
- **Species:** Mouse (*Mus musculus*)
- **Method:** RNA-seq — brain endothelial cells in health and disease models
- **Key file:** `Munji2019_S5_BBBEnrichedGenes.xlsx` — updated BBB-enriched gene list
- **Why it matters:** Modern RNA-seq confirmation and expansion of Daneman 2010; also defines a BBB dysfunction module across neurological diseases

### Winkler et al. (2022) — Human Brain Vasculature Single-Cell Atlas
- **Species:** Human (*Homo sapiens*)
- **Method:** Single-cell RNA-seq — 181,388 cells from adult cerebrovasculature
- **Key file:** `Winkler2022_S10_VascularCellState_GeneSets_UCell.xlsx` — curated gene sets per vascular cell state (artery, capillary, venule, vein, pericyte, SMC, fibroblast)
- **Why it matters:** Human reference atlas for BBB and vascular cell-type gene expression

### Yang et al. (2022) — Human Brain Vascular Atlas (Alzheimer's focus)
- **Species:** Human (*Homo sapiens*)
- **Method:** VINE-seq — 143,793 single-nucleus transcriptomes, hippocampus and cortex
- **Key file:** `Yang2022_S4_AverageExpression_PerCellType.xlsx` — average expression per cell type
- **Why it matters:** Provides human BBB gene data with Alzheimer's disease dimension; includes both healthy and AD samples

### Wälchli et al. (2024) — Human Brain Vasculature Single-Cell Atlas (Comprehensive)
- **Species:** Human (*Homo sapiens*)
- **Method:** Single-cell RNA-seq — 606,380 endothelial cells across development, adulthood and disease
- **Key file:** `Walchli2024_AdultBrain_TemporalLobe_SortedEndothelial_SeuratObject.rds.gz` — Seurat object, adult temporal lobe sorted ECs
- **Why it matters:** Most comprehensive human brain vascular atlas to date; covers fetal and adult across 5 CNS pathologies
- **Note:** Large file — requires R + Seurat to process

---

---

*Analysis lead: Yara*
*Supervisor: Dr. Clelland*

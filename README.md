# BBB Gene Conservation Analysis
**Clelland Lab Collaboration — Yara**

## Project Goal
How similar are the genomes of human, non-human primate (rhesus macaque), and mouse across blood-brain barrier (BBB) genes? This project compares BBB coding regions, non-coding regions, and regulatory regions across three species to assess evolutionary conservation.

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
│   └── giger_2010/
├── processed/          # Cleaned, standardized outputs (populated during analysis)
├── scripts/            # R and Python analysis scripts (populated during analysis)
└── docs/               # Research papers and project documentation
```

---

## Datasets

### Daneman et al. (2010) — Mouse BBB Gene Reference
- **Species:** Mouse (*Mus musculus*)
- **Method:** Microarray — FACS-sorted brain, liver, lung endothelial cells
- **Key file:** `Daneman2010_S3_CoreBBBGenes_BrainEC_Enriched.xls` — 213 core BBB genes enriched in brain ECs vs both liver and lung ECs
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

### Giger et al. (2010) — Primate Neuronal vs Endothelial Transcriptome Evolution
- **Species:** Human, Chimpanzee, Rhesus Macaque (*Homo sapiens, Pan troglodytes, Macaca mulatta*)
- **Method:** Microarray (Affymetrix HG-U133 Plus 2.0) — laser-capture microdissected neurons and endothelial cells
- **Key file:** `Giger2010_GSE12293_Human_NeuronEndothelial_ExpressionMatrix.txt` — GEO series matrix
- **Why it matters:** Shows endothelial genes are more conserved across primates than neuronal genes; provides baseline conservation rates

---

## What's Still Needed (Future Steps)

- Genome annotation files (GTF) for human (GRCh38), macaque (Mmul_10), mouse (GRCm39) — from Ensembl
- NCBI genome assemblies for sequence-level comparisons (very large files, downloaded when needed)
- Control gene set (to be defined after Dr. Clelland input)

---

## Status
- [x] Literature compiled
- [x] Raw datasets downloaded and organized
- [ ] BBB gene list standardization
- [ ] Wälchli RDS processing (R/Seurat)
- [ ] Cross-species gene ID mapping (BioMart)
- [ ] Sequence extraction and alignment
- [ ] Conservation analysis and statistics

---

*Analysis lead: Yara*
*Supervisor: Dr. Clelland*

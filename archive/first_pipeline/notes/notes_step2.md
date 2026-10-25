# Step 2 Notes — Process the Wälchli Seurat Object

## What Step 2 Did

The Wälchli 2024 file is a Seurat object — a specialised R data structure for single-cell RNA sequencing data. It contained expression measurements for **26,666 genes** across **76,125 individual human brain endothelial cells**. That is a 26,666 × 76,125 matrix — far too large to work with directly.

`AverageExpression()` collapsed that matrix by grouping cells according to their annotated cell type (stored in the `ECclusters` metadata column) and computing the mean expression of each gene within each group. The result is a 26,666 × 12 matrix — one column per cell type, one row per gene — saved as `processed/Walchli2024_AverageExpression.csv`.

---

## The 12 Cell Types Found

```
Angiogenic capillary  Arteriole       Artery
Capillary             EndoMT          Large artery
Large vein            Mitochondrial   Proliferating cell
Stem-to-EC            Vein            Venule
```

**Important caveat — Mitochondrial cluster:** The `Mitochondrial` cluster is almost certainly a technical artefact, not a real cell type. It represents cells with abnormally high mitochondrial gene expression, which is a standard indicator of low-quality or dying cells in single-cell sequencing data. This cluster should be excluded or treated with caution in any downstream expression analysis. Worth flagging to Dr. Clelland.

---

## Why We Did This — and Why Not for Daneman and Munji

Daneman and Munji already handed us a pre-filtered gene list — their papers ran the statistical analysis and gave us a yes/no answer per gene. There was nothing to average; we just extracted the gene symbols.

Wälchli, by contrast, gives us raw single-cell data. No pre-filtering has been done. To make it usable and comparable, we needed to collapse it to the same kind of summary format: a table of gene expression values per cell type.

Yang 2022 (S4) already provides average expression per cell type as a ready-to-use Excel file, so Step 2 only needed to be done for Wälchli.

---

## What This Data Is Used For

The Wälchli (and Yang) expression matrices are **context and validation**, not the primary comparison. They let us confirm:
- Which cell types express BBB genes most strongly (expected answer: capillary endothelial cells, since that is where the BBB lives)
- That the genes in our master list are genuinely active in the right cell types in human brain tissue

The **core scientific question** of this project is about evolutionary conservation at the **DNA sequence level** — how similar is the actual genomic sequence of each BBB gene across human, macaque, and mouse? That analysis runs on genome sequences (Steps 3–5), not on expression numbers. Expression data is supporting evidence, not the comparison itself.

---

## Technical Note — Double-Compressed File

The `.rds.gz` file turned out to be double-compressed (a gzip file inside a gzip file). Standard R approaches (`readRDS()`, `gzfile()`, `gzcon()`) all failed. The working solution was to pipe through gunzip twice:

```r
readRDS(pipe(paste("gunzip -c", shQuote(seurat_path), "| gunzip -c")))
```

This decompresses on-the-fly without creating a temporary file on disk.

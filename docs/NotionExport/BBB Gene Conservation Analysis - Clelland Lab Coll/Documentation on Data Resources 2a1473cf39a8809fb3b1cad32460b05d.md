# Documentation on Data Resources

> A technical doc for transparency
> 

---

- For the papers of Yang, Munji, and Daneman, I used the publicly available data sets which were listed under the materials section of their published research paper. The data sets were downloaded as Excel sheets and were not altered
- For the papers of Wälchi, Ethan, and Giger, datasets were accessed via the **UCSC Cell Browser** portals provided by the authors. For these datasets, the preprocessed `.rds` Seurat objects were downloaded locally (e.g., `adult-control-temporal-lobe-sorted.rds` for *Wälchli*)

Each dataset was processed in **RStudio (macOS, R version ≥ 4.2)** using the following repeated workflow:

```r
library(Seurat)

# Load preprocessed Seurat object
dataset <- readRDS("~/Downloads/<filename>.rds")

# Generate average expression per cluster (bulk-like profile)
avg <- AverageExpression(dataset)$RNA

# Save output for integration
write.csv(avg, "~/Desktop/<DatasetName>_Averages.csv")

cat("✅ Saved to Desktop: <DatasetName>_Averages.csv\n")

```

This code was executed identically for each single-cell dataset (Giger, Wälchli, and Winkler/Ethan) to produce comparable average-expression matrices suitable for cross-dataset comparison.
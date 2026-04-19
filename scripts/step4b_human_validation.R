# ============================================================
# Step 4b: Cross-reference master BBB gene list against three
#           independent human brain EC datasets
#
# Adds three TRUE/FALSE columns + a count column to master list:
#   In_Walchli_EC   — expressed in any EC cluster (Wälchli 2024)
#   In_Yang_EC      — marker gene in Yang 2022 S6 EC subtype list
#   In_Winkler_EC   — marker gene in Winkler 2022 S2 EC list
#   Human_BBB_Datasets_N — 0–3 (how many datasets confirm it)
#
# Output: processed/master_BBB_genelist_validated.csv
# ============================================================

setwd("~/Documents/Claude/Projects/BBB")

library(readr)
library(readxl)
library(dplyr)

master <- read_csv("processed/master_BBB_genelist.csv", show_col_types = FALSE)
cat("Master list loaded:", nrow(master), "genes\n")

# ── 1. Wälchli 2024 ──────────────────────────────────────────────────────────
# Full average-expression matrix across 12 EC clusters (all columns are EC types)
walchli_expr <- read_csv("processed/Walchli2024_AverageExpression.csv",
                          show_col_types = FALSE)

# First column is gene symbol (written as row names by R)
walchli_genes <- walchli_expr[[1]]

# A gene "counts" if it has any expression > 0 across all 12 clusters
expr_mat <- as.matrix(walchli_expr[, -1])
walchli_expressed <- walchli_genes[rowSums(expr_mat, na.rm = TRUE) > 0]

master <- master |>
  mutate(In_Walchli_EC = GeneSymbol_Human %in% walchli_expressed)

cat("Wälchli EC coverage:", sum(master$In_Walchli_EC), "/", nrow(master),
    sprintf("(%.1f%%)\n", 100 * mean(master$In_Walchli_EC)))

# ── 2. Yang 2022 S6 — Endothelial subtype markers ────────────────────────────
yang_ec <- read_xlsx("raw_data/yang_2022/Yang2022_S6_EndothelialSubtype_Markers.xlsx")
# All 112 rows are EC marker genes; Gene column holds human symbols
yang_genes <- unique(na.omit(yang_ec$Gene))

master <- master |>
  mutate(In_Yang_EC = GeneSymbol_Human %in% yang_genes)

cat("Yang EC coverage:    ", sum(master$In_Yang_EC), "/", nrow(master),
    sprintf("(%.1f%%)\n", 100 * mean(master$In_Yang_EC)))

# ── 3. Winkler 2022 S2 — Endothelial marker genes ────────────────────────────
# Row 1 is a title, row 2 is the real header → skip = 1
winkler_raw <- read_xlsx(
  "raw_data/winkler_2022/science.abi7377_tables_s1_to_s10/Winkler2022_S2_EndothelialMarkers.xlsx",
  skip = 1, col_names = TRUE
)
# Column 1 = gene, last column = cluster label
colnames(winkler_raw)[1] <- "gene"
colnames(winkler_raw)[ncol(winkler_raw)] <- "cluster"

winkler_ec <- winkler_raw |>
  filter(!is.na(gene), gene != "gene",  # remove any leftover header rows
         cluster == "EC")

winkler_genes <- unique(winkler_ec$gene)

master <- master |>
  mutate(In_Winkler_EC = GeneSymbol_Human %in% winkler_genes)

cat("Winkler EC coverage: ", sum(master$In_Winkler_EC), "/", nrow(master),
    sprintf("(%.1f%%)\n", 100 * mean(master$In_Winkler_EC)))

# ── 4. Composite count ───────────────────────────────────────────────────────
master <- master |>
  mutate(Human_BBB_Datasets_N = as.integer(In_Walchli_EC) +
                                 as.integer(In_Yang_EC) +
                                 as.integer(In_Winkler_EC))

cat("\nDataset confirmation counts:\n")
print(table(master$Human_BBB_Datasets_N))

cat(sprintf("\nGenes confirmed by all 3 datasets:   %d\n", sum(master$Human_BBB_Datasets_N == 3)))
cat(sprintf("Genes confirmed by 2+ datasets:      %d\n", sum(master$Human_BBB_Datasets_N >= 2)))
cat(sprintf("Genes with NO human EC evidence:     %d\n", sum(master$Human_BBB_Datasets_N == 0)))

# ── 5. Save ──────────────────────────────────────────────────────────────────
write_csv(master, "processed/master_BBB_genelist_validated.csv")
cat("\nSaved: processed/master_BBB_genelist_validated.csv\n")
cat("Run scripts/eval_step4b.R to validate.\n")

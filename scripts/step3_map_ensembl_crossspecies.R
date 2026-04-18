# ============================================================
# Step 3: Map BBB genes to Ensembl IDs across all three species
#
# Takes the master BBB gene list (human HGNC symbols) and uses
# Ensembl BioMart to find the Ensembl gene ID for each gene in:
#   - Homo sapiens (human)
#   - Mus musculus (mouse)
#   - Macaca mulatta (rhesus macaque)
#
# Output: processed/BBB_genes_crossspecies_ensembl.csv
# ============================================================

setwd("~/Documents/Claude/Projects/BBB")

if (!requireNamespace("biomaRt", quietly = TRUE)) {
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install("biomaRt")
}
if (!requireNamespace("readr", quietly = TRUE)) install.packages("readr")
if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")

library(biomaRt)
library(readr)
library(dplyr)

# ── STEP A: Load master BBB gene list ────────────────────────────────────────
master <- read_csv("processed/master_BBB_genelist.csv", show_col_types = FALSE)
human_symbols <- unique(master$GeneSymbol_Human)
cat("BBB genes to map:", length(human_symbols), "\n")

# ── STEP B: Connect to Ensembl human database ────────────────────────────────
# We query the human mart because human is our reference species.
# From there, biomaRt can look up orthologues in mouse and macaque directly.
# Trying mirrors in order — use whichever connects first.

cat("\nConnecting to Ensembl (requires internet)...\n")

human_mart <- tryCatch(
  useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl", mirror = "useast"),
  error = function(e) {
    cat("US East failed, trying US West...\n")
    tryCatch(
      useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl", mirror = "uswest"),
      error = function(e) {
        cat("US West failed, trying main...\n")
        useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl")
      }
    )
  }
)

cat("Connected.\n")

# ── STEP C: Three separate queries — biomaRt does not allow mixing gene
#            attributes and homologue attributes in a single query ──────────────

cat("\nQuery 1: human gene symbols → human Ensembl IDs...\n")
human_ids <- getBM(
  attributes = c("hgnc_symbol", "ensembl_gene_id"),
  filters    = "hgnc_symbol",
  values     = human_symbols,
  mart       = human_mart
)
cat("Human IDs returned:", nrow(human_ids), "\n")

cat("Query 2: human Ensembl IDs → mouse orthologue IDs...\n")
mouse_ids <- getBM(
  attributes = c("ensembl_gene_id", "mmusculus_homolog_ensembl_gene"),
  filters    = "ensembl_gene_id",
  values     = human_ids$ensembl_gene_id,
  mart       = human_mart
)
cat("Mouse orthologue rows returned:", nrow(mouse_ids), "\n")

cat("Query 3: human Ensembl IDs → macaque orthologue IDs...\n")
macaque_ids <- getBM(
  attributes = c("ensembl_gene_id", "mmulatta_homolog_ensembl_gene"),
  filters    = "ensembl_gene_id",
  values     = human_ids$ensembl_gene_id,
  mart       = human_mart
)
cat("Macaque orthologue rows returned:", nrow(macaque_ids), "\n")

# Join all three tables on human Ensembl ID
ensembl_table <- human_ids |>
  left_join(mouse_ids,   by = "ensembl_gene_id") |>
  left_join(macaque_ids, by = "ensembl_gene_id")

cat("Combined table rows:", nrow(ensembl_table), "\n")

# ── STEP D: Clean and rename columns ─────────────────────────────────────────
ensembl_clean <- ensembl_table |>
  rename(
    HumanSymbol    = hgnc_symbol,
    HumanEnsembl   = ensembl_gene_id,
    MouseEnsembl   = mmusculus_homolog_ensembl_gene,
    MacaqueEnsembl = mmulatta_homolog_ensembl_gene
  ) |>
  # Remove rows where we got no human Ensembl ID at all
  filter(!is.na(HumanEnsembl), HumanEnsembl != "") |>
  # Replace empty strings in orthologue columns with NA for clarity
  mutate(
    MouseEnsembl   = na_if(MouseEnsembl, ""),
    MacaqueEnsembl = na_if(MacaqueEnsembl, "")
  )

cat("\nAfter cleaning:", nrow(ensembl_clean), "rows\n")

# ── STEP E: Summarise coverage ────────────────────────────────────────────────
n_human   <- sum(!is.na(ensembl_clean$HumanEnsembl))
n_mouse   <- sum(!is.na(ensembl_clean$MouseEnsembl))
n_macaque <- sum(!is.na(ensembl_clean$MacaqueEnsembl))

cat("\n--- Ensembl ID coverage ---\n")
cat("Human Ensembl IDs found:  ", n_human, "\n")
cat("Mouse orthologues found:  ", n_mouse, "\n")
cat("Macaque orthologues found:", n_macaque, "\n")

# ── STEP F: Join back to master list (keep Source and Daneman_Filter) ─────────
output <- master |>
  left_join(ensembl_clean, by = c("GeneSymbol_Human" = "HumanSymbol")) |>
  select(GeneSymbol_Mouse, GeneSymbol_Human, Source, Daneman_Filter,
         HumanEnsembl, MouseEnsembl, MacaqueEnsembl) |>
  arrange(GeneSymbol_Human)

cat("\nFinal cross-species table:", nrow(output), "rows\n")
cat("\nFirst 5 rows:\n")
print(head(output, 5))

# ── STEP G: Save ─────────────────────────────────────────────────────────────
write_csv(output, "processed/BBB_genes_crossspecies_ensembl.csv")
cat("\nSaved: processed/BBB_genes_crossspecies_ensembl.csv\n")
cat("Now run scripts/eval_step3.R to validate.\n")

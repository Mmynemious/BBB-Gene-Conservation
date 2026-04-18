# ============================================================
# Step 1: Extract and standardise the master BBB gene list
# Run in VSCode: open this file, then Ctrl+Enter line by line
#                OR Ctrl+Shift+S to source the whole file
# ============================================================

# ── Set working directory to the BBB project folder ─────────────────────────
setwd("~/Documents/Claude/Projects/BBB")

# ── Install missing packages (only runs once) ────────────────────────────────
if (!requireNamespace("readxl",   quietly = TRUE)) install.packages("readxl")
if (!requireNamespace("dplyr",    quietly = TRUE)) install.packages("dplyr")
if (!requireNamespace("readr",    quietly = TRUE)) install.packages("readr")
if (!requireNamespace("homologene", quietly = TRUE)) install.packages("homologene")

library(readxl)
library(dplyr)
library(readr)
library(homologene)
select <- dplyr::select   # homologene masks select — force dplyr's version

# ── STEP A: Load Daneman 2010 S3 ─────────────────────────────────────────────
# This file lists 213 core BBB genes enriched in brain ECs vs BOTH liver and lung.
# The first row of the .xls is a merged title row — skip it with skip = 1.
# Real column names appear in row 2: Probe Set ID | Gene Title | Gene Symbol | ...

daneman_raw <- read_xls(
  "raw_data/daneman_2010/Daneman2010_S3_CoreBBBGenes_BrainEC_Enriched.xls",
  skip = 0,
  col_names = FALSE   # no proper header row — we'll select by position
)

cat("Daneman S3 rows:", nrow(daneman_raw), "\n")
cat("First 3 rows of column 3 (gene symbols):\n")
print(daneman_raw[1:3, 3])

# Column 3 is the gene symbol column (Probe Set | Gene Title | Gene Symbol | ...)
# Row 1 appears to be the first data row, so no skipping needed
daneman_genes <- daneman_raw |>
  pull(3) |>
  na.omit() |>
  (\(x) x[!x %in% c("", "Gene Symbol", "---")])() |>  # drop header text and unannotated probes
  unique()

cat("Daneman unique mouse gene symbols:", length(daneman_genes), "\n")

# ── STEP B: Load Munji 2019 S5 ───────────────────────────────────────────────
# This file has multiple sheets. The one we want is called "BBB-enriched".
# It lists 519 mouse BBB-enriched genes with their Ensembl IDs.

munji_raw <- read_xlsx(
  "raw_data/munji_2019/Munji2019_S5_BBBEnrichedGenes.xlsx",
  sheet = "BBB-enriched"
)

cat("Munji S5 columns:\n")
print(colnames(munji_raw))
cat("Rows:", nrow(munji_raw), "\n\n")

# Pull the gene symbol column (first column, labelled "Gene")
munji_genes <- munji_raw |>
  filter(!is.na(Gene...1), Gene...1 != "") |>
  pull(Gene...1) |>
  unique()

cat("Munji unique mouse gene symbols:", length(munji_genes), "\n")

# ── STEP C: Tag each gene with its source before converting ──────────────────
# We create one table that knows where each mouse symbol came from.
# If a gene appears in both papers, we label it "Both" — that's a strong signal.

source_table <- bind_rows(
  data.frame(GeneSymbol_Mouse = daneman_genes, Source = "Daneman",
             stringsAsFactors = FALSE),
  data.frame(GeneSymbol_Mouse = munji_genes,   Source = "Munji",
             stringsAsFactors = FALSE)
) |>
  group_by(GeneSymbol_Mouse) |>
  summarise(Source = if (n() > 1) "Both" else first(Source), .groups = "drop")

cat("\nAll unique mouse genes:", nrow(source_table), "\n")
cat("Source breakdown:\n")
print(table(source_table$Source))

# ── STEP D: Convert mouse symbols → human HGNC using biomaRt ─────────────────
# biomaRt queries the Ensembl database over the internet.
# We ask: "for each mouse gene symbol, what is the equivalent human gene symbol?"
# This is called finding orthologues (evolutionary equivalent genes across species).

cat("\nConverting mouse → human symbols using homologene (no internet needed)...\n")
# homologene uses the NCBI HomoloGene database, bundled locally in the package.
# Tax IDs: mouse = 10090, human = 9606

conversion_raw <- homologene(
  genes   = source_table$GeneSymbol_Mouse,
  inTax   = 10090,   # Mus musculus
  outTax  = 9606     # Homo sapiens
)

cat("Conversion table rows returned:", nrow(conversion_raw), "\n")

conversion_table <- data.frame(
  GeneSymbol_Mouse = conversion_raw[["10090"]],
  GeneSymbol_Human = conversion_raw[["9606"]],
  stringsAsFactors = FALSE
) |>
  filter(!is.na(GeneSymbol_Human), GeneSymbol_Human != "")

cat("Valid mouse→human pairs:", nrow(conversion_table), "\n")

# ── STEP E: Build the final master list ──────────────────────────────────────
# Join the source labels with the human gene symbols.
# Where one mouse gene maps to multiple human genes, we keep all rows (biologically valid).
# Where the same human gene symbol appears twice, keep the row with the best source label.

master <- source_table |>
  left_join(conversion_table, by = "GeneSymbol_Mouse") |>
  filter(!is.na(GeneSymbol_Human), GeneSymbol_Human != "") |>
  # If same human gene came from multiple mouse genes, prefer "Both" > "Daneman" > "Munji"
  arrange(factor(Source, levels = c("Both", "Daneman", "Munji"))) |>
  distinct(GeneSymbol_Human, .keep_all = TRUE) |>
  mutate(
    GeneSymbol_Human  = toupper(GeneSymbol_Human),
    # All Daneman genes in this list come from S3: enriched vs BOTH liver AND lung
    # S4 (liver only) and S5 (lung only) genes could be added in future
    Daneman_Filter = case_when(
      Source %in% c("Daneman", "Both") ~ "liver_and_lung",
      TRUE ~ NA_character_
    )
  ) |>
  select(GeneSymbol_Mouse, GeneSymbol_Human, Source, Daneman_Filter) |>
  arrange(GeneSymbol_Human)

cat("\nFinal master BBB gene list:", nrow(master), "genes\n")
cat("Source breakdown:\n")
print(table(master$Source))
cat("\nFirst 10 rows:\n")
print(head(master, 10)) 

# ── STEP F: Save output ───────────────────────────────────────────────────────
dir.create("processed", showWarnings = FALSE)
write_csv(master, "processed/master_BBB_genelist.csv")
cat("\nSaved: processed/master_BBB_genelist.csv\n")
cat("Now run scripts/eval_step1.R to validate.\n")



  
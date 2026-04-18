# ============================================================
# Step 1: Extract and standardise the master BBB gene list
# Run in VSCode: open this file, then Ctrl+Enter line by line
#                OR Ctrl+Shift+S to source the whole file
# ============================================================

# ── Set working directory to the BBB project folder ─────────────────────────
setwd("~/Documents/Claude/Projects/BBB")

# ── Install missing packages (only runs once) ────────────────────────────────
if (!requireNamespace("readxl",     quietly = TRUE)) install.packages("readxl")
if (!requireNamespace("dplyr",      quietly = TRUE)) install.packages("dplyr")
if (!requireNamespace("readr",      quietly = TRUE)) install.packages("readr")
if (!requireNamespace("homologene", quietly = TRUE)) install.packages("homologene")

library(readxl)
library(dplyr)
library(readr)
library(homologene)
select <- dplyr::select   # homologene masks select — force dplyr's version

# ── Helper: extract gene symbols from a Daneman .xls file ────────────────────
# All Daneman files share the same layout: col 3 = gene symbol, row 1 = header
read_daneman_genes <- function(path) {
  raw <- read_xls(path, skip = 0, col_names = FALSE)
  raw |>
    pull(3) |>
    na.omit() |>
    (\(x) x[!x %in% c("", "Gene Symbol", "---")])() |>
    unique()
}

# ── STEP A: Load all three Daneman files ─────────────────────────────────────
# S3 = enriched in brain ECs vs BOTH liver AND lung  → most stringent
# S4 = enriched in brain ECs vs liver only           → less stringent
# S5 = enriched in brain ECs vs lung only            → less stringent

cat("Loading Daneman S3 (liver AND lung)...\n")
daneman_s3 <- read_daneman_genes(
  "raw_data/daneman_2010/Daneman2010_S3_CoreBBBGenes_BrainEC_Enriched.xls"
)
cat("S3 unique genes:", length(daneman_s3), "\n")

cat("Loading Daneman S4 (liver only)...\n")
daneman_s4 <- read_daneman_genes(
  "raw_data/daneman_2010/Daneman2010_S4_BrainEC_vs_LiverEC_Enriched.xls"
)
cat("S4 unique genes:", length(daneman_s4), "\n")

cat("Loading Daneman S5 (lung only)...\n")
daneman_s5 <- read_daneman_genes(
  "raw_data/daneman_2010/Daneman2010_S5_BrainEC_vs_LungEC_Enriched.xls"
)
cat("S5 unique genes:", length(daneman_s5), "\n")

# ── STEP B: Assign Daneman_Filter to each Daneman gene ───────────────────────
# Priority: S3 (most stringent) overrides S4/S5 labels.
# A gene in both S4 and S5 (but not S3) also gets "liver_and_lung" —
# it passed both comparisons separately, equivalent to S3 logic.

all_daneman <- unique(c(daneman_s3, daneman_s4, daneman_s5))

daneman_table <- data.frame(
  GeneSymbol_Mouse = all_daneman,
  stringsAsFactors = FALSE
) |>
  mutate(
    in_s3 = GeneSymbol_Mouse %in% daneman_s3,
    in_s4 = GeneSymbol_Mouse %in% daneman_s4,
    in_s5 = GeneSymbol_Mouse %in% daneman_s5,
    Daneman_Filter = case_when(
      in_s3               ~ "liver_and_lung",   # passed both simultaneously
      in_s4 & in_s5       ~ "liver_and_lung",   # passed both separately
      in_s4               ~ "liver_only",
      in_s5               ~ "lung_only"
    )
  ) |>
  select(GeneSymbol_Mouse, Daneman_Filter)

cat("\nDaneman filter breakdown:\n")
print(table(daneman_table$Daneman_Filter))

# ── STEP C: Load Munji 2019 S5 ───────────────────────────────────────────────
cat("\nLoading Munji S5...\n")
munji_raw <- read_xlsx(
  "raw_data/munji_2019/Munji2019_S5_BBBEnrichedGenes.xlsx",
  sheet = "BBB-enriched"
)

munji_genes <- munji_raw |>
  filter(!is.na(Gene...1), Gene...1 != "") |>
  pull(Gene...1) |>
  unique()

cat("Munji unique mouse gene symbols:", length(munji_genes), "\n")

# ── STEP D: Build source table ────────────────────────────────────────────────
# Tag each gene as Daneman / Munji / Both, then attach the Daneman_Filter label.

source_table <- bind_rows(
  data.frame(GeneSymbol_Mouse = all_daneman, Source = "Daneman",
             stringsAsFactors = FALSE),
  data.frame(GeneSymbol_Mouse = munji_genes, Source = "Munji",
             stringsAsFactors = FALSE)
) |>
  group_by(GeneSymbol_Mouse) |>
  summarise(Source = if (n() > 1) "Both" else first(Source), .groups = "drop") |>
  left_join(daneman_table, by = "GeneSymbol_Mouse")

cat("\nAll unique mouse genes:", nrow(source_table), "\n")
cat("Source breakdown:\n")
print(table(source_table$Source))

# ── STEP E: Convert mouse → human HGNC using homologene ──────────────────────
cat("\nConverting mouse → human symbols (no internet needed)...\n")

conversion_raw <- homologene(
  genes  = source_table$GeneSymbol_Mouse,
  inTax  = 10090,   # Mus musculus
  outTax = 9606     # Homo sapiens
)

conversion_table <- data.frame(
  GeneSymbol_Mouse = conversion_raw[["10090"]],
  GeneSymbol_Human = conversion_raw[["9606"]],
  stringsAsFactors = FALSE
) |>
  filter(!is.na(GeneSymbol_Human), GeneSymbol_Human != "")

cat("Valid mouse→human pairs:", nrow(conversion_table), "\n")

# ── STEP F: Build the final master list ──────────────────────────────────────
master <- source_table |>
  left_join(conversion_table, by = "GeneSymbol_Mouse") |>
  filter(!is.na(GeneSymbol_Human), GeneSymbol_Human != "") |>
  # Where the same human gene appears multiple times, keep most informative row:
  # Both > Daneman > Munji; within Daneman, liver_and_lung > liver_only > lung_only
  arrange(
    factor(Source, levels = c("Both", "Daneman", "Munji")),
    factor(Daneman_Filter, levels = c("liver_and_lung", "liver_only", "lung_only"))
  ) |>
  distinct(GeneSymbol_Human, .keep_all = TRUE) |>
  mutate(GeneSymbol_Human = toupper(GeneSymbol_Human)) |>
  select(GeneSymbol_Mouse, GeneSymbol_Human, Source, Daneman_Filter) |>
  arrange(GeneSymbol_Human)

cat("\nFinal master BBB gene list:", nrow(master), "genes\n")
cat("Source breakdown:\n")
print(table(master$Source))
cat("\nDaneman_Filter breakdown:\n")
print(table(master$Daneman_Filter, useNA = "ifany"))
cat("\nFirst 10 rows:\n")
print(head(master, 10))

# ── STEP G: Save output ───────────────────────────────────────────────────────
dir.create("processed", showWarnings = FALSE)
write_csv(master, "processed/master_BBB_genelist.csv")
cat("\nSaved: processed/master_BBB_genelist.csv\n")
cat("Now run scripts/eval_step1.R to validate.\n")

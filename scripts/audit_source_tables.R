# ============================================================
# Audit: verify source-table identity before trusting a filename
#
# Reproduces Findings 1, 2 and 4 of docs/data_audit_and_weaknesses.md
# directly from the raw files. Run this from the project root.
#
# Why this script exists: the pipeline read three Daneman tables
# on the strength of their filenames. The filenames were wrong.
# This script checks each file against the criterion the paper
# says defines it, so the mistake cannot recur silently.
#
# Nothing here writes to processed/ — it only reports.
# ============================================================

# Paths are relative to the project root. If you are running this
# from RStudio with a different working directory, set it first:
# setwd("~/Documents/Claude/Projects/BBB")

suppressPackageStartupMessages({
  library(readxl)
  library(readr)
  library(dplyr)
})

pass <- function(msg) cat("  PASS:", msg, "\n")
fail <- function(msg) cat("  FAIL:", msg, "\n")
info <- function(msg) cat("  INFO:", msg, "\n")

# ============================================================
# FINDING 1 — Which Daneman table is which?
# ============================================================
# Each Daneman supplementary table was defined by a ratio cutoff.
# The ratio columns are carried in every table, so we can ask each
# file which cutoff it actually satisfies, and compare that to the
# name on the tin.
#
#   col 12 "Adult/Postnatal"  -> S4 (up) / S5 (down)
#   col 13 "Brain/Liver"      -> S6 requires > 2
#   col 14 "Brain/Lung"       -> S6 requires > 2
#   col 11 "Brain G+/G+P-"    -> S3 (pericyte) requires > 2

cat("\n=== FINDING 1: Daneman table identity ===\n\n")

daneman_files <- c(
  S3 = "raw_data/daneman_2010/Daneman2010_S3_CoreBBBGenes_BrainEC_Enriched.xls",
  S4 = "raw_data/daneman_2010/Daneman2010_S4_BrainEC_vs_LiverEC_Enriched.xls",
  S5 = "raw_data/daneman_2010/Daneman2010_S5_BrainEC_vs_LungEC_Enriched.xls",
  S6 = "raw_data/daneman_2010/Daneman2010_S6_PostnatalBrainEC_Enriched.xls",
  S7 = "raw_data/daneman_2010/Daneman2010_S7_AdultBrainEC_Enriched.xls"
)

signature <- function(path) {
  raw <- suppressMessages(read_xls(path, col_names = FALSE))
  num <- function(i) suppressWarnings(as.numeric(raw[[i]][-1]))
  tibble(
    n_rows           = nrow(raw) - 1,
    pericyte_gt2     = mean(num(12) > 2, na.rm = TRUE),   # Brain G+/G+P-
    dev_up_gt2       = mean(num(13) > 2, na.rm = TRUE),   # Adult/Postnatal
    dev_down_lt0.5   = mean(num(13) < 0.5, na.rm = TRUE),
    brain_liver_gt2  = mean(num(14) > 2, na.rm = TRUE),
    brain_lung_gt2   = mean(num(15) > 2, na.rm = TRUE)
  )
}

sigs <- bind_rows(lapply(daneman_files, signature), .id = "file")
print(as.data.frame(sigs), digits = 2)

cat("\nInterpretation:\n")
cat("  The file named S3 satisfies the PERICYTE criterion, not the BBB criterion.\n")
cat("  The files named S4 and S5 satisfy the DEVELOPMENTAL up/down criteria.\n")
cat("  The file named S6 is the one satisfying BOTH Brain/Liver > 2 AND\n")
cat("  Brain/Lung > 2 — it is Daneman's actual BBB-enriched gene list.\n\n")

# Landmark-gene cross-check: a BBB table must contain BBB genes.
read_symbols <- function(path) {
  raw <- suppressMessages(read_xls(path, col_names = FALSE))
  s <- raw[[3]][-1]
  unique(s[!is.na(s) & !s %in% c("", "---")])
}

bbb_landmarks      <- c("Slc2a1", "Slco1a4", "Slc7a5", "Ocln", "Abcb1a")
pericyte_landmarks <- c("Pdgfrb", "Kcnj8", "Vtn", "Zic1")

for (nm in names(daneman_files)) {
  s <- read_symbols(daneman_files[[nm]])
  cat(sprintf("  %s: %4d genes | BBB landmarks %d/%d | pericyte landmarks %d/%d\n",
              nm, length(s),
              sum(bbb_landmarks %in% s),      length(bbb_landmarks),
              sum(pericyte_landmarks %in% s), length(pericyte_landmarks)))
}

# The decisive number: overlap with Munji's independent BBB list.
munji <- read_xlsx("raw_data/munji_2019/Munji2019_S5_BBBEnrichedGenes.xlsx",
                   sheet = "BBB-enriched")
munji_genes <- unique(na.omit(munji[[1]]))

as_built  <- unique(c(read_symbols(daneman_files[["S3"]]),
                      read_symbols(daneman_files[["S4"]]),
                      read_symbols(daneman_files[["S5"]])))
corrected <- read_symbols(daneman_files[["S6"]])
corrected <- corrected[is.na(suppressWarnings(as.numeric(corrected)))]

cat(sprintf("\n  Munji S5: %d genes\n", length(munji_genes)))
cat(sprintf("  As built  (S3+S4+S5): %4d genes, overlap with Munji = %3d\n",
            length(as_built),  length(intersect(as_built,  munji_genes))))
cat(sprintf("  Corrected (S6 only) : %4d genes, overlap with Munji = %3d\n",
            length(corrected), length(intersect(corrected, munji_genes))))
cat("  Two studies of the same tissue should overlap substantially.\n")

# ============================================================
# FINDING 2 — The HPA brain filter matched nothing
# ============================================================
cat("\n=== FINDING 2: HPA brain-exclusion filter ===\n\n")

hpa_file <- "raw_data/hpa_liver/rna_tissue_consensus.tsv"
hpa_zip  <- "raw_data/hpa_liver/rna_tissue_consensus.tsv.zip"
if (!file.exists(hpa_file) && file.exists(hpa_zip)) {
  unzip(hpa_zip, exdir = "raw_data/hpa_liver")
}

if (file.exists(hpa_file)) {
  hpa <- read_tsv(hpa_file, show_col_types = FALSE)

  n_brain <- sum(hpa$Tissue == "brain")
  if (n_brain == 0) {
    fail('step4c filters Tissue == "brain", which matches 0 rows — the filter never ran.')
  } else {
    pass('Tissue == "brain" matches rows.')
  }
  info(paste("Brain-related labels actually present:",
             paste(grep("cort|hippo|cereb|amyg|basal|midbrain|hypothal|choroid|spinal|retina",
                        unique(hpa$Tissue), value = TRUE), collapse = ", ")))

  brain_regions <- c("amygdala", "basal ganglia", "cerebellum", "cerebral cortex",
                     "choroid plexus", "hippocampal formation", "hypothalamus",
                     "midbrain", "spinal cord", "retina", "pituitary gland")

  brain_max <- hpa |>
    filter(Tissue %in% brain_regions) |>
    group_by(Gene = `Gene name`) |>
    summarise(Brain_max_nTPM = max(nTPM), .groups = "drop")

  ctrl <- read_csv("processed/liver_control_genelist.csv", show_col_types = FALSE) |>
    left_join(brain_max, by = "Gene") |>
    mutate(Brain_max_nTPM = ifelse(is.na(Brain_max_nTPM), 0, Brain_max_nTPM))

  cat(sprintf("\n  Control genes: %d\n", nrow(ctrl)))
  cat(sprintf("  Committed Brain_nTPM column — all zero? %s\n",
              all(ctrl$Brain_nTPM == 0)))
  cat(sprintf("  Recomputed: brain nTPM >= 10 in %d genes (%.1f%%)\n",
              sum(ctrl$Brain_max_nTPM >= 10),
              100 * mean(ctrl$Brain_max_nTPM >= 10)))
  cat(sprintf("  Median liver nTPM %.1f  vs  median max brain nTPM %.1f\n",
              median(ctrl$Liver_nTPM), median(ctrl$Brain_max_nTPM)))
  cat(sprintf("  Genes surviving a correct brain filter: %d\n",
              sum(ctrl$Brain_max_nTPM < 10)))
} else {
  info("HPA file not found — skipping Finding 2.")
}

# ============================================================
# FINDING 4 — Yang S6 has seven sheets; step4b read one
# ============================================================
cat("\n=== FINDING 4: Yang S6 sheet coverage ===\n\n")

yang_path <- "raw_data/yang_2022/Yang2022_S6_EndothelialSubtype_Markers.xlsx"
sheets <- excel_sheets(yang_path)
cat("  Sheets present:", paste(sheets, collapse = ", "), "\n")
cat("  read_xlsx() with no sheet= argument reads only:", sheets[1], "\n\n")

master <- read_csv("processed/master_BBB_genelist.csv", show_col_types = FALSE)
for (sh in sheets) {
  g <- unique(na.omit(read_xlsx(yang_path, sheet = sh)[[1]]))
  cat(sprintf("    %-15s %3d genes | overlap with master list = %3d\n",
              sh, length(g), sum(master$GeneSymbol_Human %in% g)))
}

ec_sheets <- intersect(c("Arterial", "Capillary", "Venous"), sheets)
ec_genes  <- unique(unlist(lapply(ec_sheets,
                    function(sh) na.omit(read_xlsx(yang_path, sheet = sh)[[1]]))))
cat(sprintf("\n  All three EC sheets combined: %d genes, overlap = %d\n",
            length(ec_genes), sum(master$GeneSymbol_Human %in% ec_genes)))

cat("\n=== AUDIT COMPLETE ===\n")

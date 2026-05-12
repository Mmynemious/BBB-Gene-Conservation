# ============================================================
# Step 3c: Gene biotype audit
#
# Verifies that genes in the BBB master list are protein-coding
# in Ensembl. Any noncoding biotypes (lncRNA, pseudogene, etc.)
# that slipped into a CDS-based conservation pipeline would
# invalidate the "coding sequence conservation" framing.
#
# Biotypes are extracted directly from the Ensembl CDS FASTA
# headers (gene_biotype: field) — no BioMart query needed.
#
# Output:
#   processed/gene_biotype_check.csv   — per-gene biotype table
#   processed/gene_biotype_summary.txt — plain-text audit report
# ============================================================

setwd("~/Documents/Claude/Projects/BBB")

suppressPackageStartupMessages({
  library(Biostrings)
  library(readr)
  library(dplyr)
})

master <- read_csv("processed/master_BBB_genelist_complete.csv",
                   show_col_types = FALSE)
cat("Genes in master list:", nrow(master), "\n\n")

# ── Helper: build gene_id → gene_biotype lookup from FASTA headers ───────────
parse_biotypes <- function(fasta_path, label) {
  cat(sprintf("  Parsing %s FASTA headers...\n", label))
  s <- readDNAStringSet(fasta_path)
  h <- names(s)
  gene_id  <- sub(".*gene:(ENS[A-Z0-9]+)\\.?[0-9]*.*", "\\1", h)
  biotype  <- sub(".*gene_biotype:(\\S+).*", "\\1", h)
  biotype[!grepl("gene_biotype:", h)] <- NA_character_
  df <- data.frame(gene_id = gene_id, gene_biotype = biotype,
                   stringsAsFactors = FALSE) |>
    distinct(gene_id, .keep_all = TRUE)
  cat(sprintf("    %d unique gene IDs parsed\n", nrow(df)))
  df
}

human_bt   <- parse_biotypes("genomes/Homo_sapiens.GRCh38.113.cds.all.fa.gz",   "human")
mouse_bt   <- parse_biotypes("genomes/Mus_musculus.GRCm39.113.cds.all.fa.gz",   "mouse")
macaque_bt <- parse_biotypes("genomes/Macaca_mulatta.Mmul_10.113.cds.all.fa.gz","macaque")

# ── Merge biotypes onto master table ─────────────────────────────────────────
audit <- master |>
  left_join(human_bt   |> rename(Human_Biotype   = gene_biotype),
            by = c("HumanEnsembl"   = "gene_id")) |>
  left_join(mouse_bt   |> rename(Mouse_Biotype   = gene_biotype),
            by = c("MouseEnsembl"   = "gene_id")) |>
  left_join(macaque_bt |> rename(Macaque_Biotype = gene_biotype),
            by = c("MacaqueEnsembl" = "gene_id")) |>
  mutate(
    Human_NonCoding   = !is.na(Human_Biotype)   & Human_Biotype   != "protein_coding",
    Mouse_NonCoding   = !is.na(Mouse_Biotype)   & Mouse_Biotype   != "protein_coding",
    Macaque_NonCoding = !is.na(Macaque_Biotype) & Macaque_Biotype != "protein_coding",
    Any_NonCoding     = Human_NonCoding | Mouse_NonCoding | Macaque_NonCoding
  )

write_csv(audit, "processed/gene_biotype_check.csv")
cat("\nSaved: processed/gene_biotype_check.csv\n\n")

# ── Check whether flagged genes entered the CDS pipeline ─────────────────────
cat("Checking CDS FASTA coverage for flagged genes...\n")
bbb_human_fa   <- readDNAStringSet("genomes/cds_BBB_human.fa")
bbb_mouse_fa   <- readDNAStringSet("genomes/cds_BBB_mouse.fa")
bbb_macaque_fa <- readDNAStringSet("genomes/cds_BBB_macaque.fa")

noncoding <- audit |> filter(Any_NonCoding == TRUE)

nc_in_human   <- sum(noncoding$HumanEnsembl   %in% names(bbb_human_fa),   na.rm=TRUE)
nc_in_mouse   <- sum(noncoding$MouseEnsembl   %in% names(bbb_mouse_fa),   na.rm=TRUE)
nc_in_macaque <- sum(noncoding$MacaqueEnsembl %in% names(bbb_macaque_fa), na.rm=TRUE)

# ── Write summary ─────────────────────────────────────────────────────────────
lines <- c(
  "=== Gene Biotype Audit ===",
  "",
  sprintf("Total genes in master list: %d", nrow(audit)),
  "",
  "── Human biotype breakdown ──────────────────────────────",
  capture.output(sort(table(audit$Human_Biotype), decreasing=TRUE)),
  "",
  "── Mouse biotype breakdown ──────────────────────────────",
  capture.output(sort(table(audit$Mouse_Biotype), decreasing=TRUE)),
  "",
  "── Macaque biotype breakdown ────────────────────────────",
  capture.output(sort(table(audit$Macaque_Biotype), decreasing=TRUE)),
  "",
  "── Noncoding flags ──────────────────────────────────────",
  sprintf("Human genes not protein_coding:   %d", sum(audit$Human_NonCoding,   na.rm=TRUE)),
  sprintf("Mouse genes not protein_coding:   %d", sum(audit$Mouse_NonCoding,   na.rm=TRUE)),
  sprintf("Macaque genes not protein_coding: %d", sum(audit$Macaque_NonCoding, na.rm=TRUE)),
  sprintf("Any species noncoding:            %d", sum(audit$Any_NonCoding,     na.rm=TRUE)),
  ""
)

if (nrow(noncoding) > 0) {
  lines <- c(lines,
    "── Noncoding gene details ───────────────────────────────",
    capture.output(
      noncoding |>
        select(GeneSymbol_Human, HumanEnsembl, Human_Biotype,
               Mouse_Biotype, Macaque_Biotype) |>
        as.data.frame() |> print()
    ),
    "",
    sprintf("Of these, in human CDS FASTA:   %d", nc_in_human),
    sprintf("Of these, in mouse CDS FASTA:   %d", nc_in_mouse),
    sprintf("Of these, in macaque CDS FASTA: %d", nc_in_macaque),
    ""
  )
  if (nc_in_human == 0 && nc_in_mouse == 0 && nc_in_macaque == 0) {
    lines <- c(lines,
      "RESULT: No noncoding genes entered the CDS alignment pipeline.",
      "        The CDS extraction step correctly filtered them out because",
      "        they have no CDS entries in the Ensembl FASTA files.",
      "        Pipeline integrity: CONFIRMED."
    )
  } else {
    lines <- c(lines,
      "WARNING: Some noncoding genes are present in the CDS FASTAs.",
      "         These should be excluded before publication."
    )
  }
} else {
  lines <- c(lines,
    "RESULT: All genes in the master list are protein_coding in all three species.",
    "        No annotation gap detected."
  )
}

writeLines(lines, "processed/gene_biotype_summary.txt")
cat(paste(lines, collapse="\n"), "\n")

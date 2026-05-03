# ============================================================
# Step 5b alternative: Use MANE Select / Ensembl Canonical
# transcripts instead of the longest-transcript heuristic
#
# For each gene, Ensembl designates one transcript as
# "canonical." For human, this matches the MANE Select
# annotation (jointly curated by NCBI + EMBL-EBI). For mouse
# and macaque, Ensembl Canonical is its own designation.
#
# This script:
#   1. Queries Ensembl BioMart for canonical transcript IDs
#      per gene in each species
#   2. Filters the CDS FASTA files to retain only those
#      canonical transcripts
#   3. Saves canonical-only FASTAs that can be swapped in
#      place of the longest-transcript FASTAs in Step 5c
#
# To use: replace step5b extraction outputs with the canonical
# versions and re-run alignment scripts.
#
# Output:
#   processed/canonical_tx_human.csv
#   processed/canonical_tx_mouse.csv
#   processed/canonical_tx_macaque.csv
#   genomes/cds_canonical_human.fa
#   genomes/cds_canonical_mouse.fa
#   genomes/cds_canonical_macaque.fa
# ============================================================

setwd("~/Documents/Claude/Projects/BBB")

suppressPackageStartupMessages({
  library(biomaRt)
  library(Biostrings)
  library(readr)
  library(dplyr)
})

# ── Connect to Ensembl ───────────────────────────────────────────────────────
cat("Connecting to Ensembl BioMart (main server)...\n")
human_mart   <- useEnsembl("genes", dataset = "hsapiens_gene_ensembl",
                            host = "https://www.ensembl.org")
mouse_mart   <- useEnsembl("genes", dataset = "mmusculus_gene_ensembl",
                            host = "https://www.ensembl.org")
macaque_mart <- useEnsembl("genes", dataset = "mmulatta_gene_ensembl",
                            host = "https://www.ensembl.org")
cat("Connected.\n\n")

# ── Helper: pull canonical transcript IDs for a species ─────────────────────
get_canonical <- function(mart, label, cache_file) {
  if (file.exists(cache_file)) {
    cat("  Loading cached", label, "canonical transcripts\n")
    return(read_csv(cache_file, show_col_types = FALSE))
  }
  cat("  Querying canonical transcripts for", label, "...\n")
  res <- getBM(
    attributes = c("ensembl_gene_id", "ensembl_transcript_id",
                    "transcript_is_canonical"),
    mart = mart
  )
  out <- res |>
    filter(!is.na(transcript_is_canonical), transcript_is_canonical == 1) |>
    select(gene_id = ensembl_gene_id, canonical_tx = ensembl_transcript_id) |>
    distinct()
  write_csv(out, cache_file)
  cat("  ", nrow(out), label, "canonical transcripts saved\n")
  out
}

human_canonical   <- get_canonical(human_mart,   "human",
                                    "processed/canonical_tx_human.csv")
mouse_canonical   <- get_canonical(mouse_mart,   "mouse",
                                    "processed/canonical_tx_mouse.csv")
macaque_canonical <- get_canonical(macaque_mart, "macaque",
                                    "processed/canonical_tx_macaque.csv")

# ── Filter CDS FASTA to canonical transcripts only ────────────────────────────
filter_canonical <- function(fasta_path, canonical_table, out_path, label) {
  cat("\nFiltering", label, "CDS FASTA...\n")
  s <- readDNAStringSet(fasta_path)
  # First "word" of header is the transcript ID with version (e.g. ENST00000390473.1)
  tx_ids_raw <- sub(" .*$", "", names(s))
  tx_ids     <- sub("\\.[0-9]+$", "", tx_ids_raw)  # strip version
  gene_ids   <- sub(".*gene:(ENS[A-Z0-9]+)\\.?[0-9]*.*", "\\1", names(s))

  keep_idx <- tx_ids %in% canonical_table$canonical_tx
  out <- s[keep_idx]
  names(out) <- gene_ids[keep_idx]  # rename by gene ID for easy downstream lookup

  writeXStringSet(out, out_path)
  cat("  Kept", length(out), "/", length(s), "transcripts (canonical only)\n")
  cat("  Saved to:", out_path, "\n")
  invisible(out)
}

filter_canonical("genomes/Homo_sapiens.GRCh38.113.cds.all.fa.gz",
                  human_canonical,
                  "genomes/cds_canonical_human.fa", "Human")
filter_canonical("genomes/Mus_musculus.GRCm39.113.cds.all.fa.gz",
                  mouse_canonical,
                  "genomes/cds_canonical_mouse.fa", "Mouse")
filter_canonical("genomes/Macaca_mulatta.Mmul_10.113.cds.all.fa.gz",
                  macaque_canonical,
                  "genomes/cds_canonical_macaque.fa", "Macaque")

cat("\nStep 5b alternative complete.\n")
cat("To use these instead of the longest-transcript FASTAs:\n")
cat("  - Edit step5c/d/h to point to genomes/cds_canonical_*.fa\n")
cat("  - Re-run the alignment pipeline\n")

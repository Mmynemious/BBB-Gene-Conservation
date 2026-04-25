# ============================================================
# Step 5b: Extract CDS sequences for BBB genes and liver control
#
# Reads the Ensembl CDS FASTA files and pulls the coding sequence
# for each gene in our master list and liver control.
# For genes with multiple transcripts, keeps the longest CDS
# (most complete representation of the gene).
#
# Output: genomes/cds_BBB_human.fa, cds_BBB_mouse.fa, cds_BBB_macaque.fa
#         genomes/cds_liver_human.fa, cds_liver_mouse.fa, cds_liver_macaque.fa
# ============================================================

setwd("~/Documents/Claude/Projects/BBB")

library(Biostrings)
library(readr)
library(dplyr)

# ── Helper: parse CDS FASTA and index by gene ID ─────────────────────────────
load_cds_by_gene <- function(fasta_path) {
  cat("  Loading", basename(fasta_path), "...\n")
  seqs <- readDNAStringSet(fasta_path)

  # Extract gene ID from header: "gene:ENSG00000123456.1" → "ENSG00000123456"
  headers <- names(seqs)
  gene_ids <- sub(".*gene:(ENS[A-Z0-9]+)\\.?[0-9]*.*", "\\1", headers)

  # For each gene, keep the longest CDS transcript
  df <- data.frame(
    gene_id = gene_ids,
    length  = width(seqs),
    idx     = seq_along(seqs),
    stringsAsFactors = FALSE
  )

  best <- df |>
    group_by(gene_id) |>
    slice_max(length, n = 1, with_ties = FALSE) |>
    ungroup()

  seqs_by_gene <- seqs[best$idx]
  names(seqs_by_gene) <- best$gene_id

  cat("  ", length(seqs_by_gene), "unique genes loaded\n")
  return(seqs_by_gene)
}

# ── Load all three species CDS ────────────────────────────────────────────────
cat("Loading CDS FASTA files...\n")
human_cds   <- load_cds_by_gene("genomes/Homo_sapiens.GRCh38.113.cds.all.fa.gz")
mouse_cds   <- load_cds_by_gene("genomes/Mus_musculus.GRCm39.113.cds.all.fa.gz")
macaque_cds <- load_cds_by_gene("genomes/Macaca_mulatta.Mmul_10.113.cds.all.fa.gz")

# ── Load gene lists ───────────────────────────────────────────────────────────
master <- read_csv("processed/master_BBB_genelist_complete.csv",
                   show_col_types = FALSE) |>
  filter(!is.na(HumanEnsembl))

liver <- read_csv("processed/liver_control_genelist.csv",
                  show_col_types = FALSE)

# ── Function: extract sequences for a gene set ────────────────────────────────
extract_sequences <- function(gene_table, human_cds, mouse_cds, macaque_cds,
                              label = "gene set") {
  cat("\nExtracting sequences for", label, "...\n")

  human_ids   <- unique(na.omit(gene_table$HumanEnsembl))
  mouse_ids   <- unique(na.omit(gene_table$MouseEnsembl))
  macaque_ids <- unique(na.omit(gene_table$MacaqueEnsembl))

  h <- human_cds[names(human_cds) %in% human_ids]
  m <- mouse_cds[names(mouse_cds) %in% mouse_ids]
  q <- macaque_cds[names(macaque_cds) %in% macaque_ids]

  cat(sprintf("  Human:   %d / %d genes with CDS\n", length(h), length(human_ids)))
  cat(sprintf("  Mouse:   %d / %d genes with CDS\n", length(m), length(mouse_ids)))
  cat(sprintf("  Macaque: %d / %d genes with CDS\n", length(q), length(macaque_ids)))

  list(human = h, mouse = m, macaque = q)
}

# ── BBB gene sequences ────────────────────────────────────────────────────────
bbb_seqs <- extract_sequences(master, human_cds, mouse_cds, macaque_cds,
                               label = "BBB genes")

# For liver we need Ensembl IDs — get them from biomaRt or use gene symbols
# For now, match liver genes by human gene symbol via the master list mapping
# We'll get liver Ensembl IDs from the human CDS headers directly
cat("\nMatching liver gene symbols to Ensembl IDs via CDS headers...\n")

# Parse gene symbols from human CDS FASTA headers
human_fasta_raw <- readDNAStringSet("genomes/Homo_sapiens.GRCh38.113.cds.all.fa.gz")
headers <- names(human_fasta_raw)
symbol_map <- data.frame(
  gene_id = sub(".*gene:(ENS[A-Z0-9]+)\\.?[0-9]*.*", "\\1", headers),
  symbol  = sub(".*gene_symbol:([^ ]+).*", "\\1", headers),
  length  = width(human_fasta_raw),
  stringsAsFactors = FALSE
) |>
  filter(grepl("^ENS", gene_id), !grepl("gene_symbol:", symbol)) |>
  group_by(symbol) |>
  slice_max(length, n = 1, with_ties = FALSE) |>
  ungroup()

liver_matched <- liver |>
  left_join(symbol_map |> select(Gene = symbol, HumanEnsembl = gene_id),
            by = "Gene") |>
  filter(!is.na(HumanEnsembl))

cat(sprintf("Liver genes with Ensembl IDs: %d / %d\n",
    nrow(liver_matched), nrow(liver)))

# Extract liver human sequences (liver control is human-only for now —
# we compare the SAME liver genes across species using orthologue IDs from biomaRt)
liver_human_ids <- unique(liver_matched$HumanEnsembl)
liver_human_seqs <- human_cds[names(human_cds) %in% liver_human_ids]
cat(sprintf("Liver human CDS extracted: %d genes\n", length(liver_human_seqs)))

# ── Save FASTA files ──────────────────────────────────────────────────────────
cat("\nSaving FASTA files...\n")
writeXStringSet(bbb_seqs$human,   "genomes/cds_BBB_human.fa")
writeXStringSet(bbb_seqs$mouse,   "genomes/cds_BBB_mouse.fa")
writeXStringSet(bbb_seqs$macaque, "genomes/cds_BBB_macaque.fa")
writeXStringSet(liver_human_seqs, "genomes/cds_liver_human.fa")

cat("Saved:\n")
cat("  genomes/cds_BBB_human.fa   —", length(bbb_seqs$human), "genes\n")
cat("  genomes/cds_BBB_mouse.fa   —", length(bbb_seqs$mouse), "genes\n")
cat("  genomes/cds_BBB_macaque.fa —", length(bbb_seqs$macaque), "genes\n")
cat("  genomes/cds_liver_human.fa —", length(liver_human_seqs), "genes\n")
cat("\nStep 5b complete. Run step5c_align_and_score.R next.\n")

# ============================================================
# Step 5f: Recompute dN/dS using proper codon-level alignment
#
# The original step5c/d passed unaligned trimmed sequences to
# seqinr::kaks(), which assumes its input is already codon-aligned.
# For divergent species (mouse-human), this can produce biased
# dN/dS values.
#
# This script recomputes dN/dS using a codon-aware alignment:
#   1. Translate both nucleotide sequences to protein
#   2. Align the proteins (handles divergence + indels properly)
#   3. Back-map the protein alignment to codons
#   4. Strip codon columns containing gaps
#   5. Pass the cleaned codon-aligned pair to seqinr::kaks()
#
# Output: processed/conservation_scores_combined_v2.csv
# ============================================================

setwd("~/Documents/Claude/Projects/BBB")

library(Biostrings)
library(pwalign)
library(seqinr)
library(readr)
library(dplyr)

# ── Load all CDS FASTAs ──────────────────────────────────────────────────────
cat("Loading CDS FASTA files...\n")
load_indexed <- function(path) {
  s <- readDNAStringSet(path)
  if (grepl("BBB|liver_human", path)) return(s)  # already indexed
  gene_ids <- sub(".*gene:(ENS[A-Z0-9]+)\\.?[0-9]*.*", "\\1", names(s))
  df <- data.frame(gene_id = gene_ids, length = width(s),
                   idx = seq_along(s), stringsAsFactors = FALSE)
  best <- df |> group_by(gene_id) |>
    slice_max(length, n = 1, with_ties = FALSE) |> ungroup()
  out <- s[best$idx]
  names(out) <- best$gene_id
  out
}

bbb_human   <- readDNAStringSet("genomes/cds_BBB_human.fa")
bbb_mouse   <- readDNAStringSet("genomes/cds_BBB_mouse.fa")
bbb_macaque <- readDNAStringSet("genomes/cds_BBB_macaque.fa")
liver_human <- readDNAStringSet("genomes/cds_liver_human.fa")
mouse_all   <- load_indexed("genomes/Mus_musculus.GRCm39.113.cds.all.fa.gz")
macaque_all <- load_indexed("genomes/Macaca_mulatta.Mmul_10.113.cds.all.fa.gz")
cat("Loaded.\n")

# ── Helper: codon-aware alignment via protein alignment back-translation ────
codon_align <- function(nt1, nt2) {
  nt1c <- as.character(nt1); nt2c <- as.character(nt2)
  len1 <- floor(nchar(nt1c) / 3) * 3
  len2 <- floor(nchar(nt2c) / 3) * 3
  if (len1 < 30 || len2 < 30) return(NULL)
  nt1c <- substr(nt1c, 1, len1)
  nt2c <- substr(nt2c, 1, len2)

  aa1 <- suppressWarnings(as.character(Biostrings::translate(DNAString(nt1c),
                                                   if.fuzzy.codon = "X")))
  aa2 <- suppressWarnings(as.character(Biostrings::translate(DNAString(nt2c),
                                                   if.fuzzy.codon = "X")))
  # strip stop codons (* in Biostrings translate output) at the end
  aa1 <- sub("\\*+$", "", aa1)
  aa2 <- sub("\\*+$", "", aa2)
  if (nchar(aa1) < 10 || nchar(aa2) < 10) return(NULL)

  aln <- pwalign::pairwiseAlignment(aa1, aa2, type = "global",
            substitutionMatrix = "BLOSUM62",
            gapOpening = 10, gapExtension = 4)
  aa1_aln <- as.character(pwalign::alignedPattern(aln))
  aa2_aln <- as.character(pwalign::alignedSubject(aln))

  # Build codon arrays (one codon per AA position)
  codons1 <- substring(nt1c, seq(1, nchar(aa1)*3, 3), seq(3, nchar(aa1)*3, 3))
  codons2 <- substring(nt2c, seq(1, nchar(aa2)*3, 3), seq(3, nchar(aa2)*3, 3))

  c1 <- character(0); c2 <- character(0)
  i1 <- 1; i2 <- 1
  aa1_chars <- strsplit(aa1_aln, "")[[1]]
  aa2_chars <- strsplit(aa2_aln, "")[[1]]
  for (k in seq_along(aa1_chars)) {
    g1 <- aa1_chars[k] == "-"
    g2 <- aa2_chars[k] == "-"
    if (!g1 && !g2 && i1 <= length(codons1) && i2 <= length(codons2)) {
      cd1 <- codons1[i1]; cd2 <- codons2[i2]
      # Skip codons with N or non-ACGT characters
      if (!grepl("[^ACGTacgt]", cd1) && !grepl("[^ACGTacgt]", cd2)) {
        c1 <- c(c1, cd1); c2 <- c(c2, cd2)
      }
    }
    if (!g1) i1 <- i1 + 1
    if (!g2) i2 <- i2 + 1
  }

  if (length(c1) < 10) return(NULL)
  list(seq1 = paste(c1, collapse = ""), seq2 = paste(c2, collapse = ""))
}

calc_dnds_v2 <- function(seq1, seq2) {
  tryCatch({
    if (is.null(seq1) || is.null(seq2)) return(list(dN=NA,dS=NA,dNdS=NA))
    al <- codon_align(seq1, seq2)
    if (is.null(al)) return(list(dN=NA, dS=NA, dNdS=NA))
    aln <- list(nb = 2, nam = c("a","b"),
                seq = c(tolower(al$seq1), tolower(al$seq2)),
                com = c(NA, NA))
    class(aln) <- "alignment"
    res <- seqinr::kaks(aln)
    dN <- as.numeric(res$ka); dS <- as.numeric(res$ks)
    list(dN = dN, dS = dS,
         dNdS = ifelse(!is.na(dS) && dS > 0 && is.finite(dS), dN/dS, NA_real_))
  }, error = function(e) list(dN = NA, dS = NA, dNdS = NA))
}

# ── Recompute for combined scores ────────────────────────────────────────────
df <- read_csv("processed/conservation_scores_combined.csv",
               show_col_types = FALSE)
cat(sprintf("Loaded %d rows.\n", nrow(df)))

cat("Recomputing dN/dS with codon-aware alignment...\n")
n <- nrow(df)
df$dN_Human_Mouse_v2     <- NA_real_
df$dS_Human_Mouse_v2     <- NA_real_
df$dNdS_Human_Mouse_v2   <- NA_real_
df$dN_Human_Macaque_v2   <- NA_real_
df$dS_Human_Macaque_v2   <- NA_real_
df$dNdS_Human_Macaque_v2 <- NA_real_

for (i in seq_len(n)) {
  if (i %% 500 == 0) cat(sprintf("  %d / %d\n", i, n))
  h_id <- df$HumanEnsembl[i]
  m_id <- df$MouseEnsembl[i]
  q_id <- df$MacaqueEnsembl[i]
  is_bbb <- df$Gene_Set[i] == "BBB"

  h_seq <- if (!is.na(h_id)) {
    if (is_bbb && h_id %in% names(bbb_human))   bbb_human[[h_id]]
    else if (h_id %in% names(liver_human))      liver_human[[h_id]]
    else NULL
  } else NULL
  m_seq <- if (!is.na(m_id)) {
    if (is_bbb && m_id %in% names(bbb_mouse))   bbb_mouse[[m_id]]
    else if (m_id %in% names(mouse_all))        mouse_all[[m_id]]
    else NULL
  } else NULL
  q_seq <- if (!is.na(q_id)) {
    if (is_bbb && q_id %in% names(bbb_macaque)) bbb_macaque[[q_id]]
    else if (q_id %in% names(macaque_all))      macaque_all[[q_id]]
    else NULL
  } else NULL

  rm <- calc_dnds_v2(h_seq, m_seq)
  rq <- calc_dnds_v2(h_seq, q_seq)
  df$dN_Human_Mouse_v2[i]     <- rm$dN
  df$dS_Human_Mouse_v2[i]     <- rm$dS
  df$dNdS_Human_Mouse_v2[i]   <- rm$dNdS
  df$dN_Human_Macaque_v2[i]   <- rq$dN
  df$dS_Human_Macaque_v2[i]   <- rq$dS
  df$dNdS_Human_Macaque_v2[i] <- rq$dNdS
}

write_csv(df, "processed/conservation_scores_combined_v2.csv")
cat(sprintf("\nSaved: processed/conservation_scores_combined_v2.csv (%d rows)\n",
    nrow(df)))

# ── Summary comparison: old vs new ───────────────────────────────────────────
cat("\n=== dN/dS: old (unaligned) vs new (codon-aligned) ===\n\n")
for (gs in c("BBB", "Liver_Control")) {
  d <- df[df$Gene_Set == gs, ]
  cat(sprintf("%s (n=%d)\n", gs, nrow(d)))
  cat(sprintf("  Mouse   dN/dS  old=%.3f  new=%.3f\n",
      median(d$dNdS_Human_Mouse,    na.rm=TRUE),
      median(d$dNdS_Human_Mouse_v2, na.rm=TRUE)))
  cat(sprintf("  Macaque dN/dS  old=%.3f  new=%.3f\n",
      median(d$dNdS_Human_Macaque,    na.rm=TRUE),
      median(d$dNdS_Human_Macaque_v2, na.rm=TRUE)))
  cat("\n")
}

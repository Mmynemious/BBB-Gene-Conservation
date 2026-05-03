# ============================================================
# Step 5h: Direct mouse vs macaque conservation comparison
#
# Closes the triangle of pairwise comparisons. Uses the same
# pipeline as step5c/5f but compares mouse CDS directly against
# macaque CDS (no human reference).
#
# Output: processed/conservation_scores_mouse_macaque.csv
# ============================================================

setwd("~/Documents/Claude/Projects/BBB")

suppressPackageStartupMessages({
  library(Biostrings)
  library(pwalign)
  library(seqinr)
  library(readr)
  library(dplyr)
})

# ── Load CDS FASTAs (full mouse + macaque, plus BBB-specific subsets) ────────
cat("Loading CDS FASTA files...\n")
load_indexed <- function(path) {
  s <- readDNAStringSet(path)
  gene_ids <- sub(".*gene:(ENS[A-Z0-9]+)\\.?[0-9]*.*", "\\1", names(s))
  df <- data.frame(gene_id = gene_ids, length = width(s),
                   idx = seq_along(s), stringsAsFactors = FALSE)
  best <- df |> group_by(gene_id) |>
    slice_max(length, n = 1, with_ties = FALSE) |> ungroup()
  out <- s[best$idx]
  names(out) <- best$gene_id
  out
}

bbb_mouse   <- readDNAStringSet("genomes/cds_BBB_mouse.fa")
bbb_macaque <- readDNAStringSet("genomes/cds_BBB_macaque.fa")
mouse_all   <- load_indexed("genomes/Mus_musculus.GRCm39.113.cds.all.fa.gz")
macaque_all <- load_indexed("genomes/Macaca_mulatta.Mmul_10.113.cds.all.fa.gz")
cat("Loaded.\n")

# ── Helpers (same as step5f) ──────────────────────────────────────────────────
calc_pct_identity <- function(seq1, seq2) {
  tryCatch({
    aln <- pwalign::pairwiseAlignment(seq1, seq2, type = "global",
      substitutionMatrix = pwalign::nucleotideSubstitutionMatrix(
        match = 1, mismatch = -1, baseOnly = FALSE),
      gapOpening = -5, gapExtension = -2)
    pwalign::pid(aln, type = "PID1")
  }, error = function(e) NA_real_)
}

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
  aa1 <- sub("\\*+$", "", aa1)
  aa2 <- sub("\\*+$", "", aa2)
  if (nchar(aa1) < 10 || nchar(aa2) < 10) return(NULL)
  aln <- pwalign::pairwiseAlignment(aa1, aa2, type = "global",
            substitutionMatrix = "BLOSUM62", gapOpening = 10, gapExtension = 4)
  aa1_aln <- as.character(pwalign::alignedPattern(aln))
  aa2_aln <- as.character(pwalign::alignedSubject(aln))
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

calc_dnds <- function(seq1, seq2) {
  tryCatch({
    if (is.null(seq1) || is.null(seq2)) return(list(dN=NA,dS=NA,dNdS=NA))
    al <- codon_align(seq1, seq2)
    if (is.null(al)) return(list(dN=NA,dS=NA,dNdS=NA))
    aln <- list(nb=2, nam=c("a","b"),
                seq=c(tolower(al$seq1), tolower(al$seq2)), com=c(NA,NA))
    class(aln) <- "alignment"
    res <- seqinr::kaks(aln)
    dN <- as.numeric(res$ka); dS <- as.numeric(res$ks)
    list(dN = dN, dS = dS,
         dNdS = ifelse(!is.na(dS) && dS > 0 && is.finite(dS), dN/dS, NA_real_))
  }, error = function(e) list(dN=NA,dS=NA,dNdS=NA))
}

# ── Load existing scores to grab gene pairs ──────────────────────────────────
combined <- read_csv("processed/conservation_scores_combined_v2.csv",
                     show_col_types = FALSE)

# Only keep rows with BOTH mouse and macaque orthologues
pairs <- combined |>
  filter(!is.na(MouseEnsembl), !is.na(MacaqueEnsembl)) |>
  select(GeneSymbol_Human, HumanEnsembl, MouseEnsembl, MacaqueEnsembl,
         Source, Daneman_Filter, Gene_Set)

cat(sprintf("Gene pairs to score: %d (BBB=%d, Liver=%d)\n",
    nrow(pairs),
    sum(pairs$Gene_Set == "BBB"),
    sum(pairs$Gene_Set == "Liver_Control")))

# ── Score mouse vs macaque ───────────────────────────────────────────────────
cat("\nScoring mouse vs macaque...\n")
n <- nrow(pairs)
results <- list()

for (i in seq_len(n)) {
  if (i %% 500 == 0) cat(sprintf("  %d / %d\n", i, n))

  m_id <- pairs$MouseEnsembl[i]
  q_id <- pairs$MacaqueEnsembl[i]
  is_bbb <- pairs$Gene_Set[i] == "BBB"

  m_seq <- if (is_bbb && m_id %in% names(bbb_mouse))   bbb_mouse[[m_id]]
           else if (m_id %in% names(mouse_all))        mouse_all[[m_id]]
           else NULL
  q_seq <- if (is_bbb && q_id %in% names(bbb_macaque)) bbb_macaque[[q_id]]
           else if (q_id %in% names(macaque_all))      macaque_all[[q_id]]
           else NULL

  if (!is.null(m_seq) && !is.null(q_seq)) {
    pct  <- calc_pct_identity(m_seq, q_seq)
    dnds <- calc_dnds(m_seq, q_seq)
  } else {
    pct  <- NA; dnds <- list(dN=NA, dS=NA, dNdS=NA)
  }

  results[[i]] <- data.frame(
    GeneSymbol_Human  = pairs$GeneSymbol_Human[i],
    MouseEnsembl      = m_id,
    MacaqueEnsembl    = q_id,
    Source            = pairs$Source[i],
    Daneman_Filter    = pairs$Daneman_Filter[i],
    Gene_Set          = pairs$Gene_Set[i],
    PctId_Mouse_Macaque = pct,
    dN_Mouse_Macaque    = dnds$dN,
    dS_Mouse_Macaque    = dnds$dS,
    dNdS_Mouse_Macaque  = dnds$dNdS,
    stringsAsFactors  = FALSE
  )
}

mm_scores <- bind_rows(results)
write_csv(mm_scores, "processed/conservation_scores_mouse_macaque.csv")
cat(sprintf("\nSaved: processed/conservation_scores_mouse_macaque.csv (%d rows)\n",
    nrow(mm_scores)))

# ── Summary ──────────────────────────────────────────────────────────────────
cat("\n=== Mouse vs Macaque conservation ===\n\n")
for (gs in c("BBB", "Liver_Control")) {
  d <- mm_scores[mm_scores$Gene_Set == gs, ]
  cat(sprintf("%s (n=%d)\n", gs, nrow(d)))
  cat(sprintf("  Mouse vs Macaque — %% identity: median=%.1f%%\n",
      median(d$PctId_Mouse_Macaque, na.rm=TRUE)))
  cat(sprintf("  Mouse vs Macaque — dN/dS:      median=%.3f\n",
      median(d$dNdS_Mouse_Macaque, na.rm=TRUE)))
  cat("\n")
}

cat("\n=== Three-way summary (median % identity, BBB only) ===\n")
bbb_v2 <- combined[combined$Gene_Set == "BBB", ]
mm_bbb <- mm_scores[mm_scores$Gene_Set == "BBB", ]
cat(sprintf("  Human  vs Mouse:    %.1f%%\n",
    median(bbb_v2$PctId_Human_Mouse, na.rm=TRUE)))
cat(sprintf("  Human  vs Macaque:  %.1f%%\n",
    median(bbb_v2$PctId_Human_Macaque, na.rm=TRUE)))
cat(sprintf("  Mouse  vs Macaque:  %.1f%%\n",
    median(mm_bbb$PctId_Mouse_Macaque, na.rm=TRUE)))

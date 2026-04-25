# ============================================================
# Step 5c: Pairwise alignment and conservation scoring
#
# For each BBB gene with orthologues in mouse and/or macaque:
#   1. Align human CDS vs mouse CDS (pairwise, global)
#   2. Align human CDS vs macaque CDS
#   3. Calculate % nucleotide identity from each alignment
#   4. Calculate dN/dS using the ape package
#
# Then repeats for a random sample of liver control genes
# matched in number to the BBB set.
#
# Output: processed/conservation_scores.csv
# ============================================================

setwd("~/Documents/Claude/Projects/BBB")

library(Biostrings)
library(pwalign)
library(readr)
library(dplyr)
library(seqinr)

# ── Load sequences ────────────────────────────────────────────────────────────
cat("Loading CDS sequences...\n")
bbb_human   <- readDNAStringSet("genomes/cds_BBB_human.fa")
bbb_mouse   <- readDNAStringSet("genomes/cds_BBB_mouse.fa")
bbb_macaque <- readDNAStringSet("genomes/cds_BBB_macaque.fa")

cat("BBB human:  ", length(bbb_human), "genes\n")
cat("BBB mouse:  ", length(bbb_mouse), "genes\n")
cat("BBB macaque:", length(bbb_macaque), "genes\n")

# ── Load master list for gene pairing ─────────────────────────────────────────
master <- read_csv("processed/master_BBB_genelist_complete.csv",
                   show_col_types = FALSE) |>
  filter(!is.na(HumanEnsembl))

# ── Helper: % identity from pairwise alignment ────────────────────────────────
calc_pct_identity <- function(seq1, seq2) {
  tryCatch({
    aln <- pwalign::pairwiseAlignment(seq1, seq2,
                              type = "global",
                              substitutionMatrix = pwalign::nucleotideSubstitutionMatrix(
                                match = 1, mismatch = -1, baseOnly = FALSE),
                              gapOpening = -5, gapExtension = -2)
    pwalign::pid(aln, type = "PID1")  # matches / alignment length
  }, error = function(e) NA_real_)
}

# ── Helper: dN/dS via seqinr ─────────────────────────────────────────────────
calc_dnds <- function(seq1, seq2) {
  tryCatch({
    len <- floor(min(nchar(as.character(seq1)), nchar(as.character(seq2))) / 3) * 3
    if (len < 30) return(list(dN = NA, dS = NA, dNdS = NA))

    s1 <- tolower(substr(as.character(seq1), 1, len))
    s2 <- tolower(substr(as.character(seq2), 1, len))

    aln <- list(nb = 2, nam = c("seq1", "seq2"), seq = c(s1, s2), com = c(NA, NA))
    class(aln) <- "alignment"

    res <- seqinr::kaks(aln)
    dN <- as.numeric(res$ka)
    dS <- as.numeric(res$ks)
    list(dN   = dN,
         dS   = dS,
         dNdS = ifelse(!is.na(dS) && dS > 0, dN / dS, NA_real_))
  }, error = function(e) list(dN = NA, dS = NA, dNdS = NA))
}

# ── Score all BBB genes ───────────────────────────────────────────────────────
cat("\nScoring BBB genes...\n")

results <- list()
n <- nrow(master)

for (i in seq_len(n)) {
  if (i %% 100 == 0) cat(sprintf("  %d / %d\n", i, n))

  h_id <- master$HumanEnsembl[i]
  m_id <- master$MouseEnsembl[i]
  q_id <- master$MacaqueEnsembl[i]
  gene <- master$GeneSymbol_Human[i]

  h_seq <- if (!is.na(h_id) && h_id %in% names(bbb_human))
              bbb_human[[h_id]] else NULL
  m_seq <- if (!is.na(m_id) && m_id %in% names(bbb_mouse))
              bbb_mouse[[m_id]] else NULL
  q_seq <- if (!is.na(q_id) && q_id %in% names(bbb_macaque))
              bbb_macaque[[q_id]] else NULL

  # Human vs Mouse
  if (!is.null(h_seq) && !is.null(m_seq)) {
    pct_hm  <- calc_pct_identity(h_seq, m_seq)
    dnds_hm <- calc_dnds(h_seq, m_seq)
  } else {
    pct_hm  <- NA
    dnds_hm <- list(dN = NA, dS = NA, dNdS = NA)
  }

  # Human vs Macaque
  if (!is.null(h_seq) && !is.null(q_seq)) {
    pct_hq  <- calc_pct_identity(h_seq, q_seq)
    dnds_hq <- calc_dnds(h_seq, q_seq)
  } else {
    pct_hq  <- NA
    dnds_hq <- list(dN = NA, dS = NA, dNdS = NA)
  }

  results[[i]] <- data.frame(
    GeneSymbol_Human   = gene,
    HumanEnsembl       = h_id,
    MouseEnsembl       = m_id,
    MacaqueEnsembl     = q_id,
    Source             = master$Source[i],
    Daneman_Filter     = master$Daneman_Filter[i],
    Mouse_OrthoType    = master$Mouse_OrthoType[i],
    Macaque_OrthoType  = master$Macaque_OrthoType[i],
    Human_BBB_Datasets_N = master$Human_BBB_Datasets_N[i],
    PctId_Human_Mouse    = pct_hm,
    PctId_Human_Macaque  = pct_hq,
    dN_Human_Mouse       = dnds_hm$dN,
    dS_Human_Mouse       = dnds_hm$dS,
    dNdS_Human_Mouse     = dnds_hm$dNdS,
    dN_Human_Macaque     = dnds_hq$dN,
    dS_Human_Macaque     = dnds_hq$dS,
    dNdS_Human_Macaque   = dnds_hq$dNdS,
    Gene_Set             = "BBB",
    stringsAsFactors     = FALSE
  )
}

bbb_scores <- bind_rows(results)

# Save immediately after scoring so a later error doesn't lose the data
write_csv(bbb_scores, "processed/conservation_scores_BBB.csv")
cat(sprintf("Saved: processed/conservation_scores_BBB.csv (%d rows)\n", nrow(bbb_scores)))

cat(sprintf("BBB genes scored: %d rows\n", nrow(bbb_scores)))
cat("  With mouse pct identity:  ", sum(!is.na(bbb_scores$PctId_Human_Mouse)), "\n")
cat("  With macaque pct identity:", sum(!is.na(bbb_scores$PctId_Human_Macaque)), "\n")

# ── Score liver control genes (human CDS only — use BioMart % id for others) ─
# For the liver control we use the % identity already in the master list approach:
# sample ~same number as BBB set from liver, then score human vs mouse/macaque
# using the same method but we only have human sequences for now.
# We report liver human CDS length distribution as a baseline check,
# and note that full liver cross-species scoring requires a biomaRt orthologue
# lookup for the liver genes (planned as Step 5d extension).

liver_human <- readDNAStringSet("genomes/cds_liver_human.fa")
cat(sprintf("\nLiver human CDS: %d genes loaded\n", length(liver_human)))
cat(sprintf("Liver CDS length: median=%d bp, mean=%d bp\n",
    as.integer(median(width(liver_human))),
    as.integer(mean(width(liver_human)))))
cat(sprintf("BBB CDS length:   median=%d bp, mean=%d bp\n",
    as.integer(median(width(bbb_human))),
    as.integer(mean(width(bbb_human)))))

# (already saved above)

# ── Summary statistics ────────────────────────────────────────────────────────
cat("\n=== Conservation Summary ===\n")
cat(sprintf("Human vs Mouse   — median %%id: %.1f%%  mean: %.1f%%\n",
    median(bbb_scores$PctId_Human_Mouse, na.rm=TRUE),
    mean(bbb_scores$PctId_Human_Mouse, na.rm=TRUE)))
cat(sprintf("Human vs Macaque — median %%id: %.1f%%  mean: %.1f%%\n",
    median(bbb_scores$PctId_Human_Macaque, na.rm=TRUE),
    mean(bbb_scores$PctId_Human_Macaque, na.rm=TRUE)))
cat(sprintf("Human vs Mouse   — median dN/dS: %.3f\n",
    median(bbb_scores$dNdS_Human_Mouse, na.rm=TRUE)))
cat(sprintf("Human vs Macaque — median dN/dS: %.3f\n",
    median(bbb_scores$dNdS_Human_Macaque, na.rm=TRUE)))

cat("\nStep 5c complete. Run step5d_liver_scoring.R for control group.\n")

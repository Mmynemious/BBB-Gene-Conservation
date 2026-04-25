# ============================================================
# Step 5d: Score liver control genes for conservation
#
# Gets mouse + macaque orthologues for liver control genes via
# biomaRt, extracts their CDS, then runs the same pairwise
# alignment pipeline as step5c.
#
# Output: processed/conservation_scores_liver.csv
#         processed/conservation_scores_combined.csv
# ============================================================

setwd("~/Documents/Claude/Projects/BBB")

library(Biostrings)
library(pwalign)
library(biomaRt)
library(readr)
library(dplyr)
library(seqinr)

# ── Load liver human sequences ────────────────────────────────────────────────
cat("Loading liver human CDS sequences...\n")
liver_human <- readDNAStringSet("genomes/cds_liver_human.fa")
liver_human_ids <- names(liver_human)
cat("Liver human genes:", length(liver_human_ids), "\n")

# ── Query biomaRt for mouse + macaque orthologues ─────────────────────────────
cat("\nConnecting to Ensembl BioMart...\n")
human_mart <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl",
                          host = "https://www.ensembl.org")
cat("Connected.\n")

# Query in batches of 500 to avoid server timeouts
query_orthologues <- function(ids, species_attr_prefix, mart, label, cache_file) {
  if (file.exists(cache_file)) {
    cat("  Loading cached", label, "orthologues from", cache_file, "\n")
    return(read_csv(cache_file, show_col_types = FALSE))
  }

  attrs <- c("ensembl_gene_id",
             paste0(species_attr_prefix, "_homolog_ensembl_gene"),
             paste0(species_attr_prefix, "_homolog_orthology_type"),
             paste0(species_attr_prefix, "_homolog_perc_id"))

  batch_size <- 500
  batches <- split(ids, ceiling(seq_along(ids) / batch_size))
  results <- list()

  for (i in seq_along(batches)) {
    cat(sprintf("  Batch %d / %d ...\n", i, length(batches)))
    res <- tryCatch(
      getBM(attributes = attrs, filters = "ensembl_gene_id",
            values = batches[[i]], mart = mart),
      error = function(e) { cat("  ERROR in batch", i, ":", conditionMessage(e), "\n"); NULL }
    )
    if (!is.null(res)) results[[i]] <- res
  }

  out <- bind_rows(results) |>
    filter(.data[[attrs[2]]] != "") |>
    rename(HumanEnsembl = ensembl_gene_id,
           Orthologue   = !!attrs[2],
           OrthoType    = !!attrs[3],
           PercId       = !!attrs[4]) |>
    group_by(HumanEnsembl) |>
    slice_max(PercId, n = 1, with_ties = FALSE) |>
    ungroup()

  write_csv(out, cache_file)
  cat(" ", nrow(out), label, "orthologues found and cached.\n")
  out
}

cat("Querying mouse orthologues for liver genes...\n")
liver_mouse_raw <- query_orthologues(liver_human_ids, "mmusculus", human_mart,
                                      "mouse", "processed/liver_mouse_ortho_cache.csv") |>
  rename(MouseEnsembl = Orthologue, Mouse_OrthoType = OrthoType, Mouse_PercId = PercId)

cat("Querying macaque orthologues for liver genes...\n")
liver_macaque_raw <- query_orthologues(liver_human_ids, "mmulatta", human_mart,
                                        "macaque", "processed/liver_macaque_ortho_cache.csv") |>
  rename(MacaqueEnsembl = Orthologue, Macaque_OrthoType = OrthoType, Macaque_PercId = PercId)

liver_mouse_ortho   <- liver_mouse_raw
liver_macaque_ortho <- liver_macaque_raw

# ── Load mouse + macaque CDS FASTAs ──────────────────────────────────────────
cat("\nLoading mouse and macaque CDS FASTAs...\n")
mouse_cds_all   <- readDNAStringSet("genomes/Mus_musculus.GRCm39.113.cds.all.fa.gz")
macaque_cds_all <- readDNAStringSet("genomes/Macaca_mulatta.Mmul_10.113.cds.all.fa.gz")

# Index by gene ID (longest transcript per gene)
index_by_gene <- function(seqs) {
  gene_ids <- sub(".*gene:(ENS[A-Z0-9]+)\\.?[0-9]*.*", "\\1", names(seqs))
  df <- data.frame(gene_id = gene_ids, length = width(seqs),
                   idx = seq_along(seqs), stringsAsFactors = FALSE)
  best <- df |> group_by(gene_id) |>
    slice_max(length, n = 1, with_ties = FALSE) |> ungroup()
  out <- seqs[best$idx]
  names(out) <- best$gene_id
  out
}

mouse_cds   <- index_by_gene(mouse_cds_all)
macaque_cds <- index_by_gene(macaque_cds_all)
rm(mouse_cds_all, macaque_cds_all)

# ── Alignment helpers (same as step5c) ───────────────────────────────────────
calc_pct_identity <- function(seq1, seq2) {
  tryCatch({
    aln <- pwalign::pairwiseAlignment(seq1, seq2, type = "global",
      substitutionMatrix = pwalign::nucleotideSubstitutionMatrix(
        match = 1, mismatch = -1, baseOnly = FALSE),
      gapOpening = -5, gapExtension = -2)
    pwalign::pid(aln, type = "PID1")
  }, error = function(e) NA_real_)
}

calc_dnds <- function(seq1, seq2) {
  tryCatch({
    len <- floor(min(nchar(as.character(seq1)),
                     nchar(as.character(seq2))) / 3) * 3
    if (len < 30) return(list(dN = NA, dS = NA, dNdS = NA))
    s1 <- tolower(substr(as.character(seq1), 1, len))
    s2 <- tolower(substr(as.character(seq2), 1, len))
    aln <- list(nb = 2, nam = c("seq1","seq2"), seq = c(s1,s2), com = c(NA,NA))
    class(aln) <- "alignment"
    res <- seqinr::kaks(aln)
    dN <- as.numeric(res$ka); dS <- as.numeric(res$ks)
    list(dN = dN, dS = dS,
         dNdS = ifelse(!is.na(dS) && dS > 0, dN / dS, NA_real_))
  }, error = function(e) list(dN = NA, dS = NA, dNdS = NA))
}

# ── Build liver orthologue table ──────────────────────────────────────────────
liver_table <- data.frame(HumanEnsembl = liver_human_ids,
                           stringsAsFactors = FALSE) |>
  left_join(liver_mouse_ortho,   by = "HumanEnsembl") |>
  left_join(liver_macaque_ortho, by = "HumanEnsembl")

cat("\nLiver orthologue table:", nrow(liver_table), "rows\n")
cat("  With mouse orthologue:   ", sum(!is.na(liver_table$MouseEnsembl)), "\n")
cat("  With macaque orthologue: ", sum(!is.na(liver_table$MacaqueEnsembl)), "\n")

# ── Score liver genes ─────────────────────────────────────────────────────────
cat("\nScoring liver control genes...\n")
results <- list()
n <- nrow(liver_table)

for (i in seq_len(n)) {
  if (i %% 500 == 0) cat(sprintf("  %d / %d\n", i, n))

  h_id <- liver_table$HumanEnsembl[i]
  m_id <- liver_table$MouseEnsembl[i]
  q_id <- liver_table$MacaqueEnsembl[i]

  h_seq <- if (!is.na(h_id) && h_id %in% names(liver_human))
              liver_human[[h_id]] else NULL
  m_seq <- if (!is.na(m_id) && m_id %in% names(mouse_cds))
              mouse_cds[[m_id]] else NULL
  q_seq <- if (!is.na(q_id) && q_id %in% names(macaque_cds))
              macaque_cds[[q_id]] else NULL

  if (!is.null(h_seq) && !is.null(m_seq)) {
    pct_hm  <- calc_pct_identity(h_seq, m_seq)
    dnds_hm <- calc_dnds(h_seq, m_seq)
  } else {
    pct_hm  <- NA; dnds_hm <- list(dN = NA, dS = NA, dNdS = NA)
  }

  if (!is.null(h_seq) && !is.null(q_seq)) {
    pct_hq  <- calc_pct_identity(h_seq, q_seq)
    dnds_hq <- calc_dnds(h_seq, q_seq)
  } else {
    pct_hq  <- NA; dnds_hq <- list(dN = NA, dS = NA, dNdS = NA)
  }

  results[[i]] <- data.frame(
    GeneSymbol_Human     = NA_character_,
    HumanEnsembl         = h_id,
    MouseEnsembl         = m_id,
    MacaqueEnsembl       = q_id,
    Source               = "Liver_Control",
    Daneman_Filter       = NA_character_,
    Mouse_OrthoType      = liver_table$Mouse_OrthoType[i],
    Macaque_OrthoType    = liver_table$Macaque_OrthoType[i],
    Human_BBB_Datasets_N = NA_integer_,
    PctId_Human_Mouse    = pct_hm,
    PctId_Human_Macaque  = pct_hq,
    dN_Human_Mouse       = dnds_hm$dN,
    dS_Human_Mouse       = dnds_hm$dS,
    dNdS_Human_Mouse     = dnds_hm$dNdS,
    dN_Human_Macaque     = dnds_hq$dN,
    dS_Human_Macaque     = dnds_hq$dS,
    dNdS_Human_Macaque   = dnds_hq$dNdS,
    Gene_Set             = "Liver_Control",
    stringsAsFactors     = FALSE
  )
}

liver_scores <- bind_rows(results)
write_csv(liver_scores, "processed/conservation_scores_liver.csv")
cat(sprintf("\nSaved: processed/conservation_scores_liver.csv (%d rows)\n",
    nrow(liver_scores)))

# ── Combine and summarise ─────────────────────────────────────────────────────
bbb_scores <- read_csv("processed/conservation_scores_BBB.csv",
                        show_col_types = FALSE)
combined   <- bind_rows(bbb_scores, liver_scores)
write_csv(combined, "processed/conservation_scores_combined.csv")

cat("\n=== Conservation Comparison: BBB vs Liver ===\n\n")
for (gs in c("BBB", "Liver_Control")) {
  d <- combined[combined$Gene_Set == gs, ]
  cat(sprintf("%s (n=%d)\n", gs, nrow(d)))
  cat(sprintf("  Human vs Mouse   — median %%id: %.1f%%  dN/dS: %.3f\n",
      median(d$PctId_Human_Mouse, na.rm=TRUE),
      median(d$dNdS_Human_Mouse,  na.rm=TRUE)))
  cat(sprintf("  Human vs Macaque — median %%id: %.1f%%  dN/dS: %.3f\n",
      median(d$PctId_Human_Macaque, na.rm=TRUE),
      median(d$dNdS_Human_Macaque,  na.rm=TRUE)))
  cat("\n")
}

cat("Saved: processed/conservation_scores_combined.csv\n")
cat("\nStep 5d complete.\n")

# eval_step4.R — Validate the downloaded GTF files

eval_gtf_downloads <- function() {

  cat("=== EVAL: GTF File Downloads ===\n\n")

  files <- list(
    list(path = "genomes/Homo_sapiens.GRCh38.113.gtf.gz",    label = "Human",   min_mb = 40),
    list(path = "genomes/Mus_musculus.GRCm39.113.gtf.gz",    label = "Mouse",   min_mb = 25),
    list(path = "genomes/Macaca_mulatta.Mmul_10.113.gtf.gz", label = "Macaque", min_mb = 10)
  )

  all_ok <- TRUE

  for (f in files) {
    if (!file.exists(f$path)) {
      cat("FAIL:", f$label, "— file not found at", f$path, "\n")
      all_ok <- FALSE
      next
    }

    size_mb <- round(file.size(f$path) / 1e6, 1)

    if (size_mb < f$min_mb) {
      cat(sprintf("FAIL: %s — file too small (%.1f MB, expected >%d MB) — possible failed download\n",
                  f$label, size_mb, f$min_mb))
      all_ok <- FALSE
    } else {
      cat(sprintf("PASS: %s — %.1f MB at %s\n", f$label, size_mb, f$path))
    }
  }

  if (all_ok) {
    cat("\nAll GTF files present and correct size.\n")
    cat("These files are NOT committed to git (too large) — kept local only.\n")
  }

  cat("\n=== EVAL COMPLETE ===\n")
}

eval_gtf_downloads()

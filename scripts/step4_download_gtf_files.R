# ============================================================
# Step 4: Download Ensembl GTF annotation files
#
# Downloads genome annotation files (GTF format) for all three
# species from Ensembl release 113. GTF files tell us the exact
# chromosome, start, end, exon structure of every gene.
#
# Files are saved to genomes/ — NOT committed to git (too large).
# ============================================================

setwd("~/Documents/Claude/Projects/BBB")

dir.create("genomes", showWarnings = FALSE)

# Ensembl release 113 — matches the biomaRt version used in Step 3
RELEASE <- "113"

files <- list(
  human = list(
    url  = paste0("https://ftp.ensembl.org/pub/release-", RELEASE,
                  "/gtf/homo_sapiens/Homo_sapiens.GRCh38.", RELEASE, ".gtf.gz"),
    dest = "genomes/Homo_sapiens.GRCh38.113.gtf.gz",
    label = "Human (GRCh38)"
  ),
  mouse = list(
    url  = paste0("https://ftp.ensembl.org/pub/release-", RELEASE,
                  "/gtf/mus_musculus/Mus_musculus.GRCm39.", RELEASE, ".gtf.gz"),
    dest = "genomes/Mus_musculus.GRCm39.113.gtf.gz",
    label = "Mouse (GRCm39)"
  ),
  macaque = list(
    url  = paste0("https://ftp.ensembl.org/pub/release-", RELEASE,
                  "/gtf/macaca_mulatta/Macaca_mulatta.Mmul_10.", RELEASE, ".gtf.gz"),
    dest = "genomes/Macaca_mulatta.Mmul_10.113.gtf.gz",
    label = "Macaque (Mmul_10)"
  )
)

# ── Download each file ───────────────────────────────────────────────────────
for (sp in names(files)) {
  f <- files[[sp]]

  if (file.exists(f$dest)) {
    size_mb <- round(file.size(f$dest) / 1e6, 1)
    cat(f$label, "— already downloaded (", size_mb, "MB ), skipping.\n")
    next
  }

  cat("Downloading", f$label, "...\n")
  cat("URL:", f$url, "\n")

  tryCatch({
    download.file(f$url, destfile = f$dest, mode = "wb", method = "curl",
                  extra = "--progress-bar")
    size_mb <- round(file.size(f$dest) / 1e6, 1)
    cat("Done —", size_mb, "MB\n\n")
  }, error = function(e) {
    cat("ERROR downloading", f$label, ":", conditionMessage(e), "\n")
  })
}

# ── Verify downloads ─────────────────────────────────────────────────────────
cat("\n--- Download summary ---\n")
for (sp in names(files)) {
  f <- files[[sp]]
  if (file.exists(f$dest)) {
    size_mb <- round(file.size(f$dest) / 1e6, 1)
    cat("OK ", f$label, "—", size_mb, "MB at", f$dest, "\n")
  } else {
    cat("MISSING:", f$label, "\n")
  }
}

cat("\nStep 4 complete. Run scripts/eval_step4.R to validate.\n")

# ============================================================
# Step 5a: Download CDS FASTA files from Ensembl release 113
#
# These files contain the coding sequences (CDS) of every
# annotated transcript for each species — pre-extracted by
# Ensembl, so we don't need the full genome FASTA files.
#
# Files saved to genomes/ — NOT committed to git (too large).
# ============================================================

setwd("~/Documents/Claude/Projects/BBB")

RELEASE <- "113"

files <- list(
  human = list(
    url  = paste0("https://ftp.ensembl.org/pub/release-", RELEASE,
                  "/fasta/homo_sapiens/cds/Homo_sapiens.GRCh38.cds.all.fa.gz"),
    dest  = "genomes/Homo_sapiens.GRCh38.113.cds.all.fa.gz",
    label = "Human CDS"
  ),
  mouse = list(
    url  = paste0("https://ftp.ensembl.org/pub/release-", RELEASE,
                  "/fasta/mus_musculus/cds/Mus_musculus.GRCm39.cds.all.fa.gz"),
    dest  = "genomes/Mus_musculus.GRCm39.113.cds.all.fa.gz",
    label = "Mouse CDS"
  ),
  macaque = list(
    url  = paste0("https://ftp.ensembl.org/pub/release-", RELEASE,
                  "/fasta/macaca_mulatta/cds/Macaca_mulatta.Mmul_10.cds.all.fa.gz"),
    dest  = "genomes/Macaca_mulatta.Mmul_10.113.cds.all.fa.gz",
    label = "Macaque CDS"
  )
)

for (sp in names(files)) {
  f <- files[[sp]]

  if (file.exists(f$dest)) {
    size_mb <- round(file.size(f$dest) / 1e6, 1)
    cat(f$label, "— already downloaded (", size_mb, "MB), skipping.\n")
    next
  }

  cat("Downloading", f$label, "...\n")
  tryCatch({
    download.file(f$url, destfile = f$dest, mode = "wb", method = "curl",
                  extra = "--progress-bar")
    size_mb <- round(file.size(f$dest) / 1e6, 1)
    cat("Done —", size_mb, "MB\n\n")
  }, error = function(e) {
    cat("ERROR:", conditionMessage(e), "\n")
  })
}

cat("\n--- Download summary ---\n")
for (sp in names(files)) {
  f <- files[[sp]]
  if (file.exists(f$dest)) {
    cat("OK:", f$label, "—", round(file.size(f$dest)/1e6, 1), "MB\n")
  } else {
    cat("MISSING:", f$label, "\n")
  }
}

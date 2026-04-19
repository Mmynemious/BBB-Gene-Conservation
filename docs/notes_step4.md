# Step 4 Notes — Download GTF Annotation Files

## What Step 4 Did

Downloaded genome annotation files (GTF format) from Ensembl release 113 for all three species:

| Species | Genome build | File size | Local path |
|---------|-------------|-----------|------------|
| Human | GRCh38 | 64.1 MB | `genomes/Homo_sapiens.GRCh38.113.gtf.gz` |
| Mouse | GRCm39 | 40.5 MB | `genomes/Mus_musculus.GRCm39.113.gtf.gz` |
| Macaque | Mmul_10 | 19.8 MB | `genomes/Macaca_mulatta.Mmul_10.113.gtf.gz` |

These files are **not committed to git** — they are too large and can be re-downloaded at any time from Ensembl FTP.

---

## What a GTF File Is

A GTF (Gene Transfer Format) file is the genome's postal directory. For every gene in a species it records:
- Which chromosome it lives on
- Exact start and end position in base pairs
- Where each exon begins and ends
- Which strand of DNA it is on

This is what allows Step 5 to look up the precise DNA sequence of each BBB gene in each species and compare them.

---

## Why Ensembl Release 113

Release 113 was chosen to match the biomaRt version used in Step 3, ensuring the Ensembl IDs in `processed/BBB_genes_crossspecies_ensembl.csv` correspond exactly to the gene coordinates in these GTF files. Using different releases could cause ID mismatches where a gene ID from Step 3 does not appear in the GTF file.

---

## Why the Macaque File is Smaller

The macaque (19.8 MB) GTF is smaller than human (64.1 MB) and mouse (40.5 MB) because the macaque genome is less thoroughly annotated. Fewer genes, transcripts, and exon boundaries have been confirmed experimentally in macaque compared to the two more extensively studied species. This is consistent with the 83% macaque orthologue coverage seen in Step 3.

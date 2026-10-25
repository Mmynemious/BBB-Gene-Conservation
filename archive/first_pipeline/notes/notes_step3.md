# Step 3 Notes — Map BBB Genes to Ensembl IDs Across Three Species

## What Step 3 Did

Took the 1,526 human gene symbols from the master BBB gene list and queried Ensembl BioMart to find the permanent genomic identifier for each gene in all three species:

- **Human** (*Homo sapiens*) → `ENSG...` IDs
- **Mouse** (*Mus musculus*) → `ENSMUSG...` IDs
- **Macaque** (*Macaca mulatta*) → `ENSMMUG...` IDs

**Output:** `processed/BBB_genes_crossspecies_ensembl.csv` — 1,868 rows, 7 columns including the original master list columns plus the three Ensembl ID columns.

---

## Why Ensembl IDs and Not Just Gene Symbols

Gene symbols (like SLC2A1) are convenient for humans but unreliable for computational work. The same gene can have different names in different databases, and symbols sometimes change as annotation improves. Ensembl IDs are permanent, species-specific identifiers that point to an exact location in a specific genome build. Every tool used in Steps 4 and 5 (GTF files, sequence databases, alignment tools) speaks Ensembl, not gene symbols.

---

## Coverage Results

| Species | Ensembl IDs found | Coverage |
|---------|------------------|----------|
| Human | 1,776 | 95.1% |
| Mouse | 1,622 | 86.8% |
| Macaque | 1,552 | 83.1% |

The ~5% of human genes with no Ensembl ID are likely genes that have been retired, renamed, or are not in the current Ensembl release. The lower mouse and macaque coverage reflects missing orthologue annotations — not all human genes have a confirmed equivalent in every species.

---

## Why Three Separate Queries Were Needed

BioMart organises its attributes into "pages" — gene attributes (like Ensembl ID and gene symbol) are on one page, and homologue attributes (like mouse/macaque orthologue IDs) are on another. BioMart does not allow mixing attributes from different pages in a single query.

The solution was three queries joined on the human Ensembl ID:
1. Human gene symbol → human Ensembl ID
2. Human Ensembl ID → mouse orthologue Ensembl ID
3. Human Ensembl ID → macaque orthologue Ensembl ID

---

## Many-to-Many Relationships

The final table has 1,868 rows despite starting from 1,526 genes. This is because some human genes have multiple orthologues in mouse or macaque — a consequence of gene duplication events during evolution where one ancestral gene became two or more related genes in one lineage. These are biologically valid and are kept in the table. Downstream steps will need to decide how to handle genes with multiple orthologues (e.g. take the best-matching one, or analyse all).

---

## Missing Macaque Orthologues — Known Limitation

ABCB1 (P-glycoprotein — a major BBB drug transporter) and ABCG2 have no macaque orthologue recorded in Ensembl. This does not necessarily mean the gene is absent from the macaque genome — it may simply not be annotated as an orthologue in the current Ensembl release. The macaque (*Macaca mulatta*) genome is less thoroughly annotated than human or mouse, and some orthologue relationships are missing or not yet confirmed.

This is a known limitation of macaque genomic resources and is why macaque coverage (83%) is lower than human or mouse.

---

## Macaque Enters the Pipeline Here

Steps 1 and 2 worked entirely with human and mouse data — macaque had no published BBB gene list to contribute. Starting from Step 3, macaque enters as a **genome**, not a gene list. The orthologue lookup answers: does each human BBB gene have an equivalent in the macaque genome? If yes, that equivalent gene's DNA sequence will be compared in Step 5. This is the standard approach when no species-specific experimental data exists — use the genome annotation to find the equivalent gene computationally.

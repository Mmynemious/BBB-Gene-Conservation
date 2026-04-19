# Step 3b Notes — Orthology Confidence Scoring

## What Step 3b Did

Re-queried Ensembl BioMart (release 115, main server) to:
1. Recover genes that HomoloGene could not convert in Step 1
2. Add `orthology_type` for each mouse and macaque orthologue
3. Add `perc_id` (% sequence identity between human gene and its orthologue)

**Output:** `processed/BBB_genes_orthology_confidence.csv` — 1,724 rows (expanded from 1,526 due to paralogue rows), 12 columns including all validation columns from Step 4b.

---

## New Columns Added

| Column | Meaning |
|--------|---------|
| `Mouse_OrthoType` | `ortholog_one2one`, `ortholog_one2many`, `ortholog_many2many` |
| `Mouse_PercId` | % sequence identity between human gene and mouse orthologue |
| `Macaque_OrthoType` | Same as above for macaque |
| `Macaque_PercId` | % sequence identity between human gene and macaque orthologue |

---

## Key Findings

### % Sequence Identity

| Species | Median | Mean | Min |
|---------|--------|------|-----|
| Mouse | 88.7% | 82.7% | 12.7% |
| Macaque | 98.1% | 95.8% | 24.4% |

Macaque BBB genes are substantially more conserved relative to human than mouse BBB genes. This is expected — macaque and human diverged ~25 million years ago; mouse and human ~87 million years ago. The high macaque conservation supports its validity as a model.

### Orthologue Type Distribution

**Mouse:**
- `ortholog_one2one`: 1,327 (clean 1-to-1 correspondence)
- `ortholog_one2many`: 196 (gene duplicated in rodent lineage)
- `ortholog_many2many`: 18 (complex duplication events in both lineages)

**Macaque:**
- `ortholog_one2one`: 1,401
- `ortholog_one2many`: 65
- `ortholog_many2many`: 8

### Notable Gene Examples

- **SLC2A1** (GLUT1, the main BBB glucose transporter): one-to-one orthologue in both species, 97.2% identical to mouse and 98.6% to macaque — extremely conserved
- **ABCB1** (P-gp drug transporter): two mouse orthologues (`Abcb1a` 87.3% identity, `Abcb1b` 80.5% identity) — the classic BBB paralogue case
- **ABCG2** (BCRP): 54.5% identity to its mouse orthologue despite being a one-to-one match — lower conservation than most

---

## Why the Row Count Expanded (1,526 → 1,724)

Genes with `ortholog_one2many` relationships produce multiple rows — one per orthologue. For example, human ABCB1 generates two rows: one for Abcb1a and one for Abcb1b. These are biologically real (gene duplication in the rodent lineage) and are kept with their individual % identity scores.

**Decision for downstream analysis:** both paralogue rows are retained and flagged by their `orthology_type`. Step 5 will report both % identities for these genes and flag them as duplicated in the rodent lineage. This is more informative than arbitrarily choosing one.

---

## Why HomoloGene Was Not Sufficient

The original Step 1 used the HomoloGene R package (offline, database frozen ~2014) for mouse→human gene symbol conversion. BioMart uses the current Ensembl annotation, which:
- Captures more orthologue relationships
- Provides confidence scores (% identity, orthologue type)
- Distinguishes clean one-to-one from ambiguous one-to-many relationships

The BioMart approach is now the authoritative source for Ensembl IDs and orthology data in this pipeline.

---

## Note on Ensembl Release Version

Step 3 used release 113 (via the archive server); Step 3b used the current release (115, main server) because the archive server was unavailable. Release 115 is a newer annotation and may have slightly different gene models. This version difference is noted here as a limitation — IDs from release 113 (in the GTF files downloaded in Step 4) and release 115 (in this table) should be cross-checked if any ID mismatches appear in Step 5.

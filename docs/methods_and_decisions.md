# Methods, Decisions, and Assumptions
*This document explains the scientific reasoning behind every major decision in this pipeline. It is not a code guide — the scripts and step notes cover that. This is about why we did things the way we did, where the uncertainties lie, and what still needs external validation.*

---

## The Central Scientific Question

**How similar are the genomes of human, rhesus macaque, and mouse across blood-brain barrier (BBB) genes?**

Specifically: are BBB genes more conserved across these three species than randomly selected genes? If yes, that supports the validity of mouse models for BBB research. If not, it raises questions about how well mouse BBB biology translates to humans and macaques.

---

## Pipeline Overview

| Step | What it does | Output |
|------|-------------|--------|
| Step 1 | Build unified BBB master gene list from Daneman + Munji | `master_BBB_genelist.csv` (1,526 genes) |
| Step 2 | Process Wälchli Seurat object — average expression per EC cell type | `Walchli2024_AverageExpression.csv` |
| Step 3 | Map BBB genes to Ensembl IDs across human, mouse, macaque | `BBB_genes_crossspecies_ensembl.csv` |
| Step 3b | Re-query BioMart for orthology type + % sequence identity | `master_BBB_genelist_complete.csv` |
| Step 4b | Cross-reference master list against 3 human EC datasets | Adds `Human_BBB_Datasets_N` to master list |
| Step 4c | Build liver control gene set from Human Protein Atlas | `liver_control_genelist.csv` (7,131 genes) |
| Step 5a | Download Ensembl CDS FASTA files (coding sequences only) | `genomes/*.cds.all.fa.gz` |
| Step 5b | Extract longest CDS per gene for BBB and liver gene sets | `cds_BBB_*.fa`, `cds_liver_human.fa` |
| Step 5c | Pairwise alignment + % identity + dN/dS for BBB genes | `conservation_scores_BBB.csv` |
| Step 5d | Same alignment pipeline for liver control genes | `conservation_scores_liver.csv` |
| Step 5e | Statistical comparison (Wilcoxon) BBB vs liver | `conservation_summary.txt` |
| Step 5f | Recompute dN/dS with proper codon-aware alignment | `conservation_scores_combined_v2.csv` |

**Note on the original Step 4:** An earlier draft of this pipeline included a Step 4 that downloaded GTF (genome annotation) files. The plan was to use those coordinates plus full genome FASTA files to extract coding sequences manually. Step 5a replaced this with a simpler approach — Ensembl provides pre-extracted CDS FASTA files (~20 MB per species vs ~800 MB for full genomes), so the GTF approach was unused and has been removed from the pipeline.

---

## Scientific Assumptions

### Assumption 1: Mouse-derived BBB genes are valid proxies for human BBB genes

**What we did:** The master gene list was built from two mouse studies (Daneman 2010, Munji 2019). We converted mouse gene symbols to human HGNC format using HomoloGene (orthologue mapping) and treated the resulting list as the foundation for all downstream analysis.

**Why this is reasonable:** HomoloGene identifies orthologues — genes in different species that descended from the same ancestral gene and retained equivalent function. When it maps mouse *Slc2a1* to human *SLC2A1*, these genes share ~90 million years of evolutionary history and make functionally equivalent proteins. Using mouse BBB data as a starting point is standard practice in the field — it is the basis of nearly all BBB drug development research.

**Where the assumption could break down:** Some genes are BBB-enriched in mice but not in humans, or the expression pattern differs between species. We are not claiming the master list is a confirmed human BBB gene list — it is a candidate list based on mouse evidence.

**How we address it:** Step 4b cross-referenced the full master list against three independent human brain EC datasets (Wälchli 2024, Yang 2022, Winkler 2022). Each gene now has a `Human_BBB_Datasets_N` score (0–3) indicating how many human datasets confirm it. This allows downstream analysis to be stratified by human validation confidence.

---

### Assumption 2: Orthologues represent genuine functional equivalents

**What we did:** Used the HomoloGene R package for the initial mouse→human conversion, then biomaRt (Ensembl) for all three species in Steps 3 and 3b.

**Why this is reasonable:** Orthologue databases are built from genome-wide alignments and phylogenetic analysis. They are the standard tool for cross-species gene mapping in genomics.

**Where the assumption could break down:** HomoloGene has not been updated since ~2014. 138 mouse genes had no human orthologue in HomoloGene and were initially dropped.

**How we address it:** Step 3b re-ran the full orthologue lookup using biomaRt (Ensembl release 115), recovering previously dropped genes and adding two confidence metrics per gene: `orthology_type` (`ortholog_one2one`, `ortholog_one2many`, `ortholog_many2many`) and `perc_id` (% sequence identity). The final table has 1,724 rows. Mouse BBB genes are 88.7% identical to human on average; macaque genes are 98.1% identical.

---

### Assumption 3: The Daneman_Filter column reflects meaningful biological stringency

**What we did:** Added a `Daneman_Filter` column with values `liver_and_lung`, `liver_only`, and `lung_only`, based on which Daneman supplementary table each gene came from (S3, S4, or S5 respectively).

**Why this is reasonable:** Daneman S3 required enrichment in brain ECs vs both liver AND lung simultaneously — a very strict double filter. S4 and S5 each tested only one peripheral tissue.

**Scientific value:** The conservation analysis will run on S3 only first (strictest), then on the full list, to test whether the most stringently defined BBB genes are also the most evolutionarily conserved.

---

### Assumption 4: Liver-expressed genes are a biologically appropriate control

**What we did:** Downloaded HPA RNA tissue consensus expression data and selected 7,131 genes with liver nTPM ≥ 10 that are not in the BBB master list as the conservation control group.

**Why liver:** Liver is metabolically active and highly conserved across mammals — it provides a meaningful conservation baseline. If BBB genes are more conserved than liver genes, that is a strong result specific to the BBB. A random gene set would not provide this biological context.

**Source:** Human Protein Atlas (proteinatlas.org), RNA tissue consensus data. Landmark liver genes ALB (198,523 nTPM), CYP3A4 (3,367 nTPM), APOB (1,114 nTPM) are all confirmed present.

---

## Key Decisions Made

### Decision 1: Use Daneman S3, S4, and S5 — not just S3

**Initial approach:** Started with only S3 (most stringent), giving 545 genes with only 1 overlap with Munji.

**Problem identified:** Canonical BBB genes like *Slc2a1* (GLUT1) and *Cldn5* are absent from S3 because they are also expressed in liver endothelium and fail the double-comparison filter.

**Decision:** Include S3, S4, and S5 with filter labels. This expanded the list to 1,526 genes and increased Daneman/Munji overlap from 1 to 80 genes — a far more biologically plausible result.

---

### Decision 2: HomoloGene for Step 1, biomaRt for Steps 3 and 3b

**Reason:** Ensembl BioMart was unreachable at the time of Step 1 due to server downtime. HomoloGene (offline, bundled in R) was used as a fallback.

**Resolution:** Step 3b re-ran all orthologue lookups using biomaRt once servers were stable. The HomoloGene conversion is retained in the master list for traceability but biomaRt is now the authoritative source for Ensembl IDs and orthology data.

---

### Decision 3: Group Wälchli cells by ECclusters, not by patient age

**Problem:** The default Seurat identity labels in the Wälchli object were patient ages (15, 19 years etc.), not cell types.

**Decision:** Used the `ECclusters` metadata column, giving a 26,666 × 12 average expression matrix across biologically meaningful EC subtypes (Capillary, Arteriole, Artery, Venule, etc.).

---

### Decision 4: Exclude the Mitochondrial cluster from downstream analysis

**Reason:** The `Mitochondrial` cluster is a technical artefact — cells flagged for high mitochondrial gene expression are dying or low-quality cells. It does not represent a real endothelial cell type.

---

### Decision 5: Report both paralogue rows, flag orthology type

**Context:** Some human genes have multiple orthologues in mouse due to rodent-lineage gene duplications (e.g. human ABCB1 → mouse Abcb1a + Abcb1b). These produce multiple rows in the master table.

**Decision:** Retain all paralogue rows. Each has its own `perc_id` score. Flag them via the `orthology_type` column (`ortholog_one2many`). Step 5 will report both % identity scores for these genes rather than arbitrarily choosing one.

---

### Decision 6: CLDN5 is intentionally absent from the master list

**Finding:** CLDN5 (Claudin-5), the canonical BBB tight junction protein, does not appear in any of the Daneman or Munji source tables and is therefore absent from the master list.

**Reason:** Both papers built their lists by comparing brain ECs to other endothelia. CLDN5 is expressed in all endothelia (just most highly at the BBB) and failed the enrichment cutoff. This is a known limitation of enrichment-based gene lists — they capture what makes brain ECs *different*, not everything that is important to them.

**Implication:** The 74 genes with `Human_BBB_Datasets_N = 0` have no human EC evidence; CLDN5 represents the opposite problem — important but not captured. Both limitations should be stated in the methods section of any manuscript.

---

## The Role of Each Dataset in the Pipeline

| Dataset | Species | Role |
|---------|---------|------|
| Daneman 2010 | Mouse | Source of BBB gene list (S3, S4, S5 supplementary tables) |
| Munji 2019 | Mouse | Second source of BBB gene list (S5 — top CNS endothelial genes) |
| Wälchli 2024 | Human | Human EC validation (full expression matrix, Step 4b) |
| Yang 2022 | Human | Human EC validation (endothelial subtype markers, Step 4b) |
| Winkler 2022 | Human | Human EC validation (EC marker genes, Step 4b) |
| Human Protein Atlas | Human | Liver control gene set (Step 4c) |
| Giger 2010 | Human/Macaque | Future extension — expression in macaque endothelium |

Macaque does not contribute a gene list — it contributes a genome. The macaque orthologue of each BBB gene is identified via BioMart and its DNA sequence will be compared in Step 5.

---

## Project Milestones and Results

### Milestones

| # | Milestone | Status |
|---|-----------|--------|
| 1 | Master BBB gene list built (Daneman + Munji, 1,526 genes, filter-tagged) | ✅ |
| 2 | Wälchli single-cell EC expression processed (26,666 × 12 matrix) | ✅ |
| 3 | Cross-species Ensembl IDs mapped (human, mouse, macaque) | ✅ |
| 3b | Orthology confidence (type + % identity) added via BioMart | ✅ |
| 4b | Human EC validation against 3 datasets (Wälchli, Yang, Winkler) | ✅ |
| 4c | Liver control set built from HPA (7,131 genes) | ✅ |
| 5a-b | CDS sequences downloaded and extracted for all gene sets | ✅ |
| 5c-d | Pairwise alignment + % identity + dN/dS for BBB and liver | ✅ |
| 5e | Statistical comparison (Wilcoxon, stratified) | ✅ |
| 5f-g | dN/dS recomputed with codon-aware alignment + multiple-testing correction | ✅ |
| 5h | Direct mouse vs macaque comparison (closing the triangle) | ✅ |

### Headline Results

**Three-species conservation triangle (BBB genes, median % nucleotide identity):**

```
                Human
              /        \
          84.1%        96.8%
            /              \
        Mouse  ── 82.6% ── Macaque
```

- Macaque is the closest model to human BBB at the DNA level (96.8% identical).
- Mouse is more divergent (84.1%) but still substantially conserved.
- Mouse vs Macaque (82.6%) reflects their similar evolutionary distance from each other.

**BBB vs Liver — purifying selection (codon-aligned dN/dS):**

| Comparison | BBB dN/dS | Liver dN/dS | p (Bonferroni) | Interpretation |
|------------|-----------|-------------|----------------|----------------|
| Human vs Macaque | **0.152** | **0.196** | 1.97 × 10⁻⁶ | BBB under stronger purifying selection |
| Human vs Mouse | 0.117 | 0.119 | 0.504 | No detectable difference |
| Mouse vs Macaque | 0.128 | 0.129 | n.s. | No detectable difference |

The BBB-vs-liver conservation differential is real and statistically robust in the human-vs-macaque comparison. The signal is not visible in deeper evolutionary comparisons (mouse-human, mouse-macaque), likely because all genes have accumulated extensive change at that distance, washing out subtle BBB-specific signal.

**Paralogue effect (one-to-many orthologues):**

| Group | Median % identity | Median dN/dS |
|-------|-------------------|--------------|
| One-to-one orthologues | 85.0% | 0.107 |
| One-to-many (paralogues) | 47.3% | 0.354 |
| Effect size r | -0.72 (huge) | 0.59 (large) |

Gene duplications in the rodent lineage (e.g. ABCB1 → Abcb1a + Abcb1b) lead to substantial divergence between paralogues. Both paralogue rows are retained in the master table with their orthologue type flagged.

---

## Decisions Confirmed by Dr. Clelland

1. **Control gene set:** ✅ Liver-expressed genes (HPA) — more biologically interpretable than a random matched set.

2. **Human validation timing:** ✅ Done before Step 5 (Step 4b). Full list enters Step 5 with `Human_BBB_Datasets_N` as a stratification variable.

3. **Daneman filter:** ✅ Run Step 5 on S3 only first, then repeat on full list and compare.

4. **Alignment metrics:** ✅ Compute both % nucleotide identity and dN/dS. Report % identity as primary metric; dN/dS as supplementary.

5. **Regulatory regions:** ✅ Coding sequence only for now. Noncoding regions are under less evolutionary pressure and will differ substantially — flagged as future work.

6. **Extended analyses:** ✅ Start simple, validate results manually, then extend. Do not overextend before the primary analysis is confirmed.

---

## Things to Consider — Caveats and Future Work

These are honest limitations of the current pipeline. They do not invalidate the results but should be stated alongside any presentation.

### Methodological caveats

- **Pairwise dN/dS, not phylogenetic.** Current dN/dS values use the Nei-Gojobori method on pairs of sequences. A more rigorous analysis would use a multi-species phylogenetic model (PAML codeml) that accounts for the tree topology and lineage-specific rates. Pairwise gives a reasonable first approximation but is sensitive to local sequence ambiguities.
- **Longest-transcript heuristic.** When a gene has multiple isoforms, we used the longest CDS. A cleaner choice would be the MANE Select transcript (a curated canonical isoform per gene), which requires an extra lookup but is more standard.
- **Alignment parameters are defaults.** We used standard scoring (match +1, mismatch -1, gap -5/-2 for nucleotide; BLOSUM62 with gap 10/4 for protein). These are fine for closely related species but could be tuned for sensitivity in divergent regions.
- **Multiple testing was applied to primary tests only.** The four BBB-vs-liver primary comparisons were corrected (Bonferroni + FDR). Stratified analyses (by Daneman filter, validation score, orthologue type) were not corrected and should be treated as exploratory.

### Biological caveats

- **The starting gene list is enrichment-based, not exhaustive.** Daneman and Munji built their lists by comparing brain ECs to other endothelia. CLDN5 — the canonical BBB tight junction protein — is missing because it failed the differential cutoff. Other broadly endothelial genes important to the BBB may be similarly absent.
- **Liver control is one specific choice.** Liver was chosen for biological interpretability (a metabolically active, highly conserved tissue), but it has functional differences from BBB endothelium that could confound the comparison. A complementary control would be brain-expressed non-BBB genes, or genes matched on protein domain composition.
- **Coding sequence only.** Regulatory regions (promoters, enhancers, UTRs) are deliberately excluded. Many BBB-relevant cross-species differences may live in regulation rather than protein sequence.
- **HPA "brain" is whole brain, not endothelium.** When excluding brain-high genes from the liver control set, no genes were removed because HPA's brain sample is dominated by neuronal expression. A finer EC-specific filter would be cleaner.
- **The BBB-vs-liver signal is small in absolute terms.** Effect sizes (rank-biserial r) for the primary comparisons are under 0.1. Statistical significance comes largely from the large sample size. The biological story is real but subtle.

### Future extensions (in roughly increasing scope)

1. **Run PAML codeml** on the same alignments for proper phylogenetic dN/dS.
2. **Add promoter sequence comparison** (2 kb upstream of TSS) to test whether regulatory regions diverge more than coding regions, as expected.
3. **Re-run with MANE Select transcripts** instead of longest-transcript heuristic.
4. **Extend to additional control sets** — brain non-BBB genes; genes matched on domain composition or expression level.
5. **Bring the Giger macaque dataset back in** to ask whether the most-conserved BBB genes are also the most expressed in macaque endothelium.
6. **Compare different BBB cell types** (capillary vs arteriole vs venule) using the Wälchli per-cluster expression to ask whether conservation correlates with cell-type expression specificity.

---

## AI Transparency

I have decided to use **Claude** (specifically Sonnet 4.6) to assist with R scripting and pipeline development. Claude was not involved in:

- The project structure and deliverables
- Dataset search and selection
- Setting up evaluation metrics to ensure meaningful and correct R script results
- Which species to include, whether to use S3 only vs S3+S4+S5, reporting both paralogues, choosing DNA over protein alignment, adding human validation before Step 5
- Reading and confirming results against the source papers
- Catching anomalies (such as the CLDN5 absence)
- Biological interpretation of data

**Claude** was used for:

- Writing R scripts based on the approach I defined
- Debugging errors in the scripts
- Advising on table and file structure
- Markdown formatting and grammar

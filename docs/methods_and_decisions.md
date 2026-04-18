# Methods, Decisions, and Assumptions
*This document explains the scientific reasoning behind every major decision in this pipeline. It is not a code guide — the scripts and step notes cover that. This is about why we did things the way we did, where the uncertainties lie, and what still needs external validation.*

---

## The Central Scientific Question

**How similar are the genomes of human, rhesus macaque, and mouse across blood-brain barrier (BBB) genes?**

Specifically: are BBB genes more conserved across these three species than randomly selected genes? If yes, that supports the validity of mouse models for BBB research. If not, it raises questions about how well mouse BBB biology translates to humans and macaques.

---

## Scientific Assumptions

### Assumption 1: Mouse-derived BBB genes are valid proxies for human BBB genes

**What we did:** The master gene list was built from two mouse studies (Daneman 2010, Munji 2019). We converted mouse gene symbols to human HGNC format using HomoloGene (orthologue mapping) and treated the resulting list as the foundation for all downstream analysis.

**Why this is reasonable:** HomoloGene identifies orthologues — genes in different species that descended from the same ancestral gene and retained equivalent function. When it maps mouse *Slc2a1* to human *SLC2A1*, these genes share ~90 million years of evolutionary history and make functionally equivalent proteins. Using mouse BBB data as a starting point is standard practice in the field — it is the basis of nearly all BBB drug development research.

**Where the assumption could break down:** Some genes are BBB-enriched in mice but not in humans, or the expression pattern differs between species. We are not claiming the master list is a confirmed human BBB gene list — it is a candidate list based on mouse evidence.

**How we test it:** The Winkler 2022 and Yang 2022 human single-cell datasets are in the pipeline precisely to validate which genes in the master list are also expressed in human brain endothelial cells. Cross-referencing the master list against these human datasets is how we move from "mouse-derived candidate" to "human-validated BBB gene."

**Flag for Dr. Clelland:** Is this framing acceptable for the methods section? Should human validation against Winkler/Yang be done before or after the conservation analysis?

---

### Assumption 2: Orthologues found by HomoloGene represent genuine functional equivalents

**What we did:** Used the HomoloGene R package to convert mouse gene symbols to human HGNC symbols, and will use biomaRt in Step 3 to find macaque orthologues.

**Why this is reasonable:** Orthologue databases are the standard tool for cross-species gene mapping in genomics. They are built from genome-wide alignments and phylogenetic analysis, not guesswork.

**Where the assumption could break down:** HomoloGene has not been updated since ~2014. Some gene annotations, especially for macaque (*Macaca mulatta*), have changed substantially since then. This means some orthologues may be missing, outdated, or incorrectly mapped. 138 mouse genes from the master list had no human orthologue in HomoloGene and were silently dropped.

**How we address it:** Step 3 uses biomaRt (Ensembl), which is actively maintained and more current. Once Ensembl was back online, biomaRt replaced HomoloGene for the macaque orthologue lookup. A future action item is to re-run the mouse→human conversion using biomaRt and compare against the HomoloGene output to recover any dropped genes.

---

### Assumption 3: The Daneman_Filter column reflects meaningful biological stringency

**What we did:** Added a `Daneman_Filter` column to the master list with values `liver_and_lung`, `liver_only`, and `lung_only`, based on which Daneman supplementary table each gene came from (S3, S4, or S5 respectively).

**Why this is reasonable:** Daneman S3 required enrichment in brain ECs vs both liver AND lung simultaneously — a very strict double filter. S4 and S5 each tested only one peripheral tissue. A gene that passes both filters independently is a stronger BBB candidate than one that passes only one.

**Scientific value:** In the downstream conservation analysis, this column allows us to ask: are the most stringently defined BBB genes (liver_and_lung) also the most evolutionarily conserved? That is a more interesting question than treating all BBB genes as equal.

---

## Key Decisions Made

### Decision 1: Use Daneman S3, S4, and S5 — not just S3

**Initial approach:** We started with only S3 (most stringent), giving 545 genes with only 1 overlap with Munji.

**Problem identified:** After reading both papers and inspecting the files, we found that canonical BBB genes like *Slc2a1* (GLUT1) and *Cldn5* are absent from S3 because they are also expressed in liver endothelium and fail the double-comparison filter. Using S3 alone would exclude well-established BBB genes.

**Decision:** Include S3, S4, and S5 with filter labels. This expanded the list to 1,526 genes and increased the Daneman/Munji overlap from 1 to 80 genes — a far more biologically plausible result.

---

### Decision 2: Use HomoloGene instead of biomaRt for Step 1

**Reason:** Ensembl BioMart servers were unreachable at the time of analysis due to server downtime. HomoloGene is bundled locally in the R package and works offline.

**Trade-off:** HomoloGene is less current (last updated ~2014) but was the only option available. The macaque orthologue lookup in Step 3 uses biomaRt, which is more up-to-date.

**Action item:** Re-run mouse→human conversion with biomaRt when Ensembl is stable and compare outputs.

---

### Decision 3: Group Wälchli cells by ECclusters, not by patient age

**Problem:** The default Seurat identity labels in the Wälchli object were patient ages (15 years, 19 years, etc.), not cell types. Averaging by age would give us nothing biologically useful.

**Decision:** Used the `ECclusters` metadata column, which contains the actual cell type annotations (Capillary, Arteriole, Artery, Venule, etc.), giving a 26,666 × 12 expression matrix.

---

### Decision 4: Exclude the Mitochondrial cluster from downstream analysis

**Reason:** The `Mitochondrial` cluster in the Wälchli data is a technical artefact — cells with abnormally high mitochondrial gene expression, which is the standard indicator of dying or low-quality cells in single-cell sequencing. It does not represent a real endothelial cell type and should not be used in expression comparisons.

---

## The Role of Each Dataset in the Pipeline

| Dataset | Species | Role |
|---------|---------|------|
| Daneman 2010 | Mouse | Source of BBB gene list (master list foundation) |
| Munji 2019 | Mouse | Second source of BBB gene list (expands and validates Daneman) |
| Winkler 2022 | Human | Validates which master list genes are expressed in human BBB |
| Yang 2022 | Human | Second human validation source; average expression already computed |
| Wälchli 2024 | Human | Third human validation source; average expression computed in Step 2 |
| Giger 2010 | Human/Macaque | Background context; supports hypothesis that BBB genes are primate-conserved |

Macaque does not contribute a gene list — it contributes a genome. The macaque orthologue of each BBB gene is identified in Step 3 and its DNA sequence is compared in Step 5.

---

## Open Questions — Flag for Dr. Clelland

1. **Control gene set strategy (Step 5):** A randomly sampled set of ~3x genes matched for length and chromosome distribution is proposed as the control group. Dr. Clelland must approve this before the conservation comparison runs.

2. **Human validation timing:** Should we cross-reference the master list against Winkler/Yang expression data before or after the conservation analysis? Doing it before would give a cleaner, human-confirmed gene list going into Step 5. Doing it after keeps the analysis broader.

3. **Daneman S4/S5 inclusion:** Is it scientifically appropriate to include less stringent Daneman genes (liver_only, lung_only) in the conservation analysis, or should we restrict to S3 (liver_and_lung) only? The Daneman_Filter column makes it easy to filter either way.

4. **HomoloGene dropped genes:** 138 mouse genes had no human orthologue in HomoloGene. Is it worth recovering these using biomaRt before proceeding, or is the current 1,526-gene list sufficient?

5. **Giger macaque dataset:** After Step 5 is complete, the Giger data could be used to ask whether highly conserved BBB genes are also the most actively expressed in macaque endothelium. Should this be included in the final analysis?

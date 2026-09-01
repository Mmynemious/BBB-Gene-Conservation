# Data Audit and Weaknesses Review

*An independent audit of what data entered the pipeline at each step, what each step
actually did to it, and where the pipeline is weak. Every claim below was checked
against the files on disk and against the source papers, not against the project's
own documentation. Where this document contradicts `CONTEXT.md` or
`methods_and_decisions.md`, this document is the corrected version.*

**Date:** 2026-09-01
**Scope:** Steps 1 through 5h, all six raw datasets, all `processed/` outputs.

---

## Executive summary

Two findings change the results and need fixing before any further analysis:

1. **The Daneman gene list was built from the wrong supplementary tables.** The files
   named S3/S4/S5 in this repo are, per the paper's own legends, *pericyte-specific
   genes*, *developmentally up-regulated vascular genes*, and *developmentally
   down-regulated vascular genes*. Daneman's actual BBB-enriched gene list is the file
   named `S6_PostnatalBrainEC_Enriched.xls`, and it was never used. 78% of the master
   gene list traces to the three wrong tables.

2. **The liver control's brain-exclusion filter never ran.** It filtered on a tissue
   label (`"brain"`) that does not exist in the HPA consensus file, so it silently
   matched zero rows and excluded nothing. 86.6% of the "liver control" genes are also
   expressed in brain at nTPM >= 10. The control set is not liver-specific; it is a
   broad, housekeeping-dominated gene set.

Four further issues materially affect the numbers (orthologue-handling asymmetry
between the BBB and control sets, which flips the sign of the human-vs-mouse result;
a sheet-reading bug in the Yang validation; a near-vacuous Wälchli validation
criterion; and 92 genes lost to gene-symbol formatting). The remainder are ordinary
methodological limitations.

The codon-aware dN/dS implementation (Steps 5f/5h), the statistical framework, and
the macaque conservation result all survive scrutiny.

---

## Part 1 — Data provenance: what went in at each step

| Step | Script | Input files actually read | What the code did | Output |
|------|--------|---------------------------|-------------------|--------|
| 1 | `step1_extract_BBB_genelist.R` | `daneman_2010/` S3, S4, S5 (column 3 = Gene Symbol); `munji_2019/S5` sheet `BBB-enriched` (column `Gene`) | Union of mouse symbols; `Daneman_Filter` tag; mouse→human via **homologene** (offline, ~2014 vintage); dedup to one row per human symbol; `toupper()` | `master_BBB_genelist.csv` — 1,526 genes (Daneman 1,109 / Munji 337 / Both 80) |
| 2 | `step2_process_walchli_seurat.R` | `walchli_2024/*.rds.gz` — **not present in the repo** (gitignored) | Seurat `AverageExpression(slot = "data", group.by = "ECclusters")` | `Walchli2024_AverageExpression.csv` — 26,666 × 12 |
| 3 | `step3_map_ensembl_crossspecies.R` | `master_BBB_genelist.csv` + BioMart | Symbol → human Ensembl; human Ensembl → mouse and macaque orthologue IDs | `BBB_genes_crossspecies_ensembl.csv` — 1,868 rows |
| 4b | `step4b_human_validation.R` | master list; `Walchli2024_AverageExpression.csv`; `yang_2022/S6`; `winkler_2022/S2` | Three TRUE/FALSE flags plus `Human_BBB_Datasets_N` (0–3) | `master_BBB_genelist_validated.csv` — 1,526 |
| 3b | `step3b_biomart_orthology_confidence.R` | validated list + BioMart (release 113) | Re-queried orthologues; added `orthology_type` and `perc_id`; one row per orthologue pair | `master_BBB_genelist_complete.csv` — **1,724 rows** |
| 3c | `step3c_biotype_check.R` | master_complete + the three CDS FASTAs | Parsed `gene_biotype:` out of FASTA headers | `gene_biotype_check.csv`, `gene_biotype_summary.txt` |
| 4c | `step4c_liver_control_geneset.R` | `hpa_liver/rna_tissue_consensus.tsv.zip` + master_complete | liver nTPM >= 10; minus BBB genes; minus "brain-high" genes | `liver_control_genelist.csv` — 7,131 genes |
| 5a | `step5a_download_cds_fastas.R` | Ensembl **release 113** FTP | Downloaded CDS FASTAs: human GRCh38, mouse GRCm39, macaque Mmul_10 | `genomes/*.cds.all.fa.gz` (gitignored) |
| 5b | `step5b_extract_cds_sequences.R` | CDS FASTAs + master_complete + liver list | Longest CDS per gene; liver symbols matched to Ensembl IDs via FASTA headers | `genomes/cds_*.fa` |
| 5c | `step5c_align_and_score.R` | BBB CDS FASTAs + master_complete | Global nucleotide `pairwiseAlignment` → `pid(type = "PID1")`; dN/dS via `seqinr::kaks` on **truncated, unaligned** sequences | `conservation_scores_BBB.csv` — 1,632 rows |
| 5d | `step5d_liver_scoring.R` | liver CDS + BioMart orthologues + full mouse/macaque FASTAs | Same scoring, applied to liver. **Kept only the best-`perc_id` orthologue per gene** | `conservation_scores_liver.csv` — 7,091 rows |
| 5e | `step5e_statistics.R` | combined scores | Wilcoxon rank-sum + rank-biserial r | `conservation_summary.txt` |
| 5f | `step5f_recompute_dnds.R` | combined scores + FASTAs | Recomputed **dN/dS only** with codon-aware alignment (translate → BLOSUM62 align → back-map → strip gapped codons) | `conservation_scores_combined_v2.csv` — 8,723 rows |
| 5g | `step5g_stats_v2.R` | v2 scores | Wilcoxon + Bonferroni/FDR on the four primary tests | `conservation_summary_v2.txt` |
| 5h | `step5h_mouse_vs_macaque.R` | mouse/macaque FASTAs + v2 table | Closed the three-species triangle | `conservation_scores_mouse_macaque.csv` — 4,578 rows |

### Data collected but never used

- **Daneman:** S1, S2, **S6 (the real BBB-enriched list)**, S7
- **Munji:** S3 (present only as a dangling Git LFS pointer, 134 bytes), S7
- **Winkler:** S1, S3–S10 (including S10, listed in `CONTEXT.md` as a key file)
- **Yang:** S1–S5, S7, S8 (including S4, listed in `CONTEXT.md` as a key file)
- **Giger 2010:** the entire dataset

---

## Part 2 — Findings

### FINDING 1 (critical) — The Daneman gene list comes from the wrong tables

The supplementary files are numbered correctly (S1–S7 match the paper). The
**descriptive suffixes added to the filenames are wrong**, and the pipeline trusted
the suffixes rather than the paper.

Legends extracted from `docs/Daneman2010_Paper.pdf`, checked against the numeric
signature of each file:

| File in this repo | Paper's actual legend | Signature in the file (proof) |
|---|---|---|
| `S3_CoreBBBGenes_BrainEC_Enriched.xls` | **S3 — Pericyte specific genes** | `BrainG+/G+P-` > 2 in **95%** of rows |
| `S4_BrainEC_vs_LiverEC_Enriched.xls` | **S4 — Developmentally UP-regulated CNS vascular genes** | `Adult/Postnatal` > 2 in **95%** of rows |
| `S5_BrainEC_vs_LungEC_Enriched.xls` | **S5 — Developmentally DOWN-regulated CNS vascular genes** | `Adult/Postnatal` < 0.5 in **88%**; median `Brain/Lung` = **0.78** |
| `S6_PostnatalBrainEC_Enriched.xls` | **S6 — BBB enriched genes** *(the list the project needs)* | `Brain/Liver` > 2 in **99%** AND `Brain/Lung` > 2 in **99%** |
| `S7_AdultBrainEC_Enriched.xls` | **S7 — Peripheral endothelial enriched genes** *(the inverse of BBB)* | `Brain/Liver` < 0.5 in **97%**; `Brain/Lung` < 0.5 in **97%** |

Landmark genes confirm this independently:

- The file named S3 contains **Pdgfrb, Kcnj8, Vtn, Zic1** — canonical *pericyte*
  markers — and contains no BBB genes.
- The file named S6 contains **Slc2a1 (GLUT1), Slco1a4, Slc7a5, Ocln, Abcb1a, Rgs5**
  — canonical *BBB* genes.

**Consequences**

- 1,189 of 1,526 master-list genes (78%) derive from files S3/S4/S5.
- The 130 genes tagged `Daneman_Filter = liver_and_lung` — the "strictest BBB" stratum
  used for the stratified analysis in Steps 5e and 5g — are **pericyte genes**.
- The 567 genes tagged `lung_only` are genes CNS vessels *down-regulate* with age.
- `Daneman_Filter` as a stringency ordering (`liver_and_lung` > `liver_only` >
  `lung_only`) has no basis in the data. Assumption 3 in `methods_and_decisions.md`
  should be withdrawn.

**This explains an anomaly the project already recorded and mis-diagnosed.** Decision 1
in `methods_and_decisions.md` notes that S3-only gave a **1-gene** overlap with Munji,
concluded the filter was too stringent, and added S4 and S5 to raise the overlap to 80.
The real cause was the wrong table. Using Daneman's actual BBB list:

| | Genes | Overlap with Munji S5 (518 genes) |
|---|---|---|
| As built (S3 ∪ S4 ∪ S5) | 1,699 mouse symbols | 80 |
| **Corrected (S6)** | **411 mouse symbols** | **151 (37%)** |

A 37% overlap between two independent BBB transcriptome studies is what the biology
predicts. The corrected master list is **778 unique mouse genes** (S6 ∪ Munji S5),
not 1,526. Only 204 of the 1,699 genes currently in use appear in the real BBB list.

**Not affected by this finding:** CLDN5 is genuinely absent from Daneman S6 *and* from
Munji S5. Decision 6 in `methods_and_decisions.md` stands as written.

**Fix:** point Step 1 at `Daneman2010_S6_PostnatalBrainEC_Enriched.xls`; drop S3/S4/S5;
rename the files to match the paper's legends; correct the inventory in `CONTEXT.md`;
retire the `Daneman_Filter` column.

---

### FINDING 2 (critical) — The liver control's brain filter silently did nothing

`step4c_liver_control_geneset.R` filters `Tissue == "brain"`. **There is no tissue
called `brain` in the HPA RNA consensus file.** The labels are anatomical:
`cerebral cortex`, `hippocampal formation`, `basal ganglia`, `cerebellum`, `amygdala`,
`midbrain`, `hypothalamus`, `choroid plexus`, `spinal cord`, `retina`,
`pituitary gland`.

The filter matched zero rows. The subsequent `left_join` produced all-`NA`, and
`ifelse(is.na(Brain_nTPM), 0, Brain_nTPM)` converted every value to `0`. The
`Brain_nTPM < 10` condition then passed all 7,131 genes.

This is visible in the committed output: **every** `Brain_nTPM` value in
`liver_control_genelist.csv` is exactly `0.0`.

`methods_and_decisions.md` records this as *"no genes were removed because HPA's brain
sample is dominated by neuronal expression."* That explanation is incorrect — the
filter never executed.

**Recomputed with the real brain-region labels (max nTPM across regions):**

| Measure | Value |
|---|---|
| Control genes also expressed in brain at nTPM >= 10 | **6,176 / 7,129 (86.6%)** |
| Control genes also expressed in brain at nTPM >= 1 | 6,907 / 7,129 (96.9%) |
| Median max brain-region nTPM among control genes | **31.3** |
| Median liver nTPM among control genes | 31.2 |

The control group is therefore not "liver genes". It is ~7,100 broadly expressed genes
dominated by housekeeping genes — which are among the most conserved in the genome.
This undermines the framing of the primary BBB-vs-liver comparison and helps explain
why every effect size is below 0.1. Applying a correct filter leaves **953** genuinely
liver-biased genes.

**Fix:** filter on `max(nTPM)` across the real brain-region labels. Then decide with
Dr. Clelland whether the control should be the ~953 liver-biased genes, or stay broad
and be renamed honestly (e.g. "broadly expressed control").

---

### FINDING 3 (material) — BBB and liver genes were scored under different orthologue rules

The two gene sets were not treated the same way:

- **BBB** (`step5c`): retained *every* orthologue row returned by BioMart, including
  one-to-many. 196 of 1,541 scored rows (12.7%) are `ortholog_one2many`, and BioMart
  itself rates those pairs at a median of **25.8% protein identity** — they are distant
  family members, not orthologue pairs.
- **Liver** (`step5d`): applied `slice_max(PercId, n = 1)` — only the single
  best-matching orthologue per gene.

The BBB set therefore carries low-quality pairings that the control set does not.
Restricting both sides to one-to-one orthologues reverses the published human–mouse
result:

| Human vs Mouse, % identity | BBB | Liver | p | Direction |
|---|---|---|---|---|
| As published (all rows) | 84.12 (n=1,539) | 84.32 (n=5,265) | 0.022 | BBB **less** conserved |
| **One-to-one only** | **85.04** (n=1,327) | **84.47** (n=5,046) | **0.0075** | BBB **more** conserved |

The human–macaque result is robust to this (96.88 vs 96.22, p = 5.2e-05 — essentially
unchanged from the published 96.82 vs 96.15).

This also means the "paralogue effect" table in `methods_and_decisions.md`
(one2one 85.0% vs one2many 47.3%, r = -0.72) is partly measuring an artefact of
retaining non-orthologous pairs rather than a biological property of gene duplication.

**Fix:** apply the same orthologue-selection rule to both sets, and report the
one-to-one-only sensitivity analysis alongside the primary result either way.

---

### FINDING 4 (material) — The Yang validation read only one of seven sheets

`Yang2022_S6_EndothelialSubtype_Markers.xlsx` contains seven sheets: `Arterial`,
`Capillary`, `Venous`, `SMC`, `Pericyte`, `P. Fibroblast`, `Astrocyte`.
`read_xlsx(path)` with no `sheet =` argument reads **only the first**.

`In_Yang_EC` therefore means "is an *arterial* marker gene". The Capillary sheet — the
most BBB-relevant of the three endothelial sheets — was never read.

| Sheet | Genes | Overlap with master list |
|---|---|---|
| Arterial *(the only one read)* | 112 | **36** |
| Capillary | 82 | 28 |
| Venous | 97 | 34 |
| All seven sheets combined | 463 | 134 |

`notes_step4b.md` attributes the 2.4% coverage to Yang S6 being "a very stringent list";
the actual cause is that six sheets were skipped.

**Fix:** loop over `Arterial`, `Capillary`, `Venous` for an EC flag.

---

### FINDING 5 (material) — `In_Walchli_EC` is a detection test, not a validation test

The criterion is "average expression > 0 in any of the 12 EC clusters".
**88.5% of all 23,600 detected genes in the Wälchli matrix pass it** (23,600 of 26,666
genes have non-zero expression somewhere). It therefore carries almost no information
about BBB identity.

It is nonetheless the dominant contributor to the stratification variable: 1,131 of
1,526 genes score `Human_BBB_Datasets_N = 1` on the strength of this flag alone.

The three flags are also incommensurate and should not be summed as equal votes:

| Flag | Nature of the criterion | Pass rate on the master list |
|---|---|---|
| `In_Walchli_EC` | detection (expression > 0) | 95.2% |
| `In_Yang_EC` | membership of a 112-gene marker list (one sheet of seven) | 2.4% |
| `In_Winkler_EC` | membership of a 1,386-gene EC marker list | 20.6% |

**Fix:** use an expression threshold (or cluster-specificity) for Wälchli rather than
`> 0`, and either weight the three flags or report them separately instead of as a
0–3 count.

---

### FINDING 6 (material) — 92 genes lost to symbol formatting; Step 3b did not recover them

92 genes in `master_BBB_genelist_complete.csv` have `HumanEnsembl = NA` and are absent
from all conservation scoring. Two distinct causes:

**(a) 18 genes broken by the blanket `toupper()` in Step 1.** HGNC writes these symbols
with a lowercase `orf` (`C10orf10`). `toupper()` produced `C10ORF10`, which BioMart does
not recognise:

`C10ORF10`, `C11ORF31`, `C11ORF83`, `C12ORF57`, `C14ORF2`, `C16ORF52`, `C16ORF54`,
`C16ORF89`, `C17ORF49`, `C19ORF70`, `C1ORF115`, `C1ORF131`, `C1ORF198`, `C1ORF64`,
`C1ORF74`, `C5ORF45`, `C7ORF50`, `C8ORF4`

**(b) 74 genes are retired symbols** from homologene's ~2014 vocabulary, never
reconciled against current HGNC. All are recognisable renames — `CTGF`→CCN2,
`CYR61`→CCN1, `GPR56`→ADGRG1, `LPHN1`→ADGRL1, `GUCY1A3`→GUCY1A1, `SDPR`→CAVIN2,
`ELTD1`→ADGRL4 — plus the entire `ATP5*`, `H2AF*`, `HIST1*` and `KIAA*` families.

The header of `step3b_biomart_orthology_confidence.R` states that it "replaces the
HomoloGene-based Step 1 conversion for genes that were dropped". It does not: it
re-queries the same stale symbols against the same `hgnc_symbol` filter and recovers
none of them. Assumption 2 in `methods_and_decisions.md` overstates what Step 3b
achieved.

**Related:** Munji S5 ships mouse **Ensembl gene IDs** in column 2. Step 1 ignored them
and routed through symbol → homologene instead, which is where these losses originate.

**Fix:** drop the blanket `toupper()`; add an `external_synonym` fallback query; use
Munji's Ensembl IDs directly where available.

---

### FINDING 7 (methodological) — Smaller issues

- **`PID1` is length-sensitive.** `pid(type = "PID1")` is matches ÷ *alignment length*.
  Combined with the longest-transcript heuristic, a gene whose mouse longest isoform is
  twice the length of the human one is capped near 50% no matter how conserved it is.
  The one2many median of 47.3% (IQR 46.3–76.7, clustering tightly just above 47%) has
  that shape. `PID2` would separate real divergence from isoform-length mismatch.
  Note that **% identity was never recomputed in v2** — Step 5f corrected dN/dS only,
  so all reported % identity values still come from the original Step 5c/5d method.
- **Differential macaque coverage.** BBB genes have a macaque orthologue in 90% of rows
  (1,474 / 1,632); liver genes in only 52% (3,701 / 7,091). The two groups are missing
  data at very different — and probably non-random — rates.
- **The biotype audit (Step 3c) is circular.** Biotypes are read from the CDS FASTA,
  which by construction contains only protein-coding genes. A non-coding gene has no CDS
  entry, returns `NA`, and `NA` is not flagged as non-coding. The conclusion
  ("no annotation gap detected") is guaranteed by the method. 250 rows have no macaque
  biotype at all and are invisible to the check.
- **Ensembl release drift.** Step 3 used whichever release was live; `step3b`'s comments
  say release 113 while `methods_and_decisions.md` says 115; Steps 5a–5h use release 113
  FASTAs. Orthology assignments and sequences should come from a single pinned release.
- **`GeneSymbol_Mouse` is unreliable.** Step 1's
  `distinct(GeneSymbol_Human, .keep_all = TRUE)` after an `arrange()` keeps whichever
  mouse paralogue happens to sort first, with no record of the choice. For ABCB1 you get
  either `Abcb1a` or `Abcb1b` arbitrarily.
- **Two duplicate symbols** in `liver_control_genelist.csv`: `NPIPA9`, `MKKS`.
- **Reproducibility.** Every script hardcodes
  `setwd("~/Documents/Claude/Projects/BBB")`. The Wälchli `.rds.gz` and the `genomes/`
  FASTAs are gitignored, and `Munji2019_S3` is a dangling LFS pointer — Steps 2 and 5
  cannot be re-run from a fresh clone.

---

## Part 3 — What holds up

- **The codon-aware dN/dS (Steps 5f, 5h) is correctly implemented** — translate,
  align on protein with BLOSUM62, back-map to codons, strip gapped and ambiguous codons,
  then call `seqinr::kaks`. The v1 → v2 correction was necessary and real: it flipped
  the one2many mouse dN/dS from 0.14 to 0.354.
- **The statistical framework is appropriate** — Wilcoxon rank-sum for skewed
  distributions, rank-biserial effect sizes, Bonferroni and FDR on the four primary
  tests, with stratified analyses correctly labelled exploratory.
- **Munji 2019 S5 is genuinely a BBB-enriched list** and is parsed correctly. It is the
  only sound BBB gene source currently in the pipeline.
- **The Winkler S2 parse is correct** — `cluster == "EC"` selects 1,386 of 19,472 rows
  across 15 cell types, which is right.
- **The macaque conservation result is robust.** BBB 96.8% vs liver 96.2%
  (p = 1.8e-05), and dN/dS 0.152 vs 0.196 (Bonferroni p = 1.97e-06), both hold under
  every sensitivity check run here, including the one-to-one-only restriction.
- **The original anomaly detection was right.** The 1-gene Daneman/Munji overlap was
  correctly identified as implausible; only the diagnosis was wrong.

---

## Part 4 — Recommended order of repair

| # | Action | Changes conclusions? |
|---|---|---|
| 1 | Repoint Step 1 at Daneman S6; drop S3/S4/S5; rename files to match the paper; fix `CONTEXT.md`; retire `Daneman_Filter` | **Yes** |
| 2 | Fix the HPA brain filter; agree the control definition with Dr. Clelland | **Yes** |
| 3 | Make the orthologue-selection rule symmetric across BBB and control | **Yes** — flips the human–mouse direction |
| 4 | Yang: read all three EC sheets. Wälchli: use a real expression threshold | Numbers, not direction |
| 5 | Remove blanket `toupper()`; add synonym fallback; use Munji Ensembl IDs | Numbers, not direction |
| 6 | Switch to `PID2`, pin one Ensembl release, use MANE Select transcripts (`step5b_alt` already exists) | Numbers, not direction |
| 7 | Re-run Steps 5c–5h; rewrite the results section | — |

Items 1–3 must be done before any result is presented or written up.

---

## Part 5 — Process note

Both critical findings were preceded by a visible warning sign that was explained away
rather than investigated:

- The 1-gene Daneman/Munji overlap was treated as excessive filter stringency and
  "fixed" by adding two more wrong tables, which raised the overlap to 80 and made the
  problem look solved.
- The brain filter reporting `removed 0 brain-high genes` was rationalised in the
  methods document as a property of HPA's brain sample, rather than checked.

The `methods_and_decisions.md` AI Transparency section states that Claude was not
involved in anomaly-catching or biological interpretation. In practice, both of the
above rationalisations are recorded in that document as scientific reasoning, and both
are wrong. Anomalies that a proposed explanation makes disappear should receive more
scrutiny than others, not less.

**Suggested standing check for future data intake:** before any supplementary file is
read by a script, confirm its identity against the paper's own table legend *and*
against a numeric signature in the data (a defining ratio column, or the presence of
2–3 landmark genes the table must contain). Filenames and inventory documents are not
evidence.

---

## Appendix — How to reproduce these checks

`scripts/audit_source_tables.R` reproduces Findings 1, 2 and 4 from the raw files.

The diagnostic criteria for Finding 1, for manual verification in Excel:

| File | Column to check | Expected if the filename is right | Actual |
|---|---|---|---|
| S3 | `Brain/Liver` and `Brain/Lung` | both > 2 for most rows | 56% / 39%; `BrainG+/G+P-` > 2 for 95% |
| S4 | `Brain/Liver` | > 2 for most rows | 44%; `Adult/Postnatal` > 2 for 95% |
| S5 | `Brain/Lung` | > 2 for most rows | 5%; `Adult/Postnatal` < 0.5 for 88% |
| S6 | `Adult/Postnatal` | < 0.5 for most rows (if "postnatal-enriched") | 5%; `Brain/Liver` and `Brain/Lung` both > 2 for 99% |

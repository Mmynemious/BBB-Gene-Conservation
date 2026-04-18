# Step 1 Notes — Extract and Standardise the Master BBB Gene List

## What Step 1 Did

The two source papers (Daneman 2010 and Munji 2019) both studied BBB genes in **mice**. Mouse gene symbols are written differently from human ones — for example, the mouse gene is *Slc2a1* but the human equivalent is *SLC2A1*. Every downstream step in this project compares across species, so everything needs to be in one standard format: human HGNC (all uppercase).

Step 1 did three things:
1. Extracted the gene symbol columns from both Excel files (Daneman S3 and Munji S5)
2. Converted mouse symbols → human HGNC symbols using the HomoloGene database (which knows which mouse gene corresponds to which human gene — these are called *orthologues*)
3. Merged the two lists into one master file, tagging each gene as Daneman / Munji / Both

**Output:** `processed/master_BBB_genelist.csv` — 545 genes in human HGNC format, ready to be looked up in any genome database in downstream steps.

---

## Why Only Daneman and Munji?

These are the only two datasets in this project that provide a clean **mouse BBB gene list** derived from controlled enrichment experiments. Winkler 2022 and Yang 2022 provide human single-cell data, which enters the pipeline in Steps 3+. Daneman and Munji are the starting point because they explicitly compared brain endothelial cells to peripheral endothelial cells to identify BBB-specific genes.

---

## Confidence Tiers: Different Statistical Thresholds

The two studies use different statistical criteria, meaning the master list contains genes at **two different confidence levels**:

| Study | Technology | Statistical threshold | Enrichment comparison |
|-------|-----------|----------------------|----------------------|
| Daneman 2010 | Affymetrix microarray | FDR < 1% (SAM) | Brain EC vs liver AND lung EC |
| Munji 2019 | RNA-seq (DESeq2) | p < 0.05 | Brain EC vs heart, kidney, lung, AND liver EC |

Daneman S3 is the more stringent list — FDR < 1% is a harder cutoff than p < 0.05. When interpreting downstream conservation results, genes sourced from Daneman should be treated as higher-confidence BBB markers.

---

## Daneman_Filter Column

All Daneman genes currently in the master list come from **S3**, which requires enrichment in brain ECs vs **both** liver and lung simultaneously. This is marked in the `Daneman_Filter` column as `liver_and_lung`.

Future enhancement: Daneman S4 (brain vs liver only) and S5 (brain vs lung only) could be added with filter labels `liver_only` and `lung_only` respectively. This would allow downstream analyses to ask whether the most stringent BBB genes (both filters) are more evolutionarily conserved than the less stringent ones — a more nuanced scientific question than a flat gene list.

---

## Anomalies and Known Issues

### 1. Only 1 gene overlaps between Daneman and Munji

We found only 1 gene tagged as "Both" out of 685 unique mouse gene symbols. This is lower than expected given that both papers study mouse BBB genes using similar enrichment logic.

**Resolution after investigation:** This is likely real, not a processing error. Two reasons:

- **Different stringency:** Daneman S3 requires enrichment vs *both* liver and lung simultaneously (FDR < 1%), which is a very high bar. Canonical BBB genes like *Slc2a1* (GLUT1) and *Cldn5* are absent from Daneman S3 because they are also expressed in liver endothelium and therefore do not pass the double-comparison filter. Munji uses a less stringent threshold and more comparison organs.

- **Confirmed by file inspection:** A search of all 19 columns of the Daneman S3 file found no trace of *Slc2a1*, *Cldn5*, or *Abcb1a* in any capitalisation format. Their absence is not a read error — they genuinely did not pass Daneman's enrichment criteria.

**To flag to Dr. Clelland:** Whether the low overlap is scientifically acceptable for the conservation analysis, or whether a less stringent Daneman list (S4 or S5) should be incorporated.

---

### 2. HomoloGene database not updated since ~2014

The HomoloGene database used for mouse→human symbol conversion was quietly deprioritised by NCBI around 2014 — never officially retired but no longer maintained. This means:

- 138 mouse genes (685 input → 547 converted) had no human orthologue found and were silently dropped
- Some gene symbols or orthologue mappings from recent genome annotation updates may be missing or incorrect
- This is particularly relevant for macaque (*Macaca mulatta*), whose genome annotation has changed substantially since 2014

**Workaround used:** HomoloGene was used because Ensembl BioMart (the preferred, up-to-date alternative) was unreachable at the time of analysis due to server downtime.

**Action item:** Once Ensembl is accessible, re-run the mouse→human conversion using `biomaRt::getLDS()` and compare outputs. Any genes recovered by biomaRt but not HomoloGene should be added to the master list.

---

## Eval Checks (all passed)

1. All three required columns present (`GeneSymbol_Mouse`, `GeneSymbol_Human`, `Source`)
2. Gene count in expected range (545, within 300–800)
3. No duplicate human gene symbols
4. All human symbols are uppercase HGNC format
5. Source column contains only valid values (Daneman / Munji / Both)
6. Landmark gene check — WARN: CLDN5, TJP1, PECAM1, VWF absent (expected — these are pan-endothelial, not BBB-specific under Daneman's criteria)

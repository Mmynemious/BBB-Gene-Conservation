# Step 4b Notes — Human BBB Validation Cross-Reference

## What Step 4b Did

Cross-referenced the 1,526-gene master BBB list against three independent human brain endothelial cell (EC) datasets to add human expression evidence for each gene.

Three new columns were added to `processed/master_BBB_genelist_validated.csv`:

| Column | Source | Criterion |
|--------|--------|-----------|
| `In_Walchli_EC` | Wälchli 2024 | Expression > 0 in any EC cluster (full expression matrix) |
| `In_Yang_EC` | Yang 2022 S6 | Appears in the 112-gene endothelial subtype marker list |
| `In_Winkler_EC` | Winkler 2022 S2 | Appears in the EC marker gene list |
| `Human_BBB_Datasets_N` | — | Count of the three above (0–3) |

---

## Results

| N datasets | Gene count | Interpretation |
|-----------|-----------|----------------|
| 3 | 29 | Rock-solid human BBB genes |
| 2 | 292 | Strong human evidence |
| 1 | 1131 | Detected in at least one human EC dataset |
| 0 | 74 | No human EC evidence — came from mouse data only |

- **Wälchli coverage (95.2%)**: High because this is a full average-expression matrix across all 12 EC clusters — any detectable gene is counted.
- **Yang coverage (2.4%)**: Low because Yang S6 only contains 112 curated endothelial subtype marker genes — a very stringent list.
- **Winkler coverage (20.6%)**: EC-specific markers from a large human brain vascular atlas.

---

## Why These Three Datasets?

All three were chosen because they profile **human brain endothelial cells** specifically:

- **Wälchli 2024** — human adult brain sorted endothelial cells, the full expression landscape
- **Yang 2022** — human cerebrovascular atlas including endothelial subtypes
- **Winkler 2022** — large human brain vascular atlas identifying cell-type markers

Using human data as validation is important because the master BBB gene list was originally built from **mouse** data (Daneman 2010, Munji 2019). The validation columns answer: does each mouse-derived BBB gene also appear in human brain endothelium?

---

## Notable Finding: CLDN5 Is Absent From the Master List

CLDN5 (Claudin-5) — the most well-known BBB tight junction protein — is **not present** in the master gene list at all, and this is not a bug.

Daneman and Munji built their lists by comparing brain ECs to other endothelia (liver, lung, peripheral). CLDN5 is expressed in all endothelia, just more highly in the BBB. Because the enrichment cutoff required *differential* expression relative to other EC types, CLDN5 was excluded from both source lists.

This reveals a known limitation of enrichment-based BBB gene lists: they capture what makes brain ECs *different* from other ECs, but miss genes that are important to brain ECs yet shared with all endothelia. The 74 genes with `Human_BBB_Datasets_N = 0` represent the opposite problem — mouse-derived genes with no human evidence.

This limitation should be stated in the methods section of any manuscript.

---

## Downstream Use of Human_BBB_Datasets_N

The `Human_BBB_Datasets_N` column will allow downstream analyses to:
1. Report conservation statistics stratified by human validation confidence
2. Flag genes with `N = 0` as "mouse-only" with no human EC evidence
3. Provide a tier structure: N=3 (high confidence), N=2 (moderate), N=1 (low), N=0 (unvalidated in human)

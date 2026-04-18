# Analysis Methods & Pipeline

> **Status:** Draft methodology - to be refined after dataset approval
> 

---

## 📦 Data Processing Pipeline

### Step 1: Gene List Standardization

**Input:** BBB protein/gene lists from approved sources

**Process:**

- Map protein names to official gene symbols (HGNC, MGI)
- Cross-reference with Ensembl/NCBI gene databases
- Handle alternative gene names and synonyms
- Create unified gene identifier list

**Tools:**

- UniProt ID mapping
- Ensembl BioMart
- HGNC/MGI databases

### Step 2: Genomic Coordinate Extraction

**For each species:**

- Extract gene coordinates (chromosome, start, end, strand)
- Identify all transcript isoforms
- Annotate exons, introns, UTRs
- Define regulatory regions (promoter: -2kb to +500bp from TSS)

**Databases:**

- Human: Ensembl GRCh38
- Macaque: Ensembl Mmul_10
- Mouse: Ensembl GRCm39

### Step 3: Sequence Retrieval

**Extract sequences for:**

- Coding sequences (CDS)
- Intronic regions
- 5' and 3' UTRs
- Promoter regions (-2kb to +500bp)
- Enhancer regions (if annotated)

---

## 🔍 Comparative Analysis Methods

### Sequence Alignment Strategy

**Pairwise alignments:**

- Human vs. Macaque
- Human vs. Mouse
- Macaque vs. Mouse

**Tools:**

- BLAST+ for initial similarity
- MUSCLE/ClustalW for multiple sequence alignment
- MAFFT for large-scale alignments

### Conservation Scoring

**Metrics to calculate:**

1. **Sequence identity %** - Exact nucleotide matches
2. **Sequence similarity %** - Including conservative substitutions
3. **Alignment coverage %** - Proportion of sequence aligned
4. **Gap frequency** - Insertions/deletions

**Advanced metrics:**

- PhyloP conservation scores
- GERP++ constraint scores (if available)
- Synonymous vs. non-synonymous substitution rates

---

## 🎯 Statistical Analysis Plan

### Conservation Comparison

**Primary analysis:**

- Compare BBB gene conservation vs. control genes
- Statistical test: Wilcoxon rank-sum test
- Effect size: Cohen's d
- Multiple testing correction: Benjamini-Hochberg FDR

### Regional Analysis

**Compare conservation across:**

- Coding vs. non-coding regions
- Exons vs. introns
- Promoters vs. gene bodies
- Different functional gene categories

### Species-Specific Patterns

**Analyze:**

- Which species pair shows highest similarity?
- Are there BBB genes unique to primates?
- Mouse-specific BBB adaptations

---

## 📏 Control Group Design

### Control Strategy (Pending Dr. Clelland Input)

**Control Group 1: Random Genomic Background**

- Random sampling of genes matched for:
    - Chromosome distribution
    - Gene length
    - Exon number
    - Expression level (if data available)

**Control Group 2: Functional Controls**

- Other barrier genes (gut-blood, blood-retinal)
- General endothelial genes
- Tight junction genes (non-BBB)

**Control Group 3: Brain-Expressed Genes**

- Genes expressed in brain but not BBB-specific
- Matched for expression level and cell type

**Sample size:** 3x BBB gene number for robust statistics

---

## 📋 Output & Visualization Plan

### Summary Statistics Tables

1. **Overall conservation by region type**
2. **Species pairwise conservation matrices**
3. **Top/bottom conserved BBB genes**
4. **Functional category analysis**

### Visualizations

1. **Conservation heatmaps** - Genes vs. species pairs
2. **Box plots** - BBB vs. control conservation
3. **Scatter plots** - Coding vs. regulatory conservation
4. **Phylogenetic trees** - Based on BBB gene conservation
5. **Genome browser tracks** - Highlighting conserved regions

### Statistical Outputs

- P-values and effect sizes
- Confidence intervals
- Power analysis results
- Multiple testing corrections

---

## 🛠️ Computational Requirements

### Software Stack

- **Python 3.9+** with BioPython, pandas, numpy
- **R 4.0+** for statistical analysis and visualization
- **BLAST+ suite** for sequence alignment
- **Bedtools** for genomic interval operations
- **UCSC tools** for genome data handling

### Hardware Needs

- **RAM:** 16GB minimum (32GB preferred)
- **Storage:** 100GB for genome files and results
- **CPU:** Multi-core for parallel processing

### Expected Runtime

- Data collection: 2-3 days
- Sequence extraction: 1 day
- Alignments: 3-5 days (depending on gene number)
- Analysis: 2-3 days
- **Total:** ~2 weeks computational time

---

## 📚 Quality Control Measures

### Data Validation

- Cross-check gene IDs across databases
- Verify sequence lengths and coordinates
- Validate alignment quality scores
- Check for batch effects

### Analysis Validation

- Reproduce key findings with alternative methods
- Sensitivity analysis with different parameters
- Comparison with published conservation studies
- Peer review of analysis code

---

*Pipeline designed by: Yara Shamshoum*

*Subject to revision based on Dr. Clelland feedback*

*Last updated: [Current Date]*
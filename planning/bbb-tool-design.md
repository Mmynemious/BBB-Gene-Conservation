# Design proposal: turning the pipeline into a tool

**Status:** proposal, not committed work. Nothing here is built yet.

**Goal.** Take two or more gene sets as input, plus a choice of species, and return a
step-by-step PDF of the conservation analysis. Usable by the lab without writing code.

## Decisions taken

| Question | Decision |
|---|---|
| Who runs it | Yara and the lab, no coding required. Implies a simple upload-and-download interface, most likely Shiny |
| Species | Any Ensembl species, chosen at run time, rather than fixed human, mouse and macaque |
| Report | Step by step, including the diagnostics, so the report surfaces its own problems |

## The finding that shapes the design

The [audit](https://bbb-gene-conservation.vercel.app/audit.html) splits the 15 steps in two,
and the split is not where you would expect.

**Bespoke ingestion** — `step1`, `step2`, `step4b`, `step4c`. These parse specific
supplementary files from specific papers. **Two of the three breaking defects live here**: the
wrong Daneman source tables, and the liver filter that matched nothing.

**Generic engine** — `step3`, `step3b`, `step3c`, `step5a` through `step5h`. Orthologue
mapping, CDS retrieval, alignment, dN/dS, statistics. **Almost all of this is clean.**

So the clean steps are exactly the generalisable ones, and the buggy steps are exactly the
ones that cannot be generalised anyway.

### What that means for the three breaking defects

- **`step1` (wrong source tables)** and **`step4c` (filter never ran)** are input selection.
  In a tool the user supplies the gene set and the control set, so this code stops existing.
  These defects are not fixed, they are designed out.
- **`step5d` (different orthologue rules on the two sides)** is fixed **by construction**. If
  there is one scoring function called twice, once per gene set, the asymmetry cannot be
  written. This is the strongest single argument for refactoring rather than patching.

### Two engine defects that must be fixed, because they travel

- **`step5b`, longest-transcript heuristic.** Must default to MANE Select or Ensembl Canonical.
  `step5b_alt_canonical_transcripts.R` already implements this and was never switched in.
- **`step5c`, % identity never recomputed.** Identity must come from the same codon-aware
  alignment that produces dN/dS, not from a separate uncorrected pass.

## Proposed architecture

```
bbbtool/                      R package, the engine
  R/orthologues.R             symbols to Ensembl IDs to orthologues, any species pair
  R/sequences.R               CDS retrieval and caching, MANE Select selection
  R/align.R                   codon-aware alignment, % identity and dN/dS from one alignment
  R/score.R                   score_gene_set(): ONE function, called once per input set
  R/stats.R                   Wilcoxon, rank-biserial, BH and Bonferroni
  R/diagnostics.R             what went in, what came out, what was lost and why
  inst/report/report.qmd      parameterised step-by-step PDF
app/                          Shiny front end: upload, submit job, collect PDF
```

Command-line form, which is also what the app calls:

```bash
quarto render report.qmd \
  -P sets=focal.csv,control.csv \
  -P species=homo_sapiens,mus_musculus,macaca_mulatta
```

### Design rules carried over from the audit

1. **One scoring path for every input set.** No special casing of a control set.
2. **Every step returns data plus a diagnostics record**, and the report prints it. Genes lost
   at a step are counted and named, never dropped silently. This is the direct answer to the
   92 genes that vanished without an error.
3. **Checks that can fail.** `step3c` passed because a CDS file contains only coding genes by
   definition. Any check that cannot fail does not go in.
4. **Saturation is reported, not hidden.** The proportion of genes at the dS ceiling is a
   headline diagnostic, since it was 38% before the codon-aware fix and 2% after.
5. **Orthology type is carried through**, so one-to-many pairs can be included or excluded
   explicitly rather than by accident.

## The main risk: runtime and data volume

This is the constraint that shapes everything, and it needs solving before the UI is built.

- `step5c` and `step5d` take **hours** on a few thousand genes. A lab member cannot wait on a
  loading page for that.
- CDS FASTAs are **hundreds of megabytes per species**, so they cannot be fetched per request.

Implications:

- The app must **submit a job and return later**, with a link or an email when the PDF is
  ready. It cannot be a synchronous request.
- A **persistent genome cache** is needed, keyed by species and Ensembl release, shared across
  runs.
- Free Shiny hosting will not carry this. Hosting needs a real answer before the app stage:
  a lab machine, a university VM, or a paid container with a mounted cache volume.
- Alignment is embarrassingly parallel, so the per-gene loop should be parallelised. That is
  likely the difference between hours and minutes.

## Suggested stages

1. **Extract the engine.** Port `step3`, `step3b`, `step5a`, `step5b`, `step5f`, `step5g` into
   functions with explicit arguments. Remove `setwd()` and every hardcoded path and species.
   Fix the MANE default and the % identity source while porting.
2. **Parameterised report.** One Quarto document producing the step-by-step PDF from the
   engine's outputs and diagnostics. At this point the tool is usable from the command line.
3. **Validate against the current project.** Run the corrected engine on the existing BBB and
   liver gene lists and compare to the published numbers. This doubles as the pipeline fix the
   audit asks for, so the two pieces of work merge here rather than competing.
4. **Parallelise and cache**, until runtime is tolerable.
5. **Shiny front end**, only once 1 to 4 hold.

Stages 1 to 3 are worth doing regardless of whether the app is ever built, because they are
the same work as fixing the pipeline.

## Open questions

- Where would this be hosted, and who pays for it?
- Two gene sets, or arbitrarily many? The statistics are pairwise comparisons against a
  nominated control, so more than two needs a decision about which comparisons are run.
- Does the tool fetch from Ensembl live, which is reproducible only if the release is pinned,
  or ship with pinned release data?
- Is this a project deliverable, or a side build? It is a substantial piece of software
  engineering and it competes for time with finishing the analysis itself.

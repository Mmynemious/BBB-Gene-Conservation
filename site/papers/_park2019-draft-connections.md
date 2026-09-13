# Draft: how Park 2019 relates to this project

Written by Claude before Yara worked through the connections herself. Parked here
deliberately so it does not shape her own reading. Compare afterwards if useful.

## How it relates to this project

This is the most **different** paper in the review so far. It has no cross-species comparison,
no sequence, and no conservation measure. Its relevance is indirect but real.

1. **It is an argument against the pipeline this project is premised on.** The opening frame
   of this project is mouse, then macaque, then human, and the question of whether sequence
   conservation makes that chain trustworthy. Park's work attacks the same problem from the
   other end: build a **human** system in vitro and skip the species jump entirely. The paper
   states plainly that neither animal models nor existing in vitro cultures effectively mimic
   the human BBB. If human chips mature, cross-species conservation becomes less load-bearing
   as a translational argument, and more of a question about evolutionary biology. Worth
   saying somewhere in the project's framing.

2. **It independently nominates several genes on this project's list as functionally
   important**, and does so by *function* rather than by expression: P-gp (`ABCB1`), BCRP
   (`ABCG2`), MRP1 and MRP4, GLUT-1 (`SLC2A1`), plus LRP1 and the transferrin receptor as
   shuttling routes. These are not just detected here, they are shown to pump, and knocking
   each one out with a specific inhibitor changes what crosses.

3. **Claudin-5 again, and this time it matters.** Park treats claudin-5 as a defining marker
   of a working barrier. This project's gene list **does not contain CLDN5**, because it failed
   the differential-expression cutoff in the source studies: it is expressed in all endothelia,
   not only brain. That is a defensible reason for a *enrichment-based* list, but it means the
   list omits a gene the functional literature treats as central. See the
   [limitations](../results.qmd).

4. **It connects directly to [Kumabe 2025](kumabe2025.qmd) on in vitro model fidelity.**
   Kumabe found claudin-5 **absent** from hCMEC/D3, an immortalised line, and found
   immortalised lines drifting towards generic vein endothelium. Park's chip **does** express
   claudin-5 and reaches barrier values two orders of magnitude above primary-cell chips. Read
   together, the two papers say the same thing from opposite directions: how the endothelium is
   produced determines whether an in vitro BBB is a BBB at all.

5. **Function is not presence.** `step4b` scores a gene on whether it was *detected* in human
   datasets. Park's whole method rests on the difference between a transporter being present
   and a transporter working. That is a fair challenge to how this project's validation score
   is constructed, quite apart from the parsing bugs the [audit](../audit.qmd) records.

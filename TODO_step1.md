# Step 1 to-do list

Open items for the transporter target list (`scripts/v2/step1_transporter_targets.py`).
Written 25 October 2026.

## Verify

- [ ] **ABCG2 mouse partner.** Ensembl (release 113 and 116) pairs human ABCG2 with mouse
      Abcg3 (ENSMUSG00000029299, 54.9% identity, low confidence). NCBI pairs it with Abcg2
      (ENSMUSG00000029802). Check which is right, then pair by hand if needed, like ABCB1.

## Fix before Step 2 (scoring)

- [ ] **SLC10A4.** Ensembl pairs it with the mouse pseudogene Slc10a4-ps (0.33x human length).
      The real gene Slc10a4 (ENSMUSG00000029219) is protein-coding and NCBI pairs it correctly.
      Decide: pair with Slc10a4 by hand, or leave out until checked. Must not be scored as is.
- [ ] **Merged macaque gene models.** Add an alignment coverage check for every gene: align
      each macaque protein to its human protein and flag large stretches that match nothing.
      Known suspects: TAP1 (PSMB8), TAP2 (MAMU-DOB), ABCB6 (ATG9A), SLC10A3 (UBL4A),
      SLC25A44 (BGLAP), SLC39A3 (DIRAS1), SLC71A1 (SLC35A3).
      Fix order: another macaque isoform, then NCBI's macaque sequence, else drop from the
      three-species set only.
- [ ] **Incomplete CDS (12 genes).** Keep and flag the 11 with a missing start or stop codon.
      TAP2 has a broken reading frame (length not divisible by 3): treat it like a merged model.
      - Macaque: ABCC2, TAP2, OCA2, SLC13A5, SLC25A4, SLC26A7, SLC2A9, SLC37A2, SLC8A1
      - Mouse: SLC35D1, SLC6A1, SLC22A17 (also flagged in human, probably a genuine non-ATG start)

## Decide

- [ ] **The 49 genes with review flags:** keep all, drop all, or check one by one? And why?
- [ ] **9 genes still outside 0.8 to 1.25x human length**
      - Mouse: MPC1, TMEM184C, SLC10A4
      - Macaque: NIPAL4, SLC9A2, SLC16A11, SLC25A1, SLC25A41, SLC49A3
- [ ] **6 low-confidence orthologue calls:** ABCB6, MPC1, SLC10A4, SLC16A10, SLC22A16, SLC6A16

## Add to the write-up

- [ ] ABCB1 rule: rescued as the main BBB efflux pump; Abcb1a chosen because it is the copy in
      the mouse BBB datasets (Daneman 2010, Munji 2019). NCBI's own pick would be Abcb1b.
- [ ] 21 genes lost to mouse one-to-many (listed in `processed/v2/step1_exclusions.csv`); 2
      more lost on the macaque side.
- [ ] 18 genes with no one-to-one macaque partner (incl. ABCG2, SLC29A1); they stay in the
      human and mouse set.
- [ ] Macaque isoform rule: 14 transcripts swapped for the isoform closest to human length
      (length, not sequence similarity).
- [ ] NCBI macaque check is weak: 163 macaque genes have no NCBI ID in Ensembl.
- [ ] ABCE1, ABCF1, ABCF2, ABCF3 removed: no membrane domain, so not transporters.
- [ ] Versions: Ensembl 113; HGNC and NCBI downloaded 25 October 2026, copies in
      `raw_data/v2_step1/`.

## Housekeeping

- [ ] Commit and push the step card and step page changes (Method-Driven Pipeline list,
      `site/steps/`).
- [ ] Delete the stash of old home page edits once sure they are not needed (`git stash list`).
- [ ] Send the update email: the draft is in Gmail as a reply in the "schedule meeting" thread,
      to Dr. Pickard and Dr. Graf with Dr. Clelland in cc. Push the site first.

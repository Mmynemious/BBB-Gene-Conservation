# /// script
# requires-python = ">=3.10"
# dependencies = ["pandas>=2.0"]
# ///
"""
Step 1 (pipeline v2): SLC and ABC transporter target lists with one-to-one
orthologues, for two analyses.

Run from the project root:
    uv run scripts/v2/step1_transporter_targets.py

Decisions this script implements
  Gene set         All human genes in the HGNC "Solute carriers" (group 752) and
                   "ATP binding cassette transporters" (group 417) gene groups,
                   protein-coding only, minus the ABC members with no membrane
                   domain (ABCE1, ABCF1, ABCF2, ABCF3). Membership comes from
                   HGNC, not from a symbol regex, so SLC members with other names
                   (MFSD2A, UCP2, SPNS2, TAP1, CFTR ...) are included.
  Expression set   Human-mouse ortholog_one2one. This is the set calibrated
                   against Miao 2026, which compares human and mouse only.
  Mouse ABCB1      Exception to one-to-one: human ABCB1 is paired with mouse
                   Abcb1a. Mouse has two copies (Abcb1a, Abcb1b) from a rodent
                   duplication; this is documented in Methods.
  3-species set    The expression-set genes that are also human-macaque
                   ortholog_one2one. Used for sequence conservation across all
                   three species.
  Transcripts      Human MANE Select; mouse and macaque Ensembl Canonical.
                   Macaque override: when the canonical CDS is outside
                   0.8-1.25x the human MANE CDS length, use the macaque
                   protein-coding transcript whose CDS length is closest to human
                   MANE (complete CDSs preferred).

With one-to-one there is exactly one candidate gene per species, so no
"best orthologue" is ever chosen; % identity and dN/dS are measured later.
Every gene that drops out is written to the exclusions file with the stage,
the reason, and which set it is missing from.

Inputs  (downloaded on first run, then cached in raw_data/v2_step1/)
  HGNC gene group downloads (date stamped)
  Ensembl 113 BioMart (oct2024 archive): genes, transcripts, homologues
  genomes/*.113.cds.all.fa.gz
Outputs
  processed/v2/step1_expression_set.csv      human-mouse, one row per gene
  processed/v2/step1_three_species_set.csv   human-mouse-macaque subset
  processed/v2/step1_exclusions.csv
  processed/v2/step1_summary.txt
"""

import datetime
import gzip
import hashlib
import io
import urllib.parse
import urllib.request
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
RAW = ROOT / "raw_data" / "v2_step1"
OUT = ROOT / "processed" / "v2"
RAW.mkdir(parents=True, exist_ok=True)
OUT.mkdir(parents=True, exist_ok=True)

BIOMART = "https://oct2024.archive.ensembl.org/biomart/martservice"  # Ensembl 113
HGNC_GROUPS = {"SLC": 752, "ABC": 417}
NON_MEMBRANE_ABC = {"ABCE1", "ABCF1", "ABCF2", "ABCF3"}
# human gene -> (chosen mouse gene, why); used where one-to-one does not hold
MOUSE_OVERRIDES = {
    "ENSG00000085563": ("ENSMUSG00000040584",
                        "ABCB1 paired with Abcb1a; mouse has Abcb1a and Abcb1b (rodent duplication)"),
}
LENGTH_RANGE = (0.8, 1.25)
CDS_FASTA = {
    "human": ROOT / "genomes" / "Homo_sapiens.GRCh38.113.cds.all.fa.gz",
    "mouse": ROOT / "genomes" / "Mus_musculus.GRCm39.113.cds.all.fa.gz",
    "macaque": ROOT / "genomes" / "Macaca_mulatta.Mmul_10.113.cds.all.fa.gz",
}
STOPS = {"TAA", "TAG", "TGA"}


# ── Download helpers (cached so reruns are reproducible) ─────────────────────

def fetch(url: str, dest: Path, data: bytes | None = None) -> Path:
    if not dest.exists():
        print(f"  downloading -> {dest.relative_to(ROOT)}")
        with urllib.request.urlopen(url, data=data, timeout=300) as r:
            dest.write_bytes(r.read())
    return dest


def biomart(name: str, dataset: str, attributes: list[str], gene_ids: list[str]) -> pd.DataFrame:
    """One BioMart query, filtered on Ensembl gene IDs, cached as TSV.
    The cache name carries a hash of the query, so a changed gene list refetches."""
    attrs = "".join(f'<Attribute name="{a}"/>' for a in attributes)
    xml = (
        '<?xml version="1.0" encoding="UTF-8"?><!DOCTYPE Query>'
        '<Query virtualSchemaName="default" formatter="TSV" header="1" uniqueRows="1">'
        f'<Dataset name="{dataset}" interface="default">'
        f'<Filter name="ensembl_gene_id" value="{",".join(sorted(set(gene_ids)))}"/>'
        f"{attrs}</Dataset></Query>"
    )
    key = hashlib.sha1(xml.encode()).hexdigest()[:8]
    path = fetch(BIOMART, RAW / f"biomart113_{name}_{key}.tsv",
                 urllib.parse.urlencode({"query": xml}).encode())
    text = path.read_text()
    if "ERROR" in text[:200]:
        path.unlink()
        raise RuntimeError(f"BioMart error for {name}: {text[:300]}")
    df = pd.read_csv(io.StringIO(text), sep="\t", dtype=str).fillna("")
    df.columns = attributes  # BioMart returns display names; use the attribute names
    return df


def read_cds(path: Path, wanted: set[str]) -> dict[str, str]:
    """Return {transcript_id_without_version: CDS} for the wanted transcripts."""
    seqs, tid, chunks = {}, None, []
    with gzip.open(path, "rt") as f:
        for line in f:
            if line.startswith(">"):
                if tid in wanted:
                    seqs[tid] = "".join(chunks)
                tid, chunks = line[1:].split()[0].split(".")[0], []
            else:
                chunks.append(line.strip())
    if tid in wanted:
        seqs[tid] = "".join(chunks)
    return seqs


def cds_problem(seq: str | None) -> str:
    """Empty string if the CDS is complete and clean, otherwise what is wrong."""
    if seq is None:
        return "not in CDS FASTA"
    seq = seq.upper()
    problems = []
    if len(seq) % 3:
        problems.append("length not divisible by 3")
    if not seq.startswith("ATG"):
        problems.append("no ATG start")
    if seq[-3:] not in STOPS:
        problems.append("no terminal stop")
    codons = [seq[i:i + 3] for i in range(0, len(seq) - 3, 3)]
    if any(c in STOPS for c in codons):
        problems.append("internal stop")
    return "; ".join(problems)


exclusions = []


def exclude(rows: pd.DataFrame, stage: str, missing_from: str, reason) -> None:
    for _, r in rows.iterrows():
        exclusions.append({
            "symbol": r["symbol"], "family": r["family"], "hgnc_subfamily": r["hgnc_subfamily"],
            "human_gene_id": r["human_gene_id"], "stage": stage, "missing_from": missing_from,
            "reason": reason(r) if callable(reason) else reason,
        })


def orthologues(h: pd.DataFrame, prefix: str) -> pd.DataFrame:
    """Per human gene: number of orthologues, their types, and the partner if unique."""
    h = h[h[f"{prefix}_homolog_ensembl_gene"] != ""]
    return h.groupby("ensembl_gene_id").agg(
        n=(f"{prefix}_homolog_ensembl_gene", "size"),
        types=(f"{prefix}_homolog_orthology_type", lambda s: ",".join(sorted(set(s)))),
        candidates=(f"{prefix}_homolog_ensembl_gene", lambda s: ",".join(sorted(s))),
        partner=(f"{prefix}_homolog_ensembl_gene", "first"),
        confidence=(f"{prefix}_homolog_orthology_confidence", "first"),
    )


def orthology_reason(sp: str):
    return lambda r: (f"no {sp} orthologue" if pd.isna(r.n)
                      else f"{sp}: {r.types} ({int(r.n)} {sp} genes)")


def canonical(tx: pd.DataFrame, sp: str) -> pd.DataFrame:
    return tx[tx.transcript_is_canonical == "1"].groupby("ensembl_gene_id").agg(
        n=("ensembl_transcript_id", "size"),
        **{f"{sp}_canonical_transcript_id": ("ensembl_transcript_id", "first"),
           f"{sp}_symbol": ("external_gene_name", "first")},
    )


# ── Stage 0: gene set from HGNC gene groups ──────────────────────────────────

print("Stage 0: HGNC gene groups")
today = datetime.date.today().isoformat()
frames = []
for family, gid in HGNC_GROUPS.items():
    cached = sorted(RAW.glob(f"hgnc_group{gid}_*.tsv"))
    path = cached[-1] if cached else fetch(
        f"https://www.genenames.org/cgi-bin/genegroup/download?id={gid}&type=branch",
        RAW / f"hgnc_group{gid}_{today}.tsv")
    g = pd.read_csv(path, sep="\t", dtype=str).fillna("")
    g["family"] = family
    frames.append(g)
hgnc = pd.concat(frames, ignore_index=True).rename(columns={
    "Approved symbol": "symbol", "Locus type": "locus_type",
    "Ensembl gene ID": "human_gene_id", "Group name": "hgnc_subfamily",
    "HGNC ID": "hgnc_id", "NCBI Gene ID": "human_ncbi_id",
})[["hgnc_id", "symbol", "family", "hgnc_subfamily", "locus_type", "human_gene_id", "human_ncbi_id"]]
hgnc = hgnc.drop_duplicates("hgnc_id")
print(f"  {len(hgnc)} HGNC entries ({(hgnc.family == 'SLC').sum()} SLC, {(hgnc.family == 'ABC').sum()} ABC)")

genes = hgnc
bad = genes[genes.locus_type != "gene with protein product"]
exclude(bad, "0 gene set", "both", lambda r: f"not protein-coding in HGNC ({r.locus_type})")
genes = genes[genes.locus_type == "gene with protein product"]
bad = genes[genes.symbol.isin(NON_MEMBRANE_ABC)]
exclude(bad, "0 gene set", "both", "ABC family member with no membrane domain (not a transporter)")
genes = genes[~genes.symbol.isin(NON_MEMBRANE_ABC)]
bad = genes[genes.human_gene_id == ""]
exclude(bad, "0 gene set", "both", "no Ensembl gene ID in HGNC")
genes = genes[genes.human_gene_id != ""].copy()
print(f"  {len(genes)} protein-coding transporter genes with an Ensembl ID")

# ── Stage 1a: human-mouse one-to-one (expression set) ────────────────────────

print("Stage 1: orthology (Ensembl 113)")
ids = genes.human_gene_id.tolist()
human_tx = biomart("human_transcripts", "hsapiens_gene_ensembl",
                   ["ensembl_gene_id", "external_gene_name", "chromosome_name",
                    "ensembl_transcript_id", "transcript_mane_select",
                    "transcript_is_canonical"], ids)
homolog_attrs = lambda p: ["ensembl_gene_id", f"{p}_homolog_ensembl_gene",
                           f"{p}_homolog_orthology_type", f"{p}_homolog_orthology_confidence"]
mouse_h = biomart("human_mouse_homologs", "hsapiens_gene_ensembl", homolog_attrs("mmusculus"), ids)
macaque_h = biomart("human_macaque_homologs", "hsapiens_gene_ensembl", homolog_attrs("mmulatta"), ids)

bad = genes[~genes.human_gene_id.isin(set(human_tx.ensembl_gene_id))]
exclude(bad, "1 orthology", "both", "HGNC Ensembl ID not found in Ensembl 113")
genes = genes[genes.human_gene_id.isin(set(human_tx.ensembl_gene_id))]

genes = genes.join(orthologues(mouse_h, "mmusculus"), on="human_gene_id")
genes["mouse_orthology_note"] = ""
for human_id, (mouse_id, note) in MOUSE_OVERRIDES.items():
    row = genes.human_gene_id == human_id
    assert mouse_id in genes.loc[row, "candidates"].iloc[0].split(","), \
        f"override {mouse_id} is not an Ensembl orthologue of {human_id}"
    genes.loc[row, ["n", "partner", "mouse_orthology_note"]] = [1, mouse_id, note]
    genes.loc[row, "types"] = "ortholog_one2one"  # treated as one-to-one from here on
ok = (genes.n == 1) & (genes.types == "ortholog_one2one")
exclude(genes[~ok], "1 orthology", "both", orthology_reason("mouse"))
genes = genes[ok].rename(columns={"partner": "mouse_gene_id",
                                  "confidence": "mouse_orthology_confidence"}) \
    .drop(columns=["n", "types", "candidates"])
print(f"  {len(genes)} genes one-to-one human-mouse")

# ── Stage 2a: human MANE Select, mouse Ensembl Canonical ─────────────────────

print("Stage 2: transcripts")
mane = human_tx[human_tx.transcript_mane_select != ""].groupby("ensembl_gene_id").agg(
    n=("ensembl_transcript_id", "size"),
    human_transcript_id=("ensembl_transcript_id", "first"),
    human_mane_refseq=("transcript_mane_select", "first"),
)
genes = genes.join(mane, on="human_gene_id")
bad = genes[genes.n != 1]
exclude(bad, "2 transcripts", "both", lambda r: "no human MANE Select transcript"
        if pd.isna(r.n) else f"{int(r.n)} human MANE Select transcripts")
genes = genes[genes.n == 1].drop(columns="n")

mouse_tx = biomart("mouse_transcripts", "mmusculus_gene_ensembl",
                   ["ensembl_gene_id", "external_gene_name", "ensembl_transcript_id",
                    "transcript_is_canonical"], genes.mouse_gene_id.tolist())
genes = genes.join(canonical(mouse_tx, "mouse"), on="mouse_gene_id")
bad = genes[genes.n != 1]
exclude(bad, "2 transcripts", "both", lambda r: "no mouse Ensembl Canonical transcript"
        if pd.isna(r.n) else f"{int(r.n)} mouse canonical transcripts")
genes = genes[genes.n == 1].drop(columns="n") \
    .rename(columns={"mouse_canonical_transcript_id": "mouse_transcript_id"})
print(f"  expression set: {len(genes)} genes (human MANE + mouse canonical)")

# ── Stage 1b: human-macaque one-to-one (three-species subset) ────────────────

genes = genes.join(orthologues(macaque_h, "mmulatta"), on="human_gene_id")
ok = (genes.n == 1) & (genes.types == "ortholog_one2one")
exclude(genes[~ok], "1b macaque orthology", "three-species set only", orthology_reason("macaque"))
genes["in_three_species_set"] = ok
genes.loc[~ok, ["partner", "confidence"]] = ""
genes = genes.rename(columns={"partner": "macaque_gene_id",
                              "confidence": "macaque_orthology_confidence"}) \
    .drop(columns=["n", "types", "candidates"])

# ── Stage 2b: macaque Ensembl Canonical, with the length-based override ──────

macaque_tx = biomart("macaque_transcripts", "mmulatta_gene_ensembl",
                     ["ensembl_gene_id", "external_gene_name", "ensembl_transcript_id",
                      "transcript_is_canonical", "transcript_biotype"],
                     genes.loc[genes.in_three_species_set, "macaque_gene_id"].tolist())
genes = genes.join(canonical(macaque_tx, "macaque"), on="macaque_gene_id")
bad = genes[genes.in_three_species_set & (genes.n != 1)]
exclude(bad, "2b macaque transcripts", "three-species set only",
        lambda r: "no macaque Ensembl Canonical transcript"
        if pd.isna(r.n) else f"{int(r.n)} macaque canonical transcripts")
genes.loc[genes.in_three_species_set & (genes.n != 1), "in_three_species_set"] = False
genes = genes.drop(columns="n")

# CDS for every chosen transcript, plus every macaque protein-coding transcript
# of the three-species genes (candidates for the override)
macaque_pc = macaque_tx[macaque_tx.transcript_biotype == "protein_coding"]
seqs = {
    "human": read_cds(CDS_FASTA["human"], set(genes.human_transcript_id)),
    "mouse": read_cds(CDS_FASTA["mouse"], set(genes.mouse_transcript_id)),
    "macaque": read_cds(CDS_FASTA["macaque"], set(macaque_pc.ensembl_transcript_id)),
}
genes["human_cds_length"] = genes.human_transcript_id.map(lambda t: len(seqs["human"].get(t, "")))
macaque_by_gene = macaque_pc.groupby("ensembl_gene_id").ensembl_transcript_id.apply(list)


def pick_macaque(r) -> pd.Series:
    """Canonical, unless its CDS length is outside LENGTH_RANGE x human MANE;
    then the protein-coding transcript closest to human MANE length."""
    if not r.in_three_species_set:
        return pd.Series(["", ""])
    canon = r.macaque_canonical_transcript_id
    ratio = len(seqs["macaque"].get(canon, "")) / r.human_cds_length
    if LENGTH_RANGE[0] <= ratio <= LENGTH_RANGE[1]:
        return pd.Series([canon, "canonical"])
    options = [t for t in macaque_by_gene.get(r.macaque_gene_id, []) if t in seqs["macaque"]]
    best = min(options, key=lambda t: (cds_problem(seqs["macaque"][t]) != "",
                                       abs(len(seqs["macaque"][t]) - r.human_cds_length)))
    if best == canon:
        return pd.Series([canon, f"canonical (no closer isoform; canonical {ratio:.2f}x human)"])
    return pd.Series([best, f"closest to human MANE length (canonical was {ratio:.2f}x human)"])


genes[["macaque_transcript_id", "macaque_transcript_rule"]] = genes.apply(pick_macaque, axis=1)

for sp in CDS_FASTA:
    tid = genes[f"{sp}_transcript_id"]
    genes[f"{sp}_cds_length"] = tid.map(lambda t: len(seqs[sp].get(t, "")) if t else float("nan"))
    genes[f"{sp}_cds_problem"] = tid.map(lambda t: cds_problem(seqs[sp].get(t)) if t else "")
for sp in ["mouse", "macaque"]:
    genes[f"{sp}_length_ratio"] = (genes[f"{sp}_cds_length"].astype(float)
                                   / genes.human_cds_length).round(2)
genes["cds_complete_human_mouse"] = (genes.human_cds_problem == "") & (genes.mouse_cds_problem == "")
genes["cds_complete_all"] = genes.cds_complete_human_mouse & (genes.macaque_cds_problem == "") \
    & genes.in_three_species_set

# ── Cross-check against NCBI orthologues (flag only) ─────────────────────────
# Ensembl gene trees occasionally pair the wrong paralogue (Ensembl 113 makes
# human ABCG2 one-to-one with mouse Abcg3 and gives mouse Abcg2 no orthologue).
# Miao 2026 cross-checked Ensembl against NCBI Gene; we do the same and flag
# every disagreement rather than silently switching partners.

print("Cross-check: NCBI gene_orthologs")
TAXA = {"human": "9606", "mouse": "10090", "macaque": "9544"}
ncbi_path = RAW / f"ncbi_gene_orthologs_human_mouse_macaque_{today}.tsv"
cached = sorted(RAW.glob("ncbi_gene_orthologs_human_mouse_macaque_*.tsv"))
if cached:
    ncbi_path = cached[-1]
else:
    print("  streaming ftp.ncbi.nlm.nih.gov/gene/DATA/gene_orthologs.gz (keeps human-mouse/macaque rows)")
    keep = {(TAXA["human"], TAXA["mouse"]), (TAXA["human"], TAXA["macaque"]),
            (TAXA["mouse"], TAXA["human"]), (TAXA["macaque"], TAXA["human"])}
    with urllib.request.urlopen("https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene_orthologs.gz",
                                timeout=600) as r, gzip.open(r, "rt") as f, open(ncbi_path, "w") as out:
        out.write("human_ncbi_id\tother_taxon\tother_ncbi_id\n")
        for line in f:
            p = line.rstrip("\n").split("\t")
            if (p[0], p[3]) in keep:
                h, o = (p[1], p[4]) if p[0] == TAXA["human"] else (p[4], p[1])
                out.write(f"{h}\t{p[3] if p[0] == TAXA['human'] else p[0]}\t{o}\n")
ncbi = pd.read_csv(ncbi_path, sep="\t", dtype=str).drop_duplicates()

for sp, dataset in [("mouse", "mmusculus_gene_ensembl"), ("macaque", "mmulatta_gene_ensembl")]:
    partners = genes.loc[genes[f"{sp}_gene_id"] != "", f"{sp}_gene_id"].tolist()
    xref = biomart(f"{sp}_ncbi_xref", dataset, ["ensembl_gene_id", "entrezgene_id"], partners)
    ens_to_ncbi = xref[xref.entrezgene_id != ""].groupby("ensembl_gene_id").entrezgene_id.apply(set)
    sp_ncbi = ncbi[ncbi.other_taxon == TAXA[sp]].groupby("human_ncbi_id").other_ncbi_id.apply(set)

    def agrees(r, sp=sp, ens_to_ncbi=ens_to_ncbi, sp_ncbi=sp_ncbi) -> str:
        if not r[f"{sp}_gene_id"]:
            return ""
        ncbi_partners = sp_ncbi.get(r.human_ncbi_id, set())
        ens_partner = ens_to_ncbi.get(r[f"{sp}_gene_id"], set())
        if not ncbi_partners:
            return "no NCBI orthologue"
        if not ens_partner:
            return "Ensembl partner has no NCBI ID"
        return "agrees" if ens_partner & ncbi_partners else \
            f"disagrees (NCBI partner GeneID {','.join(sorted(ncbi_partners))})"

    genes[f"{sp}_ncbi_check"] = genes.apply(agrees, axis=1)


# Review flags (do not filter). Canonical transcripts outside human MANE are less
# curated: a gene model merged with a neighbouring gene shows up as the
# orthologue carrying the neighbour's symbol (e.g. TAP1 -> PSMB8 in macaque).
def review_flags(r) -> str:
    flags = []
    species = ["mouse", "macaque"] if r.in_three_species_set else ["mouse"]
    for sp in species:
        ratio = r[f"{sp}_length_ratio"]
        if not LENGTH_RANGE[0] <= ratio <= LENGTH_RANGE[1]:
            flags.append(f"{sp} CDS {ratio}x human")
        if str(r[f"{sp}_symbol"]).upper() != r.symbol.upper():
            flags.append(f"{sp} symbol {r[f'{sp}_symbol']}")
        if r[f"{sp}_orthology_confidence"] == "0":
            flags.append(f"{sp} low-confidence orthology")
        if r[f"{sp}_ncbi_check"].startswith("disagrees"):
            flags.append(f"{sp} NCBI {r[f'{sp}_ncbi_check']}")
    for sp in ["human"] + species:
        if r[f"{sp}_cds_problem"]:
            flags.append(f"{sp} CDS: {r[f'{sp}_cds_problem']}")
    return "; ".join(flags)


genes["review_flags"] = genes.apply(review_flags, axis=1)

# ── Write outputs ────────────────────────────────────────────────────────────

three = genes[genes.in_three_species_set]
assert genes.human_gene_id.is_unique and genes.mouse_gene_id.is_unique \
    and three.macaque_gene_id.is_unique, "a gene ID appears twice; one-to-one is broken"

cols = ["symbol", "family", "hgnc_subfamily", "hgnc_id",
        "human_gene_id", "human_transcript_id", "human_mane_refseq",
        "mouse_symbol", "mouse_gene_id", "mouse_transcript_id",
        "mouse_orthology_confidence", "mouse_orthology_note", "mouse_ncbi_check",
        "in_three_species_set",
        "macaque_symbol", "macaque_gene_id", "macaque_canonical_transcript_id",
        "macaque_transcript_id", "macaque_transcript_rule", "macaque_orthology_confidence",
        "macaque_ncbi_check",
        "human_cds_length", "mouse_cds_length", "macaque_cds_length",
        "mouse_length_ratio", "macaque_length_ratio",
        "human_cds_problem", "mouse_cds_problem", "macaque_cds_problem",
        "cds_complete_human_mouse", "cds_complete_all", "review_flags"]
genes = genes[cols].sort_values(["family", "symbol"])
genes.to_csv(OUT / "step1_expression_set.csv", index=False)
genes[genes.in_three_species_set].drop(columns="in_three_species_set") \
    .to_csv(OUT / "step1_three_species_set.csv", index=False)
excl = pd.DataFrame(exclusions).sort_values(["family", "stage", "symbol"])
excl.to_csv(OUT / "step1_exclusions.csv", index=False)
for old in ["step1_transporter_targets.csv", "step1_transporter_exclusions.csv"]:
    (OUT / old).unlink(missing_ok=True)  # superseded outputs of the first version

three = genes[genes.in_three_species_set]
lines = [f"Step 1 transporter targets  (run {today}, Ensembl 113)", ""]
for fam in ["SLC", "ABC"]:
    e = excl[excl.family == fam]
    lines.append(f"{fam}: {(hgnc.family == fam).sum()} HGNC entries -> "
                 f"{(genes.family == fam).sum()} expression set -> "
                 f"{(three.family == fam).sum()} three-species set")
    for stage, n in e.groupby("stage").size().items():
        lines.append(f"    excluded at {stage}: {n}")
lines += ["", "Orthology exclusions by reason:"]
reasons = excl[excl.stage.str.contains("orthology")].reason \
    .str.replace(r" \(\d+ \w+ genes\)", "", regex=True)
lines += [f"    {n:4d}  {r}" for r, n in reasons.value_counts().items()]
lines += ["", "Mouse orthology overrides:"]
lines += [f"    {r.symbol}: {r.mouse_orthology_note}"
          for r in genes[genes.mouse_orthology_note != ""].itertuples()]
overridden = three[three.macaque_transcript_rule.str.startswith("closest")]
lines += ["", f"Macaque transcripts overriding canonical: {len(overridden)}"]
lines += [f"    {r.symbol}: {r.macaque_transcript_rule}, now {r.macaque_length_ratio}x"
          for r in overridden.itertuples()]
lines += ["", "NCBI cross-check (Ensembl partner vs NCBI gene_orthologs):"]
for sp, subset in [("mouse", genes), ("macaque", three)]:
    counts = subset[f"{sp}_ncbi_check"].str.replace(r" \(.*", "", regex=True).value_counts()
    lines.append(f"    {sp}: " + ", ".join(f"{k} {v}" for k, v in counts.items()))
    lines += [f"        {r.symbol}: {r[f'{sp}_symbol']} {r[f'{sp}_ncbi_check']}"
              for _, r in subset[subset[f"{sp}_ncbi_check"].str.startswith("disagrees")].iterrows()]
lines += ["",
          f"Expression set, genes with any review flag: {(genes.review_flags != '').sum()}",
          f"Three-species set, genes with any review flag: {(three.review_flags != '').sum()}",
          f"Three-species set, CDS outside {LENGTH_RANGE[0]}-{LENGTH_RANGE[1]}x human after override: "
          f"{(~three.mouse_length_ratio.between(*LENGTH_RANGE) | ~three.macaque_length_ratio.between(*LENGTH_RANGE)).sum()}"]
(OUT / "step1_summary.txt").write_text("\n".join(lines) + "\n")
print("\n" + "\n".join(lines))

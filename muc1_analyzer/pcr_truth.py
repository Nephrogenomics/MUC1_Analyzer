"""Score the LR-PCR lot against its clinical truth — on FOUR axes, not one.

Why this set changes what can be claimed. Everything measured so far rested on 12 labelled carriers of a
single variant family, or on panels of non-carriers where specificity is the only computable number. The
180-PCR lot carries, per sample: BOTH allele lengths (27-114 repeats), the variant, its position, which
allele bears it, and the clinical conclusion. That is five times the carriers of anything before it and,
for the first time, real data on the axes the score is actually built from.

⚠ WHICH WORKBOOK, and counted HOW — both had drifted. The runs use
`Clinical-DATA/Synthese_MUC1-corrected-180PCR.xlsx`, measured 2026-08-04:
**58 carrier ROWS / 119 negatives / 5 unknown, over 55 DISTINCT carrier subjects** (three carriers appear
in both substrates, and 6 subjects in all are paired). The "57 positives" this docstring used to quote
came from the earlier `Synthese_MUC1-lot-180PCR.xlsx` and was read as if it described the corrected file —
which is where the unexplained 58-vs-57 came from. Rows and subjects are different denominators; say
which one a figure uses.

The four axes, and what each one would otherwise never be tested on:
  1. CARRIER — 58 rows / 55 subjects against 119 negatives, where before it was 12/8.
  2. LENGTH, both alleles — the first set that populates the top of the distribution (A2 up to 114).
  3. VARIANT IDENTITY — 12 families, against dupC alone.
  4. POSITION + which allele carries it — 60 rows with a position; this is the ONSET axis's only real test.

⚠ TWO SUBSTRATES are mixed in the lot: raw `*.fastq.gz` and `*_minimap2_ont_T2T-clean.fq`, the latter
already aligned to T2T and cleaned. Score them SEPARATELY or a preprocessing effect will be read as a
method effect — the mistake logged on 2026-07-29 in another form.

⚠ What this lot already says about our own rules, before a single run: `60dupA` appears as a POSITIVE
variant at repeat 39, while the same motif is a recurrent artifact at repeat 10 in unrelated negatives.
`ins18` alone is NEGATIVE twice (in-frame, 18 bp), and `58_59delCC` comes out once positive and once
negative. The motif NAME cannot be the discriminator.

    python3 -m muc1_analyzer.pcr_truth --xlsx Clinical-DATA/Synthese_MUC1-corrected-180PCR.xlsx --out truth.tsv
    python3 -m muc1_analyzer.pcr_truth --xlsx ... --outdir pcr_lot180 --report

⚠ The truth table is PATIENT DATA: keep it on /scratch, never commit it.
"""
from __future__ import annotations
import argparse
import json
import os
import re
import sys

# Families we can recognise in a caller's motif label, mapped to the workbook's spelling. `59insG` and
# `58_59insG` are the same event written two ways in the lot; `ins18` is IN-FRAME and is NOT pathogenic
# there (2/2 negative), which is why it must map to a family of its own rather than to "a frameshift".
VARIANT_FAMILIES = {
    "59dupc": "59dupC", "58_59delcc": "58_59delCC", "delcc": "58_59delCC",
    "52dupg": "52dupG", "dupg": "52dupG",
    "58_59insg": "58_59insG", "59insg": "58_59insG", "insg": "58_59insG",
    "8_27del": "del8_27", "del8_27": "del8_27",
    "11dupt": "11dupT", "60dupa": "60dupA", "2_3delinsgta": "2_3delinsGTA",
    "52_33del": "52_33del", "ins18": "ins18",
}


# `ins18` is IN-FRAME (18 bp) and NOT pathogenic in this lot — 2/2 negative when it is the only finding.
# It matters here because a subject can carry BOTH: the workbook writes `Variant='ins18 | 59dupC'` with
# `VNTR_muté_repeats='80 | 82'`, one entry per ALLELE. Picking the frameshift is what makes the row
# scoreable; picking `ins18` would file an ADTKD carrier under a benign in-frame insertion.
IN_FRAME = {"ins18"}


def split_compound(cell) -> list:
    """A workbook cell holding one value per ALLELE (`'80 | 82'`, `'ins18 | 59dupC'`) → its parts. Pure.

    Found 2026-08-04: two carriers reported "no mut_len" for three days because `_num('80 | 82')` is None.
    The data was never missing — the cell holds a pair, one entry per allele, and the parser read only
    scalars."""
    if cell is None:
        return []
    s = str(cell).strip()
    if not s or s.lower() == "nan":
        return []
    return [p.strip() for p in s.split("|") if p.strip()]


def variant_family(label) -> "str | None":
    """A caller's motif label or a workbook cell → a canonical family. Pure."""
    s = re.sub(r"[^0-9a-zA-Z_]", "", str(label or "")).lower()
    if not s:
        return None
    for key in sorted(VARIANT_FAMILIES, key=len, reverse=True):
        if key in s:
            return VARIANT_FAMILIES[key]
    return None


def is_data_row(fichier) -> bool:
    """Workbook rows include SECTION HEADERS ('POSITIFS — 59dupC (45 fichiers)') and blanks. A data row
    names a file. Pure."""
    s = str(fichier or "").strip()
    return bool(s) and s.lower() != "nan" and bool(re.search(r"\.(fastq|fq)(\.gz)?$", s, re.I))


def subject_key(sample) -> str:
    """The SUBJECT behind a file name. Pure.

    The lot ships the same subject twice — as raw fastq and as a `_minimap2`-aligned-and-cleaned file —
    so a per-file table has two rows per paired subject. Collapsing on this key is what defines `paired`,
    and it is the only way to compare the two substrates WITHIN a subject instead of between arms."""
    return re.sub(r"_minimap2.*$", "", str(sample or ""))


def substrate(fichier: str) -> str:
    """Raw reads vs an already-aligned-and-cleaned file. Pure — the lot mixes both and they must not be
    pooled: one measures the pipeline, the other measures the pipeline after someone else's preprocessing."""
    return "t2t_clean" if re.search(r"T2T-clean", str(fichier), re.I) else "raw_fastq"


def _compound(variant_cell, mutlen_cell, _num) -> dict:
    """Pair each variant with ITS allele length, and report the FRAMESHIFT one. Pure.

    Both cells are per-allele lists in the same order, so entry i of one belongs with entry i of the
    other. When a subject carries an in-frame `ins18` alongside a frameshift, the frameshift is what the
    row must be scored on — reporting `ins18` would file an ADTKD carrier under a benign insertion, and
    reporting no length at all (what happened until 2026-08-04) drops the row out of the length axis."""
    fams = [variant_family(v) for v in split_compound(variant_cell)]
    lens = [_num(x) for x in split_compound(mutlen_cell)]
    pairs = [(f, lens[i] if i < len(lens) else None) for i, f in enumerate(fams) if f]
    if not pairs:
        return {"variant": variant_family(variant_cell), "mut_len": _num(mutlen_cell),
                "variant_other": None}
    frameshift = [p for p in pairs if p[0] not in IN_FRAME] or pairs
    fam, mlen = frameshift[0]
    others = [f for f, _ in pairs if f != fam]
    return {"variant": fam,
            "mut_len": mlen if mlen is not None else _num(mutlen_cell),
            "variant_other": ("|".join(others) or None)}


def parse_truth(xlsx: str) -> list:
    """Workbook → one dict per DATA row. Impure (pandas)."""
    import pandas as pd
    d = pd.read_excel(xlsx)
    out = []
    for _, r in d.iterrows():
        fich = r.get("Fichier")
        if not is_data_row(fich):
            continue
        concl = str(r.get("Conclusion") or "").strip().lower()
        status = "pos" if concl.startswith("posit") else ("neg" if concl.startswith("nég") or
                                                          concl.startswith("neg") else "unknown")
        def _num(x):
            try:
                return int(float(x))
            except (TypeError, ValueError):
                return None
        out.append({
            "file": str(fich).strip(),
            "sample": re.sub(r"\.(fastq|fq)(\.gz)?$", "", str(fich).strip(), flags=re.I),
            "substrate": substrate(fich),
            "status": status,
            "a1": _num(r.get("A1_repeats")), "a2": _num(r.get("A2_repeats")),
            **_compound(r.get("Variant"), r.get("VNTR_muté_repeats"), _num),
            "variant_raw": (None if str(r.get("Variant") or "").lower() in ("nan", "") else
                            str(r.get("Variant")).strip()),
            "variant_pos": _num(r.get("Variant_pos")),
            # The sequencing run. Cross-contamination cannot cross a run, so this column is the CONTROL
            # of the within/between-run test in `contamination.run_enrichment` — without it that test
            # has no null and the whole detector rests on assuming some panel is clean.
            "run": re.sub(r"^\w+\s*/\s*", "", str(r.get("Source (feuille / run)") or "")).replace(" ", "")
                   or None,
        })
    # PAIRED subjects — the same person present as BOTH raw reads and a T2T-cleaned file. Six of them in
    # this lot, and they are the only way to measure what the third-party cleaning COSTS instead of
    # assuming it: same subject, same pipeline, two substrates. If the cleaning kept only reads that
    # mapped, it removed long-allele reads preferentially, and the 150 independent cleaned files would
    # carry that bias into every length measurement.
    seen = {}
    for r in out:
        seen.setdefault(subject_key(r["sample"]), set()).add(r["substrate"])
    for r in out:
        r["paired"] = len(seen[subject_key(r["sample"])]) > 1
    return out


# ── the four axes ─────────────────────────────────────────────────────────────────────────────────
# The truth vocabulary, in ONE place. It is written here (`pos`/`neg`/`unknown`) and read by every tool
# that scores against this table; a consumer that invents its own spelling silently matches nothing and
# reports an empty cohort as if it were a clean result (measured 2026-08-03: `miss_report` filtered on
# "carrier"/"positive" and printed "no labelled carrier matched" on a table holding 44 of them).
STATUS_POS, STATUS_NEG, STATUS_UNKNOWN = "pos", "neg", "unknown"


def is_carrier(status) -> bool:
    """True for a labelled CARRIER row, tolerant of the spellings the sheets have used. Pure."""
    return str(status or "").strip().lower() in (STATUS_POS, "positive", "positif", "carrier")


def is_negative(status) -> bool:
    """True for a labelled NON-CARRIER row. Pure."""
    return str(status or "").strip().lower() in (STATUS_NEG, "negative", "negatif", "négatif")


def check_row(r) -> list:
    """Everything IMPOSSIBLE about one truth row. Pure. Empty when the row is self-consistent.

    A truth table is not a measurement, it is the ruler — so a row that cannot be true must be visible
    before it silently sets a denominator. Three kinds were found on 2026-08-04, all of which had been
    passing through the pipeline unremarked:

    · `mut_len` matching NEITHER allele. Two carriers read 78 with alleles 43|44 and 95 with 42|43. The
      mutant allele is one of the two alleles by definition, so one of the three numbers is wrong — and
      the signature (a long mutant against two short alleles) is what a clinical method that under-calls
      long alleles produces, which points at `a1`/`a2`. Not decidable here: it needs the source report.
    · a CARRIER whose variant family did not parse. `variant_family` returns None when nothing matches,
      and a boolean `False` cell survives as the string "False" — a carrier with no usable family is
      excluded from every per-family denominator without anyone noticing.
    · a carrier with no `mut_len` at all: the length axis cannot score it, and the length-bucket analysis
      silently drops it into `unknown`."""
    out = []
    if not is_carrier(r.get("status")):
        return out
    v = str(r.get("variant") or "").strip()
    if not v or v.lower() in ("false", "true", "nan", "none"):
        # `variant_raw` is what the WORKBOOK holds; without it the message says a family is missing but
        # not what to fix. A boolean cell arrives as "False" and a typo arrives as itself.
        out.append(f"carrier with no usable variant family — workbook cell "
                   f"{r.get('variant_raw')!r} did not match any family in VARIANT_FAMILIES")
    m, a1, a2 = r.get("mut_len"), r.get("a1"), r.get("a2")
    if m is None:
        out.append("carrier with no mut_len — cannot be scored on the length axis")
    elif a1 is not None and a2 is not None:
        if min(abs(int(m) - int(a1)), abs(int(m) - int(a2))) > 3:
            out.append(f"mut_len {m} matches NEITHER allele ({a1}|{a2}) — one of the three is wrong")
    return out


def check_truth(rows) -> list:
    """[(sample, [problems])] over a whole truth table. Pure."""
    return [(r.get("sample"), p) for r in rows for p in [check_row(r)] if p]


def score_carrier(pred, truth_status) -> "str | None":
    """TP/FP/FN/TN, or None when the truth is unknown (never folded into the negatives). Pure."""
    if truth_status not in ("pos", "neg"):
        return None
    if truth_status == "pos":
        return "TP" if pred else "FN"
    return "FP" if pred else "TN"


def score_lengths(pred_alleles, a1, a2, *, tol: int = 1) -> dict:
    """Both alleles against the clinical pair, unordered, with a tolerance. Pure.

    `tol=1` by default because the two methods have differed by exactly 1 on samples both call correct
    (45 vs 44 on one subject, and the VNTRPipeline concordance table already treated delta<=1 as agreement). The
    LONG allele is reported separately: it is the one that goes missing, and the one assessability
    collapses on."""
    truth = sorted(x for x in (a1, a2) if x is not None)
    pred = sorted(int(x) for x in (pred_alleles or []) if x is not None)
    out = {"truth": truth, "pred": pred, "n_truth": len(truth), "n_pred": len(pred),
           "exact": None, "within_tol": None, "long_within_tol": None, "delta_long": None}
    if not truth or not pred:
        return out
    out["exact"] = (pred == truth)
    if len(pred) == len(truth):
        out["within_tol"] = all(abs(p - t) <= tol for p, t in zip(pred, truth))
    out["delta_long"] = pred[-1] - truth[-1]
    out["long_within_tol"] = abs(out["delta_long"]) <= tol
    return out


def score_variant(pred_label, truth_variant) -> "str | None":
    """Identity of the called variant vs the clinical one. Pure.
    None when the sample is a negative or the truth carries no variant."""
    if not truth_variant:
        return None
    p = variant_family(pred_label)
    if p is None:
        return "missed"
    return "match" if p == truth_variant else "wrong_family"


def score_position(pred_repeat, truth_pos, *, tol: int = 2) -> dict:
    """Repeat index against the clinical one. Pure.

    ⚠ Both must be in the PHYSICAL frame. The C-run rank published as a repeat number cost 0.11 on the
    onset axis (2026-07-29); with 60 positions in this lot, that class of error becomes measurable rather
    than anecdotal."""
    if truth_pos is None or pred_repeat is None:
        return {"delta": None, "within_tol": None}
    delta = int(pred_repeat) - int(truth_pos)
    return {"delta": delta, "within_tol": abs(delta) <= tol}


def summarize(rows: list) -> dict:
    """Per-substrate confusion on the four axes. Pure."""
    out = {}
    for sub in sorted({r.get("substrate", "raw_fastq") for r in rows}) + ["ALL"]:
        grp = rows if sub == "ALL" else [r for r in rows if r.get("substrate") == sub]
        c = {"TP": 0, "FP": 0, "FN": 0, "TN": 0}
        for r in grp:
            v = r.get("carrier_verdict")
            if v:
                c[v] += 1
        lens = [r for r in grp if r.get("length") and r["length"].get("within_tol") is not None]
        longs = [r for r in grp if r.get("length") and r["length"].get("long_within_tol") is not None]
        vars_ = [r.get("variant_verdict") for r in grp if r.get("variant_verdict")]
        poss = [r for r in grp if (r.get("position") or {}).get("within_tol") is not None]
        out[sub] = {
            "n": len(grp), "carrier": c,
            "sensitivity": c["TP"] / max(1, c["TP"] + c["FN"]),
            "specificity": c["TN"] / max(1, c["TN"] + c["FP"]),
            "length_pairs_scored": len(lens),
            "length_within_tol": sum(1 for r in lens if r["length"]["within_tol"]),
            "length_exact": sum(1 for r in lens if r["length"]["exact"]),
            "long_allele_scored": len(longs),
            "long_allele_within_tol": sum(1 for r in longs if r["length"]["long_within_tol"]),
            "variant_match": vars_.count("match"), "variant_wrong": vars_.count("wrong_family"),
            "variant_missed": vars_.count("missed"),
            "position_scored": len(poss),
            "position_within_tol": sum(1 for r in poss if r["position"]["within_tol"]),
        }
    return out


def render(summ: dict) -> str:
    L = []
    for sub, s in summ.items():
        c = s["carrier"]
        L.append(f"── {sub}  (n={s['n']})")
        L.append(f"   CARRIER   TP={c['TP']} FN={c['FN']} sens={s['sensitivity']:.2f} | "
                 f"FP={c['FP']} TN={c['TN']} spec={s['specificity']:.2f}")
        if s["length_pairs_scored"]:
            L.append(f"   LENGTH    {s['length_within_tol']}/{s['length_pairs_scored']} both alleles "
                     f"within tol ({s['length_exact']} exact) · LONG allele "
                     f"{s['long_allele_within_tol']}/{s['long_allele_scored']}")
        if s["variant_match"] + s["variant_wrong"] + s["variant_missed"]:
            L.append(f"   VARIANT   match={s['variant_match']} wrong_family={s['variant_wrong']} "
                     f"missed={s['variant_missed']}")
        if s["position_scored"]:
            L.append(f"   POSITION  {s['position_within_tol']}/{s['position_scored']} within tol")
        L.append("")
    L.append("⚠ Substrates are scored apart on purpose: `t2t_clean` files went through someone else's")
    L.append("  alignment and cleaning, so pooling them would read a preprocessing effect as a method one.")
    return "\n".join(L) + "\n"


def collect(outdir: str, sample: str) -> dict:
    """What `run` left for one sample → the predictions the four axes need. Pure-ish (fs)."""
    base = os.path.join(outdir, sample)
    pred = {"carrier": None, "alleles": None, "variant_label": None, "repeat": None}
    try:
        s = json.load(open(base + ".score.json"))
        pred["carrier"] = s.get("carrier")
        pred["alleles"] = [x for x in (s.get("mut_len"), s.get("healthy_len")) if x is not None]
        pred["variant_label"] = s.get("variant")
        oi, ml = s.get("onset_index"), s.get("mut_len")
        if oi is not None and ml:
            pred["repeat"] = round((1 - float(oi)) * int(ml))
    except Exception:
        pass
    try:
        d = json.load(open(base + ".dupc.json"))
        v = ((d.get("result") or {}).get("variant") or {})
        pred["repeat"] = v.get("repeat") if v.get("repeat") is not None else pred["repeat"]
        pred["variant_label"] = pred["variant_label"] or v.get("label")
    except Exception:
        pass
    return pred


def main(argv=None):
    ap = argparse.ArgumentParser(prog="muc1_analyzer.pcr_truth", description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--xlsx", required=True)
    ap.add_argument("--out", help="write the parsed truth as TSV (PATIENT DATA — keep on /scratch)")
    ap.add_argument("--outdir", help="directory of `run` outputs → score the four axes")
    ap.add_argument("--tol", type=int, default=1, help="length tolerance in repeats (default 1)")
    ap.add_argument("--pos-tol", dest="pos_tol", type=int, default=2)
    ap.add_argument("--json", dest="out_json", default=None)
    a = ap.parse_args(argv)

    truth = parse_truth(a.xlsx)
    print(f"[pcr_truth] {len(truth)} data rows · "
          f"pos={sum(1 for r in truth if r['status'] == 'pos')} "
          f"neg={sum(1 for r in truth if r['status'] == 'neg')} "
          f"unknown={sum(1 for r in truth if r['status'] == 'unknown')}", file=sys.stderr)
    for sub in sorted({r["substrate"] for r in truth}):
        print(f"[pcr_truth]   substrate {sub}: {sum(1 for r in truth if r['substrate'] == sub)}",
              file=sys.stderr)
    problems = check_truth(truth)
    if problems:
        print(f"[pcr_truth] ⚠ {len(problems)} INCONSISTENT truth row(s) — fix at the SOURCE workbook, "
              f"not here; every one of these sets a denominator:", file=sys.stderr)
        for sample, probs in problems:
            for msg in probs:
                print(f"[pcr_truth]     {sample}: {msg}", file=sys.stderr)
    print(f"[pcr_truth]   paired (both substrates): "
          f"{sum(1 for r in truth if r.get('paired'))} rows = "
          f"{sum(1 for r in truth if r.get('paired')) // 2} subjects", file=sys.stderr)

    if a.out:
        cols = ["file", "sample", "substrate", "status", "a1", "a2", "mut_len", "variant",
                "variant_other", "variant_pos", "paired", "run"]
        with open(a.out, "w") as fh:
            fh.write("\t".join(cols) + "\n")
            for r in truth:
                fh.write("\t".join("" if r.get(c) is None else str(r.get(c)) for c in cols) + "\n")
        print(f"[pcr_truth] truth → {a.out}  ⚠ patient data, keep on /scratch", file=sys.stderr)

    if not a.outdir:
        return 0
    scored = []
    for r in truth:
        pred = collect(a.outdir, r["sample"])
        if pred["carrier"] is None:
            continue                                    # not run (yet) — never scored as a miss
        row = dict(r)
        row["carrier_verdict"] = score_carrier(bool(pred["carrier"]), r["status"])
        row["length"] = score_lengths(pred["alleles"], r["a1"], r["a2"], tol=a.tol)
        row["variant_verdict"] = score_variant(pred["variant_label"], r["variant"])
        row["position"] = score_position(pred["repeat"], r["variant_pos"], tol=a.pos_tol)
        scored.append(row)
    print(f"[pcr_truth] {len(scored)}/{len(truth)} samples have a `run` output", file=sys.stderr)
    summ = summarize(scored)
    sys.stdout.write(render(summ))
    if a.out_json:
        with open(a.out_json, "w") as fh:
            json.dump({"summary": summ, "rows": scored}, fh, indent=1, default=str)
        print(f"[pcr_truth] → {a.out_json}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

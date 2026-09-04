#!/usr/bin/env python3
"""miss_report — why a labelled carrier was NOT called, gate by gate.

A cohort sweep prints one sensitivity figure and hides four different failures behind it. The 152-sample
run scored 20/44 all-families and 0.62 dupC among the adequately powered; "12 dupC missed AT ADEQUATE
POWER" is the dominant term and a single number cannot say whether those 12 are a threshold problem or a
substrate problem. Those have opposite fixes: a gate can be moved, a signal that is not in the reads
cannot.

So this reads the `.dupc.json` a `run` already wrote and splits each miss into ONE exit:

  SIGNAL_ABSENT        no read anywhere on either allele renders the mutant C-tract length (k_max = 0).
                       Nothing statistical is happening — the variant is not in the basecalls. Moving a
                       gate cannot recover this; it is a substrate/basecalling question.
  AT_BACKGROUND        k > 0, but the sample's OWN null already explains it (uncorrected p > 0.05). Also
                       unrecoverable by any threshold. This exit exists because `k > 0` was first treated
                       as evidence: the one paired carrier read k 11/541 at f_null 0.02, where the null
                       predicts 10.8 — reported as GATE_REJECTED, i.e. as if a gate change would find it.
  NO_POSITION_MEASURED no C-run reached `min_reads` — coverage, not detection.
  NULL_NOT_ESTIMABLE   fewer than `min_null_obs` observations to calibrate against (leave-one-out null).
  p_bonf / min_k / min_frac   the signal IS there and a named gate rejected it (possibly several).

`k_max` is taken over EVERY scored position of EVERY allele, not just the best-p candidate: the candidate
is the most significant position, which need not be the variant's. `candidate.position` is a C-run RANK
and does not map to a repeat index without the reference, so this tool deliberately does not try — a
`k_max` of 0 is conclusive on its own and needs no position mapping.

⚠ JSON only, no pysam: this must stay runnable on a login node.

    python3 -m muc1_analyzer.miss_report --outdir pcr_lot180_t2t_clean --truth truth.tsv \
        [--family 59dupC] [--out miss.tsv]
"""
from __future__ import annotations
import argparse
import json
import os
import sys

# Mirrors runlen_shift.DEFAULTS. Imported, not copied — a gate that drifts here would explain a miss
# with a threshold the caller never applied. Same reasoning for the status vocabulary: `pcr_truth` writes
# `pos`/`neg`, and a consumer that spells it "carrier" matches nothing and prints an empty cohort as
# though it were a result. Neither module imports pysam; this stays login-node safe.
from .detectors.runlen_shift import DEFAULTS
from .pcr_truth import is_carrier, subject_key

EXITS = ("CALLED", "SIGNAL_ABSENT", "AT_BACKGROUND", "NO_POSITION_MEASURED", "NULL_NOT_ESTIMABLE",
         "GATE_REJECTED", "NO_OUTPUT")

# The line between "not signal" and "signal a gate rejected". Deliberately the conventional 0.05 on the
# caller's own UNCORRECTED p: a position that cannot clear that before any multiple-testing correction is
# not separable from the sample's own null, and no threshold change makes it separable.
BACKGROUND_ALPHA = 0.05

# Fixed schema. Every row carries every column, and the TSV header is written from this — NOT from
# `rows[0].keys()`, which silently takes its columns from whichever sample happened to come first.
ROW_COLS = ("sample", "variant", "substrate", "run", "called", "interpretation", "exit", "gates_failed",
            "arbiter_lengths", "mut_len_clinical", "a1", "a2",
            "k_max", "n", "frac", "f_null", "p", "p_bonf", "n_null_obs", "ctract_rank", "contig",
            "mut_contig", "mut_k", "mut_n", "mut_frac", "mut_p", "mut_p_bonf", "mut_n_tested",
            "n_cover", "adequately_powered", "families_tested", "alarms")


def read_truth(path: str) -> list:
    with open(path) as fh:
        head = fh.readline().rstrip("\n").split("\t")
        return [dict(zip(head, ln.rstrip("\n").split("\t"))) for ln in fh if ln.strip()]


def alleles_of(doc: dict) -> dict:
    """{contig: runlen result} from a `.dupc.json`, whatever wrapper it arrived in. Pure."""
    res = (doc or {}).get("result") or doc or {}
    return res.get("per_allele") or {}


def best_evidence(per_allele: dict) -> dict:
    """The strongest MUTANT evidence anywhere on the sample: max k over every scored position. Pure.

    Not the most significant position — the most POSITIVE one. A p-value ranks positions against the
    sample's own null; `k` asks the prior question of whether any read carries the variant at all."""
    best, k_max, n_tested = None, 0, 0
    for contig, r in (per_allele or {}).items():
        if not isinstance(r, dict) or r.get("error"):
            continue
        n_tested += r.get("n_tested") or 0
        for s in r.get("scored") or []:
            if (s.get("k_ge_mut") or 0) >= k_max:
                k_max = s.get("k_ge_mut") or 0
                best = {**s, "contig": contig}
    return {"k_max": k_max, "row": best, "n_tested_total": n_tested}


def mutant_side(per_allele: dict, mut_len) -> dict:
    """The candidate on the scaffold nearest the CLINICAL mutant allele. Pure. {} when unresolvable.

    `k_max` answers "is the variant anywhere", and for that a raw count is right. It is the WRONG number
    for asking why a carrier was missed: on a deep short-allele scaffold the maximum k is simply the
    background — k=22 of 918 at f_null 0.0153 is exactly what noise predicts — so `k_max` lands on the
    best-covered contig rather than on the one carrying the variant. Reading a long-allele carrier off
    that row shows `p_bonf 0.105` next to `called=True`, which is not a contradiction but a wrong column.
    This reports what was seen ON THE ALLELE WHERE THE VARIANT IS EXPECTED."""
    import re
    try:
        L = int(float(mut_len))
    except (TypeError, ValueError):
        return {}
    best, dist = None, None
    for contig, r in (per_allele or {}).items():
        if not isinstance(r, dict) or r.get("error"):
            continue
        m = re.search(r"MUC1_VNTR_(\d+)repeats", str(contig))
        if not m:
            continue
        d = abs(int(m.group(1)) - L)
        if dist is None or d < dist:
            scored = [x for x in (r.get("scored") or []) if x.get("p_bonf") is not None]
            if not scored:
                continue
            cand = sorted(scored, key=lambda x: (x["p_bonf"], -(x.get("k_ge_mut") or 0)))[0]
            # The RAW p and the number of positions tested, so the Bonferroni burden can be re-priced
            # without a re-run: p_bonf = p x n_tested, and n_tested scales with the allele's length —
            # 81 positions on a 103-copy contig against ~40 on a 40-copy one. Whether shrinking the
            # search space could recover a miss is then arithmetic rather than a rebuild.
            best, dist = {"mut_contig": contig, "mut_k": cand.get("k_ge_mut"), "mut_n": cand.get("n"),
                          "mut_frac": cand.get("frac"), "mut_p": cand.get("p"),
                          "mut_p_bonf": cand.get("p_bonf"), "mut_n_tested": r.get("n_tested")}, d
    return best or {}


def classify_miss(per_allele: dict, gates: dict = None) -> dict:
    """Which single exit explains a non-call. Pure — this is the whole point of the module."""
    g = {**DEFAULTS, **(gates or {})}
    ev = best_evidence(per_allele)
    row = ev["row"] or {}
    if not ev["n_tested_total"]:
        return {"exit": "NO_POSITION_MEASURED", "gates_failed": [], **ev}
    if (row.get("n_null_obs") or 0) < g["min_null_obs"]:
        return {"exit": "NULL_NOT_ESTIMABLE", "gates_failed": [], **ev}
    if ev["k_max"] == 0:
        return {"exit": "SIGNAL_ABSENT", "gates_failed": [], **ev}
    # k > 0 is NOT signal. The paired carrier read k 11/541 at f_null 0.02 — the null predicts 10.8, so
    # its "evidence" was the background itself, and calling that GATE_REJECTED advertises it as
    # recoverable by a threshold. Use the caller's own UNCORRECTED p at that position: above 0.05 the
    # position is not separable from the sample's own null, and no gate change makes it so.
    if (row.get("p") if row.get("p") is not None else 1.0) > BACKGROUND_ALPHA:
        return {"exit": "AT_BACKGROUND", "gates_failed": [], **ev}
    failed = []
    # `min_cover` belongs in this list. It was missing, so a carrier refused by the power floor showed
    # only the OTHER gates it failed — and "would lowering the floor recover anyone?" could not be read
    # off the table at all. It is named FIRST because it is the gate that fires before the evidence is
    # even weighed.
    if (row.get("n") or 0) < g["min_cover"]:
        failed.append("min_cover")
    if (row.get("p_bonf") if row.get("p_bonf") is not None else 1.0) > g["alpha"]:
        failed.append("p_bonf")
    if ev["k_max"] < g["min_k"]:
        failed.append("min_k")
    if (row.get("frac") or 0.0) < g["min_frac"]:
        failed.append("min_frac")
    return {"exit": "GATE_REJECTED" if failed else "CALLED", "gates_failed": failed, **ev}


def arbiter_lengths(outdir: str, sample: str) -> "str | None":
    """The allele lengths the run SCAFFOLDED on, from `<sample>.score.json`. Impure (reads a file).

    `dispatch_dupc_vntr` gets `lengths=` from the arbiter when it measured them, and `allele_contigs`
    picks contigs on its own when it did not. That second regime is what produced the 5C debacle: reads
    from both alleles pile onto one contig, the mutant fraction halves, and a copy-number mismatch plants
    an indel artifact. A carrier whose mutant fraction reads 0.05-0.23 where a heterozygote owes ~0.50 is
    exactly what that looks like, so whether the arbiter spoke belongs next to every miss."""
    path = os.path.join(outdir, f"{sample}.score.json")
    if not os.path.exists(path):
        return None
    try:
        with open(path) as fh:
            d = json.load(fh)
    except Exception:
        return None
    a = arbiter_alleles_raw(outdir, sample)
    if a is not None and len(a) == 2:
        return f"{a[0]},{a[1]}"
    return "none"


def arbiter_alleles_raw(outdir: str, sample: str) -> "list | None":
    """EVERY allele length the arbiter measured, however many, or None when it measured nothing. Impure.

    `arbiter_lengths` above answers "did the arbiter hand the caller two lengths", so it collapses ONE
    measured allele and NO score.json into the same "none" — which is the right answer for scaffolding and
    the wrong one for evaluability, where "the arbiter saw two alleles and the caller reconstructed one" is
    exactly the case to isolate. Same file, same key, so this stays the single place to be wrong about
    where the arbiter's lengths live."""
    path = os.path.join(outdir, f"{sample}.score.json")
    if not os.path.exists(path):
        return None
    try:
        with open(path) as fh:
            d = json.load(fh)
    except Exception:
        return None
    a = (d.get("arbiter_vs_caller_length") or {}).get("arbiter")
    if not isinstance(a, (list, tuple)) or not a:
        return None
    if not all(isinstance(x, (int, float)) and not isinstance(x, bool) for x in a):
        return None
    return [int(x) for x in a]


def power_profile(doc: dict) -> dict:
    """How many alleles were TESTED and how many of those cleared the power floor. Pure.

    Distinct from `power_of` below, which answers "was there any power at all" with `any`/`max` — the
    best case. Evaluability needs the worst case AND the alleles that were never scaffolded: the
    2026-08-07 note in `dupc_dispatch` records a subject reported `fully_powered: True` with half its
    genotype missing, because `underpowered` can only speak about alleles that were tested. So this
    returns the two counts separately and lets the caller decide."""
    res = (doc or {}).get("result") or doc or {}
    pw = res.get("power_by_allele") or {}
    per = res.get("per_allele") or {}
    tested = [c for c, r in per.items() if (r or {}).get("n_tested")]
    if not tested:                       # fall back to the power map when per_allele is absent
        tested = list(pw)
    return {"n_tested": len(tested),
            "n_powered": sum(1 for c in tested if (pw.get(c) or {}).get("adequately_powered")),
            "has_record": bool(pw or per)}


def power_of(doc: dict) -> dict:
    """Worst-case power over the sample's alleles — a carrier is missed on ONE allele, not on average."""
    res = (doc or {}).get("result") or doc or {}
    pw = res.get("power_by_allele") or {}
    if not pw:
        return {"n_cover": None, "adequately_powered": None}
    cov = [(v or {}).get("n_cover") or 0 for v in pw.values()]
    return {"n_cover": max(cov) if cov else 0,
            "adequately_powered": any((v or {}).get("adequately_powered") for v in pw.values())}


def diagnose(outdir: str, truth: list, family: str = None) -> list:
    """One row per labelled CARRIER that has a `.dupc.json`. Impure (reads files)."""
    rows = []
    for t in truth:
        if not is_carrier(t.get("status")):
            continue
        if family and (t.get("variant") or "") != family:
            continue
        base = {c: None for c in ROW_COLS}
        # The lot mixes RAW fastq with fastq made from an already-aligned-and-cleaned T2T bam, and
        # `pcr_lot180.sbatch` runs ONE substrate per directory. Carrying the label makes a single-substrate
        # arm say so itself, instead of letting "13 misses" read as a property of the caller.
        base.update(sample=t["sample"], variant=t.get("variant"), substrate=t.get("substrate"),
                    run=t.get("run"), mut_len_clinical=t.get("mut_len"), a1=t.get("a1"), a2=t.get("a2"))
        path = os.path.join(outdir, f"{t['sample']}.dupc.json")
        if not os.path.exists(path):
            # A sample the sweep never produced is NOT a caller miss — it never reached the caller. It
            # keeps the full schema so it cannot break a consumer, and `called` is False, never absent.
            rows.append({**base, "exit": "NO_OUTPUT", "called": False, "gates_failed": ""})
            continue
        with open(path) as fh:
            doc = json.load(fh)
        res = (doc or {}).get("result") or doc or {}
        cls = classify_miss(alleles_of(doc))
        if res.get("called"):
            cls["exit"], cls["gates_failed"] = "CALLED", []
        pw = power_of(doc)
        row = {**cls["row"]} if cls.get("row") else {}
        base["arbiter_lengths"] = arbiter_lengths(outdir, t["sample"])
        base.update(mutant_side(alleles_of(doc), t.get("mut_len")))
        rows.append({**base,
                     "called": bool(res.get("called")), "interpretation": res.get("interpretation"),
                     "exit": cls["exit"], "gates_failed": ",".join(cls["gates_failed"]),
                     "k_max": cls["k_max"], "n": row.get("n"), "frac": row.get("frac"),
                     "f_null": row.get("f_null"), "p": row.get("p"), "p_bonf": row.get("p_bonf"),
                     "n_null_obs": row.get("n_null_obs"), "ctract_rank": row.get("position"),
                     "contig": row.get("contig"), "n_cover": pw["n_cover"],
                     "adequately_powered": pw["adequately_powered"],
                     "families_tested": ",".join(str(x) for x in (res.get("ctract_families_tested") or [])),
                     "alarms": ",".join(res.get("frameshift_alarms") or [])})
    return rows


def pair_up(rows: list) -> list:
    """[(subject, raw_row, clean_row)] for subjects present on BOTH substrates. Pure.

    Between-arm rates compare different patients, so a difference can be the patients. Within a subject
    everything but the substrate is held constant, which is the only design that can attribute a miss to
    the input material."""
    by = {}
    for r in rows:
        by.setdefault(subject_key(r["sample"]), {})[r["substrate"]] = r
    return [(s, d["raw_fastq"], d["t2t_clean"]) for s, d in sorted(by.items())
            if "raw_fastq" in d and "t2t_clean" in d]


def paired_block(rows: list) -> list:
    """The 2x2 a paired design exists to produce: who is called on which substrate. Pure."""
    pairs = pair_up(rows)
    if not pairs:
        return ["no subject has BOTH substrates here — the paired comparison is not available", ""]
    cells = {(True, True): [], (True, False): [], (False, True): [], (False, False): []}
    for s, raw, clean in pairs:
        cells[(bool(raw["called"]), bool(clean["called"]))].append((s, raw, clean))
    lost = cells[(True, False)]
    gained = cells[(False, True)]
    out = [f"PAIRED (same subject, both substrates) — {len(pairs)} subjects:"]
    if len(pairs) < 5:
        # `substrate()` keys on "T2T-clean" in the file name while `subject_key()` collapses "_minimap2".
        # When the two markers do not co-occur, the substrate is still detected but nothing pairs, and the
        # arm quietly shrinks to a handful of subjects instead of reporting that it could not pair.
        out += ["  ⚠ TOO FEW TO CONCLUDE. Before reading anything into this, check that the pairing key is",
                "    finding the subjects: `substrate()` matches 'T2T-clean' in the file name, `subject_key()`",
                "    collapses '_minimap2...'. If those markers do not co-occur, pairing fails silently."]
    out += [
           f"  called on both                     : {len(cells[(True, True)])}",
           f"  raw_fastq ONLY (lost by cleaning)  : {len(lost)}   ← the substrate hypothesis",
           f"  t2t_clean ONLY (gained by cleaning): {len(gained)}",
           f"  missed on both                     : {len(cells[(False, False)])}"]
    # Counts for EVERY paired subject, not just the discordant ones. A subject missed on both substrates
    # still says whether cleaning removed reads — the verdict is thresholded, the counts are not, and with
    # a handful of subjects the counts are the only thing carrying information.
    out.append("  reads at the best position, raw → clean (verdict in brackets):")
    for s, raw, clean in pairs[:20]:
        v = f"{'call' if raw['called'] else 'miss'}→{'call' if clean['called'] else 'miss'}"
        out.append(f"    [{v}]  k {raw['k_max']}/{raw['n']} (frac {raw['frac']})"
                   f"  →  k {clean['k_max']}/{clean['n']} (frac {clean['frac']})")
    if len(pairs) > 20:
        out.append(f"    … {len(pairs) - 20} more, see the TSV")
    if len(lost) > len(gained):
        out += ["", "  ⚠ Discordance is one-directional. That is consistent with the cleaning having kept",
                "    only reads that mapped — which removes long-allele reads preferentially, and the",
                "    mutant sits on the long allele. NOT proof: confirm on the read counts above, and",
                "    remember a between-arm rate difference alone would not have shown this."]
    return out + [""]


def scaffold_block(rows: list) -> list:
    """Called vs missed against the mutant FRACTION and the scaffold's provenance. Pure.

    A heterozygote owes ~0.50 of the reads on its mutant allele's scaffold. Measured on the cohort, the
    CALLED carriers sit at 0.29-0.77 and the MISSED at 0.05-0.23 — the mutant is diluted, not absent, and
    that dilution is what saturates p_bonf (k=3 of 24 gives a raw p near 0.02, which Bonferroni over ~45
    C-runs takes to 1.0). So the ceiling is not the threshold, it is the numerator.

    Whether the arbiter measured the lengths is the first thing to rule in or out: without them
    `allele_contigs` picks contigs on its own, both alleles land on one scaffold, and the mutant fraction
    halves by construction. That regime is exactly what the 5C route turned out to be made of."""
    def frac_of(r):
        return r["frac"] if isinstance(r.get("frac"), (int, float)) else None

    called = [r for r in rows if r["called"]]
    missed = [r for r in rows if not r["called"]]
    out = ["MUTANT FRACTION — a heterozygote owes ~0.50 on its own allele:"]
    for label, grp in (("called", called), ("missed", missed)):
        fr = sorted(x for x in (frac_of(r) for r in grp) if x is not None)
        if fr:
            mid = fr[len(fr) // 2]
            out.append(f"  {label:<7} n={len(fr):<3} median {mid:.2f}   range {fr[0]:.2f}-{fr[-1]:.2f}")
        else:
            out.append(f"  {label:<7} n=0")
    # LENGTH of the mutant allele. `dupc_dispatch`s own header records the mechanism — "the PCR-depleted
    # mutant dilutes to the ~5 % ONT homopolymer floor" — and LR-PCR depletes the LONGER allele, which is
    # usually the mutant one. If the miss rate rises with the mutant allele's length, the dilution is the
    # amplification and no threshold reaches it.
    buckets = {}
    for r in rows:
        try:
            L = int(float(r["mut_len_clinical"]))
        except (TypeError, ValueError):
            L = None
        key = "unknown" if L is None else ("<50" if L < 50 else "50-69" if L < 70 else ">=70")
        b = buckets.setdefault(key, [0, 0])
        b[0] += 1
        b[1] += int(bool(r["called"]))
    order = {"<50": 0, "50-69": 1, ">=70": 2, "unknown": 3}
    out.append("")
    out.append("MUTANT ALLELE LENGTH (called / carriers with an output):")
    for key in sorted(buckets, key=lambda k: order.get(k, 9)):
        tot, ok = buckets[key]
        out.append(f"  {key:<8} {ok}/{tot}" + (f"   {ok / tot:.2f}" if tot else ""))
    out.append("  → a rate that falls as the allele grows is LR-PCR depleting the mutant, not a gate.")

    by = {}
    for r in rows:
        key = "arbiter measured" if (r["arbiter_lengths"] or "none") != "none" else "NO arbiter lengths"
        s = by.setdefault(key, [0, 0])
        s[0] += 1
        s[1] += int(bool(r["called"]))
    out.append("")
    out.append("SCAFFOLD provenance (called / carriers with an output):")
    for key, (tot, ok) in sorted(by.items()):
        out.append(f"  {key:<20} {ok}/{tot}")
    if len(by) > 1:
        out += ["  → if the miss rate concentrates on `NO arbiter lengths`, the scaffold is the cause and",
                "    the fix is upstream of the caller: both alleles are landing on one contig, which",
                "    halves the mutant fraction before any threshold is applied."]
    return out + [""]


def summarize(rows: list) -> str:
    # A carrier the sweep never produced is not a sensitivity denominator — counting it as a miss would
    # blame the caller for a sample it never saw. Reported on its own line instead.
    absent = [r for r in rows if r["exit"] == "NO_OUTPUT"]
    seen = [r for r in rows if r["exit"] != "NO_OUTPUT"]
    miss = [r for r in seen if not r["called"]]
    powered = [r for r in miss if r["adequately_powered"]]
    out = [f"labelled carriers       : {len(rows)}",
           f"  no output from the run: {len(absent)}   ← not scored; the sweep never produced these",
           f"carriers with an output : {len(seen)}",
           f"  called                : {sum(1 for r in seen if r['called'])}",
           f"  missed                : {len(miss)}  (adequately powered: {len(powered)})", ""]
    if not seen:
        return "\n".join(out) + "\nNothing to decompose — no carrier in this table has a .dupc.json.\n"
    subs = {}
    for r in seen:
        s = subs.setdefault(r["substrate"] or "?", [0, 0])
        s[0] += 1
        s[1] += int(bool(r["called"]))
    if len(subs) == 1:
        only = next(iter(subs))
        out += [f"⚠ SINGLE SUBSTRATE ({only}): `pcr_lot180.sbatch` runs one substrate per directory, so",
                "  every miss below is a miss ON THIS SUBSTRATE. Nothing here separates the caller from the",
                "  input material — that needs the other arm, or better SUB=paired (same patients, both).", ""]
    else:
        out.append("by substrate (called / with output):")
        for s, (tot, ok) in sorted(subs.items()):
            out.append(f"  {s:<12} {ok}/{tot}")
        out.append("")
        out += paired_block(seen)
    out += scaffold_block(seen)
    out.append("MISSES AT ADEQUATE POWER, by exit — this is the decomposition that matters:")
    by = {}
    for r in powered:
        key = r["exit"] + (f" [{r['gates_failed']}]" if r["gates_failed"] else "")
        by.setdefault(key, []).append(r)
    for key, grp in sorted(by.items(), key=lambda kv: -len(kv[1])):
        fams = {}
        for r in grp:
            fams[r["variant"]] = fams.get(r["variant"], 0) + 1
        out.append(f"  {len(grp):>3}  {key:<28} " + " ".join(f"{k}×{v}" for k, v in sorted(fams.items())))
    out += ["", "Reading it: SIGNAL_ABSENT is a substrate/basecalling result — no gate change recovers it.",
            "            AT_BACKGROUND is k>0 that the sample's own null already explains — also unrecoverable.",
            "            GATE_REJECTED names a threshold, and a threshold can be re-measured on the negatives.",
            "            Anything else is coverage, not detection."]
    if any(r["exit"] == "GATE_REJECTED" for r in powered):
        out += ["", "⚠ Do NOT move a gate on this table: it is the carrier arm. Price any change on the",
                "  independent negatives first (`lists/neg_validation.txt`), per decisions.md 2026-07-28."]
    return "\n".join(out) + "\n"


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="why a labelled carrier was not called, gate by gate")
    ap.add_argument("--outdir", required=True, help="directory holding <sample>.dupc.json")
    ap.add_argument("--truth", required=True, help="truth.tsv from pcr_truth")
    ap.add_argument("--family", default=None, help="restrict to one clinical variant family (e.g. 59dupC)")
    ap.add_argument("--out", default=None, help="write the per-sample table here (TSV)")
    a = ap.parse_args(argv)

    truth = read_truth(a.truth)
    rows = diagnose(a.outdir, truth, family=a.family)
    if not rows:
        # An empty cohort must never be silent: it reads exactly like a clean sweep. Print what the table
        # actually holds, so a vocabulary mismatch is visible in one line instead of one round-trip.
        import collections
        st = collections.Counter((t.get("status") or "").strip() for t in truth)
        fam = collections.Counter((t.get("variant") or "").strip()
                                  for t in truth if is_carrier(t.get("status")))
        print(f"no labelled carrier matched in {a.truth} ({len(truth)} rows)", file=sys.stderr)
        print(f"  status values present : {dict(st)}", file=sys.stderr)
        print(f"  carrier families      : {dict(fam) or '(none — no row passed is_carrier)'}",
              file=sys.stderr)
        if a.family:
            print(f"  --family {a.family!r} matched none of those families", file=sys.stderr)
        return 1
    if a.out:
        with open(a.out, "w") as fh:
            fh.write("\t".join(ROW_COLS) + "\n")
            for r in rows:
                fh.write("\t".join("" if r.get(c) is None else str(r.get(c)) for c in ROW_COLS) + "\n")
        print(f"[miss_report] per-sample table → {a.out}  ⚠ patient data, keep on /scratch", file=sys.stderr)
    sys.stdout.write(summarize(rows))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

def arbiter_lengths_of(outdir: str, sample: str) -> "list | None":
    """The arbiter's two allele lengths as ints, or None. Impure.

    Sits beside `arbiter_lengths`, which owns the read — the `score.json` key is
    `arbiter_vs_caller_length.arbiter`, which was guessed wrong three times on 2026-08-04 before anyone
    opened a real file. One reader, one place to be wrong."""
    s = arbiter_lengths(outdir, sample)
    if not s or s == "none":
        return None
    try:
        return [int(x) for x in str(s).split(",")]
    except ValueError:
        return None

def find_bam(sample: str, outdirs) -> "str | None":
    """The `.vntr.bam` `run` wrote for `sample`, from the FIRST directory holding it. Impure."""
    dirs = [outdirs] if isinstance(outdirs, str) else list(outdirs)
    for d in dirs:
        p = os.path.join(d, f"{sample}.vntr.bam")
        if os.path.exists(p):
            return p
    return None

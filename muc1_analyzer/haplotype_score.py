#!/usr/bin/env python3
"""haplotype_score — per-carrier MUC1 haplotype score, formalizing the two orthogonal prognostic axes.

Turns `plot_two_axis` into numbers. Two SEPARATE sub-scores (they are ORTHOGONAL — do NOT sum them into one
clinical number without the caveat):

  ONSET sub-score      = w_onset * onset_index                       (onset_index = N_m/L_m = frameshift-tail
                         fraction of the MUTANT allele; early frameshift -> long tail -> more MUC1fs neoprotein
                         -> earlier onset).  Validated: Spearman(onset_index, age_onset) rho ~ -0.89 (n=11).

  SEVERITY sub-score   = splice_term + buffer_term
                         splice_term = +w_splice  if the MUTANT allele is VNTR-conserved (rs4072037-C, MUC1-TR
                                                   -> keeps VNTR+frameshift -> full neoprotein -> SEVERE)
                                       -w_splice  if the mutant allele is splice-out (rs4072037-T, MUC1-Y
                                                   -> VNTR+frameshift excluded -> PROTECTED).
                                       Validated (direction): Mann-Whitney severity p ~ 0.02, ESRD Fisher p ~ 0.07.
                         buffer_term = -w_buffer  if the HEALTHY allele is LONGER than the mutant allele
                                                   (more WT MUC1 buffers the dominant-negative).  ⚠ EXPLORATORY,
                                                   OFF BY DEFAULT (w_buffer=0): when tested it DEGRADES the
                                                   severity direction (Spearman +0.63 p 0.03 -> +0.49 p 0.11;
                                                   P6 is a splice-T ESRD case the buffer pushes the wrong way) —
                                                   the mechanistic intuition is not borne out at n=12, same as
                                                   the raw ratio being ns for onset. Kept as a `--w-buffer` knob
                                                   for exploration, not in the default score.

⚠ n=12, single cohort. The weights are MECHANISTIC PRIORS (set from biology + signal strength), NOT fit to the
data — fitting on n=12 would overfit. The outcome correlations below are a DIRECTION sanity check only
(hypothesis-generating), not a validated predictive model. The enhancer block (rs2070803/rs4971101) is NOT a
term: it is r^2=1-collinear with VNTR length, which the mutant/healthy lengths already carry directly.

    python -m muc1_analyzer.haplotype_score --matrix Clinical-DATA/muc1_matrix.xlsx
"""
from __future__ import annotations
import argparse

_SEVRANK = {"Normal": 1, "CKD III-IV": 2, "CKD III/IV": 2, "ESRD": 3}

# rs4072037 splice status of the MUTANT allele, DIRECTLY PHASED single-molecule on LR-PCR (phase_fs_snp,
# 2026-07-17) for the 4 rs4072037-HET carriers — the only ones that previously needed the cis inference.
# Value = mut_T (True = T/MUC1-Y/protected, False = C/MUC1-TR/severe). All 4 MATCH the cis inference; this only
# upgrades the provenance (dosage-homozygous carriers are already direct). See decisions.md 2026-07-17.
RS4072037_PHASED_HET = {"P3": True, "P4": True, "P10": False, "P11": False}


def _num(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def _mut_splice_is_T(rd, cd, mut_long):
    """Does the MUTANT allele carry rs4072037-T (splice-out, MUC1-Y)? True/False/None.

    dosage 2 -> both alleles T -> mutant is T. dosage 0 -> mutant is C. het (1) -> resolve by cis phase:
    rs4072037-T is cis with the SHORT VNTR (~94% rule) and, when a DEL anchor exists, with the ALT/REF as
    recorded in `rs4072037_cis_del`. Falls back to the short/long rule when the DEL anchor is absent."""
    rdi = int(rd) if str(rd).isdigit() else None
    if rdi == 2:
        return True
    if rdi == 0:
        return False
    if rdi == 1:
        cd = str(cd)
        if cd == "DEL_cis_ALT":          # T cis with DEL/ALT -> T on that allele
            return not mut_long           # ALT/T is cis-short -> mutant is T iff mutant is the short allele
        if cd == "DEL_cis_REF":
            return mut_long
        # no DEL anchor: T tags the SHORT VNTR -> mutant carries T iff the mutant allele is the SHORT one
        return not mut_long
    return None


def splice_call(sample, dosage, mut_T_inf):
    """PURE: resolve the mutant-allele splice status + its provenance.
    Precedence: single-molecule long-read phasing (het carriers) > dosage-homozygous (direct) > cis inference."""
    if sample in RS4072037_PHASED_HET:
        return RS4072037_PHASED_HET[sample], "phased(long-read)"
    if str(dosage) in ("0", "2"):
        return mut_T_inf, "dosage(direct)"
    return mut_T_inf, "inferred(cis)"


def load(path):
    import openpyxl
    wb = openpyxl.load_workbook(path, data_only=True)
    ws = wb[wb.sheetnames[0]]
    rows = list(ws.iter_rows(values_only=True))
    hdr = [str(c) for c in rows[0]]
    gi = {h: i for i, h in enumerate(hdr)}
    out = []
    for r in rows[1:]:
        if str(r[gi["carrier"]]) != "True":
            continue
        h1, h2, ratio, pos = (_num(r[gi[k]]) for k in ("hap1_len", "hap2_len", "ratio_vntr", "variant_position"))
        if None in (h1, h2, ratio, pos):
            continue
        l_mut = max(h1, h2) if ratio >= 1 else min(h1, h2)
        l_healthy = min(h1, h2) if ratio >= 1 else max(h1, h2)
        mut_long = l_mut == max(h1, h2)
        dosage = r[gi["rs4072037_dosage"]]
        mut_T_inf = _mut_splice_is_T(dosage, r[gi["rs4072037_cis_del"]], mut_long)
        sample = r[gi["Sample"]]
        mut_T, source = splice_call(sample, dosage, mut_T_inf)
        rf = str(r[gi["renal_function"]]).strip()
        out.append({
            "s": sample,
            "l_mut": l_mut, "l_healthy": l_healthy, "pos": pos,
            "onset_index": 1 - pos / l_mut if l_mut else None,
            "age_onset": _num(r[gi["age_onset"]]),
            "mut_T": mut_T,                      # True=protected(T), False=severe(C), None=unknown
            "mut_T_inferred": mut_T_inf,         # the cis-rule call (for the observed-vs-inferred check)
            "splice_source": source,             # phased(long-read) | dosage(direct) | inferred(cis)
            "buffer": l_healthy > l_mut,         # healthy allele longer -> protective
            "rf": rf, "sev": _SEVRANK.get(rf),
        })
    return out


def score(rows, *, w_onset=2.0, w_splice=1.0, w_buffer=0.5):
    """Attach onset_score / severity_score (two orthogonal sub-scores). Pure; no I/O."""
    for d in rows:
        d["onset_score"] = round(w_onset * d["onset_index"], 3) if d["onset_index"] is not None else None
        splice = None if d["mut_T"] is None else (-w_splice if d["mut_T"] else +w_splice)
        buffer = -w_buffer if d["buffer"] else 0.0
        d["splice_term"] = splice
        d["buffer_term"] = buffer
        d["severity_score"] = round(splice + buffer, 3) if splice is not None else None
        # optional composite (secondary, reconflates the two axes) — reported, not primary
        d["composite"] = (round(d["onset_score"] + d["severity_score"], 3)
                          if None not in (d["onset_score"], d["severity_score"]) else None)
    return rows


def severity_mwu(rows, flips=None):
    """PURE: renal-rank Mann-Whitney (C-severe > T-protected, one-sided 'greater'), optionally flipping the
    mutant-allele splice (mut_T) for the samples in `flips` ({sample: new_bool}). Returns {p, nC, nT} or None
    if either group is empty. Powers the het-carrier sensitivity check: flip one carrier's splice and re-run
    to bound the axis. (For P10 the flip is a ROBUSTNESS check, not a real ambiguity — its mutant allele is
    the 80-copy one, del8_27 @ repeat 37, confirmed by LR-PCR + AS + visual review; the automated caller
    misses del8_27 but the allele assignment is grounded — see DATA-Patients-MUC1-vntr.xlsx.)"""
    from scipy.stats import mannwhitneyu
    flips = flips or {}
    sv = [d for d in rows if d.get("sev") is not None and d.get("mut_T") is not None]

    def mt(d):
        return flips[d["s"]] if d["s"] in flips else d["mut_T"]

    C = [d["sev"] for d in sv if mt(d) is False]
    T = [d["sev"] for d in sv if mt(d) is True]
    if not (C and T):
        return None
    _, p = mannwhitneyu(C, T, alternative="greater")
    return {"p": round(float(p), 4), "nC": len(C), "nT": len(T)}


def _report(rows):
    import numpy as np
    from scipy.stats import spearmanr, mannwhitneyu, fisher_exact

    hdr = f"{'P':<5}{'L_mut':>6}{'L_heal':>7}{'pos':>5}{'onset_i':>8}{'onsetSc':>8}{'mutSpl':>7}" \
          f"{'splice_src':>18}{'sevSc':>7}{'comp':>6}{'age_on':>7}  renal"
    print(hdr)
    print("-" * len(hdr))
    for d in sorted(rows, key=lambda x: (x["onset_score"] is None, -(x["onset_score"] or 0))):
        spl = {True: "T(-)", False: "C(+)", None: "?"}[d["mut_T"]]
        print(f"{str(d['s']):<5}{d['l_mut']:>6.0f}{d['l_healthy']:>7.0f}{d['pos']:>5.0f}"
              f"{d['onset_index']:>8.2f}{(d['onset_score'] if d['onset_score'] is not None else float('nan')):>8.2f}"
              f"{spl:>7}{d.get('splice_source', '?'):>18}"
              f"{(d['severity_score'] if d['severity_score'] is not None else float('nan')):>7.2f}"
              f"{(d['composite'] if d['composite'] is not None else float('nan')):>6.1f}"
              f"{(d['age_onset'] if d['age_onset'] is not None else float('nan')):>7.0f}  {d['rf']}")

    # provenance summary + observed-vs-inferred concordance on the het carriers
    from collections import Counter
    prov = Counter(d.get("splice_source") for d in rows if d["mut_T"] is not None)
    het = [d for d in rows if d.get("splice_source") == "phased(long-read)"]
    concord = sum(1 for d in het if d["mut_T"] == d.get("mut_T_inferred"))
    print(f"\nsplice provenance: {dict(prov)}")
    if het:
        disc = [(d["s"], d.get("mut_T_inferred"), d["mut_T"]) for d in het
                if d["mut_T"] != d.get("mut_T_inferred")]
        print(f"observed (long-read) vs inferred (cis) on het carriers: {concord}/{len(het)} concordant")
        for s, i, o in disc:
            tag = {True: "T", False: "C", None: "?"}
            print(f"  CORRECTED by observation: {s}  inferred={tag[i]} -> observed={tag[o]} "
                  f"(the cis anchor mis-phased this carrier)")

    # ── DIRECTION sanity checks (confirmatory only, underpowered) ──────────────
    print("\nDirection checks (n small, hypothesis-generating):")
    on = [d for d in rows if d["onset_score"] is not None and d["age_onset"] is not None]
    if len(on) >= 3:
        rho, p = spearmanr([d["onset_score"] for d in on], [d["age_onset"] for d in on])
        print(f"  ONSET   : onset_score vs age_onset   Spearman rho={rho:+.2f} p={p:.1e} (n={len(on)}) "
              f"[expect NEGATIVE: higher score -> earlier onset]")
    sv = [d for d in rows if d["severity_score"] is not None and d["sev"] is not None]
    if len(sv) >= 3:
        rho, p = spearmanr([d["severity_score"] for d in sv], [d["sev"] for d in sv])
        print(f"  SEVERITY: severity_score vs renal-rank Spearman rho={rho:+.2f} p={p:.1e} (n={len(sv)}) "
              f"[expect POSITIVE: higher score -> worse]")
        C = [d for d in sv if d["mut_T"] is False]
        T = [d for d in sv if d["mut_T"] is True]
        if C and T:
            eC, eT = sum(d["sev"] == 3 for d in C), sum(d["sev"] == 3 for d in T)
            _, pf = fisher_exact([[eC, len(C) - eC], [eT, len(T) - eT]])
            _, pm = mannwhitneyu([d["sev"] for d in C], [d["sev"] for d in T], alternative="greater")
            print(f"            splice split: ESRD {eC}/{len(C)} (conserved-C) vs {eT}/{len(T)} (splice-T)  "
                  f"Fisher p={pf:.2f}; severity MWU p={pm:.3f}")

    # ── SEVERITY-axis SENSITIVITY to the 4 phased-het splice calls (robustness) ──────────
    # Flip each het carrier one at a time and re-run the MWU to see which assignments the axis rests on.
    # P10: its mutant allele is the 80-copy one (del8_27 @ repeat 37, confirmed LR-PCR + AS +
    # visual) → C is grounded; its flip is a robustness check, not a real ambiguity. The axis is also robust
    # to it (significant either way). The pivotal carriers (P11, P3) are directly long-read-phased.
    base = severity_mwu(rows)
    if base:
        print("\nSeverity-axis sensitivity (renal-rank MWU, C-severe > T-protected; flip one het carrier):")
        print(f"  as-called (phased)          : p={base['p']}   (nC={base['nC']} nT={base['nT']})")
        for s in sorted(RS4072037_PHASED_HET):
            cur = next((d for d in rows if str(d["s"]) == s), None)
            if cur is None or cur["mut_T"] is None:
                continue
            fl = severity_mwu(rows, {s: not cur["mut_T"]})
            if not fl:
                continue
            frm, to = ("T", "C") if cur["mut_T"] else ("C", "T")
            tag = "   <- del8_27 on the 80-copy allele (confirmed LR-PCR+AS+visual) -> robustness only" if s == "P10" else ""
            print(f"  if {s} splice flipped {frm}->{to}   : p={fl['p']}   "
                  f"(L_mut={cur['l_mut']:.0f} / L_heal={cur['l_healthy']:.0f}){tag}")


_SEVCOL = {"Normal": "#2e7d32", "CKD III-IV": "#ef6c00", "CKD III/IV": "#ef6c00", "ESRD": "#b71c1c"}


def render_2d(rows, out_prefix):
    """2-D positioning figure: onset_score (x) vs severity_score (y), coloured by renal function.

    The two axes are ORTHOGONAL by construction; this plots WHERE each carrier sits. Bottom-right = early
    onset but splice-protected (e.g. the P6 ESRD exception surfaces as a low-severity/mid-onset point)."""
    import os
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    pts = [d for d in rows if d["onset_score"] is not None and d["severity_score"] is not None]
    fig, ax = plt.subplots(figsize=(7.6, 6.2))
    ax.axhline(0, color="#bbb", lw=1, zorder=1)
    for d in pts:
        ax.scatter(d["onset_score"], d["severity_score"], s=150,
                   color=_SEVCOL.get(d["rf"], "#888"), edgecolor="white", lw=1.3, zorder=3)
        ax.annotate(str(d["s"]), (d["onset_score"], d["severity_score"]), textcoords="offset points",
                    xytext=(8, 4), fontsize=8, color="#333")
    ax.set_xlabel("ONSET sub-score  (= w·onset_index, frameshift-tail fraction) →  earlier onset", fontsize=9.5)
    ax.set_ylabel("SEVERITY sub-score  (mutant-allele splice: +C conserved / −T protected) →  worse", fontsize=9.5)
    ax.set_title("MUC1 haplotype score — two orthogonal axes, per carrier", fontsize=11.5, fontweight="bold")
    ax.text(0.98, 0.02, "colour = renal function\n(green Normal · orange CKD · red ESRD)",
            transform=ax.transAxes, ha="right", va="bottom", fontsize=8,
            bbox=dict(boxstyle="round,pad=0.3", fc="#f6f6f6", ec="#ccc"))
    ax.grid(True, alpha=0.22)
    plt.tight_layout()
    os.makedirs(os.path.dirname(out_prefix) or ".", exist_ok=True)
    for ext in ("png", "pdf"):
        fig.savefig(f"{out_prefix}.{ext}", dpi=200, bbox_inches="tight")
    print(f"\nwrote {out_prefix}.png / .pdf  (n={len(pts)})")


def main(argv=None):
    ap = argparse.ArgumentParser(prog="muc1_analyzer.haplotype_score",
                                 description="Per-carrier MUC1 haplotype score (two orthogonal sub-scores)")
    ap.add_argument("--matrix", default="Clinical-DATA/muc1_matrix.xlsx")
    ap.add_argument("--w-onset", type=float, default=2.0, help="weight of onset_index in the ONSET sub-score")
    ap.add_argument("--w-splice", type=float, default=1.0, help="weight of the mutant-allele splice term (SEVERITY)")
    ap.add_argument("--w-buffer", type=float, default=0.0,
                    help="weight of the contralateral healthy-length buffer (EXPLORATORY, OFF by default: "
                         "degrades severity direction at n=12; try 0.5 to explore)")
    ap.add_argument("--tsv", default=None, help="optional: write the per-carrier table to this TSV")
    ap.add_argument("--fig", default=None,
                    help="optional: write the 2-D onset×severity positioning figure to this prefix "
                         "(e.g. figures/muc1_haplotype_axes)")
    a = ap.parse_args(argv)
    rows = score(load(a.matrix), w_onset=a.w_onset, w_splice=a.w_splice, w_buffer=a.w_buffer)
    _report(rows)
    if a.fig:
        render_2d(rows, a.fig)
    if a.tsv:
        import csv
        cols = ["s", "l_mut", "l_healthy", "pos", "onset_index", "onset_score", "mut_T",
                "splice_term", "buffer", "buffer_term", "severity_score", "composite", "age_onset", "rf", "sev"]
        with open(a.tsv, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", extrasaction="ignore")
            w.writeheader()
            w.writerows(rows)
        print(f"\n-> {a.tsv}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

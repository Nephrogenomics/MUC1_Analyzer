"""Two-axis clinical summary appended to the VNTR caller report when a GENOMIC (chr1) BAM is given.

The VNTR-reference caller reports motif structure + frameshift (the ONSET axis) but works on the
VNTR-only reference, which (a) on LR-PCR competitively DEPLETES the long allele → the top-2 haplotype
lengths under-call it, and (b) cannot see rs4072037 (the splice SNP sits ~370 bp 5' of the tandem,
outside the reference). This module reads the GENOMIC BAM and adds the two signals the score needs:

  * ONSET  ← the frameshift position (already in the caller) on the ARBITER length (flank-to-flank,
    alignment-free `vntr_raw_length.call_alleles`) that recovers the long allele the reference caller drops.
  * SEVERITY ← rs4072037 phased to each allele (`phase_fs_snp`): which base (C = MUC1-TR, retained →
    severe / T = MUC1-Y, spliced out → protected) sits on hap1, hap2, or both, and specifically on the
    MUTANT allele.

Pure assembly (`build_summary`) + renderers (text / PDF flowables). Gated behind `--genomic-bam`: the
default caller path is untouched, so the non-regression bit-identity holds.
"""
from __future__ import annotations

import re

# C = MUC1-TR (VNTR retained → severe) · T = MUC1-Y (VNTR spliced out → protected)
SPLICE = {"C": ("MUC1-TR", "VNTR retained → severe"),
          "T": ("MUC1-Y", "VNTR spliced out → protected")}


def _genotype(pairs, min_af=0.15):
    """rs4072037 genotype from the raw (base, copies) pairs, by minor-allele fraction."""
    n = {"C": 0, "T": 0}
    for b, _ in pairs:
        if b in n:
            n[b] += 1
    tot = n["C"] + n["T"]
    if tot == 0:
        return None
    minor = min(n["C"], n["T"]) / tot
    if minor >= min_af:
        return "C/T"
    return "C/C" if n["C"] >= n["T"] else "T/T"


def _nearest_hap(length, haps, tol=8):
    """The hap (rank) whose length is closest to `length` within tol, else None."""
    if length is None:
        return None
    cand = sorted(((abs(h["length"] - length), h["rank"]) for h in haps))
    if cand and cand[0][0] <= tol:
        return cand[0][1]
    return None


def _parse_indel(indel):
    """('59dupC', 16) from a summary like '59dupC ~ repeat 16' (or 'X-59dupC ~ repeat 16')."""
    if not indel:
        return None, None
    m = re.search(r"~\s*repeat\s*(\d+)", indel)
    repeat = int(m.group(1)) if m else None
    variant = re.sub(r"\s*~\s*repeat.*$", "", indel).lstrip("X-").strip()
    return (variant or None), repeat


def _coverage_gate(cov, min_reads, min_strand):
    """Is a frameshift decoded on this haplotype's consensus reliable?

    A frameshift is only trusted when the consensus it was decoded from rests on enough
    reads AND on BOTH strands (a single-strand or thin consensus fabricates homopolymer-slip
    artefacts — e.g. the false 60dupA on a low-coverage allele-specific pass). Returns
    (pass, cov_dict). `cov` is {'n_reads','n_fwd','n_rev'} or None; None → no support info was
    supplied → do NOT gate (backward compatible: the gate is inert unless coverage is provided).
    """
    if not cov:
        return True, None
    n = cov.get("n_reads", 0)
    f = cov.get("n_fwd", 0)
    r = cov.get("n_rev", 0)
    ok = (n >= min_reads) and (f >= min_strand) and (r >= min_strand)
    return ok, {"n_reads": n, "n_fwd": f, "n_rev": r}


_HG38_CHR1_LEN = {248_956_422, 249_250_621}
_T2T_CHR1_LEN = {248_387_328}


def _genomic_build(bam, ref=None):
    """(build, chr1_contig_name) from the genomic BAM header. build ∈ {'hg38','t2t',None}. rs4072037 is an
    hg38 COORDINATE, so on a T2T BAM the pileup must move to the CHM13 coordinate AT THE T2T CONTIG NAME
    (chr1 / NC_060925.1 / CP068277.2) — reading the hg38 position on T2T is a different locus (a spurious
    genotype). Returns the contig name so both the SNP pileup and the length fetch can target it. Cheap."""
    try:
        import pysam
        mode = "rc" if str(bam).endswith(".cram") else "rb"
        kw = {"reference_filename": ref} if (mode == "rc" and ref) else {}
        with pysam.AlignmentFile(bam, mode, **kw) as af:
            pairs = list(zip(af.references, af.lengths))
    except Exception:
        return None, None
    lens = {ln for _, ln in pairs}
    if lens & _HG38_CHR1_LEN:
        build, target, prefer = "hg38", _HG38_CHR1_LEN, ("chr1", "1", "CM000663.2", "NC_000001.11")
    elif lens & _T2T_CHR1_LEN:
        build, target, prefer = "t2t", _T2T_CHR1_LEN, ("chr1", "NC_060925.1", "CP068277.2", "1")
    else:
        return None, None
    names = {n for n, _ in pairs}
    c1 = next((c for c in prefer if c in names), None) or next((n for n, ln in pairs if ln in target), None)
    return build, c1


def build_summary(genomic_bam, haps, *, chrom="chr1", snp_pos=None, copy_offset=4,
                  maxmm=9, ref=None, start=None, end=None, mutation=None,
                  gate_min_reads=5, gate_min_strand=2, smear=False, position_confident=None):
    """Assemble the two-axis summary from a genomic BAM + the caller's haplotypes.

    `haps` = [{'rank': int, 'length': int, 'indel': str, 'coverage': {...}?}]  (indel = '' for a WT
    allele). Optional per-hap `coverage` = {'n_reads','n_fwd','n_rev'} of the VNTR-reference read group
    that built that haplotype's consensus; when present, a plain-call (caller_decomposition) frameshift
    onset is COVERAGE-GATED — flagged unreliable when the consensus rests on < `gate_min_reads` reads or
    on fewer than `gate_min_strand` reads per strand (kills low-coverage decode artefacts). Absent →
    the coverage arm is inert. `smear=True` (the caller's multi-contig-smear / thin-top length warning)
    gates ANY plain-call onset wholesale — the per-contig consensus is untrustworthy, so a frameshift
    decoded on it is an artefact (the false 60dupA measured on a negative). A SUPPLIED mutation is authoritative, never gated.
    `mutation` (optional) = a VALIDATED call supplied from outside — the clinical/lab nomenclature or a
    careful per-allele caller run — as {'variant': str, 'repeat': int, 'mutant_length': int}. When the
    plain-call decomposition found no frameshift (e.g. a PCR-DEPLETED long mutant the top-N misses), this
    is AUTHORITATIVE for the onset: it attaches to the allele nearest `mutant_length`, promotes it out of
    low-confidence (the mutation confirms it is real), and drives the severity. The report only DISPLAYS
    what it is given — it never re-runs a caller here, so it cannot introduce a false carrier.
    Returns a JSON-friendly dict; `available` is False when the BAM has no locus-spanning reads.
    Lazy-imports the two top-level tools so the default caller path never loads them.
    """
    import phase_fs_snp as P
    import vntr_raw_length as V

    explicit_snp = snp_pos is not None
    snp_pos = snp_pos if explicit_snp else P.SNP_POS
    build, gchr1 = _genomic_build(genomic_bam, ref)

    # ── ARBITER length (flank-to-flank, recovers the long allele) ──
    # alignment-free (scans read sequence). RESTRICT the fetch to the MUC1 locus: on a WHOLE-GENOME CRAM
    # (adaptive-sampling/WGS) an until-eof scan decodes every read → ~18 min/sample; the locus reads are all
    # we need. The default window is hg38 (snp_pos-based). If it yields NOTHING — a contig-name mismatch, OR
    # a genomic BAM on a DIFFERENT BUILD (T2T puts the VNTR at ~154.3 Mb, off the hg38 window) — resolve the
    # locus from the BAM header (build-aware), then fall back to a full until-eof scan as a last resort.
    a_start, a_end = (start, end) if start is not None else (snp_pos - 12000, snp_pos + 3000)

    def _data(c, s, e):
        return V.copies_from_bam(genomic_bam, c, s, e, copy_offset, maxmm, ref)   # [(copies, hp), …]

    try:
        data = _data(chrom, a_start, a_end)
    except Exception:
        data = []
    if not data and start is None:                            # off-locus / wrong build / misnamed contig
        from .detectors import vntr as _VN
        region = _VN._resolve_region(genomic_bam, ref)        # 'contig:start-end' per build, or None
        if region and ":" in region:
            rc, rng = region.split(":", 1)
            try:
                rs, re_ = (int(v) for v in rng.split("-", 1))
                data = _data(rc, rs, re_)
            except Exception:
                data = []
    if not data:
        try:
            data = _data(chrom, None, None)                   # last resort: whole-BAM scan (flank-filtered)
        except Exception:
            data = []
    copies = [x for x, _ in data]
    # HP-aware allele call: a haplotagged genomic BAM → ONE allele per HP (consistent with the reliable
    # block), else the pooled peak-caller. Fixes the two-axis printing HOM while the reliable block prints HET.
    arb = (V._alleles_by_hp(data) or V.call_alleles(sorted(copies))) if data \
        else {"alleles": [], "counts": [], "long_low_conf": False}

    # ── rs4072037 phased to each allele length (needs the SNP contig by name; isolate its failure so the
    #    arbiter length still shows if the pileup contig is missing/misnamed). SNP_POS is an hg38 coordinate:
    #    on a T2T genomic BAM read it at the CHM13v2.0 coordinate (config.T2T, chr1:154,330,839) AT THE T2T
    #    CONTIG NAME — the hg38 position on T2T is a different locus (a spurious genotype). An explicit
    #    snp_pos from the caller always wins. ──
    if explicit_snp:
        snp_chrom, snp_at = chrom, snp_pos
    elif build == "t2t":
        from . import config as _CFG
        snp_chrom, snp_at = (gchr1 or chrom), _CFG.T2T["SPLICE_SNP_EXON2"].start
    else:
        snp_chrom, snp_at = chrom, snp_pos
    snp_pos = snp_at                                            # the coordinate actually read (for the report)
    try:
        pairs = P.collect_reads(genomic_bam, snp_chrom, snp_at, copy_offset, maxmm, ref)
    except Exception:
        pairs = []
    phase = P.phase_snp_to_length(pairs)
    genotype = _genotype(pairs)

    if not copies and not pairs:
        return {"available": False,
                "note": "no reads span the MUC1 flank anchors in the genomic BAM "
                        "(wrong locus / coordinates / too shallow)"}

    # per-base rs4072037 → allele length → hap
    per_base = []
    for base in ("C", "T"):
        v = phase["phase"].get(base)
        if not v:
            continue
        name, sev = SPLICE[base]
        per_base.append({"base": base, "length": v["length"], "n": v["n"],
                         "splice": f"{name} ({sev})",
                         "hap": _nearest_hap(v["length"], haps)})

    # LOW-CONFIDENCE length: `_allele_length` deliberately picks the HIGHEST substantial peak (to catch a
    # PCR-depleted true long allele), but a handful of chimeric/miscalled long reads can fabricate one on a
    # length-homozygous sample (the CTRL false-79). Length alone cannot tell a depleted TRUE allele from
    # noise (only the mutation call can), so we FLAG a base whose phased peak rests on FEW reads rather than
    # assert a hard HET/severity. Support = reads of that base within tol of its phased length.
    import collections as _c
    cs_by_base = _c.defaultdict(list)
    for b, cc in pairs:
        if b in SPLICE and cc is not None:
            cs_by_base[b].append(cc)
    for pb in per_base:
        cs = cs_by_base[pb["base"]]
        support = sum(1 for c in cs if pb["length"] is not None and abs(c - pb["length"]) <= 8)
        pb["support_at_length"] = support
        pb["low_conf"] = bool(pb["length"] is not None and cs
                              and (support < 15 or support < 0.15 * len(cs)))

    # base → hap by nearest length, then BY ELIMINATION (a het 2-hap case where the reference caller
    # under-called the long allele → its length is far from the true one, so nearest-length leaves it
    # unassigned; the non-mutant hap must carry the OTHER base).
    base_by_hap = {pb["hap"]: pb["base"] for pb in per_base if pb["hap"] is not None}
    if genotype == "C/T" and len(haps) == 2:
        free_haps = [h["rank"] for h in haps if h["rank"] not in base_by_hap]
        free_bases = [pb["base"] for pb in per_base if pb["base"] not in base_by_hap.values()]
        if len(free_haps) == 1 and len(free_bases) == 1:
            base_by_hap[free_haps[0]] = free_bases[0]
            for pb in per_base:
                if pb["base"] == free_bases[0] and pb["hap"] is None:
                    pb["hap"], pb["by_elimination"] = free_haps[0], True

    snp_len = {pb["base"]: pb["length"] for pb in per_base}   # SNP-phased (authoritative) allele length
    lc_by_base = {pb["base"]: pb["low_conf"] for pb in per_base}
    hom_base = genotype[0] if genotype in ("C/C", "T/T") else None
    out_haps = []
    for h in haps:
        variant, repeat = _parse_indel(h.get("indel", ""))
        base = hom_base or base_by_hap.get(h["rank"])
        out_haps.append({"rank": h["rank"], "length": h["length"], "rs4072037": base,
                         "arbiter_length": snp_len.get(base) if base else None,
                         "rs4072037_low_conf": bool(base and lc_by_base.get(base)),
                         "variant": variant, "repeat": repeat, "mutant": bool(variant)})

    # ── SUPPLIED mutation (authoritative): attach a validated call the plain-call decomposition missed
    # (e.g. a PCR-depleted long mutant, supplied from the lab nomenclature) to the nearest allele and
    # promote it out of low-confidence. AUTO-detecting it here was tried and dropped — a standalone
    # per-read/​per-allele re-run over-calls a homopolymer dupC on the abundant WT allele; robust auto
    # onset belongs in the validated detector→reporter pipeline, not this report. ──
    mut_source = None
    if mutation and mutation.get("variant"):
        # SUPPLIED = the validated lab nomenclature → AUTHORITATIVE. Clear ANY plain-call decomposition
        # mutant first: on noisy low-coverage AS reads the caller over-calls a spurious frameshift (a false
        # 60dupA on a NEGATIVE control, a wrong delCC where the truth is dupC), which must not win over the
        # validated call. Then attach the supplied mutation to the ALLELE whose SNP-PHASED length is nearest
        # mut_len (the true allele size; the caller under-calls the long allele and can split the short into
        # two haps that elimination crosses base↔hap), and promote it out of low-confidence.
        for h in out_haps:
            h["mutant"], h["variant"], h["repeat"] = False, None, None
        mlen = mutation.get("mutant_length")
        tgt = None
        if mlen is not None and per_base:
            pb = min(per_base, key=lambda p: abs((p["length"] if p["length"] is not None else 10 ** 9) - mlen))
            tgt = next((h for h in out_haps if h["rank"] == pb.get("hap")), None)
        if tgt is None and mlen is not None:
            tgt = min(out_haps, key=lambda h: abs((h.get("arbiter_length") or h["length"]) - mlen),
                      default=None)
        if tgt is not None:
            tgt["mutant"] = True
            tgt["variant"], tgt["repeat"] = mutation["variant"], mutation.get("repeat")
            tgt["rs4072037_low_conf"] = False
            arb["long_low_conf"] = False          # the validated mutation confirms the het structure
            for pb in per_base:
                if pb["base"] == tgt["rs4072037"]:
                    pb["low_conf"] = False
            mut_source = "supplied"

    # ── onset + severity from the MUTANT allele ──
    onset = severity = None
    mut = next((h for h in out_haps if h["mutant"]), None)
    if mut:
        # Mutant allele length: a SUPPLIED value wins (the caller under-calls the long allele on PCR AND can
        # split one allele into two haps whose base↔hap the elimination then crosses); else the phased/arbiter
        # length; else the caller length.
        mut_len = (mutation or {}).get("mutant_length") or mut.get("arbiter_length") or mut["length"]
        # Mutant allele rs4072037: a HOMOZYGOUS SNP wins (both alleles carry it — a minority miscalled base
        # near the mutant length must not flip it, e.g. a length-homozygous T/T carrier); else the base whose
        # SNP-PHASED length is nearest mut_len (authoritative and robust to the base↔hap crossing above);
        # else the hap base. Keep the hap display consistent with the severity.
        mb = hom_base or P.mutant_snp(phase, mut_len) or mut["rs4072037"]
        mut["rs4072037"], mut["rs4072037_low_conf"] = mb, False
        src = mut_source or "caller_decomposition"
        # a plain-call frameshift (no supplied mutation) can be a low-coverage decode artifact → flag it
        onset = {"mutant_length": mut_len, "repeat": mut["repeat"], "variant": mut["variant"],
                 "source": src, "unvalidated": src != "supplied"}
        # COVERAGE-GATE: only for a plain-call frameshift (a supplied mutation is authoritative). Two
        # triggers: (1) SMEAR — the caller flagged the length/nomenclature unreliable (multi-contig smear
        # or a thin top contig), so ANY frameshift decoded on that consensus is an artefact; (2) LOW
        # COVERAGE — the mutant haplotype's own consensus rests on too few or single-strand reads.
        if src == "caller_decomposition":
            cov_by_rank = {h["rank"]: h.get("coverage") for h in haps}
            ok, cov = _coverage_gate(cov_by_rank.get(mut["rank"]), gate_min_reads, gate_min_strand)
            if cov is not None:
                onset["coverage"] = cov
            low_cov = cov is not None and not ok
            if smear or low_cov:
                onset["coverage_gated"] = True
                onset["gate_reason"] = "multi_contig_smear" if smear else "low_coverage"
            elif cov is not None:
                onset["coverage_gated"] = False
        if mb in SPLICE:
            name, sev = SPLICE[mb]
            severity = {"mutant_length": mut_len, "rs4072037": mb,
                        "interpretation": f"{name} ({sev})", "observed": True, "low_conf": False}

    # SINGLE-AXIS MUC1_Score, folded in HERE so the report/PDF is the SINGLE deliverable (no separate
    # score.json to reconcile). The score is the FRAMESHIFT TAIL = repeats downstream of the variant on the
    # mutant allele (mut_len − position), placed on the calibration Gaussian → a 4-category verdict. onset is
    # abandoned; rs4072037 is reported for information only (no prognostic weight).
    prognosis = None
    if onset and onset.get("repeat") is not None and onset.get("mutant_length"):
        from .clinical_call import score_from_fields
        healthy = next((a for a in arb["alleles"] if a != onset["mutant_length"]), None)
        mut_base = severity.get("rs4072037") if severity else None
        # `position_confident=False` (the reporting layer's `repeat_confident`) does not withhold the score
        # but flags the category as provisional: the tail is only as reliable as the position it is measured
        # from (a prognosis must not silently ride on a repeat number the carrier reads never agreed on).
        sc = score_from_fields(onset["mutant_length"], healthy, onset["repeat"], rs4072037_mut=mut_base,
                               position_confident=position_confident)
        prognosis = {"fs_tail": sc.get("fs_tail"), "fs_percentile": sc.get("fs_percentile"),
                     "fs_z": sc.get("fs_z"), "fs_category": sc.get("fs_category"),
                     "fs_label": sc.get("fs_label"), "fs_color": sc.get("fs_color"),
                     "onset_index": sc.get("onset_index"),   # kept for back-compat / audit only
                     "coverage_gated": bool(onset.get("coverage_gated")),
                     "position_confident": position_confident,
                     "fs_confident": sc.get("fs_confident")}
        if sc.get("fs_caveat"):
            prognosis["fs_caveat"] = sc["fs_caveat"]

    snp_note = (f"rs4072037 read at the CHM13v2.0 coordinate {snp_chrom}:{snp_at}"
                if (not explicit_snp and build == "t2t") else None)
    return {"available": True, "snp_pos": snp_pos, "copy_offset": copy_offset,
            "genomic_build": build,
            "arbiter": {"alleles": arb["alleles"], "counts": arb.get("counts", []),
                        "long_low_conf": arb.get("long_low_conf", False),
                        "n_spanning": len(copies)},
            "rs4072037": {"genotype": genotype, "per_base": per_base, "n_spanning": len(pairs),
                          "note": snp_note},
            "haplotypes": out_haps, "onset": onset, "severity": severity, "prognosis": prognosis}


def _len_str(arb):
    al, ct = arb["alleles"], arb.get("counts", [])
    if not al:
        return "no clear peak"
    if len(al) == 1:
        return f"HOM ~{al[0]} copies (arbiter, n={arb['n_spanning']})"
    if arb.get("long_low_conf"):
        # do NOT assert a clean HET: the long allele is thinly supported
        lcn = ct[1] if len(ct) > 1 else "?"
        return (f"predominant ~{al[0]} copies; a LONG-allele signal ~{al[1]} (n={lcn}) is present but "
                f"LOW-CONFIDENCE — a depleted true allele or chimeric noise (confirm with the mutation call) "
                f"(arbiter, n={arb['n_spanning']})")
    return f"HET {al[0]} & {al[1]} copies (arbiter, n={arb['n_spanning']})"


def render_text(s):
    """Plain-text block for the caller report."""
    if not s.get("available"):
        return "\n  ── Two-axis clinical summary ──\n  (unavailable) %s\n" % s.get("note", "")
    # The VNTR length and the rs4072037 GENOTYPE are owned by the report's "reliable numbers" block above —
    # printing them again here duplicated (and, before the HP-aware fix, CONTRADICTED) it. This section is the
    # PROGNOSTIC layer only: per-allele phasing, onset, severity, composite. The length is re-stated ONLY when
    # it carries a caveat, because a low-confidence length undermines the allele assignment below.
    L = ["", "  " + "─" * 76,
         "  MUC1_SCORE PROGNOSIS  (frameshift tail — from the genomic BAM)",
         "  " + "─" * 76]
    if s["arbiter"].get("long_low_conf") or not s["arbiter"].get("alleles"):
        L.append(f"  ⚠ length caveat : {_len_str(s['arbiter'])}")
    # per-allele phasing detail is only meaningful when the SNP is HET (homozygous → both alleles same base)
    if s["rs4072037"]["genotype"] == "C/T":
        for pb in s["rs4072037"]["per_base"]:
            hp = f"Haplotype {pb['hap']}" if pb["hap"] else "unassigned hap"
            elim = " (by elimination)" if pb.get("by_elimination") else ""
            lc = "  [LOW-CONFIDENCE length]" if pb.get("low_conf") else ""
            L.append(f"      rs4072037-{pb['base']} ({pb['splice']})  ↔  allele {pb['length']} copies "
                     f"(n={pb['n']})  →  {hp}{elim}{lc}")
    L.append("  Haplotypes  :")
    for h in s["haplotypes"]:
        base = f"rs4072037-{h['rs4072037']}" if h["rs4072037"] else "rs4072037 unresolved"
        # the reference-caller length under-calls the long allele on LR-PCR → show the arbiter length too,
        # unless it is the low-confidence long signal (then keep it tentative)
        arb = ("" if h.get("rs4072037_low_conf") else
               (f" · ~{h['arbiter_length']} (arbiter)"
                if h.get("arbiter_length") and abs(h["arbiter_length"] - h["length"]) > 6 else ""))
        lc = " · long-allele length LOW-CONFIDENCE" if h.get("rs4072037_low_conf") else ""
        mut = f"  ◄ MUTANT: {h['variant']} ~ repeat {h['repeat']}" if h["mutant"] else ""
        L.append(f"      Hap {h['rank']}: {h['length']} copies (caller){arb} | {base}{lc}{mut}")
    if s.get("onset"):
        o = s["onset"]
        src = {"supplied": "  [supplied — validated mutation]",
               "pcr_report": "  [per-read decomposition]"}.get(o.get("source"), "")
        rep = f" at repeat {o['repeat']}" if o.get("repeat") is not None else ""
        L.append(f"  → FRAMESHIFT       : {o['variant']}{rep} "
                 f"(mutant allele = {o['mutant_length']} copies){src}")
        if o.get("coverage_gated"):
            c = o.get("coverage", {})
            why = ("multi-contig smear (unreliable consensus)" if o.get("gate_reason") == "multi_contig_smear"
                   else "thin/single-strand consensus")
            L.append(f"      ⚠ COVERAGE-GATED — frameshift decoded on a {why} "
                     f"(n={c.get('n_reads', '?')} reads, fwd/rev {c.get('n_fwd', '?')}/{c.get('n_rev', '?')}): "
                     f"UNRELIABLE — treat as NO confident onset (confirm with the validated mutation).")
    if s.get("severity"):
        v = s["severity"]
        L.append(f"  → rs4072037 (info)  : mutant allele carries rs4072037-{v['rs4072037']} "
                 f"= {v['interpretation']}  [for information only — NOT used in the score]")
    if not s.get("onset"):
        L.append("  → no frameshift on either allele (non-carrier / caller-invisible variant)")
    if s.get("prognosis"):
        p = s["prognosis"]
        gate = "  [⚠ coverage-gated — position unreliable]" if p.get("coverage_gated") else ""
        if p.get("fs_category"):
            pct = p.get("fs_percentile")
            pct_txt = f"{pct * 100:.0f}%" if pct is not None else "n/a"
            L.append(f"  → MUC1_Score (frameshift tail): {p.get('fs_tail')} repeats → "
                     f"{p.get('fs_label')} (percentile {pct_txt}, z {p.get('fs_z')}){gate}")
        else:
            L.append(f"  → MUC1_Score (frameshift tail): n/a — no measured tail{gate}")
        if p.get("fs_confident") is False:
            L.append("     ⚠ frameshift POSITION not reliable (the carrier reads did not agree on it): "
                     "the tail length, and therefore the category, is PROVISIONAL.")
    L.append("  " + "─" * 76)
    return "\n".join(L) + "\n"


def render_pdf_blocks(s, styles):
    """reportlab flowables for the PDF (styles = the dict built in write_pdf_report)."""
    from reportlab.platypus import Paragraph, Spacer
    from reportlab.lib.units import cm
    out = [Spacer(1, 0.3 * cm)]
    if not s.get("available"):
        out.append(Paragraph("Two-axis clinical summary: unavailable (%s)" % s.get("note", ""),
                             styles["meta"]))
        return out
    # length + rs4072037 genotype live in the "reliable numbers" block — not repeated here (see render_text)
    out.append(Paragraph("MUC1_Score prognosis — frameshift tail (genomic BAM)", styles["result"]))
    rows = []
    if s["arbiter"].get("long_low_conf") or not s["arbiter"].get("alleles"):
        rows.append(f"<b>⚠ length caveat:</b> {_len_str(s['arbiter'])}")
    if s["rs4072037"]["genotype"] == "C/T":
        for pb in s["rs4072037"]["per_base"]:
            hp = f"Haplotype {pb['hap']}" if pb["hap"] else "unassigned"
            rows.append(f"&nbsp;&nbsp;rs4072037-{pb['base']} ({pb['splice']}) ↔ {pb['length']} copies "
                        f"(n={pb['n']}) → {hp}")
    if s.get("onset"):
        o = s["onset"]
        gated = ("  <font color='#b00'>[⚠ COVERAGE-GATED — thin/single-strand consensus, unreliable]</font>"
                 if o.get("coverage_gated") else "")
        rows.append(f"<b>FRAMESHIFT:</b> {o['variant']} ~ repeat {o['repeat']} "
                    f"(mutant allele {o['mutant_length']} copies){gated}")
    if s.get("severity"):
        v = s["severity"]
        rows.append(f"<b>rs4072037 (info):</b> mutant allele carries rs4072037-{v['rs4072037']} "
                    f"= {v['interpretation']} — for information only, NOT used in the score")
    if s.get("prognosis"):
        p = s["prognosis"]
        gate = ("  <font color='#b00'>[⚠ coverage-gated]</font>" if p.get("coverage_gated") else "")
        if p.get("fs_category"):
            pct = p.get("fs_percentile")
            pct_txt = f"{pct * 100:.0f}%" if pct is not None else "n/a"
            col = p.get("fs_color") or "#333"
            rows.append(f"<b>MUC1_Score (frameshift tail):</b> {p.get('fs_tail')} repeats → "
                        f"<font color='{col}'><b>{p.get('fs_label')}</b></font> "
                        f"(percentile {pct_txt}, z {p.get('fs_z')}){gate}")
        else:
            rows.append(f"<b>MUC1_Score (frameshift tail):</b> n/a — no measured tail{gate}")
        if p.get("fs_confident") is False:
            rows.append("<font color='#b00'>⚠ frameshift POSITION not reliable — the tail length, and "
                        "therefore the category, is PROVISIONAL.</font>")
    for r in rows:
        out.append(Paragraph(r, styles["meta"]))
    return out

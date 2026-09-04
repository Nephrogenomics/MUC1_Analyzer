#!/usr/bin/env python3
"""allele_balance — balance allélique d'un locus MUC1 (PCR vs AS) pour départager
« biais d'amplification PCR » de « signal nanopore intrinsèque ».

Pour un BAM/CRAM, par read spanning le VNTR :
  1. mesure la **longueur d'array** (nb de C-tracts `_CTX`) → assigne l'allèle **muté vs sain** par
     la longueur attendue (si les 2 allèles diffèrent d'au moins `min_gap` répétitions) ;
  2. calcule le **taux de 8C PAR ALLÈLE** (dé-poolé) : parmi les reads de l'allèle muté, la fraction
     portant `mut_ctract` (8C = dupC) à la meilleure position — le VRAI taux muté, sans dilution par le sain ;
  3. lit **rs4072037** (base par read, parcours CIGAR) → balance allélique du SNP + **phasage intra-read**
     (quelle base SNP co-occurre avec l'allèle-muté-par-longueur) → repli quand les longueurs sont ÉGALES
     (homozygote-longueur) et intérêt propre pour le score matriciel (④ rs4072037).

Lecture : si `frac_mut` (PCR) ≪ `frac_mut` (AS) → l'allèle muté est **déplété en PCR** = biais d'amplification.
Si `frac_mut` PCR ≈ AS ≈ 0,5 mais `mut_allele_8C.frac` reste bas dans LES DEUX → **signal intrinsèque**.

⚠ pysam (nœud de calcul). **DEUX modes** :
- `--mode hp` (données HAPLOTAGGÉES = AS) : balance par tags whatshap → **FIABLE**, l'allèle muté
  s'auto-identifie par son 8C. À PRÉFÉRER dès que les reads ont un tag HP.
- `--mode length` (PCR non taggée) : proxy par le SPAN. ⚠ **DISTORDU en échelle absolue sur GRCh38
  collapsé** (le VNTR = insertion mal résolue → l'histogramme pique bien plus bas que le vrai nombre de
  répétitions ; cf. gotcha « le span sous-estime le VNTR collapsé »). Le RATIO court/long reste
  indicatif (il matche les profondeurs par-allèle du manifeste PCR), mais l'assignation par les
  longueurs ABSOLUES échoue → pour une profondeur par-allèle PCR exacte, réaligner sur la réf VNTR
  multi-contigs (voie `detect_analyzer`). `n_rep_hist` est émis pour juger la bimodalité à l'œil.

Résultat clé (2026-07-10, un porteur) : PCR déplète l'allèle LONG muté (~21× / 999× total) ; l'AS (HP) est
~2:1 → le biais est **spécifique à la PCR**. cf. docs/ROADMAP_MUC1_score.md.
"""
from __future__ import annotations
import argparse
import collections
import json
import os
import sys

from .config import GRCh38


def assign_allele(n_rep, mut_len, healthy_len, min_gap: int = 5) -> str:
    """'mut'/'healthy' par plus-proche longueur attendue ; 'ambig' si les 2 allèles sont trop
    proches (< min_gap répétitions) pour être séparés par la longueur → utiliser le SNP."""
    if abs(mut_len - healthy_len) < min_gap:
        return "ambig"
    return "mut" if abs(n_rep - mut_len) <= abs(n_rep - healthy_len) else "healthy"


def best_8c_index(idxmap: dict, mut_ctract: int = 8, min_tot: int = 8):
    """Meilleure position (index de répétition) pour le C-tract `mut_ctract` : {index,n8C,tot,frac}
    sur les positions couvertes ≥ min_tot. None si aucune. `idxmap` = {index: Counter{ctract_len}}."""
    best = None
    for i, c in idxmap.items():
        tot = sum(c.values())
        if tot >= min_tot:
            fr = c.get(mut_ctract, 0) / tot
            if best is None or fr > best["frac"]:
                best = {"index": i, "n8C": c.get(mut_ctract, 0), "tot": tot, "frac": round(fr, 3)}
    return best


def _minor_frac(counts: dict):
    """Fraction de l'allèle mineur d'un Counter de bases (None si < 2 allèles)."""
    vals = sorted(counts.values(), reverse=True)
    tot = sum(vals)
    return round(vals[1] / tot, 3) if (len(vals) >= 2 and tot) else None


def allele_balance(bam, *, mut_len_array, healthy_len_array, genome_ref=None,
                   snp_pos=None, mut_ctract=8, min_gap=5, max_span=12000, min_tot=8) -> dict:
    import pysam
    from .detectors.vntr import _qpos_at_ref, _VNTR_ANCHORS
    from .detectors.vntr_dupc import _CTX, gene_oriented
    from .detectors.splice_snp import _query_index_at

    a = _VNTR_ANCHORS
    chrom, l_anchor, r_anchor = a["chrom"], a["l_anchor"], a["r_anchor"]
    flank_const, unit = a["flank_const"], a["unit"]      # span = flank_const + unit*repeats (VNTR = insertion)
    snp_pos = snp_pos if snp_pos is not None else GRCh38["SPLICE_SNP_EXON2"].start
    mode = "rc" if str(bam).endswith(".cram") else "rb"
    kw = {"reference_filename": genome_ref} if (mode == "rc" and genome_ref) else {}
    separable = abs(mut_len_array - healthy_len_array) >= min_gap

    counts = {"mut": 0, "healthy": 0, "ambig": 0}
    idx8c = {"mut": collections.defaultdict(collections.Counter),
             "healthy": collections.defaultdict(collections.Counter)}
    n_rep_hist = collections.Counter()
    snp_counts = collections.Counter()
    phase = collections.Counter()           # (allele, snp_base) -> n
    n_span = 0

    with pysam.AlignmentFile(bam, mode, **kw) as af:
        for r in af.fetch(chrom, l_anchor, r_anchor):
            if r.is_secondary or r.is_supplementary or r.is_unmapped or r.query_sequence is None:
                continue
            snp_base = None
            qi = _query_index_at(r, snp_pos)
            if qi is not None and qi < len(r.query_sequence):
                snp_base = r.query_sequence[qi].upper()
                snp_counts[snp_base] += 1
            qL, qR = _qpos_at_ref(r, l_anchor), _qpos_at_ref(r, r_anchor)
            if qL is None or qR is None or not (0 < qR - qL <= max_span):
                continue
            # LONGUEUR d'array par le SPAN (scale linéairement) : span = flank_const + unit*repeats.
            # ⚠ NE PAS binner par le nombre de matches `_CTX` : `_CTX` ne matche que les unités à contexte
            # dupC (pas toutes les répétitions) → ne discrimine pas 44 vs 70. Les matches servent au 8C seul.
            est_rep = round((qR - qL - flank_const) / unit)
            matches = list(_CTX.finditer(gene_oriented(r.query_sequence[qL:qR])))
            if not matches:
                continue
            n_span += 1
            n_rep_hist[est_rep] += 1
            allele = assign_allele(est_rep, mut_len_array, healthy_len_array, min_gap) if separable else "ambig"
            counts[allele] += 1
            if allele in ("mut", "healthy"):
                for i, m in enumerate(matches):
                    idx8c[allele][i][len(m.group(1))] += 1
                if snp_base:
                    phase[(allele, snp_base)] += 1

    tot_bin = counts["mut"] + counts["healthy"]
    return {
        "n_span": n_span,
        "length_separable": separable,
        "expected": {"mut_len_array": mut_len_array, "healthy_len_array": healthy_len_array},
        "by_length": {**counts, "frac_mut": round(counts["mut"] / tot_bin, 3) if tot_bin else None},
        "mut_allele_8C": best_8c_index(idx8c["mut"], mut_ctract, min_tot),      # VRAI taux muté (dé-poolé)
        "healthy_allele_8C": best_8c_index(idx8c["healthy"], mut_ctract, min_tot),  # contrôle (~erreur ONT)
        "rs4072037": {"counts": dict(snp_counts), "frac_minor": _minor_frac(snp_counts),
                      "n": sum(snp_counts.values())},
        "phase_snp_x_length": {f"{al}:{b}": n for (al, b), n in sorted(phase.items())},
        "n_rep_hist": dict(sorted(n_rep_hist.items())),
    }


def _read_hp(r):
    try:
        return str(r.get_tag("HP"))
    except KeyError:
        return None


def allele_balance_hp(bam, *, genome_ref=None, mut_ctract=8, min_tot=8, max_span=12000) -> dict:
    """Mode HP-TAGS (données haplotaggées = AS) : balance allélique FIABLE par tag whatshap, sans
    deviner la longueur (le proxy span est distordu sur GRCh38 collapsé). Rapporte, PAR HP, le nombre
    de reads spanning + le taux 8C dé-poolé (HP-isolé) → l'allèle **muté s'auto-identifie** par son 8C élevé."""
    import pysam
    from .detectors.vntr import _qpos_at_ref, _VNTR_ANCHORS
    from .detectors.vntr_dupc import _CTX, gene_oriented
    a = _VNTR_ANCHORS
    chrom, l_anchor, r_anchor = a["chrom"], a["l_anchor"], a["r_anchor"]
    mode = "rc" if str(bam).endswith(".cram") else "rb"
    kw = {"reference_filename": genome_ref} if (mode == "rc" and genome_ref) else {}
    counts = collections.Counter()
    idx8c = {"1": collections.defaultdict(collections.Counter),
             "2": collections.defaultdict(collections.Counter)}
    n_span = 0
    with pysam.AlignmentFile(bam, mode, **kw) as af:
        for r in af.fetch(chrom, l_anchor, r_anchor):
            if r.is_secondary or r.is_supplementary or r.is_unmapped or r.query_sequence is None:
                continue
            hp = _read_hp(r)
            qL, qR = _qpos_at_ref(r, l_anchor), _qpos_at_ref(r, r_anchor)
            if qL is None or qR is None or not (0 < qR - qL <= max_span):
                continue
            matches = list(_CTX.finditer(gene_oriented(r.query_sequence[qL:qR])))
            if not matches:
                continue
            n_span += 1
            counts[hp if hp in ("1", "2") else "untagged"] += 1
            if hp in ("1", "2"):
                for i, m in enumerate(matches):
                    idx8c[hp][i][len(m.group(1))] += 1
    n1, n2 = counts.get("1", 0), counts.get("2", 0)
    tot = n1 + n2
    return {
        "mode": "hp_tag", "n_span": n_span,
        "by_hp": {"HP1": n1, "HP2": n2, "untagged": counts.get("untagged", 0),
                  "frac_hp1": round(n1 / tot, 3) if tot else None},
        "hp1_8C": best_8c_index(idx8c["1"], mut_ctract, min_tot),   # l'allèle muté = HP au 8C élevé
        "hp2_8C": best_8c_index(idx8c["2"], mut_ctract, min_tot),
    }


def _lengths_from_json(path):
    """(mut_len_array, healthy_len_array) depuis un *.score.json (comme matrix_score.from_phased)."""
    d = json.load(open(path))
    if isinstance(d, dict):
        d = d.get("phased", d)
    lm, lh, hp = d.get("vntr_len_mut"), d.get("vntr_len_healthy"), d.get("hp_mut")
    if lm is None and hp in ("1", "2"):
        other = "2" if hp == "1" else "1"
        lm, lh = d.get(f"vntr_len_hp{hp}"), d.get(f"vntr_len_hp{other}")
    return lm, lh


def main(argv=None):
    ap = argparse.ArgumentParser(prog="muc1_analyzer.allele_balance",
                                 description="Balance allélique VNTR (longueur) + rs4072037 — PCR vs AS")
    ap.add_argument("-b", "--bam", required=True)
    ap.add_argument("--mode", choices=["length", "hp"], default="length",
                    help="'hp' = balance par tags whatshap (FIABLE, données haplotaggées/AS) ; "
                         "'length' = proxy span (⚠ distordu sur GRCh38 collapsé, PCR non taggée)")
    ap.add_argument("--genome-ref", default=None, help="requis pour un CRAM")
    ap.add_argument("--from-json", default=None, help="*.score.json → longueurs muté/sain")
    ap.add_argument("--mut-len-array", type=int, default=None, help="longueur (répétitions) allèle muté")
    ap.add_argument("--healthy-len-array", type=int, default=None)
    ap.add_argument("--snp-pos", type=int, default=None, help="défaut rs4072037 (config)")
    ap.add_argument("--mut-ctract", type=int, default=8, help="8 = dupC (59dupC), 5 = delCC")
    ap.add_argument("--min-gap", type=int, default=5, help="écart min de longueur pour séparer par array")
    a = ap.parse_args(argv)
    if a.mode == "hp":
        print(json.dumps(allele_balance_hp(a.bam, genome_ref=a.genome_ref, mut_ctract=a.mut_ctract), indent=2))
        return 0
    lm, lh = (a.mut_len_array, a.healthy_len_array)
    if a.from_json:
        jm, jh = _lengths_from_json(a.from_json)
        lm, lh = (lm if lm is not None else jm), (lh if lh is not None else jh)
    if lm is None or lh is None:
        ap.error("fournir --from-json OU (--mut-len-array --healthy-len-array)")
    out = allele_balance(a.bam, mut_len_array=lm, healthy_len_array=lh, genome_ref=a.genome_ref,
                         snp_pos=a.snp_pos, mut_ctract=a.mut_ctract, min_gap=a.min_gap)
    print(json.dumps(out, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

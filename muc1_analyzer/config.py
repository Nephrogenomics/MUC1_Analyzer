"""Configuration du module MUC1_Score : coordonnées génomiques, seuils, pondérations.

Toutes les coordonnées GRCh38 proviennent des scripts 01/03/04 du cluster (validées).
Les coordonnées T2T-CHM13 restent à renseigner (référence disponible sur le cluster).
Les pondérations du score sont des HYPOTHÈSES par défaut, à calibrer en Phase 3.
"""
from __future__ import annotations
from dataclasses import dataclass, field


@dataclass(frozen=True)
class Locus:
    """Un intervalle génomique (0/1-based selon usage documenté au point d'emploi)."""
    chrom: str
    start: int
    end: int
    note: str = ""


# ── Coordonnées GRCh38 (« chr » naming, réf /path/to/GRCh38_chr.fa) ──
GRCh38 = {
    # DEL régulatrice SV_8533 (~365 pb, MAF~0.49, cis-eQTL ↑ expression MUC1)
    "DEL_SV8533": Locus("chr1", 155_188_853, 155_189_217, "DEL ~365bp, intron6-exon7"),
    # INS SV_8535 (~55 pb, intron2)
    "INS_SV8535": Locus("chr1", 155_189_147, 155_189_147, "INS ~55bp, intron2"),
    # Régions fonctionnelles MUC1 (gène sur brin −, coordonnées décroissantes en 5'→3')
    "PROMOTER":   Locus("chr1", 155_193_000, 155_195_000, "promoteur"),
    "UTR5":       Locus("chr1", 155_192_300, 155_192_915, "5'UTR (hypométhylé en cis par la DEL)"),
    "UTR5_CPG":   Locus("chr1", 155_192_448, 155_192_913, "fenêtre CpG 5'UTR pour l'ASM"),
    "CDS":        Locus("chr1", 155_186_500, 155_192_300, "CDS (contient le VNTR)"),
    "UTR3":       Locus("chr1", 155_185_824, 155_186_500, "3'UTR"),
    "LOCUS":      Locus("chr1", 155_185_000, 155_192_000, "fenêtre locus pour fetch"),
    # SNP splice exon 2 : rs4072037, A/G, exon2 nt8, splice-régulateur fonctionnel.
    # ~3.06 kb du DEL SV_8533 ; ~170 pb de la fenêtre CpG 5'UTR -> LD/collinéarité à démêler.
    "SPLICE_SNP_EXON2": Locus("chr1", 155_192_276, 155_192_276, "rs4072037 A/G exon2 nt8"),
}

# rs4072037 : A/G est l'écriture BRIN-GÈNE (MUC1 sur brin −). Sur le brin FORWARD
# de la référence (ce que renvoie pysam), les allèles sont le complément = T/C.
# Le détecteur est strand-agnostique : REF lu depuis la FASTA, ALT = allèle ségrégeant.
SNP_RS4072037 = {"rsid": "rs4072037", "gene_strand": "A/G", "fwd_alleles": ("C", "T")}

# Région VNTR en T2T-CHM13 (réf du README MUC1_Analyzer)
T2T = {
    "VNTR": Locus("NC_060925.1", 154_324_904, 154_335_104, "VNTR MUC1 (Bam_cleaner default)"),
    # rs4072037 in CHM13v2.0 (dbSNP maps hg38/hg37 only). Derived by EXACT anchor match: the 41 bp
    # detectors/splice_snp.RS4072037_VNTR_ANCHOR (unique 5' flank, SNP = middle base) searched in
    # chm13v2.0.fa → chr1:154,330,839, forward ref base C (anchor matches the − strand). 278 bp 5' of the
    # tandem top (154,330,561); ~195 bp above the arbiter's upper flank anchor (aH end 154,330,644).
    # Contig also seen as NC_060925.1 / CP068277.2 depending on the fasta (chm13v2.0.fa uses `chr1`).
    "SPLICE_SNP_EXON2": Locus("chr1", 154_330_839, 154_330_839, "rs4072037 A/G exon2 nt8 (CHM13v2.0)"),
    # ⚠ TODO: T2T equivalent of DEL_SV8533 if a CRAM is ever aligned to T2T (not needed for the score).
}

# ── Seuils de détection ────────────────────────────────────────────────────────
MIN_MQ = 20
MIN_BQ = 20
MIN_SPAN_PAD = 150          # bp d'ancrage requis de part et d'autre d'un SV (génotypage spanning-read)
DEL_SIZE = 364              # GRCh38 end-start ; modkit/paper ~365 (inclusif)
DEL_FRAC_HET = 0.20         # VAF min pour appeler un porteur (0/1)
DEL_FRAC_HOM = 0.80         # VAF min pour 1/1

# ── Alignement LR-PCR sur GRCh38 chr1 (DEL SV8533 / rs4072037 / dupC flanc) ─────
# Params minimap2 FIGÉS après un sweep de paramètres (verdict 2026-07-08,
# cf. MUC1_log.md) : le jeu B garde le grand indel (DEL ~365 pb + gros indels VNTR) DANS
# UNE SEULE alignement (CIGAR D unique → genotype_del le voit) sans splitter en
# supplementaires. C (`-g 20000`) n'ajoute rien ; D (`-w 5`) DÉGRADE. Le défaut nu
# `-ax map-ont` splitte les grands indels → sous-appel. NB : ceci concerne l'espace
# chr1 (flancs) ; l'espace VNTR multi-contig (longueur/copies) N'utilise PAS ces flags.
# ⚠ La PCR ne génotype PAS la DEL de façon fiable (VNTR-collapse GRCh38, cf. journal) →
# ces flags servent surtout rs4072037 + dupC ; la DEL reste appelée sur l'AS/WGS.
# ── VNTR length: flank → physical unit conversion ─────────────────────────────
# The alignment-free arbiter has two anchor schemes and they do NOT count the same units:
#   · FLANK (AL→AH): the two unique anchors sit ~4 units INSIDE what the array actually spans, so this
#     scheme systematically reports 4 units SHORT. Verified on the reference: every contig measures
#     N − 4.33 (44→39.67, 70→65.67, 150→145.67).
#   · CASSETTE (motifs 1→9): anchored on the true first/last units → reports the physical count exactly.
# So `cassette = flank + 4`, as a matter of geometry, on the reference AND on real reads.
# This is a property of the MEASUREMENT METHOD, not of the assay: the offset must be applied to every
# flank measurement (PCR, adaptive sampling, WGS alike). Treating it as a PCR-only quirk left AS/WGS
# lengths 4 short and made AS and PCR incomparable. Calibrated against 7 clinical ground truths
# (14 alleles): flank+4 = 1.07 copies mean absolute error vs flank+0 = 4.21 (a constant −4 bias).
FLANK_TO_PHYSICAL_OFFSET = 4

# Strand balance: below this minority-strand fraction the library/PCR is skewed enough that a variant
# carried by the weak strand becomes under-powered — even at amplicon depth — because our variant gates
# REQUIRE support on both strands (`pcr_dupc.strand_ok`, `two_axis._coverage_gate`). 15 % is deliberately
# generous: at 15 % of 6000 capped reads the weak strand still holds ~900 reads, so warning below it flags
# genuinely lopsided libraries without crying wolf on ordinary jitter.
STRAND_SKEW_WARN = 0.15
# MEASURED on our 23 LR-PCR amplicons: minority strand 35.5 % ± 2.9 (range 26-41), i.e. the assay carries a
# SYSTEMATIC but moderate ~65/35 skew — reproducible, not sample noise, and worth stating in the methods.
# Consequence: the absolute floor above almost never fires here; the informative signal is the DEVIATION from
# this norm (the one sample at 26 % is z=-3.3, and is also the one with the thinnest long allele — strand skew
# tracks overall library quality). Also reassuring for the --pcr cap: at 36 % minority, capping to 6000 still
# leaves ~2160 weak-strand reads, so the cap does not endanger bi-strand calling on this chemistry.
STRAND_MINOR_ASSAY_NORM = 0.355

PCR_MINIMAP_PRESET = "map-ont"
PCR_MINIMAP_EXTRA = ["-z", "600,200", "-r", "2000,20000"]   # jeu B (figé)
PCR_MINIMAP_FLAGS = ["-ax", PCR_MINIMAP_PRESET, "--secondary=no", *PCR_MINIMAP_EXTRA]

# `prepare` auto-detects an LR-PCR amplicon by read count: WGS/AS at the MUC1 locus is a few HUNDRED
# reads, LR-PCR is THOUSANDS+ (amplicon depth). At/above this many reads, prepare assumes amplicon and
# applies the frozen chaining flags (with a warning) unless --pcr / --no-pcr-autodetect is given.
PCR_AUTODETECT_MIN_READS = 2000

# A DEEP amplicon (tens of thousands of reads) vs the ~150 near-identical VNTR contigs makes minimap2 stall
# (chaining × ~150 contigs per read). In --pcr mode `prepare` subsamples the reads to this cap BEFORE the
# alignment — a random sample that PRESERVES the allele ratio, so the PCR-depleted long allele survives at a
# representative fraction while the aligner stays bounded. Tunable via --pcr-cap. Only amplicons are capped.
PCR_ALIGN_CAP = 6000


# ── Garde focal de l'appel de mutation VNTR ────────────────────────────────────
# Un vrai frameshift ADTKD touche UNE unité. Si la même Δ récurre sur une grande fraction
# de l'array = artefact de décomposition (motif de réf mal aligné à basse couverture), pas
# une mutation. (Ici, pas de dépendance lourde -> importable par refocus_vntr sans pysam.)
ARTIFACT_MIN_UNITS = 3      # jamais considéré array-wide en-dessous (tolérance bruit)
FOCAL_MAX_FRACTION = 0.15   # au-delà de cette fraction d'unités frameshiftées -> artefact


# ── Pondérations du MUC1_Score (v0 — À CALIBRER en Phase 3) ────────────────────
@dataclass(frozen=True)
class ScoreWeights:
    w_mut:  float = 40.0     # présence de la mutation dupC (déterminant principal)
    w_del:  float = 0.0      # DÉPRÉCIÉ : la « DEL » (allèle court / ex-SV8533) est MÉCANISTIQUE,
                             # non pronostique (2026-07-09) → poids 0 ; non utilisé par compute_score
    w_snp:  float = 0.0      # DÉPRÉCIÉ : rs4072037 est MÉCANISTIQUE, non pronostique (rho −0.00,
                             # 2026-07-07) → poids 0 ; non utilisé par compute_score
    w_vntr: float = 20.0     # composante longueur/ratio VNTR
    # total nominal = 100 ; le score final est clampé [0,100]


DEFAULT_WEIGHTS = ScoreWeights()


# ── Calibration de la composante pronostique VNTR COMBINÉE (ratio + position mutation) ──────
# Modèle `onset ~ intercept + w_frac·z(frac_to_mut) + w_ratio·z(ratio)` ajusté sur le survey n=27
# (2026-07-09). frac_to_mut = position mutation / longueur allèle muté ; ratio = muté/sain.
# CV leave-one-out : combiné rho +0.515 (p 0.006) vs frac seul 0.28 (ns) / ratio seul 0.31 (ns).
# ⚠ PROVISOIRE (n=27, in-sample rho 0.59 → LOO 0.52) — à recalibrer si la cohorte grandit.
SEVERITY_CAL = {
    "mean_frac": 0.4625, "std_frac": 0.2061,     # frac_to_mut (position/total muté)
    "mean_ratio": 1.1868, "std_ratio": 0.4904,   # ratio longueur muté/sain
    "w_frac": 5.8277, "w_ratio": -6.6143,        # poids OLS sur z-scores (prédiction d'onset)
    "dev_std": 8.2458,                            # écart-type de l'onset_dev (pour le squash logistique)
}


# ── Sélection de CRAM par LISTE (léger, pysam-free) ────────────────────────────
def read_cram_list(path: str) -> list:
    """Liste EXPLICITE de CRAM (contourne le glob quand les chemins sont hétérogènes : `whatshap`
    vs `02_whatshap`, dossiers noyés dans un fourre-tout). Accepte soit un TSV avec un en-tête
    contenant une colonne `cram` (ex. `lists/cohort_labels.tsv`) → extrait cette colonne, soit un
    fichier brut = un chemin par ligne (`#`/vides ignorés). Partagé (batch_score_phased, del_snp_matrix,
    12b ASM) ; pas de dépendance lourde."""
    with open(path) as fh:
        lines = [ln.rstrip("\n") for ln in fh]
    if not lines:
        return []
    header = lines[0].split("\t")
    if "cram" in header:                                   # TSV type cohort_labels.tsv
        ci = header.index("cram")
        out = []
        for ln in lines[1:]:
            cols = ln.split("\t")
            if len(cols) > ci and cols[ci].strip():
                out.append(cols[ci].strip())
        return out
    return [ln.strip() for ln in lines if ln.strip() and not ln.lstrip().startswith("#")]

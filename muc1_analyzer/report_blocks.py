#!/usr/bin/env python3
"""report_blocks — the report blocks shared by every renderer, so there is ONE clinical PDF.

Two writers had grown in parallel: `caller.write_pdf_report` (wired into `run`, carrying the arbiter
length, the T2T-native rs4072037, the length warning and the LAYER CONGRUENCE block) and
`MUC1_Analyzer_fromfastq.write_merged_pdf` (carrying the blue clinical-call box, the variant motif in red,
the power-aware negative wording and the cohort onset plot). A clinician could be handed two
different-looking reports for the same sample, and a field fixed on one side stayed wrong on the other.

These are the blocks the second writer had and the first lacked, extracted so both build from one source.
Everything here is presentation only — no calling, no thresholds, no decisions.
"""
from __future__ import annotations

COL_VARIANT = "#C00000"      # the pathogenic motif, and only it, is red
COL_BLUE = "#185FA5"
COL_BOX_BG = "#E6F1FB"


def cohort_ref():
    """(mean_frac, std_frac, n) of the onset calibration cohort. Pure-ish (reads config, never raises)."""
    try:
        from .config import SEVERITY_CAL as S
        return float(S["mean_frac"]), float(S["std_frac"]), 27
    except Exception:
        return 0.4625, 0.2061, 27


def onset_z(onset_index):
    """(mean_oi, z) for an onset_index against the calibration cohort, or (None, None). Pure."""
    if onset_index is None:
        return None, None
    mean_frac, std_frac, _ = cohort_ref()
    mean_oi = 1.0 - mean_frac
    return mean_oi, (onset_index - mean_oi) / std_frac


def dupc_verdict_html(rel_dupc) -> str:
    """The dupC verdict as one report line, POWER-AWARE. Pure; '' when there is nothing to say.

    A bare "not called" is not a clinical statement: the reader cannot tell an adequately-powered
    negative from a sample too thin to exclude anything. Ilias's four-way wording, kept verbatim in
    substance — CONFIRMED / NEGATIVE / NEGATIVE (borderline power) / INDETERMINATE — with the covering
    depth per allele and, when underpowered, the reads that would be needed."""
    if not rel_dupc:
        return ""
    if not rel_dupc.get("assessed", True):
        return "dupC (reliable verdict): <b>NOT ASSESSED</b> — " + str(rel_dupc.get("reason", ""))
    pw = rel_dupc.get("power_by_allele") or {}
    cov = ", ".join(f"{c}: {t.get('n_cover')} covering reads"
                    f"{' (powered)' if t.get('adequately_powered') else ' (underpowered)'}"
                    for c, t in pw.items())
    interp = rel_dupc.get("interpretation")
    if interp == "CONFIRMED":
        meth = {"clair3": "Clair3",
                "runlen_shift": "homopolymer run-length shift, calibrated specificity"}.get(
                    rel_dupc.get("method"), "context+paired, calibrated specificity")
        var = rel_dupc.get("variant") or {}
        extra = f", repeat {var.get('repeat')}" if var.get("repeat") is not None else ""
        return f"<b>dupC: CONFIRMED</b> on {rel_dupc.get('carrier_contig')}{extra} ({meth})"
    if interp == "NEGATIVE":
        return f"<b>dupC: NEGATIVE (reliable)</b> — not detected, with adequate power [{cov}]"
    if interp == "NEGATIVE_BORDERLINE":
        return (f"<b>dupC: NEGATIVE (borderline power)</b> — not detected; coverage at the edge of "
                f"detectability [{cov}]: a low-d carrier could be missed")
    if interp:
        need = max((t.get("reads_for_power0.90", {}).get("d=0.3", 12) for t in pw.values()), default=12)
        return (f"<b>dupC: INDETERMINATE</b> — not detected but coverage underpowered, a dupC cannot be "
                f"excluded [{cov}]; ~{need} reads/allele required (d=0.3)")
    return ""


def nomenclature_html(motifs, *, is_carrier=False, variant_repeat=None, variant_label=None,
                      indel_colour="#c0392b") -> str:
    """Motif nomenclature with the PATHOGENIC unit in red. Pure.

    The variant repeat is a 1-based motif index, so motif #variant_repeat is the one recoloured. The
    consensus may already name the motif with its frameshift (e.g. "X-59dupC"); the label is appended only
    when it is not already there, so the variant is never written twice."""
    if not motifs:
        return ""
    parts = []
    for idx, m in enumerate(motifs, 1):
        name = m.get("name", "?")
        if is_carrier and variant_repeat is not None and idx == variant_repeat:
            disp = name if (variant_label and variant_label in name) else f"{name}({variant_label})"
            parts.append(f'<font color="{COL_VARIANT}"><b>{disp}</b></font>')
        elif int(m.get("length", 60)) != 60:
            parts.append(f'<font color="{indel_colour}">{name}</font>')
        else:
            parts.append(name)
    return "-".join(parts)


# ── reportlab flowables (imported lazily so the pure helpers stay importable without it) ──

def clinical_call_box(call_str):
    """The clinical call as a boxed headline, not a footer — it is the result the reader came for."""
    from reportlab.lib import colors
    from reportlab.lib.styles import ParagraphStyle
    from reportlab.lib.units import cm
    from reportlab.platypus import Paragraph, Table, TableStyle
    style = ParagraphStyle("clinbox", fontName="Helvetica-Bold", fontSize=11,
                           textColor=colors.HexColor("#0d2b45"), leading=14)
    tbl = Table([[Paragraph(f"VNTR — clinical call: {call_str}", style)]], colWidths=[17 * cm])
    tbl.setStyle(TableStyle([("BACKGROUND", (0, 0), (-1, -1), colors.HexColor(COL_BOX_BG)),
                             ("BOX", (0, 0), (-1, -1), 1.2, colors.HexColor(COL_BLUE)),
                             ("LEFTPADDING", (0, 0), (-1, -1), 10), ("RIGHTPADDING", (0, 0), (-1, -1), 10),
                             ("TOPPADDING", (0, 0), (-1, -1), 9), ("BOTTOMPADDING", (0, 0), (-1, -1), 9)]))
    return tbl


def cohort_plot(onset_index):
    """Where this sample's onset_index sits in the calibration cohort, or None when there is no position.

    onset_index = the fraction of the mutant VNTR downstream of the frameshift, i.e. what is translated as
    the MUC1fs neoprotein. It predicts onset better than the absolute length, which is why it gets the
    figure rather than the length."""
    if onset_index is None:
        return None
    import math

    from reportlab.graphics.shapes import Circle, Drawing, Line, PolyLine, String
    from reportlab.lib import colors

    mean_oi, z = onset_z(onset_index)
    _, std_oi, n = cohort_ref()
    W, Hd, lx, rx, by, ty = 470, 150, 45, 440, 48, 126
    d = Drawing(W, Hd)
    X = lambda oi: lx + max(0.0, min(1.0, oi)) * (rx - lx)                       # noqa: E731
    G = lambda oi: math.exp(-0.5 * ((oi - mean_oi) / std_oi) ** 2)               # noqa: E731
    Y = lambda g: by + g * (ty - by)                                            # noqa: E731
    pts = []
    for i in range(81):
        oi = i / 80.0
        pts += [X(oi), Y(G(oi))]
    d.add(PolyLine(pts, strokeColor=colors.HexColor(COL_BLUE), strokeWidth=1.8))
    d.add(Line(lx, by, rx, by, strokeColor=colors.HexColor("#888"), strokeWidth=1))
    for oi, lab in [(0.0, "0"), (mean_oi, "cohort mean"), (1.0, "1")]:
        x = X(oi)
        d.add(Line(x, by, x, by - 4, strokeColor=colors.HexColor("#888")))
        d.add(String(x, by - 14, lab, fontName="Helvetica", fontSize=7,
                     fillColor=colors.HexColor("#666"), textAnchor="middle"))
    # direction of effect: onset_index up (right) = earlier frameshift -> more neoprotein -> earlier onset
    d.add(String(lx, by - 26, "◄ later onset", fontName="Helvetica-Oblique", fontSize=7,
                 fillColor=colors.HexColor("#888"), textAnchor="start"))
    d.add(String(rx, by - 26, "earlier onset ►", fontName="Helvetica-Oblique", fontSize=7,
                 fillColor=colors.HexColor("#888"), textAnchor="end"))
    xm = X(mean_oi)
    d.add(Line(xm, by, xm, Y(G(mean_oi)), strokeColor=colors.HexColor("#9BB6D6"),
               strokeWidth=0.8, strokeDashArray=[2, 2]))
    xs, ys = X(onset_index), Y(G(onset_index))
    d.add(Line(xs, by, xs, ys, strokeColor=colors.red, strokeWidth=0.8, strokeDashArray=[2, 2]))
    d.add(Circle(xs, ys, 4, fillColor=colors.red, strokeColor=colors.red))
    d.add(String(xs, ys + 8, f"sample (z = {z:+.2f})", fontName="Helvetica-Bold", fontSize=7.5,
                 fillColor=colors.red, textAnchor="middle"))
    d.add(String((lx + rx) / 2, Hd - 9, f"Onset position in the calibration cohort (n={n})",
                 fontName="Helvetica", fontSize=8, fillColor=colors.HexColor("#444"), textAnchor="middle"))
    d.add(String((lx + rx) / 2, 6, "onset_index = fraction of the mutant VNTR translated into neoprotein",
                 fontName="Helvetica-Oblique", fontSize=6.5, fillColor=colors.HexColor("#999"),
                 textAnchor="middle"))
    return d

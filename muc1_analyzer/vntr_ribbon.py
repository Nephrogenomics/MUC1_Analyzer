#!/usr/bin/env python3
"""vntr_ribbon.py — colour-ribbon rendering of a VNTR motif array from a MUC1_Analyzer analyzer JSON.
 
Each 60 bp motif is one coloured block (same alphabet as the caller: 1-9, A-K, V, W, X, ?). Haplotypes are
stacked and index-aligned so the two alleles compare at a glance. The variant-carrying motif (an indel:
name contains dup/del/ins) is HIGHLIGHTED (bold outline + variant label above, and its MOTIF NUMBER — the
1-based array index, i.e. the "repeat N" of the clinical nomenclature — below the block) — the visual
VNTRtools' colorcode does not give. By default the highlight follows the CLINICAL rule: only when the frameshift verdict is POSITIVE
(CALLED / POSITIVE_LOW_DEPTH); ``show_all=True`` highlights any reconstructed indel regardless of verdict.
 
Used two ways:
  * ``render(analyzer_json, out_png)`` — called by ``run_pcr`` and ``MUC1_Analyzer_fromfastq`` to emit the
    reconstruction figure as a DEFAULT output (``{sample}.ribbon.png``). ``verdict=`` overrides the JSON's
    own dupc verdict, for a merged JSON that carries the haplotypes but not the dupc block.
  * CLI: ``python -m muc1_analyzer.vntr_ribbon --analyzer-json SAMPLE.analyzer.json [-o out.png] [--show-all]``
"""
import argparse
import json
import os
 
NUMBERS = list("123456789")
LETTERS = list("ABCDEFGHIJK") + ["V", "W", "X", "?"]
POSITIVE_VERDICTS = {"CALLED", "POSITIVE_LOW_DEPTH"}
UNKNOWN = (0.72, 0.72, 0.72, 1.0)
 
 
def _palette():
    from matplotlib.colors import LinearSegmentedColormap
    warm = LinearSegmentedColormap.from_list("warm", ["orange", "red"])
    cool = LinearSegmentedColormap.from_list("cool", ["yellow", "blue", "purple", "green"])
    col = {c: warm(i / max(1, len(NUMBERS) - 1)) for i, c in enumerate(NUMBERS)}
    col.update({c: cool(i / max(1, len(LETTERS) - 1)) for i, c in enumerate(LETTERS)})
    return col
 
 
def _family(name, palette):
    if not name:
        return "?"
    ch = name[0]
    if ch in palette:
        return ch
    if ch.upper() in palette:
        return ch.upper()
    return "?"
 
 
def _is_indel(m):
    n = m.get("name", "") or ""
    return any(k in n for k in ("dup", "del", "ins")) or m.get("length", 60) != 60
 
 
def _suffix(name):
    return name.split("-", 1)[1] if "-" in (name or "") else name
 
 
def _verdict_of(d):
    dp = d.get("dupc") or {}
    res = dp.get("result", dp) if isinstance(dp, dict) else {}
    return res.get("verdict") if isinstance(res, dict) else None
 
 
def render(analyzer_json, out_path=None, *, verdict=None, show_all=False, legend=True,
           block=0.34, gap=0.03, row_gap=0.7):
    """Render the motif ribbon. ``analyzer_json`` is a path or an already-loaded dict. Returns the output
    path, or None when there is nothing to draw (no haplotypes/motifs) — callers treat that as a skip."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle
 
    d = analyzer_json if isinstance(analyzer_json, dict) else json.load(open(analyzer_json))
    src = "" if isinstance(analyzer_json, dict) else os.path.basename(analyzer_json)
    sample = d.get("sample") or (src.split(".")[0] if src else "sample")
    palette = _palette()
    if verdict is None:
        verdict = _verdict_of(d)
    positive = verdict in POSITIVE_VERDICTS
    highlight = show_all or positive
    haps = [h for h in d.get("haplotypes", []) if h.get("motifs")]
    if not haps:
        return None
    nmax = max(len(h["motifs"]) for h in haps)
    fig_w = max(6.0, nmax * (block + gap) * 0.55 + 3.5)
    fig_h = 1.7 * len(haps) + 1.2
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    row_h = 1.0
    used = set()
    any_var = False
    for row, h in enumerate(haps):
        y = -row * (row_h + row_gap)
        ms = h["motifs"]
        for i, m in enumerate(ms):
            x = i * (block + gap)
            fam = _family(m.get("name"), palette)
            used.add(fam)
            color = palette.get(fam, UNKNOWN)
            var = highlight and _is_indel(m)
            any_var = any_var or var
            ax.add_patch(Rectangle((x, y), block, row_h, facecolor=color,
                                   edgecolor="black" if var else "white",
                                   linewidth=2.6 if var else 0.4, zorder=3 if var else 1))
            if var:
                ax.plot([x + block / 2], [y + row_h + 0.05], marker="v", color="black", ms=6, zorder=4)
                ax.annotate(_suffix(m.get("name", "")), (x + block / 2, y + row_h + 0.22),
                            ha="center", va="bottom", fontsize=8.5, fontweight="bold", zorder=4)
                # motif NUMBER (1-based index in the array = the "repeat N" of the clinical nomenclature)
                # printed under the block, on the always-clear side opposite the name label.
                ax.annotate(f"motif {i + 1}", (x + block / 2, y - 0.06),
                            ha="center", va="top", fontsize=8, fontweight="bold",
                            color="black", zorder=4)
        contig = h.get("contig", "")
        ax.text(-0.35, y + row_h / 2, f"{contig}\n({len(ms)} motifs)",
                ha="right", va="center", fontsize=9)
    if legend:
        from matplotlib.patches import Patch
        _key = lambda c: (0, int(c)) if c.isdigit() else (1, c)
        handles = [Patch(facecolor=palette.get(c, UNKNOWN), edgecolor="0.5", label=c)
                   for c in sorted(used, key=_key)]
        if any_var:
            handles.append(Patch(facecolor="white", edgecolor="black", linewidth=2.2,
                                 label="variant (indel)"))
        ncol = 1 if len(handles) <= 14 else 2
        ax.legend(handles=handles, loc="center left", bbox_to_anchor=(1.005, 0.5),
                  fontsize=8, title="motifs", title_fontsize=9, ncol=ncol,
                  handlelength=1.1, handleheight=1.1, labelspacing=0.35,
                  columnspacing=1.0, frameon=False)
    ttl = f"{sample} — VNTR motif ribbon"
    if verdict:
        ttl += f"    [frameshift: {verdict}{'' if highlight else '  — variant not shown'}]"
    ax.set_title(ttl, fontsize=11, fontweight="bold", loc="left")
    ax.set_xlim(-4.0, nmax * (block + gap) + 0.6)
    ax.set_ylim(-len(haps) * (row_h + row_gap) - 0.2, row_h + 1.0)
    ax.set_aspect("equal")
    ax.axis("off")
    out_path = out_path or f"{sample}_ribbon.png"
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return out_path
 
 
def render_safe(analyzer_json, out_path, *, verdict=None, show_all=False, legend=True):
    """Never let the ribbon sink a pipeline. Returns the path on success, None on any failure (missing
    matplotlib, empty haplotypes, …); the caller logs and carries on."""
    try:
        return render(analyzer_json, out_path, verdict=verdict, show_all=show_all, legend=legend)
    except Exception:
        return None
 
 
def main(argv=None):
    ap = argparse.ArgumentParser(prog="python -m muc1_analyzer.vntr_ribbon")
    ap.add_argument("--analyzer-json", required=True)
    ap.add_argument("-o", "--out", default=None)
    ap.add_argument("--show-all", action="store_true",
                    help="highlight any reconstructed indel even if the verdict is NEEDS_IGV/negative")
    ap.add_argument("--no-legend", action="store_true", help="hide the colour->motif legend")
    a = ap.parse_args(argv)
    p = render(a.analyzer_json, a.out, show_all=a.show_all, legend=not a.no_legend)
    print("wrote", p) if p else print("nothing to draw (no haplotypes/motifs)")
 
 
if __name__ == "__main__":
    main()

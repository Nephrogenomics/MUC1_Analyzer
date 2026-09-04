"""MUC1 C-tract SEGMENT probe — reads the inter-anchor region and classifies it by length+composition.

Generalizes the C-tract probe (`vntr_dupc._CTX`, which only captures `(C+)`): the **58_59insG inserts a G
INSIDE the C-tract** (pos 58-59, within the 7C run pos 53-59) -> it INTERRUPTS the C-run -> the regex `(C+)`
breaks and the insG is INVISIBLE. Here we capture the whole segment between the anchors `GGGCTCCACCG` and
`AGCCCAC` then classify it: length + presence of a G distinguish dupC (8C) / delCC (5C) / insG (7C+G) / dupCCCC (11C).

Reads the LITERAL segment of the read (like `vntr_dupc`, not `match_motifs`) -> no decoding artifact.
The WT = 7 C (`CCCCCCC`). Warning: SEPARATE probe: it does not replace the validated C-tract probe (dupC/delCC).
"""
from __future__ import annotations
import re

from .vntr_dupc import _rc

# variable segment between the conserved anchors; [CG] covers dupC/delCC/insG/dupCCCC (not dupA at pos60)
_SEG = re.compile(r"GGGCTCCACCG([CG]{2,20})AGCCCAC")


def classify_segment(seg: str) -> dict:
    """Classify an inter-anchor segment. WT = 7C. Returns {type, delta, len, n_c, n_g, frameshift}.

    - all-C : 7=WT, 8=dupC(+1), 5=delCC(-2), 6=delC(-1), 11=dupCCCC(+4), otherwise cIndel ;
    - 7 intact C + 1 G : **dupG** if the G is at the head (52dupG), **insG** if internal (58_59insG), delta +1 ;
      any other G = 'other'.
    `frameshift` = delta != 0 and not a multiple of 3 (shifts the reading frame -> MUC1-fs)."""
    seg = seg.upper()
    n = len(seg)
    n_c = seg.count("C")
    n_g = seg.count("G")
    delta = n - 7
    if n_g == 0:
        t = {7: "WT", 8: "dupC", 6: "delC", 5: "delCC", 11: "dupCCCC"}.get(
            n, "cIns" if n > 7 else "cDel")
    elif n_c == 7 and n_g == 1:
        # 7C + 1G : the POSITION of the G separates 52dupG (duplicated anchor G, at the HEAD of the segment)
        # from 58_59insG (G inserted INSIDE the C-run, internal). Same composition, different positions.
        t = "dupG" if seg[0] == "G" else "insG"
    else:
        t = "other"
    return {"type": t, "delta": delta, "len": n, "n_c": n_c, "n_g": n_g,
            "frameshift": delta != 0 and delta % 3 != 0}


def segment_runs(seq: str) -> list:
    """Inter-anchor segments found in `seq` (both orientations), classified."""
    return [classify_segment(m.group(1)) for m in _SEG.finditer(seq)] + \
           [classify_segment(m.group(1)) for m in _SEG.finditer(_rc(seq))]

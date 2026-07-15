#!/usr/bin/env python3
"""
HTML diagnostics for low-marker species: crossmap-masking + the ANI species gap.

Combines the outputs of low_marker_masking.py and ani_gap.py into one
dependency-free HTML page. For every species that ended below the marker
threshold it shows:
  * the masking verdict split of its guard-dropped markers (can masking rescue
    any, or do they cross-map whole-gene?);
  * an inline-SVG ANI strip plot of within-species vs between-species pairwise
    ANI in its genus -- so you can see at a glance whether within is always
    tighter than between (a clean species gap) or the two overlap (re-merge
    candidate).

Inputs are the TSVs those two tools write; output is a single self-contained
HTML file (stdlib only, no matplotlib).
"""
import argparse
import html
import math
import os
from collections import defaultdict


def esc(x):
    return html.escape(str(x))


def load_tsv(path):
    rows = []
    if not path or not os.path.exists(path):
        return rows, []
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) == len(header):
                rows.append(dict(zip(header, f)))
    return rows, header


def fnum(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return float("nan")


def ani_strip_svg(within, between_by_sp, min_within, max_between, w=620, h=None):
    """Strip plot of ANI (%) for one species: within vs between pairs.

    within: [ani...]; between_by_sp: {other_species: [ani...]}. Points are drawn
    on an ANI x-axis; a shaded band marks the gap between max(between) and
    min(within) when they don't overlap.
    """
    anis = list(within) + [a for v in between_by_sp.values() for a in v]
    if not anis:
        return ("<div class='gapinfo dim'>no ANI pairs (single-genome species or "
                "no siblings in genus)</div>")
    lo = min(80.0, math.floor(min(anis)))
    lo = max(70.0, min(lo, 95.0))
    hi = 100.0
    pad_l, pad_r, pad_t = 46, 14, 16
    rows = 1 + len(between_by_sp)
    row_h = 26
    h = pad_t + rows * row_h + 34

    def x(ani):
        return pad_l + (ani - lo) / (hi - lo) * (w - pad_l - pad_r)

    parts = [f"<svg viewBox='0 0 {w} {h}' width='100%' style='max-width:{w}px' "
             f"role='img'>"]
    # axis
    for t in range(int(lo), 101, 5 if (hi - lo) > 12 else 2):
        xt = x(t)
        parts.append(f"<line x1='{xt:.1f}' y1='{pad_t}' x2='{xt:.1f}' "
                     f"y2='{h-24:.1f}' class='grid'/>")
        parts.append(f"<text x='{xt:.1f}' y='{h-8}' class='ax' "
                     f"text-anchor='middle'>{t}</text>")
    # gap band (between_max .. within_min) if a clean gap exists
    if min_within == min_within and max_between == max_between and min_within > max_between:
        gx1, gx2 = x(max_between), x(min_within)
        parts.append(f"<rect x='{gx1:.1f}' y='{pad_t}' width='{gx2-gx1:.1f}' "
                     f"height='{rows*row_h:.1f}' class='gapband'/>")
        parts.append(f"<text x='{(gx1+gx2)/2:.1f}' y='{pad_t+10}' class='gaplbl' "
                     f"text-anchor='middle'>gap {min_within-max_between:.1f}</text>")

    def dots(vals, y, cls):
        for a in vals:
            parts.append(f"<circle cx='{x(a):.1f}' cy='{y:.1f}' r='3.4' "
                         f"class='{cls}'/>")

    y = pad_t + row_h * 0.7
    dots(within, y, "win")
    parts.append(f"<text x='4' y='{y+3:.1f}' class='rowlbl'>within</text>")
    for i, (osp, vals) in enumerate(sorted(between_by_sp.items(),
                                           key=lambda kv: -max(kv[1])), start=1):
        yy = pad_t + row_h * (i + 0.7)
        dots(vals, yy, "btw")
        short = osp.replace("s__", "")
        parts.append(f"<text x='4' y='{yy+3:.1f}' class='rowlbl'>"
                     f"{esc(short[:22])}</text>")
    parts.append("</svg>")
    return "".join(parts)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--low_marker_species", required=True)
    ap.add_argument("--masking_summary", required=True)
    ap.add_argument("--ani_pairs", required=True)
    ap.add_argument("--ani_gap_summary", required=True)
    ap.add_argument("--threshold", type=int, default=50)
    ap.add_argument("--ani_threshold", type=float, default=95.0)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    species_rows, _ = load_tsv(args.low_marker_species)
    mask_rows, _ = load_tsv(args.masking_summary)
    gap_rows, _ = load_tsv(args.ani_gap_summary)
    pair_rows, _ = load_tsv(args.ani_pairs)

    flagged = [r for r in species_rows if r.get("flagged") == "yes"]
    mask = {r["species"]: r for r in mask_rows}
    gap = {r["species"]: r for r in gap_rows}

    within = defaultdict(list)
    between = defaultdict(lambda: defaultdict(list))
    for r in pair_rows:
        sp = r["focal_species"]
        if r["kind"] == "within":
            within[sp].append(fnum(r["ani"]))
        else:
            between[sp][r["other_species"]].append(fnum(r["ani"]))

    n_flagged = len(flagged)
    n_merge = sum(1 for r in gap.values() if r.get("merge_candidate") == "yes")
    n_overlap = sum(1 for r in gap.values() if r.get("overlap") == "yes")
    tot_dropped = sum(int(mask[s]["dropped"]) for s in mask if mask[s].get("dropped"))
    tot_resc = sum(int(mask[s]["rescuable"]) for s in mask if mask[s].get("rescuable"))

    # Sort flagged: merge candidates first, then by final markers.
    def sortkey(r):
        sp = r["species"]
        g = gap.get(sp, {})
        return (0 if g.get("merge_candidate") == "yes" else 1,
                0 if g.get("overlap") == "yes" else 1,
                int(r.get("final_markers", 0)))
    flagged.sort(key=sortkey)

    blocks = []
    for r in flagged:
        sp = r["species"]
        m = mask.get(sp, {})
        g = gap.get(sp, {})
        min_w = fnum(g.get("min_within", "nan"))
        max_b = fnum(g.get("max_between", "nan"))
        merge = g.get("merge_candidate", "no") == "yes"
        overlap = g.get("overlap", "NA")
        tag = ("<span class='badge bad'>merge candidate</span>" if merge else
               "<span class='badge warn'>overlap</span>" if overlap == "yes" else
               "<span class='badge ok'>clean gap</span>" if overlap == "no" else
               "<span class='badge dim'>single genome</span>")
        svg = ani_strip_svg(within.get(sp, []), between.get(sp, {}), min_w, max_b)
        mask_line = ""
        if m.get("dropped") and int(m["dropped"]) > 0:
            mask_line = (
                f"<div class='mask'>Dropped markers: <b>{esc(m['dropped'])}</b> "
                f"&middot; rescuable by masking <b>{esc(m.get('rescuable','0'))}</b> "
                f"({esc(m.get('pct_rescuable','0'))}%) &middot; whole-gene "
                f"{esc(m.get('whole_gene','0'))} &middot; ambiguous "
                f"{esc(m.get('ambiguous','0'))} &middot; unresolved "
                f"{esc(m.get('target_unavailable','0'))}</div>")
        gap_line = (
            f"<div class='gapinfo'>within ANI min "
            f"<b>{esc(g.get('min_within','NA'))}</b> "
            f"(n={esc(g.get('n_within','0'))}) &middot; nearest "
            f"<b>{esc((g.get('nearest_species','NA') or 'NA').replace('s__',''))}</b> "
            f"between-max <b>{esc(g.get('max_between','NA'))}</b> &middot; "
            f"gap <b>{esc(g.get('gap','NA'))}</b></div>"
            if g else "<div class='gapinfo dim'>no ANI data</div>")
        blocks.append(
            f"<div class='sp'><div class='sphdr'><span class='spname'>{esc(sp)}</span>"
            f"<span class='spmeta'>{esc(r.get('genus',''))} &middot; "
            f"{esc(r.get('genomes','?'))} genomes &middot; "
            f"{esc(r.get('final_markers','?'))} final markers</span>{tag}</div>"
            f"{mask_line}{gap_line}{svg}</div>")

    body = "\n".join(blocks) or "<p class='note ok'>No species fell below the threshold.</p>"

    out = f"""<!doctype html><html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>meta-mage low-marker diagnostics</title><style>
:root {{ --bg:#0f1115; --panel:#171a21; --panel2:#1d212b; --fg:#e6e9ef;
  --dim:#9aa4b2; --line:#2a2f3a; --accent:#5aa9e6; --warn:#e0b341; --bad:#e0685a;
  --ok:#5ac48a; }}
@media (prefers-color-scheme: light) {{ :root {{ --bg:#f6f7f9; --panel:#fff;
  --panel2:#f0f2f5; --fg:#1a1d23; --dim:#5a6472; --line:#e2e5ea; --accent:#1f6fb8; }} }}
* {{ box-sizing:border-box; }}
body {{ margin:0; background:var(--bg); color:var(--fg);
  font:14px/1.5 -apple-system,BlinkMacSystemFont,"Segoe UI",Roboto,sans-serif; }}
.wrap {{ max-width:900px; margin:0 auto; padding:32px 20px 80px; }}
h1 {{ font-size:23px; margin:0 0 4px; }}
.meta {{ color:var(--dim); font-size:13px; margin-bottom:18px; }}
.cards {{ display:grid; grid-template-columns:repeat(auto-fit,minmax(150px,1fr));
  gap:12px; margin:14px 0 26px; }}
.card {{ background:var(--panel); border:1px solid var(--line); border-radius:10px;
  padding:13px 15px; }}
.card .val {{ font-size:24px; font-weight:650; }}
.card .lbl {{ color:var(--dim); font-size:12px; margin-top:2px; }}
.sp {{ background:var(--panel); border:1px solid var(--line); border-radius:10px;
  padding:14px 16px; margin-bottom:14px; }}
.sphdr {{ display:flex; align-items:baseline; gap:10px; flex-wrap:wrap;
  margin-bottom:6px; }}
.spname {{ font-weight:650; font-size:15px; }}
.spmeta {{ color:var(--dim); font-size:12px; flex:1; }}
.badge {{ font-size:11px; padding:2px 8px; border-radius:20px; font-weight:600; }}
.badge.bad {{ background:color-mix(in srgb,var(--bad) 22%,transparent); color:var(--bad); }}
.badge.warn {{ background:color-mix(in srgb,var(--warn) 22%,transparent); color:var(--warn); }}
.badge.ok {{ background:color-mix(in srgb,var(--ok) 20%,transparent); color:var(--ok); }}
.badge.dim {{ background:var(--panel2); color:var(--dim); }}
.mask, .gapinfo {{ font-size:12.5px; color:var(--fg); margin:3px 0; }}
.gapinfo.dim {{ color:var(--dim); }}
.note {{ color:var(--dim); }} .note.ok {{ color:var(--ok); }}
.legend {{ color:var(--dim); font-size:12px; margin:2px 0 20px; }}
svg {{ margin-top:8px; }}
.grid {{ stroke:var(--line); stroke-width:1; }}
.ax {{ fill:var(--dim); font-size:10px; }}
.rowlbl {{ fill:var(--dim); font-size:10px; }}
.win {{ fill:var(--accent); fill-opacity:.85; }}
.btw {{ fill:var(--bad); fill-opacity:.7; }}
.gapband {{ fill:color-mix(in srgb,var(--ok) 16%,transparent); }}
.gaplbl {{ fill:var(--ok); font-size:10px; font-weight:600; }}
</style></head><body><div class="wrap">
<h1>meta-mage &mdash; low-marker diagnostics</h1>
<div class="meta">Species with &lt; {args.threshold} final markers &middot;
merge candidate = a sibling within &ge; {args.ani_threshold:.0f}% ANI</div>
<div class="cards">
<div class="card"><div class="val">{n_flagged}</div><div class="lbl">flagged species</div></div>
<div class="card"><div class="val">{n_merge}</div><div class="lbl">merge candidates (ANI)</div></div>
<div class="card"><div class="val">{n_overlap}</div><div class="lbl">within/between overlap</div></div>
<div class="card"><div class="val">{tot_resc}/{tot_dropped}</div><div class="lbl">dropped markers rescuable by masking</div></div>
</div>
<p class="legend"><span style="color:var(--accent)">&#9679;</span> within-species pairs
&nbsp; <span style="color:var(--bad)">&#9679;</span> between-species pairs &nbsp;&middot;&nbsp;
green band = clean ANI gap. Overlap or a between-ANI &ge; {args.ani_threshold:.0f}% suggests the
two species may be one taxon.</p>
{body}
</div></body></html>"""

    with open(args.out, "w") as fh:
        fh.write(out)
    print(f"low_marker_report: wrote {args.out} ({n_flagged} flagged, "
          f"{n_merge} merge candidates, {tot_resc}/{tot_dropped} markers rescuable)")


if __name__ == "__main__":
    main()

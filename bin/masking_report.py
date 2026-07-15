#!/usr/bin/env python3
"""
Standalone HTML report on cross-map masking (from low_marker_masking.py outputs).

Answers, per low-marker species, whether the markers the cross-map guard dropped
could be salvaged by masking the offending slice instead of discarding the whole
gene. Reads masking_summary.tsv (per-species verdict counts) and
masking_markers.tsv (per-marker spans) and renders:
  * summary cards (dropped, rescuable, whole-gene, ambiguous, unresolved);
  * a per-species table of verdict splits;
  * for each rescuable / ambiguous marker, an inline-SVG track of the gene with
    the cross-mapped range(s) shaded and the clade-unique remainder highlighted,
    plus the exact 1-based mask range(s).

Verdicts (set upstream): rescuable = a contiguous clade-unique window >=
min_survivor remains after masking; whole_gene = the gene cross-maps end-to-end;
ambiguous = cross-mapped but the safe window is too short; target_unavailable =
off-target hits were all to non-rep genomes (couldn't localise via reps.ffn).

stdlib only; no matplotlib.
"""
import argparse
import html
import os

VERDICT_ORDER = ["rescuable", "ambiguous", "whole_gene", "target_unavailable"]
VERDICT_LABEL = {
    "rescuable": "rescuable", "ambiguous": "ambiguous",
    "whole_gene": "whole-gene", "target_unavailable": "unresolved",
}


def esc(x):
    return html.escape(str(x))


def load_tsv(path):
    rows = []
    if not path or not os.path.exists(path):
        return rows
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) == len(header):
                rows.append(dict(zip(header, f)))
    return rows


def coverage_svg(gene_len, mask_ranges, w=560, h=22):
    """Green gene bar with red cross-mapped range(s) overlaid (1-based ranges)."""
    if gene_len <= 0:
        return ""
    parts = [f"<svg viewBox='0 0 {w} {h}' width='100%' style='max-width:{w}px'>"]
    parts.append(f"<rect x='0' y='4' width='{w}' height='{h-8}' class='safe'/>")
    if mask_ranges and mask_ranges != "-":
        for rng in mask_ranges.split(";"):
            try:
                s, e = (int(v) for v in rng.split("-"))
            except ValueError:
                continue
            x1 = (s - 1) / gene_len * w
            x2 = e / gene_len * w
            parts.append(f"<rect x='{x1:.1f}' y='4' width='{max(1.0, x2-x1):.1f}' "
                         f"height='{h-8}' class='xmap'/>")
    parts.append("</svg>")
    return "".join(parts)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--masking_summary", required=True)
    ap.add_argument("--masking_markers", required=True)
    ap.add_argument("--min_survivor", type=int, default=200)
    ap.add_argument("--max_markers_shown", type=int, default=40)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    summary = load_tsv(args.masking_summary)
    markers = load_tsv(args.masking_markers)

    def gi(d, k):
        try:
            return int(d.get(k, 0) or 0)
        except ValueError:
            return 0

    tot_dropped = sum(gi(r, "dropped") for r in summary)
    tot_resc = sum(gi(r, "rescuable") for r in summary)
    tot_whole = sum(gi(r, "whole_gene") for r in summary)
    tot_amb = sum(gi(r, "ambiguous") for r in summary)
    tot_unres = sum(gi(r, "target_unavailable") for r in summary)
    n_species = sum(1 for r in summary if gi(r, "dropped") > 0)

    # Per-species rows, worst (most rescuable) first.
    summary = [r for r in summary if gi(r, "dropped") > 0]
    summary.sort(key=lambda r: (-gi(r, "rescuable"), -gi(r, "dropped")))
    srows = []
    for r in summary:
        d = gi(r, "dropped")
        srows.append(
            f"<tr><td class='l'>{esc(r['species'])}</td>"
            f"<td>{esc(r.get('genomes','?'))}</td>"
            f"<td>{esc(r.get('final_markers','?'))}</td><td>{d}</td>"
            f"<td class='ok'>{gi(r,'rescuable')}</td>"
            f"<td>{gi(r,'ambiguous')}</td><td>{gi(r,'whole_gene')}</td>"
            f"<td class='dim'>{gi(r,'target_unavailable')}</td>"
            f"<td class='strong'>{esc(r.get('pct_rescuable','0'))}%</td></tr>")

    # Per-marker detail: show rescuable then ambiguous (the actionable ones).
    shown = [m for m in markers if m.get("verdict") in ("rescuable", "ambiguous")]
    shown.sort(key=lambda m: (m["verdict"] != "rescuable",
                              -int(m.get("longest_safe_window", 0) or 0)))
    blocks = []
    for m in shown[: args.max_markers_shown]:
        glen = int(m.get("gene_len", 0) or 0)
        v = m["verdict"]
        blocks.append(
            f"<div class='mk'><div class='mkhdr'>"
            f"<span class='mkcl'>{esc(m['species'])} &middot; {esc(m['cluster'])}</span>"
            f"<span class='badge {'ok' if v=='rescuable' else 'warn'}'>{esc(v)}</span>"
            f"</div><div class='mkmeta'>gene {glen} bp &middot; cross-mapped "
            f"{esc(m.get('covered_frac','?'))} &middot; safe window "
            f"<b>{esc(m.get('longest_safe_window','?'))} bp</b> &middot; mask "
            f"<code>{esc(m.get('mask_ranges_1based','-'))}</code></div>"
            f"{coverage_svg(glen, m.get('mask_ranges_1based','-'))}</div>")
    more = (f"<p class='note'>Showing {args.max_markers_shown} of {len(shown)} "
            f"rescuable/ambiguous markers.</p>" if len(shown) > args.max_markers_shown
            else "")
    detail = "\n".join(blocks) or ("<p class='note ok'>No markers are rescuable or "
                                   "ambiguous &mdash; every dropped marker cross-maps "
                                   "whole-gene or was unresolved.</p>")

    out = f"""<!doctype html><html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>meta-mage cross-map masking report</title><style>
:root {{ --bg:#0f1115; --panel:#171a21; --panel2:#1d212b; --fg:#e6e9ef;
  --dim:#9aa4b2; --line:#2a2f3a; --accent:#5aa9e6; --warn:#e0b341; --bad:#e0685a;
  --ok:#5ac48a; }}
@media (prefers-color-scheme: light) {{ :root {{ --bg:#f6f7f9; --panel:#fff;
  --panel2:#f0f2f5; --fg:#1a1d23; --dim:#5a6472; --line:#e2e5ea; --accent:#1f6fb8; }} }}
* {{ box-sizing:border-box; }}
body {{ margin:0; background:var(--bg); color:var(--fg);
  font:14px/1.5 -apple-system,BlinkMacSystemFont,"Segoe UI",Roboto,sans-serif; }}
.wrap {{ max-width:960px; margin:0 auto; padding:32px 20px 80px; }}
h1 {{ font-size:23px; margin:0 0 4px; }}
h2 {{ font-size:16px; margin:32px 0 10px; border-bottom:1px solid var(--line);
  padding-bottom:6px; }}
.meta {{ color:var(--dim); font-size:13px; }}
.cards {{ display:grid; grid-template-columns:repeat(auto-fit,minmax(140px,1fr));
  gap:12px; margin:16px 0 8px; }}
.card {{ background:var(--panel); border:1px solid var(--line); border-radius:10px;
  padding:13px 15px; }}
.card .val {{ font-size:24px; font-weight:650; }}
.card .lbl {{ color:var(--dim); font-size:12px; margin-top:2px; }}
.tablewrap {{ overflow-x:auto; border:1px solid var(--line); border-radius:10px;
  background:var(--panel); }}
table {{ border-collapse:collapse; width:100%; font-size:13px; }}
th,td {{ padding:7px 10px; text-align:right; white-space:nowrap;
  border-bottom:1px solid var(--line); }}
th {{ background:var(--panel2); color:var(--dim); font-weight:600; font-size:12px; }}
td.l, th:first-child {{ text-align:left; }}
td.dim {{ color:var(--dim); }} td.ok {{ color:var(--ok); }}
td.strong {{ font-weight:700; color:var(--accent); }}
.mk {{ background:var(--panel); border:1px solid var(--line); border-radius:9px;
  padding:11px 14px; margin-bottom:10px; }}
.mkhdr {{ display:flex; justify-content:space-between; align-items:baseline; gap:10px; }}
.mkcl {{ font-weight:600; font-size:13px; }}
.mkmeta {{ color:var(--dim); font-size:12px; margin:3px 0 2px; }}
.badge {{ font-size:11px; padding:2px 8px; border-radius:20px; font-weight:600; }}
.badge.ok {{ background:color-mix(in srgb,var(--ok) 20%,transparent); color:var(--ok); }}
.badge.warn {{ background:color-mix(in srgb,var(--warn) 22%,transparent); color:var(--warn); }}
.safe {{ fill:color-mix(in srgb,var(--ok) 55%,var(--panel2)); }}
.xmap {{ fill:var(--bad); fill-opacity:.85; }}
code {{ background:var(--panel2); padding:1px 5px; border-radius:4px; font-size:12px; }}
.note {{ color:var(--dim); font-size:13px; }} .note.ok {{ color:var(--ok); }}
.legend {{ color:var(--dim); font-size:12px; margin:6px 0 0; }}
</style></head><body><div class="wrap">
<h1>meta-mage &mdash; cross-map masking</h1>
<div class="meta">Can the guard-dropped markers be rescued by masking? A marker is
<b>rescuable</b> when a clade-unique stretch &ge; {args.min_survivor} bp remains
after masking the cross-mapped range.</div>
<div class="cards">
<div class="card"><div class="val">{n_species}</div><div class="lbl">species with dropped markers</div></div>
<div class="card"><div class="val">{tot_dropped}</div><div class="lbl">dropped markers</div></div>
<div class="card"><div class="val" style="color:var(--ok)">{tot_resc}</div><div class="lbl">rescuable by masking</div></div>
<div class="card"><div class="val">{tot_amb}</div><div class="lbl">ambiguous</div></div>
<div class="card"><div class="val">{tot_whole}</div><div class="lbl">whole-gene cross-map</div></div>
<div class="card"><div class="val">{tot_unres}</div><div class="lbl">unresolved (non-rep hits)</div></div>
</div>

<h2>Per species</h2>
<div class="tablewrap"><table><thead><tr><th>Species</th><th>Genomes</th>
<th>Final markers</th><th>Dropped</th><th>Rescuable</th><th>Ambiguous</th>
<th>Whole-gene</th><th>Unresolved</th><th>% rescuable</th></tr></thead>
<tbody>{''.join(srows)}</tbody></table></div>

<h2>Rescuable / ambiguous markers</h2>
<p class="legend"><span style="color:var(--ok)">&#9632;</span> clade-unique (keep)
&nbsp; <span style="color:var(--bad)">&#9632;</span> cross-mapped (mask)</p>
{more}
{detail}
</div></body></html>"""

    with open(args.out, "w") as fh:
        fh.write(out)
    print(f"masking_report: wrote {args.out} ({n_species} species, {tot_dropped} "
          f"dropped, {tot_resc} rescuable)")


if __name__ == "__main__":
    main()

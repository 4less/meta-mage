#!/usr/bin/env python3
"""
Per-species marker LEAKAGE report from the specificity guard.

The cross-map guard (specificity_guard.py) records, for every marker, whether its
nucleotide sequence also occurs OUTSIDE its home clade at read-mapping identity --
i.e. whether reads leak between taxa. specificity_report.tsv gives per marker:
    rank  clade  cluster  pass  n_offtarget  max_offtarget_id  worst_offtarget_clade
A leaking marker (pass=0) has a home clade and a single worst off-target clade it
leaks INTO, so every leak is one directed edge  home -> worst.

This report summarises that per species (default rank):
  * OUTGOING -- markers whose home is the focus species that leak into other
    species (the focus species' genes are not unique; reads it emits can be stolen).
  * INCOMING -- markers of OTHER species whose worst off-target is the focus species
    (foreign genes that would steal reads FROM the focus species).
For each direction the partner species are listed with the gene count, total
off-target hit count, and the max off-target identity.

Output is one self-contained HTML page (stdlib only) with a dropdown to pick the
focus species; all data is embedded so the page needs no server.

Only the WORST off-target clade is recorded per marker, so the directed edges are a
per-marker lower bound on the full leakage graph (a marker leaking to several
clades is attributed to its single worst one).
"""
import argparse
import html
import json
import os
from collections import defaultdict


def esc(x):
    return html.escape(str(x))


def short(clade):
    """Drop the rank prefix (s__, g__, ...) for display."""
    return clade.split("__", 1)[1] if "__" in clade else clade


def load_report(path, rank):
    """rows for the requested rank -> (home, pass, n_off, max_id, worst)."""
    rows = []
    with open(path) as fh:
        header = fh.readline()
        if not header:
            return rows
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 7:
                continue
            r, clade, cluster, passed, n_off, max_id, worst = f[:7]
            if rank and r != rank:
                continue
            rows.append((clade, passed == "1", int(n_off or 0),
                         float(max_id or 0.0), worst))
    return rows


def build_data(rows):
    """Aggregate per-species outgoing/incoming leakage partner tables."""
    species = set()
    totals = defaultdict(int)      # home -> markers considered
    leaking = defaultdict(int)     # home -> leaking markers
    # home -> partner -> [genes, hits, max_id]
    out = defaultdict(lambda: defaultdict(lambda: [0, 0, 0.0]))
    inc = defaultdict(lambda: defaultdict(lambda: [0, 0, 0.0]))

    for home, passed, n_off, max_id, worst in rows:
        species.add(home)
        totals[home] += 1
        if passed or not worst:
            continue
        species.add(worst)
        leaking[home] += 1
        o = out[home][worst]
        o[0] += 1; o[1] += n_off; o[2] = max(o[2], max_id)
        i = inc[worst][home]
        i[0] += 1; i[1] += n_off; i[2] = max(i[2], max_id)

    data = {}
    for sp in sorted(species):
        def partners(tbl):
            return sorted(
                ({"clade": short(c), "genes": v[0], "hits": v[1],
                  "max_id": round(v[2] * 100, 1) if v[2] <= 1.0 else round(v[2], 1)}
                 for c, v in tbl.items()),
                key=lambda d: (-d["genes"], -d["hits"]))
        o = partners(out.get(sp, {}))
        i = partners(inc.get(sp, {}))
        data[sp] = {
            "name": short(sp),
            "total_markers": totals.get(sp, 0),
            "leaking_markers": leaking.get(sp, 0),
            "out_genes": sum(p["genes"] for p in o),
            "in_genes": sum(p["genes"] for p in i),
            "out": o,
            "in": i,
        }
    # Overview order: worst combined leakage first.
    order = sorted(data, key=lambda s: -(data[s]["out_genes"] + data[s]["in_genes"]))
    return data, order


HTML = """<!doctype html><html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>meta-mage marker leakage</title><style>
:root {{ --bg:#0f1115; --panel:#171a21; --panel2:#1d212b; --fg:#e6e9ef;
  --dim:#9aa4b2; --line:#2a2f3a; --accent:#5aa9e6; --out:#e0685a; --in:#e0b341;
  --ok:#5ac48a; }}
@media (prefers-color-scheme: light) {{ :root {{ --bg:#f6f7f9; --panel:#fff;
  --panel2:#f0f2f5; --fg:#1a1d23; --dim:#5a6472; --line:#e2e5ea; --accent:#1f6fb8; }} }}
* {{ box-sizing:border-box; }}
body {{ margin:0; background:var(--bg); color:var(--fg);
  font:14px/1.5 -apple-system,BlinkMacSystemFont,"Segoe UI",Roboto,sans-serif; }}
.wrap {{ max-width:900px; margin:0 auto; padding:32px 20px 80px; }}
h1 {{ font-size:23px; margin:0 0 4px; }}
h2 {{ font-size:15px; margin:22px 0 8px; }}
.meta {{ color:var(--dim); font-size:13px; margin-bottom:18px; }}
.controls {{ display:flex; gap:10px; align-items:center; margin:6px 0 20px;
  flex-wrap:wrap; }}
select {{ font:inherit; padding:7px 10px; border-radius:8px; border:1px solid var(--line);
  background:var(--panel); color:var(--fg); max-width:100%; min-width:280px; }}
.cards {{ display:grid; grid-template-columns:repeat(auto-fit,minmax(150px,1fr));
  gap:12px; margin:8px 0 20px; }}
.card {{ background:var(--panel); border:1px solid var(--line); border-radius:10px;
  padding:13px 15px; }}
.card .val {{ font-size:24px; font-weight:650; }}
.card .lbl {{ color:var(--dim); font-size:12px; margin-top:2px; }}
.card.out .val {{ color:var(--out); }} .card.in .val {{ color:var(--in); }}
table {{ width:100%; border-collapse:collapse; font-size:13px; margin-bottom:8px; }}
th, td {{ text-align:left; padding:6px 10px; border-bottom:1px solid var(--line); }}
th {{ color:var(--dim); font-weight:600; font-size:12px; }}
td.n {{ text-align:right; font-variant-numeric:tabular-nums; }}
.none {{ color:var(--dim); font-size:13px; padding:4px 0 10px; }}
.tag {{ font-size:11px; padding:2px 8px; border-radius:20px; font-weight:600; }}
.tag.out {{ background:color-mix(in srgb,var(--out) 20%,transparent); color:var(--out); }}
.tag.in {{ background:color-mix(in srgb,var(--in) 22%,transparent); color:var(--in); }}
.legend {{ color:var(--dim); font-size:12px; margin:2px 0 8px; }}
.rate {{ color:var(--dim); font-size:12px; }}
</style></head><body><div class="wrap">
<h1>meta-mage &mdash; marker leakage ({rank})</h1>
<div class="meta">Cross-map leakage between {rank} from the specificity guard.
<b>Outgoing</b> = the focus species' own markers that also occur in another species
(its reads can be stolen). <b>Incoming</b> = other species' markers that also occur
in the focus species (it steals their reads).</div>
<div class="cards">
<div class="card"><div class="val">{n_species}</div><div class="lbl">species with markers</div></div>
<div class="card"><div class="val">{n_leaky}</div><div class="lbl">species that leak (either direction)</div></div>
<div class="card out"><div class="val">{n_edges}</div><div class="lbl">leaking markers total</div></div>
</div>
<div class="controls"><label for="sel">Focus species:</label>
<select id="sel" onchange="render()"></select></div>
<div id="focus"></div>
<script>
const DATA = {data_json};
const ORDER = {order_json};
const sel = document.getElementById('sel');
for (const s of ORDER) {{
  const o = document.createElement('option');
  o.value = s;
  const d = DATA[s];
  o.textContent = d.name + '  (out ' + d.out_genes + ' / in ' + d.in_genes + ')';
  sel.appendChild(o);
}}
function ptable(rows, dir) {{
  if (!rows.length) return '<div class="none">none</div>';
  let h = '<table><thead><tr><th>' + (dir==='out'?'leaks into':'leaks from') +
    '</th><th class="n">genes</th><th class="n">off-target hits</th>' +
    '<th class="n">max id %</th></tr></thead><tbody>';
  for (const r of rows)
    h += '<tr><td>' + r.clade + '</td><td class="n">' + r.genes +
      '</td><td class="n">' + r.hits + '</td><td class="n">' + r.max_id + '</td></tr>';
  return h + '</tbody></table>';
}}
function render() {{
  const d = DATA[sel.value];
  if (!d) {{ document.getElementById('focus').innerHTML =
    '<div class="none">No markers with leakage data.</div>'; return; }}
  const rate = d.total_markers ?
    (100*d.leaking_markers/d.total_markers).toFixed(0) + '%' : '&mdash;';
  document.getElementById('focus').innerHTML =
    '<div class="cards">' +
    '<div class="card out"><div class="val">' + d.out_genes + '</div>' +
      '<div class="lbl">outgoing leaking genes</div></div>' +
    '<div class="card in"><div class="val">' + d.in_genes + '</div>' +
      '<div class="lbl">incoming leaking genes</div></div>' +
    '<div class="card"><div class="val">' + d.leaking_markers + '/' + d.total_markers +
      '</div><div class="lbl">own markers leaking (' + rate + ')</div></div>' +
    '</div>' +
    '<h2><span class="tag out">outgoing</span> &nbsp;' + d.name +
      "'s markers found in other species</h2>" + ptable(d.out, 'out') +
    '<h2><span class="tag in">incoming</span> &nbsp;other species&rsquo; markers found in ' +
      d.name + '</h2>' + ptable(d.in, 'in');
}}
render();
</script>
</div></body></html>"""


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--specificity", required=True,
                    help="specificity_report.tsv (or NO_FILE)")
    ap.add_argument("--rank", default="species",
                    help="clade rank to summarise (default: species)")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    if os.path.basename(args.specificity) == "NO_FILE" or \
            not os.path.exists(args.specificity):
        with open(args.out, "w") as fh:
            fh.write("<!doctype html><meta charset='utf-8'><title>marker leakage</title>"
                     "<p style='font-family:sans-serif;padding:2em'>The nucleotide "
                     "cross-map guard was not run (<code>--specificity false</code>), "
                     "so there is no leakage data.</p>")
        print("leakage_report: guard not run -> placeholder page")
        return

    rows = load_report(args.specificity, args.rank)
    data, order = build_data(rows)
    n_leaky = sum(1 for s in data if data[s]["out_genes"] or data[s]["in_genes"])
    n_edges = sum(data[s]["out_genes"] for s in data)

    page = HTML.format(
        rank=esc(args.rank),
        n_species=len(data),
        n_leaky=n_leaky,
        n_edges=n_edges,
        data_json=json.dumps(data),
        order_json=json.dumps(order),
    )
    with open(args.out, "w") as fh:
        fh.write(page)
    print(f"leakage_report: wrote {args.out} ({len(data)} {args.rank}, "
          f"{n_leaky} leaking, {n_edges} leaking markers)")


if __name__ == "__main__":
    main()

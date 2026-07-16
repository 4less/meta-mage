#!/usr/bin/env python3
"""
Per-species gene-prevalence curves -- a standalone, self-contained HTML.

For each species, every gene family (cluster) in its pangenome is placed on the
x-axis (ranked by prevalence, most-prevalent first) and its in-species prevalence
on the y-axis. The resulting curve shows the core plateau (left) falling into the
accessory tail (right); horizontal guides mark the core-prevalence cutoff and the
relaxation floor, so you can see how much of the pangenome each admits.

Prevalence is completeness-weighted to match SCORE / the funnel:
    prevalence = (size - (score_sum - in_score)) / size
(with no completeness metadata this is just in_count / size). The full curve is
downsampled to --max_points for rendering (pangenomes reach 10^5 genes).

Reads counts.tsv + clade_sizes.tsv; writes one HTML with a species dropdown.
"""
import argparse
import html
import json


def esc(x):
    return html.escape(str(x))


def load_sizes(path, rank):
    sizes = {}
    with open(path) as fh:
        next(fh)
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if f[0] != rank:
                continue
            size = int(f[2])
            score_sum = float(f[3]) if len(f) > 3 else float(size)
            sizes[f[1]] = (size, score_sum)
    return sizes


def downsample(values, max_points):
    """Keep the shape of a (sorted) list within max_points via uniform striding,
    always retaining the first and last point."""
    n = len(values)
    if n <= max_points:
        return values
    step = n / max_points
    out = [values[min(n - 1, int(i * step))] for i in range(max_points)]
    out[-1] = values[-1]
    return out


PAGE = """<!doctype html><html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>meta-mage gene prevalence</title><style>
:root {{ --bg:#0f1115; --panel:#171a21; --panel2:#1d212b; --fg:#e6e9ef;
  --dim:#9aa4b2; --line:#2a2f3a; --accent:#5aa9e6; --warn:#e0b341; --ok:#5ac48a; }}
@media (prefers-color-scheme: light) {{ :root {{ --bg:#f6f7f9; --panel:#fff;
  --panel2:#f0f2f5; --fg:#1a1d23; --dim:#5a6472; --line:#e2e5ea; --accent:#1f6fb8; }} }}
* {{ box-sizing:border-box; }}
body {{ margin:0; background:var(--bg); color:var(--fg);
  font:14px/1.5 -apple-system,BlinkMacSystemFont,"Segoe UI",Roboto,sans-serif; }}
.wrap {{ max-width:1100px; margin:0 auto; padding:28px 20px 70px; }}
h1 {{ font-size:22px; margin:0 0 4px; }}
.meta {{ color:var(--dim); font-size:13px; }}
.legend {{ color:var(--dim); font-size:13px; margin:14px 2px; }}
.ctl {{ margin:12px 2px; color:var(--dim); font-size:13px; }}
select {{ background:var(--panel2); color:var(--fg); border:1px solid var(--line);
  border-radius:6px; padding:5px 9px; font-size:14px; }}
.card {{ background:var(--panel); border:1px solid var(--line); border-radius:10px;
  padding:10px 12px; overflow-x:auto; }}
.stat {{ color:var(--dim); font-size:12px; margin:8px 2px 0; }}
.stat b {{ color:var(--fg); }}
</style></head><body><div class="wrap">
<h1>meta-mage &mdash; gene prevalence per species</h1>
<div class="meta">{meta}</div>
<p class="legend">Every gene family in the species' pangenome, ranked by
in-species prevalence (most prevalent first). The core plateau on the left falls
into the accessory tail on the right. Guides: <span style="color:var(--accent)">
core cutoff {core_prev:.0%}</span> and <span style="color:var(--warn)">relax floor
{relax_floor:.0%}</span>.</p>
<div class="ctl">Species: <select id="sel" onchange="render()"></select></div>
<div class="card"><svg id="plot" width="1040" height="360"></svg></div>
<div class="stat" id="stat"></div>
<script>
const DATA = {data_json};
const ORDER = {order_json};
const CORE = {core_prev}, FLOOR = {relax_floor};
const W=1040, H=360, PADL=52, PADR=16, PADT=16, PADB=34;
function x(f){{ return PADL + f*(W-PADL-PADR); }}
function y(v){{ return PADT + (1-v)*(H-PADT-PADB); }}
function render(){{
  const sp = document.getElementById('sel').value;
  const d = DATA[sp]; const ys = d.y; const n = d.n;
  let s='';
  [0,0.25,0.5,0.75,1].forEach(function(v){{
    s+='<line x1="'+PADL+'" y1="'+y(v)+'" x2="'+(W-PADR)+'" y2="'+y(v)+'" stroke="var(--line)"/>';
    s+='<text x="'+(PADL-7)+'" y="'+(y(v)+3)+'" text-anchor="end" font-size="10" fill="var(--dim)">'+(v*100).toFixed(0)+'%</text>';
  }});
  [[CORE,'var(--accent)'],[FLOOR,'var(--warn)']].forEach(function(g){{
    s+='<line x1="'+PADL+'" y1="'+y(g[0])+'" x2="'+(W-PADR)+'" y2="'+y(g[0])+'" stroke="'+g[1]+'" stroke-dasharray="4 3" stroke-width="1.2"/>';
  }});
  let pts='';
  for(let i=0;i<ys.length;i++){{ const f=ys.length>1? i/(ys.length-1):0;
    pts += x(f).toFixed(1)+','+y(ys[i]).toFixed(1)+' '; }}
  s+='<polyline points="'+pts+'" fill="none" stroke="var(--fg)" stroke-width="1.5"/>';
  // x-axis labels: 0 .. pangenome size
  s+='<text x="'+PADL+'" y="'+(H-8)+'" font-size="10" fill="var(--dim)">0</text>';
  s+='<text x="'+(W-PADR)+'" y="'+(H-8)+'" text-anchor="end" font-size="10" fill="var(--dim)">'+n.toLocaleString()+' gene families</text>';
  document.getElementById('plot').innerHTML=s;
  document.getElementById('stat').innerHTML =
    'Pangenome <b>'+n.toLocaleString()+'</b> &middot; core (&ge;'+(CORE*100).toFixed(0)+'%) <b>'+d.core.toLocaleString()+'</b>'
    +' &middot; &ge;'+(FLOOR*100).toFixed(0)+'% <b>'+d.floor.toLocaleString()+'</b> &middot; genomes <b>'+d.g.toLocaleString()+'</b>';
}}
(function(){{ const sel=document.getElementById('sel');
  ORDER.forEach(function(sp){{ const o=document.createElement('option');
    o.value=sp; o.textContent=sp; sel.appendChild(o); }});
  if(ORDER.length) render(); }})();
</script></div></body></html>"""


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--counts", required=True)
    ap.add_argument("--clade_sizes", required=True)
    ap.add_argument("--rank", default="species")
    ap.add_argument("--core_prev", type=float, default=0.8)
    ap.add_argument("--relax_floor", type=float, default=0.5)
    ap.add_argument("--max_points", type=int, default=2000)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    sizes = load_sizes(args.clade_sizes, args.rank)

    prev = {sp: [] for sp in sizes}   # species -> [prevalence, ...]
    with open(args.counts) as fh:
        next(fh)
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 5 or f[1] != args.rank:
                continue
            ent = sizes.get(f[2])
            if not ent:
                continue
            size, score_sum = ent
            in_count = int(f[3])
            in_score = float(f[5]) if len(f) > 5 else float(in_count)
            p = (size - (score_sum - in_score)) / size if size else 0.0
            prev[f[2]].append(min(1.0, p))

    data = {}
    for sp, vals in prev.items():
        if not vals:
            continue
        vals.sort(reverse=True)
        size = sizes[sp][0]
        data[sp] = {
            "n": len(vals),
            "g": size,
            "core": sum(1 for v in vals if v >= args.core_prev),
            "floor": sum(1 for v in vals if v >= args.relax_floor),
            "y": [round(v, 4) for v in downsample(vals, args.max_points)],
        }
    order = sorted(data, key=lambda s: -data[s]["n"])   # biggest pangenome first

    meta = (f"{len(data)} species &middot; prevalence ranked over the pangenome "
            f"&middot; curves downsampled to {args.max_points} points")
    page = PAGE.format(
        meta=meta, core_prev=args.core_prev, relax_floor=args.relax_floor,
        data_json=json.dumps(data, separators=(",", ":")).replace("</", "<\\/"),
        order_json=json.dumps(order).replace("</", "<\\/"),
    )
    with open(args.out, "w") as fh:
        fh.write(page)
    print(f"prevalence_curves: {len(data)} species -> {args.out}")


if __name__ == "__main__":
    main()

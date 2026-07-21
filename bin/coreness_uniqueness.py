#!/usr/bin/env python3
"""
Per-species coreness-vs-uniqueness scatter report.

For every species, plot each of its gene clusters in the space the marker selector
actually works in:
    x = coreness   = completeness-weighted in_prevalence  (present in ~all of the
                     species' genomes)          -> selection needs >= min_in (0.80)
    y = uniqueness = -log10(out_prevalence)      (rare in OTHER genomes)
                     out_prev = (off-clade genomes with the gene)/(off-clade genomes)
                     -> selection needs out_prev <= max_out (0.02), i.e. y >= 1.70

The selection window is the top-right box. Each point is coloured by what happens
to it downstream:
    grey   = not a candidate (fails coreness and/or uniqueness)
    amber  = candidate (in the box) dropped BEFORE cross-map (length filter)
    red    = candidate that reached the nucleotide guard and CROSS-MAPS (fails c)
    green  = final marker (survives the cross-map guard)

So a species whose box is full of red has genes that look core+unique by the
genome-pooled prevalence test but are not nucleotide-specific -- the gate-(b) vs
gate-(c) disagreement, made visible. Reads only local run tables; output is one
self-contained HTML file with a species picker.
"""
import argparse
import html
import json
import math
from collections import defaultdict


def load_manifest(path):
    """species -> (size_genomes, score_sum completeness); N genomes total."""
    size = defaultdict(int)
    score = defaultdict(float)
    n_total = 0
    with open(path) as fh:
        col = {n: i for i, n in enumerate(fh.readline().rstrip("\n").split("\t"))}
        has_comp = "completeness" in col
        for line in fh:
            f = line.rstrip("\n").split("\t")
            n_total += 1
            sp = f[col["species"]]
            if sp == "NA":
                continue
            size[sp] += 1
            score[sp] += float(f[col["completeness"]]) if has_comp else 1.0
    return size, score, n_total


def load_guard(path):
    """(clade, cluster) -> pass(bool); plus the set of clusters that reached the
    guard per species."""
    verdict = {}
    with open(path) as fh:
        col = {n: i for i, n in enumerate(fh.readline().rstrip("\n").split("\t"))}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if f[col["rank"]] != "species":
                continue
            verdict[(f[col["clade"]], f[col["cluster"]])] = (f[col["pass"]] == "1")
    return verdict


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--counts", required=True, help="counts.tsv (species rows)")
    ap.add_argument("--manifest", required=True)
    ap.add_argument("--specificity", required=True)
    ap.add_argument("--min_in", type=float, default=0.80)
    ap.add_argument("--max_out", type=float, default=0.02)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    size, score, n_total = load_manifest(args.manifest)
    guard = load_guard(args.specificity)
    y_thresh = -math.log10(args.max_out)              # 0.02 -> 1.70
    floor = 0.5 / max(n_total, 2)                     # out_prev==0 lands here
    y_cap = -math.log10(floor)

    # EVERY gene cluster is included. Non-candidates are binned into a density
    # grid (so the full pangenome cloud shows without a huge point count);
    # candidates are kept as individual points so their cross-map colour is exact.
    bins = defaultdict(lambda: defaultdict(int))   # species -> {(bx,by): count}
    pts = defaultdict(list)                         # species -> [[x,y,status]]
    summ = defaultdict(lambda: defaultdict(int))

    with open(args.counts) as fh:
        fh.readline()
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 6 or f[1] != "species":
                continue
            cluster, clade = f[0], f[2]
            sz = size.get(clade)
            if not sz:
                continue
            in_count = int(f[3]); marker_total = int(f[4]); in_score = float(f[5])
            ss = score[clade]
            in_prev = (sz - (ss - in_score)) / sz
            if in_prev > 1.0:
                in_prev = 1.0
            out_denom = n_total - sz
            if out_denom <= 0:
                continue
            out_prev = (marker_total - in_count) / out_denom
            y = -math.log10(out_prev) if out_prev > floor else y_cap
            summ[clade]["total"] += 1
            candidate = in_prev >= args.min_in and out_prev <= args.max_out
            if candidate:
                v = guard.get((clade, cluster))
                status = 2 if v is None else (4 if v else 3)   # amber/red/green
                summ[clade]["cand"] += 1
                summ[clade]["final" if v else ("len" if v is None else "cross")] += 1
                pts[clade].append([round(in_prev, 3), round(y, 3), status])
            else:
                bx = round(in_prev, 2)                         # 0.01 coreness bins
                by = round(y, 1)                               # 0.1 uniqueness bins
                bins[clade][(bx, by)] += 1

    species = sorted(size, key=lambda s: (summ[s].get("final", 0),
                                          summ[s].get("cand", 0)))
    meta = {}
    for s in species:
        if summ[s].get("total", 0) == 0:
            continue
        bl = [[bx, by, c] for (bx, by), c in bins[s].items()]
        meta[s.replace("s__", "")] = {
            "bins": bl,
            "pts": pts.get(s, []),
            "genomes": size[s],
            "total": summ[s].get("total", 0),
            "cand": summ[s].get("cand", 0),
            "final": summ[s].get("final", 0),
            "cross": summ[s].get("cross", 0),
            "len": summ[s].get("len", 0),
        }

    payload = {"species": meta, "y_thresh": round(y_thresh, 3),
               "y_cap": round(y_cap, 3), "min_in": args.min_in,
               "max_out": args.max_out,
               "default": "Bacteroides ovatus" if "Bacteroides ovatus" in meta
                          else next(iter(meta), "")}
    out = TEMPLATE.replace("__DATA__", json.dumps(payload))
    with open(args.out, "w") as fh:
        fh.write(out)
    print(f"coreness_uniqueness: {len(meta)} species -> {args.out}")
    for s in ("Bacteroides ovatus",):
        if s in meta:
            m = meta[s]
            print(f"  {s}: {m['cand']} candidates (core+unique by gate b) -> "
                  f"{m['final']} final, {m['cross']} cross-map, {m['len']} length-dropped")


TEMPLATE = r"""<!doctype html><html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>meta-mage coreness vs uniqueness</title><style>
:root{--bg:#0f1115;--panel:#171a21;--fg:#e6e9ef;--dim:#9aa4b2;--line:#2a2f3a;
  --grid:#232833;--box:#e0b341;--bg2:#4a5262;--cross:#e0685a;--final:#5ac48a;
  --amber:#e0a020;}
@media (prefers-color-scheme:light){:root{--bg:#f6f7f9;--panel:#fff;--fg:#1a1d23;
  --dim:#5a6472;--line:#e2e5ea;--grid:#eceef2;--bg2:#c2c8d2;}}
*{box-sizing:border-box}
body{margin:0;background:var(--bg);color:var(--fg);
  font:14px/1.5 -apple-system,BlinkMacSystemFont,"Segoe UI",Roboto,sans-serif}
.wrap{display:flex;gap:20px;max-width:1180px;margin:0 auto;padding:22px 20px 60px;
  flex-wrap:wrap}
.left{flex:1 1 600px;min-width:520px}
.side{flex:0 0 300px;min-width:260px}
h1{font-size:20px;margin:0 0 2px}
.meta{color:var(--dim);font-size:13px;margin-bottom:10px}
svg{width:100%;height:auto;display:block;background:var(--panel);
  border:1px solid var(--line);border-radius:10px}
.axlbl{fill:var(--dim);font-size:12px}
.axnum{fill:var(--dim);font-size:10px}
.grid{stroke:var(--grid);stroke-width:1}
.thr{stroke:var(--box);stroke-width:1.2;stroke-dasharray:4 3}
.boxfill{fill:var(--box);fill-opacity:.07}
.panel{background:var(--panel);border:1px solid var(--line);border-radius:10px;
  padding:13px 15px;margin-bottom:14px}
.panel h2{font-size:12px;margin:0 0 8px;color:var(--dim);font-weight:600;
  text-transform:uppercase;letter-spacing:.03em}
input{width:100%;padding:7px 9px;border:1px solid var(--line);border-radius:7px;
  background:var(--bg);color:var(--fg);font-size:13px}
.list{max-height:340px;overflow:auto;margin-top:8px}
.row{display:flex;justify-content:space-between;gap:8px;padding:3px 7px;
  border-radius:6px;font-size:12.5px;cursor:pointer}
.row:hover{background:var(--bg)}.row.sel{background:var(--bg)}
.row .n{overflow:hidden;text-overflow:ellipsis;white-space:nowrap}
.row .v{color:var(--dim);font-variant-numeric:tabular-nums}
.row.zero .v{color:var(--cross)}
.legend{font-size:12px;color:var(--dim);line-height:1.9}
.sw{display:inline-block;width:10px;height:10px;border-radius:50%;
  vertical-align:middle;margin-right:5px}
.stat{font-size:13px;margin:2px 0}
b.f{color:var(--box)} .kpi{font-size:22px;font-weight:650}
</style></head><body><div class="wrap">
<div class="left">
<h1>meta-mage &mdash; coreness vs uniqueness</h1>
<div class="meta">each point = one gene cluster of the selected species &middot;
box = marker selection window (core &amp; unique by the genome-pooled prevalence
test) &middot; colour = cross-map outcome</div>
<svg id="viz" viewBox="0 0 640 470" role="img"></svg>
<p class="legend" style="margin-top:8px">
<span class="sw" style="background:var(--bg2)"></span>all other genes (density) &nbsp;
<span class="sw" style="background:var(--cross)"></span>candidate, cross-maps (fails nucleotide guard) &nbsp;
<span class="sw" style="background:var(--final)"></span>final marker &nbsp;
<span class="sw" style="background:var(--amber)"></span>dropped (length)</p>
</div>
<div class="side">
  <div class="panel">
    <h2>Species</h2>
    <div class="stat" id="summary"></div>
  </div>
  <div class="panel">
    <h2>Pick a species (sorted by final markers)</h2>
    <input id="search" placeholder="type a species...">
    <div class="list" id="picker"></div>
  </div>
</div>
</div>
<script>
const D = __DATA__;
const M = D.species;
let cur = D.default;
const NS="http://www.w3.org/2000/svg";
function el(t,a){const e=document.createElementNS(NS,t);for(const k in a)e.setAttribute(k,a[k]);return e;}
const COL=["var(--bg2)",null,"var(--amber)","var(--cross)","var(--final)"];
const R=[1.3,0,2.6,2.6,3.0], OP=[.5,0,.95,.8,.95];

// plot geometry
const W=640,H=470,ml=52,mr=16,mt=16,mb=42;
const x0=0,x1=1, y0=0,y1=D.y_cap;
function px(x){return ml+(x-x0)/(x1-x0)*(W-ml-mr);}
function py(y){return H-mb-(y-y0)/(y1-y0)*(H-mt-mb);}

function draw(){
  const m=M[cur]; if(!m){return;}
  const svg=document.getElementById('viz'); svg.innerHTML='';
  // selection box (x>=min_in, y>=y_thresh)
  svg.appendChild(el('rect',{x:px(D.min_in),y:py(y1),width:px(x1)-px(D.min_in),
    height:py(D.y_thresh)-py(y1),class:'boxfill'}));
  // gridlines + axis numbers
  for(let gx=0;gx<=1.0001;gx+=0.2){
    svg.appendChild(el('line',{x1:px(gx),y1:mt,x2:px(gx),y2:H-mb,class:'grid'}));
    const t=el('text',{x:px(gx),y:H-mb+14,class:'axnum','text-anchor':'middle'});
    t.textContent=gx.toFixed(1); svg.appendChild(t);
  }
  for(let gy=0;gy<=y1+.001;gy+=1){
    svg.appendChild(el('line',{x1:ml,y1:py(gy),x2:W-mr,y2:py(gy),class:'grid'}));
    const t=el('text',{x:ml-6,y:py(gy)+3,class:'axnum','text-anchor':'end'});
    t.textContent=gy; svg.appendChild(t);
  }
  // thresholds
  svg.appendChild(el('line',{x1:px(D.min_in),y1:mt,x2:px(D.min_in),y2:H-mb,class:'thr'}));
  svg.appendChild(el('line',{x1:ml,y1:py(D.y_thresh),x2:W-mr,y2:py(D.y_thresh),class:'thr'}));
  // axis labels
  const xl=el('text',{x:(ml+W-mr)/2,y:H-6,class:'axlbl','text-anchor':'middle'});
  xl.textContent='coreness  =  in-prevalence  (→ core at 0.80)'; svg.appendChild(xl);
  const yl=el('text',{x:14,y:(mt+H-mb)/2,class:'axlbl','text-anchor':'middle',
    transform:`rotate(-90 14 ${(mt+H-mb)/2})`});
  yl.textContent='uniqueness  =  −log₁₀(out-prevalence)  (↑ unique, 0.02 = 1.70)';
  svg.appendChild(yl);
  // ALL genes: non-candidates as a density grid, candidates as points on top.
  const cellw=(px(1)-px(0))*0.01, cellh=(py(0)-py(0.1));
  let maxc=1; for(const b of m.bins) if(b[2]>maxc) maxc=b[2];
  const lm=Math.log(maxc+1);
  for(const b of m.bins){
    const o=0.12+0.66*Math.log(b[2]+1)/lm;
    svg.appendChild(el('rect',{x:(px(b[0])-cellw/2).toFixed(1),
      y:(py(b[1])-cellh/2).toFixed(1),width:Math.max(1.4,cellw).toFixed(1),
      height:Math.max(1.4,cellh).toFixed(1),fill:'var(--bg2)','fill-opacity':o.toFixed(2)}));
  }
  const arr=m.pts.slice().sort((a,b)=>a[2]-b[2]);   // candidates on top
  for(const p of arr){
    svg.appendChild(el('circle',{cx:px(p[0]).toFixed(1),cy:py(p[1]).toFixed(1),
      r:R[p[2]],fill:COL[p[2]],'fill-opacity':OP[p[2]]}));
  }
  // summary
  const nm=cur.replace('Bacteroides','B.');
  document.getElementById('summary').innerHTML=
    `<div class="kpi f">${nm}</div>`+
    `<div class="stat">${m.genomes} genomes &middot; <b>${m.total}</b> gene clusters (all plotted)</div>`+
    `<div class="stat"><b>${m.cand}</b> candidates in the box (core &amp; unique by gate b)</div>`+
    `<div class="stat"><span class="sw" style="background:var(--final)"></span><b>${m.final}</b> final markers</div>`+
    `<div class="stat"><span class="sw" style="background:var(--cross)"></span><b>${m.cross}</b> cross-map (fail gate c)</div>`+
    (m.len?`<div class="stat"><span class="sw" style="background:var(--amber)"></span>${m.len} dropped (length)</div>`:``);
}

function buildPicker(q){
  const p=document.getElementById('picker'); p.innerHTML='';
  q=(q||'').toLowerCase();
  Object.keys(M).filter(s=>s.toLowerCase().includes(q)).forEach(s=>{
    const m=M[s];
    const row=el('div',{class:'row'+(s===cur?' sel':'')+(m.final===0&&m.cand>0?' zero':'')});
    row.innerHTML=`<span class="n">${s.replace('Bacteroides','B.')}</span>`+
      `<span class="v">${m.final}/${m.cand}</span>`;
    row.onclick=()=>{cur=s;draw();buildPicker(document.getElementById('search').value);};
    p.appendChild(row);
  });
}
document.getElementById('search').addEventListener('input',e=>buildPicker(e.target.value));
buildPicker(''); draw();
</script></body></html>"""


if __name__ == "__main__":
    main()

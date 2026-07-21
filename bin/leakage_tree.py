#!/usr/bin/env python3
"""
Radial-phylogeny leakage browser: draw a subclade's tree and, for one focus
taxon at a time, the incoming and outgoing cross-map "leakage" to other taxa --
so you can SEE whether a species' leakage follows the phylogeny (arcs stay local
on the tree) or crosses it (arcs jump to distant clades).

Inputs
  --tree      Newick of the subclade (tips = genome accessions), e.g.
              data/phylogeny/bacteroides_subclade.tree
  --species   accession<TAB>species map (data/phylogeny/bacteroides_species.tsv)
  --leakage   OPTIONAL directed edge table, one row per ordered species pair:
                  query_species<TAB>target_species<TAB>weight
              meaning: a marker of query_species cross-maps onto a CDS of
              target_species (query "leaks out" to target). Derive it from
              crossmap.m8 + query.map + manifest (the same ids low_marker_masking
              already resolves). If omitted, a clearly-labelled DEMO set is drawn.
  --focus     default focus species (default: Bacteroides ovatus)
  --out       output HTML

Output is a single self-contained HTML file (stdlib only, inline SVG + JS): a
radial cladogram with a focus selector; outgoing leakage (focus -> other) and
incoming leakage (other -> focus) are drawn as directed arcs across the interior,
width proportional to weight, with the leakage partners highlighted on the tree.
"""
import argparse
import html
import json
import math
import re


# ----------------------------- Newick parser ------------------------------
class Node:
    __slots__ = ("name", "length", "children", "parent", "idx")

    def __init__(self):
        self.name = ""
        self.length = 0.0
        self.children = []
        self.parent = None
        self.idx = -1


def parse_newick(text):
    text = text.strip()
    if text.endswith(";"):
        text = text[:-1]
    # tokenizer that respects single-quoted labels
    toks = re.findall(r"'[^']*'|[(),]|[^(),]+", text)
    root = Node()
    cur = root
    i = 0

    def read_suffix(node, tok_list, j):
        # tok_list[j] may be a label and/or :length (already split by regex into
        # one chunk for the non-structural run) -- handle "'label'" then ":len"
        return j

    stack = []
    node = root
    expect_child = False
    for t in toks:
        if t == "(":
            child = Node()
            child.parent = node
            node.children.append(child)
            stack.append(node)
            node = child
        elif t == ",":
            parent = stack[-1]
            child = Node()
            child.parent = parent
            parent.children.append(child)
            node = child
        elif t == ")":
            node = stack.pop()
        else:
            # label[:length]  (label may be quoted)
            lab, length = t, None
            m = re.match(r"^('[^']*'|[^:]*)(?::([0-9eE.+-]+))?$", t)
            if m:
                lab = m.group(1)
                length = m.group(2)
            lab = lab.strip()
            if lab.startswith("'") and lab.endswith("'"):
                lab = lab[1:-1]
            node.name = lab
            if length is not None:
                try:
                    node.length = float(length)
                except ValueError:
                    node.length = 0.0
    return root


# ----------------------------- radial layout ------------------------------
def layout(root, R):
    leaves = []

    def collect(n):
        if not n.children:
            leaves.append(n)
        for c in n.children:
            collect(c)
    collect(root)
    nleaves = len(leaves)

    # topological depth (edges from root)
    depth = {}

    def setdepth(n, d):
        depth[id(n)] = d
        for c in n.children:
            setdepth(c, d + 1)
    setdepth(root, 0)
    maxdepth = max(depth.values()) or 1

    angle = {}
    for k, lf in enumerate(leaves):
        angle[id(lf)] = (k + 0.5) / nleaves * 2 * math.pi

    def setangle(n):
        if not n.children:
            return angle[id(n)]
        a = [setangle(c) for c in n.children]
        angle[id(n)] = sum(a) / len(a)
        return angle[id(n)]
    setangle(root)

    # tips aligned on outer ring; internal nodes by depth fraction
    pos = {}
    for n in depth:
        pass

    def radius(n):
        if not n.children:
            return R
        return R * depth[id(n)] / maxdepth

    nodes = []

    def walk(n):
        r = radius(n)
        a = angle[id(n)]
        x = r * math.cos(a)
        y = r * math.sin(a)
        pos[id(n)] = (x, y, r, a)
        for c in n.children:
            walk(c)
    walk(root)
    return leaves, pos, angle, depth


# ------------------------------ svg helpers -------------------------------
def polar(r, a):
    return r * math.cos(a), r * math.sin(a)


def build_tree_svg(root, pos, angle):
    """Radial cladogram edges: per internal node an arc at its radius spanning
    its children's angles, plus a radial spoke to each child."""
    parts = []
    def walk(n):
        x, y, r, a = pos[id(n)]
        if n.children:
            cas = [angle[id(c)] for c in n.children]
            a0, a1 = min(cas), max(cas)
            large = 1 if (a1 - a0) > math.pi else 0
            x0, y0 = polar(r, a0)
            x1, y1 = polar(r, a1)
            parts.append(f"<path d='M{x0:.1f},{y0:.1f} A{r:.1f},{r:.1f} 0 "
                         f"{large} 1 {x1:.1f},{y1:.1f}' class='br'/>")
            for c in n.children:
                cx, cy, cr, ca = pos[id(c)]
                sx, sy = polar(r, ca)
                parts.append(f"<line x1='{sx:.1f}' y1='{sy:.1f}' "
                             f"x2='{cx:.1f}' y2='{cy:.1f}' class='br'/>")
        for c in n.children:
            walk(c)
    walk(root)
    return "".join(parts)


# --------------------------------- main -----------------------------------
def load_species(path):
    acc2sp = {}
    with open(path) as fh:
        header = fh.readline()
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) >= 2:
                acc2sp[f[0]] = f[1]
    return acc2sp


def load_leakage(path):
    edges = []
    with open(path) as fh:
        first = fh.readline().rstrip("\n").split("\t")
        # allow a header or not
        if not (len(first) >= 3 and re.match(r"^[0-9.]+$", first[2])):
            pass  # was header
        else:
            edges.append((first[0], first[1], float(first[2])))
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) >= 3:
                try:
                    edges.append((f[0], f[1], float(f[2])))
                except ValueError:
                    continue
    return edges


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--tree", required=True)
    ap.add_argument("--species", required=True)
    ap.add_argument("--leakage")
    ap.add_argument("--focus", default="Bacteroides ovatus")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    root = parse_newick(open(args.tree).read())
    acc2sp = load_species(args.species)
    R = 360.0
    leaves, pos, angle, depth = layout(root, R)

    # tip records (species-labelled)
    tips = []
    sp2i = {}
    for lf in leaves:
        x, y, r, a = pos[id(lf)]
        sp = acc2sp.get(lf.name, lf.name)
        sp2i[sp] = len(tips)
        tips.append({"i": len(tips), "sp": sp, "acc": lf.name,
                     "x": round(x, 1), "y": round(y, 1),
                     "a": round(a, 4)})

    # NO fabricated edges: with no --leakage the tree is drawn on its own.
    raw = load_leakage(args.leakage) if args.leakage else []
    edges = []
    dropped = 0
    for q, t, w in raw:
        if q in sp2i and t in sp2i:
            edges.append({"s": sp2i[q], "d": sp2i[t], "w": w})
        else:
            dropped += 1

    tree_svg = build_tree_svg(root, pos, angle)
    focus_i = sp2i.get(args.focus, 0)
    empty = not edges

    data = {"tips": tips, "edges": edges, "focus": focus_i, "R": R, "demo": empty}
    banner = ("<div class='demo'>No leakage data loaded &mdash; tree only. Pass "
              "<code>--leakage edges.tsv</code> (columns: query_species&#9;"
              "target_species&#9;weight, derived from crossmap.m8) to draw real "
              "cross-mapping arcs.</div>") if empty else ""

    out = TEMPLATE.replace("__TREE_SVG__", tree_svg) \
                  .replace("__DATA__", json.dumps(data)) \
                  .replace("__BANNER__", banner) \
                  .replace("__NTIPS__", str(len(tips))) \
                  .replace("__FOCUS__", html.escape(args.focus))
    with open(args.out, "w") as fh:
        fh.write(out)
    print(f"leakage_tree: {len(tips)} tips, {len(edges)} leakage edges "
          f"({dropped} dropped: species not in tree)"
          f"{' [no leakage loaded]' if empty else ''} -> {args.out}")


TEMPLATE = r"""<!doctype html><html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>meta-mage leakage vs phylogeny</title><style>
:root{--bg:#0f1115;--panel:#171a21;--fg:#e6e9ef;--dim:#9aa4b2;--line:#2a2f3a;
  --br:#3a4150;--out:#e0685a;--in:#5aa9e6;--focus:#e0b341;--ok:#5ac48a;}
@media (prefers-color-scheme:light){:root{--bg:#f6f7f9;--panel:#fff;--fg:#1a1d23;
  --dim:#5a6472;--line:#e2e5ea;--br:#c2c8d2;}}
*{box-sizing:border-box}
body{margin:0;background:var(--bg);color:var(--fg);
  font:14px/1.5 -apple-system,BlinkMacSystemFont,"Segoe UI",Roboto,sans-serif}
.wrap{display:flex;gap:18px;max-width:1200px;margin:0 auto;padding:22px 20px 60px;
  flex-wrap:wrap}
h1{font-size:20px;margin:0 0 2px}
.meta{color:var(--dim);font-size:13px;margin-bottom:12px}
.demo{background:color-mix(in srgb,var(--focus) 18%,transparent);
  border:1px solid var(--focus);color:var(--fg);border-radius:8px;padding:8px 11px;
  font-size:12.5px;margin-bottom:12px}
.left{flex:1 1 620px;min-width:520px}
.side{flex:0 0 300px;min-width:260px}
svg{width:100%;height:auto;display:block}
.br{stroke:var(--br);stroke-width:.7;fill:none}
.tipdot{fill:var(--dim)}
.tiplbl{font-size:6px;fill:var(--dim);opacity:.35}
.chord{fill:none;stroke-linecap:round}
.chord.out{stroke:var(--out)}
.chord.in{stroke:var(--in)}
.focusdot{fill:var(--focus);stroke:var(--bg);stroke-width:1.2}
.hl{opacity:1!important;font-weight:600}
.hl.o{fill:var(--out)} .hl.i{fill:var(--in)} .hl.b{fill:var(--focus)}
.panel{background:var(--panel);border:1px solid var(--line);border-radius:10px;
  padding:13px 15px;margin-bottom:14px}
.panel h2{font-size:13px;margin:0 0 8px;color:var(--dim);font-weight:600;
  text-transform:uppercase;letter-spacing:.03em}
input[type=text]{width:100%;padding:7px 9px;border:1px solid var(--line);
  border-radius:7px;background:var(--bg);color:var(--fg);font-size:13px}
.list{max-height:220px;overflow:auto;margin-top:8px}
.row{display:flex;justify-content:space-between;gap:8px;padding:3px 6px;
  border-radius:6px;font-size:12.5px;cursor:pointer}
.row:hover{background:var(--bg)}
.row.o .n{color:var(--out)} .row.i .n{color:var(--in)}
.w{color:var(--dim);font-variant-numeric:tabular-nums}
.legend{font-size:12px;color:var(--dim);margin-top:6px}
.sw{display:inline-block;width:10px;height:10px;border-radius:2px;
  vertical-align:middle;margin-right:4px}
.stat{font-size:12.5px;color:var(--fg);margin:2px 0}
b.f{color:var(--focus)}
</style></head><body><div class="wrap">
<div class="left">
<h1>meta-mage &mdash; leakage vs phylogeny</h1>
<div class="meta">g__Bacteroides &middot; __NTIPS__ species (GTDB r232 bac120) &middot;
focus: <b class="f" id="fname">__FOCUS__</b></div>
__BANNER__
<svg id="viz" viewBox="-420 -420 840 840" role="img">
  <defs>
    <marker id="arrOut" viewBox="0 0 10 10" refX="8.5" refY="5" markerWidth="4.2"
      markerHeight="4.2" orient="auto"><path d="M0,1 L9,5 L0,9 z" fill="var(--out)"/></marker>
    <marker id="arrIn" viewBox="0 0 10 10" refX="8.5" refY="5" markerWidth="4.2"
      markerHeight="4.2" orient="auto"><path d="M0,1 L9,5 L0,9 z" fill="var(--in)"/></marker>
  </defs>
  <g id="tree">__TREE_SVG__</g>
  <g id="chords"></g>
  <g id="dots"></g>
  <g id="labels"></g>
</svg>
</div>
<div class="side">
  <div class="panel">
    <h2>Focus taxon</h2>
    <input type="text" id="search" placeholder="type a species name...">
    <div class="list" id="picker"></div>
  </div>
  <div class="panel">
    <h2>Leakage for focus</h2>
    <div class="legend"><span class="sw" style="background:var(--out)"></span>
      outgoing (focus&rarr; other) &nbsp;
      <span class="sw" style="background:var(--in)"></span> incoming (other&rarr;focus)</div>
    <div id="summary" class="stat"></div>
    <div class="list" id="edges"></div>
  </div>
</div>
</div>
<script>
const D = __DATA__;
const tips = D.tips, edges = D.edges;
const byName = {}; tips.forEach(t=>byName[t.sp]=t);
let focus = D.focus;
const maxW = Math.max(1, ...edges.map(e=>e.w));
const NS="http://www.w3.org/2000/svg";
function el(t,a){const e=document.createElementNS(NS,t);for(const k in a)e.setAttribute(k,a[k]);return e;}

// static tip dots + faint labels
const dots=document.getElementById('dots'), labels=document.getElementById('labels');
tips.forEach(t=>{
  dots.appendChild(el('circle',{cx:t.x,cy:t.y,r:1.4,class:'tipdot','data-i':t.i}));
  const deg=t.a*180/Math.PI, flip=(t.a>Math.PI/2&&t.a<3*Math.PI/2);
  const lx=(D.R+6)*Math.cos(t.a), ly=(D.R+6)*Math.sin(t.a);
  const lab=el('text',{x:lx,y:ly,class:'tiplbl','data-i':t.i,
    transform:`rotate(${flip?deg+180:deg} ${lx} ${ly})`,
    'text-anchor':flip?'end':'start','dominant-baseline':'middle'});
  lab.textContent=t.sp.replace('Bacteroides ','B. ');
  labels.appendChild(lab);
});

function chordPath(a,b,dir){
  // quadratic bezier bowing toward centre; nearer angles bow less. `dir` (+1 out,
  // -1 in) shifts the control point perpendicular to the chord so the two
  // directions of the SAME pair bow to opposite sides and never overlap.
  const A=tips[a], B=tips[b];
  let da=Math.abs(A.a-B.a); if(da>Math.PI) da=2*Math.PI-da;
  const pull=0.08+0.42*(1-da/Math.PI); // distant partners bow deeper
  const ma=(A.a+B.a)/2, adj=(Math.abs(A.a-B.a)>Math.PI)?Math.PI:0;
  const cr=D.R*pull;
  let cx=cr*Math.cos(ma+adj), cy=cr*Math.sin(ma+adj);
  const dx=B.x-A.x, dy=B.y-A.y, len=Math.hypot(dx,dy)||1;
  const off=dir*(14+0.06*len);        // perpendicular split between in/out
  cx+=(-dy/len)*off; cy+=(dx/len)*off;
  return `M${A.x},${A.y} Q${cx.toFixed(1)},${cy.toFixed(1)} ${B.x},${B.y}`;
}

function render(){
  document.getElementById('fname').textContent=tips[focus].sp;
  const cg=document.getElementById('chords'); cg.innerHTML='';
  document.querySelectorAll('.tiplbl').forEach(l=>{l.setAttribute('class','tiplbl');});
  const outs=edges.filter(e=>e.s===focus), ins=edges.filter(e=>e.d===focus);
  const partners=new Set();
  function draw(list,cls,dir){
    const sign=dir==='out'?1:-1;
    list.forEach(e=>{
      const other=dir==='out'?e.d:e.s;
      partners.add(other);
      const p=el('path',{d:chordPath(e.s,e.d,sign),class:'chord '+cls,
        'stroke-width':(0.8+3.6*e.w/maxW).toFixed(2),
        'stroke-opacity':(0.30+0.55*e.w/maxW).toFixed(2),
        'marker-end':dir==='out'?'url(#arrOut)':'url(#arrIn)'});
      cg.appendChild(p);
    });
  }
  draw(outs,'out','out'); draw(ins,'in','in');
  // focus + partner dots/labels
  const fd=el('circle',{cx:tips[focus].x,cy:tips[focus].y,r:4,class:'focusdot'});
  cg.appendChild(fd);
  document.querySelectorAll('.tiplbl').forEach(l=>{
    const i=+l.getAttribute('data-i');
    if(i===focus) l.setAttribute('class','tiplbl hl b');
    else if(partners.has(i)){
      const o=outs.some(e=>e.d===i), inc=ins.some(e=>e.s===i);
      l.setAttribute('class','tiplbl hl '+(o&&inc?'b':o?'o':'i'));
    }
  });
  // side summary
  const sumW=v=>v.reduce((s,e)=>s+e.w,0);
  document.getElementById('summary').innerHTML=
    `<b class="f">${tips[focus].sp.replace('Bacteroides','B.')}</b> &middot; `+
    `out ${outs.length} edges / ${sumW(outs)} &middot; in ${ins.length} / ${sumW(ins)}`;
  const eb=document.getElementById('edges'); eb.innerHTML='';
  const rows=[...outs.map(e=>({o:e.d,w:e.w,d:'o'})),
              ...ins.map(e=>({o:e.s,w:e.w,d:'i'}))].sort((a,b)=>b.w-a.w);
  rows.forEach(r=>{
    const row=el('div',{class:'row '+r.d});
    row.innerHTML=`<span class="n">${r.d==='o'?'&rarr; ':'&larr; '}`+
      `${tips[r.o].sp.replace('Bacteroides','B.')}</span><span class="w">${r.w}</span>`;
    row.onclick=()=>{focus=r.o;render();buildPicker(document.getElementById('search').value);};
    eb.appendChild(row);
  });
}

function buildPicker(q){
  const p=document.getElementById('picker'); p.innerHTML='';
  q=(q||'').toLowerCase();
  tips.filter(t=>t.sp.toLowerCase().includes(q)).slice(0,60).forEach(t=>{
    const row=el('div',{class:'row'+(t.i===focus?' o':'')});
    row.innerHTML=`<span class="n">${t.sp.replace('Bacteroides','B.')}</span>`;
    row.onclick=()=>{focus=t.i;render();};
    p.appendChild(row);
  });
}
document.getElementById('search').addEventListener('input',e=>buildPicker(e.target.value));
buildPicker(''); render();
</script></body></html>"""


if __name__ == "__main__":
    main()

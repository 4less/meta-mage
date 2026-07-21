#!/usr/bin/env python3
"""
Phylogeny-guided recursive merge probe. Instead of greedily absorbing the highest
core-overlap neighbours (merge_gain.py -- which produces polyphyletic, taxonomy-
breaking merges), this walks the bac120 tree UPWARD from a focus species and, at
every ancestral node, treats the whole monophyletic clade under that node as one
merged taxon. Every step is a real clade, so the merge always follows the phylogeny.

At each node it reports:
  * clade size (run species merged) and the sister group newly absorbed
  * final markers recomputed for the merged clade (same filters as the pipeline:
    prevalence -> score+cap -> length+cross-map guard with masking re-applied), so
    you can see WHERE climbing the tree makes core markers reappear
  * % of the focus's outgoing / incoming leakage now internalised (partners that
    have moved inside the clade), by hits and by distinct partner species

Marker recomputation is reused verbatim from merge_gain.py; the tree parser from
leakage_tree.py. Needs the run's counts.tsv (species rows), manifest, specificity
report and mask_intervals -- the same inputs merge_gain uses.
"""
import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import merge_gain as mg          # noqa: E402
import leakage_tree as lt        # noqa: E402


def load_acc2sp(path):
    a2s = {}
    with open(path) as fh:
        fh.readline()
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) >= 2:
                a2s[f[0]] = f[1]
    return a2s


def load_leakage(path):
    """query_species, target_species, genes, hits."""
    out, inc = {}, {}
    outh, inh = {}, {}
    with open(path) as fh:
        fh.readline()
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 4:
                continue
            q, t, genes, hits = f[0], f[1], int(f[2]), int(f[3])
            out.setdefault(q, {})[t] = genes
            outh.setdefault(q, {})[t] = hits
            inc.setdefault(t, {})[q] = genes
            inh.setdefault(t, {})[q] = hits
    return out, outh, inc, inh


def clade_species_under(node, tip_sp, run_species):
    """run-species at the leaves under `node` (deduped)."""
    out = set()
    stack = [node]
    while stack:
        n = stack.pop()
        if not n.children:
            sp = tip_sp.get(id(n))
            if sp in run_species:
                out.add(sp)
        else:
            stack.extend(n.children)
    return out


def pct(part, whole):
    return 0.0 if whole <= 0 else 100.0 * part / whole


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--tree", required=True)
    ap.add_argument("--species", required=True, help="accession<TAB>species map")
    ap.add_argument("--leakage", required=True, help="leakage_edges.tsv")
    ap.add_argument("--counts", required=True, help="counts.tsv (species rows)")
    ap.add_argument("--manifest", required=True)
    ap.add_argument("--specificity", required=True)
    ap.add_argument("--mask_intervals", required=True)
    ap.add_argument("--focus", default="Bacteroides ovatus")
    ap.add_argument("--threshold", type=int, default=50)
    ap.add_argument("--min_in", type=float, default=0.80)
    ap.add_argument("--max_out", type=float, default=0.02)
    ap.add_argument("--min_clade_size", type=int, default=3)
    ap.add_argument("--max_per_clade", type=int, default=100)
    ap.add_argument("--score_out_exp", type=float, default=1.0)
    ap.add_argument("--recovery_min_clean", type=int, default=300)
    ap.add_argument("--recovery_max_masked_frac", type=float, default=0.5)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    p = dict(min_in=args.min_in, max_out=args.max_out,
             min_clade_size=args.min_clade_size, max_per_clade=args.max_per_clade,
             score_out_exp=args.score_out_exp,
             recovery_min_clean=args.recovery_min_clean,
             recovery_max_masked_frac=args.recovery_max_masked_frac)

    # ---- shared marker machinery (from merge_gain) ----
    size, score, mani_genus, n_total = mg.load_manifest(args.manifest)
    run_species = {s for s in size if s.startswith("s__")} or set(size)
    guard = mg.load_guard(args.specificity)
    masks = mg.load_mask_intervals(args.mask_intervals)
    # counts clades are the s__-prefixed species; manifest species are too.
    per, total = mg.load_counts(args.counts, set(size))

    # ---- tree ----
    root = lt.parse_newick(open(args.tree).read())
    acc2sp = load_acc2sp(args.species)
    # parent links + tip->species (species may carry no s__ in the map)
    tip_sp = {}
    def link(n, parent):
        n.parent = parent
        if not n.children:
            tip_sp[id(n)] = acc2sp.get(n.name, n.name)
        for c in n.children:
            link(c, n)
    link(root, None)

    # focus species -> both spellings (map has no s__, counts/manifest have s__)
    focus_plain = args.focus.replace("s__", "")
    focus_s = "s__" + focus_plain
    # which manifest spelling is the run using?
    run_focus = focus_s if focus_s in size else (focus_plain if focus_plain in size else None)

    # locate the focus tip
    focus_tip = next((n for n in tip_sp if tip_sp[n] == focus_plain), None)
    if focus_tip is None:
        sys.exit(f"focus {args.focus} not a tip in the tree")
    # map manifest/counts species spelling for clade membership
    def to_run(sp_plain):
        s = "s__" + sp_plain
        return s if s in size else (sp_plain if sp_plain in size else None)

    run_species = set(size)
    # walk focus -> root
    node = next(n for n in _iter_nodes(root) if id(n) == focus_tip)

    out, outh, inc, inh = load_leakage(args.leakage)
    # leakage_edges.tsv uses PLAIN species names (no s__) -- key with focus_plain.
    o_edges = out.get(focus_plain, {})
    o_h = outh.get(focus_plain, {})
    i_edges = inc.get(focus_plain, {})
    i_h = inh.get(focus_plain, {})
    # leakage keys are plain species; strip s__ for comparison
    def plain(s):
        return s[3:] if s.startswith("s__") else s
    out_tot_h = sum(o_h.values()); in_tot_h = sum(i_h.values())
    out_tot_n = len(o_edges); in_tot_n = len(i_edges)

    rows = []
    prev_clade = set()
    cur = node
    step = 0
    while cur is not None:
        clade_plain = {plain(s) for s in _clade_run(cur, tip_sp, run_species, to_run)}
        # skip ancestral nodes that add no new run-species (same merged taxon)
        if clade_plain == prev_clade and cur.parent is not None:
            cur = cur.parent
            continue
        clade_run = {to_run(s) for s in clade_plain}
        clade_run = {s for s in clade_run if s}
        markers, _ = mg.marker_set(sorted(clade_run), per, total, size, score,
                                   n_total, p, guard, masks)
        added = sorted(clade_plain - prev_clade)
        # leakage internalised: partners now inside the clade (excl. focus)
        inside = clade_plain - {focus_plain}
        out_h_in = sum(h for t, h in o_h.items() if plain(t) in inside)
        in_h_in = sum(h for q, h in i_h.items() if plain(q) in inside)
        out_n_in = sum(1 for t in o_edges if plain(t) in inside)
        in_n_in = sum(1 for q in i_edges if plain(q) in inside)
        rows.append({
            "step": step,
            "clade_species": len(clade_run),
            "added": ";".join(a.replace("Bacteroides ", "B. ") for a in added) or "-",
            "markers": len(markers),
            "out_leak_removed_pct": round(pct(out_h_in, out_tot_h), 1),
            "in_leak_removed_pct": round(pct(in_h_in, in_tot_h), 1),
            "out_partners_absorbed": f"{out_n_in}/{out_tot_n}",
            "in_partners_absorbed": f"{in_n_in}/{in_tot_n}",
        })
        prev_clade = clade_plain
        if len(clade_run) >= len(run_species):
            break
        cur = cur.parent
        step += 1

    cols = ["step", "clade_species", "markers", "out_leak_removed_pct",
            "in_leak_removed_pct", "out_partners_absorbed", "in_partners_absorbed",
            "added"]
    with open(args.out, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(str(r[c]) for c in cols) + "\n")

    print(f"phylo_merge: focus={args.focus}  ({out_tot_n} out / {in_tot_n} in "
          f"leakage partners; {out_tot_h}/{in_tot_h} hits)")
    hdr = f"{'step':>4} {'clade':>5} {'markers':>7} {'out%':>6} {'in%':>6}  climbed"
    print(hdr); print("-" * len(hdr))
    for r in rows:
        print(f"{r['step']:>4} {r['clade_species']:>5} {r['markers']:>7} "
              f"{r['out_leak_removed_pct']:>6} {r['in_leak_removed_pct']:>6}  "
              f"+{r['added'][:70]}")
    print(f"-> {args.out}")


def _iter_nodes(root):
    stack = [root]
    while stack:
        n = stack.pop()
        yield n
        stack.extend(n.children)


def _clade_run(node, tip_sp, run_species, to_run):
    out = set()
    stack = [node]
    while stack:
        n = stack.pop()
        if not n.children:
            sp = tip_sp.get(id(n))
            if sp and to_run(sp):
                out.add(sp)
        else:
            stack.extend(n.children)
    return out


if __name__ == "__main__":
    main()

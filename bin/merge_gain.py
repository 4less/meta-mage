#!/usr/bin/env python3
"""
Marker-gain probe for low-marker species: does merging a species with its most
overlapping same-genus neighbours lift it to the requested marker count?

Evaluated AFTER ALL FILTERS, so the baseline equals the pipeline's real final
marker count (not an intermediate prevalence-only selection). For a clade C the
marker set is rebuilt exactly as the pipeline builds it:
    1. prevalence   in_prev_C = (size_C - (score_sum_C - in_score_C)) / size_C >= min_in
                    out_prev_C = (marker_total - in_count_C) / (N - size_C)     <= max_out
    2. score + cap  score = in_prev * (1 - out_prev) ** score_out_exp, keep top max_per_clade
    3. length+guard keep a cluster only if it survived the length filter AND the
                    cross-map guard -- read from specificity_report.tsv: a cluster
                    is kept for C if its recorded off-target clades all fall INSIDE
                    C (merging absorbs them), OR it was mask-recovered.
Steps 1-2 are recomputed from counts.tsv (merging is only relabelling -- no re-run);
step 3 uses the guard's per-marker off-target footprint, so a species that ended at
0 real markers shows 0 here too, and merging is credited only for markers it can
actually make clade-unique.

Greedy: absorb the not-yet-merged same-genus species sharing the most of the clade's
CORE genes. Stop once the marker count reaches --threshold (the minimum requested
count) -- the reported clade is the SMALLEST merge that clears the bar; if the bar
is never reached, the peak-marker clade is reported.

Output (merge_gain.tsv, one row per flagged species):
    species genus n_genomes baseline_markers merged_markers delta reached_threshold
    n_species_added merged_clade trajectory

Conservative on one point: a cluster that becomes a prevalence marker ONLY after
merging and was never emitted for any single member (so the guard never evaluated
it) can't be affirmed clade-unique without a re-run, and is not counted -- the
merged gain is a lower bound, never an over-count.
"""
import argparse
import os
import sys
from collections import defaultdict


def load_manifest(path):
    """species -> (size, score_sum, genus); N_total. score_sum=size if the
    manifest has no completeness column (legacy / unweighted)."""
    size = defaultdict(int)
    score = defaultdict(float)
    sp_genus = {}
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
            sp_genus[sp] = f[col["genus"]]
    return size, score, sp_genus, n_total


def load_flagged(path):
    """Ordered list of flagged (below-threshold) species and species -> genus."""
    flagged, sp_genus = [], {}
    with open(path) as fh:
        col = {n: i for i, n in enumerate(fh.readline().rstrip("\n").split("\t"))}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            sp = f[col["species"]]
            sp_genus[sp] = f[col["genus"]]
            if f[col["flagged"]] == "yes":
                flagged.append(sp)
    return flagged, sp_genus


def load_counts(path, wanted_species):
    """cluster -> {species: (in_count, in_score)} and cluster -> marker_total,
    restricted to species-rank rows whose clade is in wanted_species."""
    per = defaultdict(dict)
    total = {}
    with open(path) as fh:
        next(fh)
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 5 or f[1] != "species":
                continue
            clade = f[2]
            if clade not in wanted_species:
                continue
            cluster = f[0]
            in_count = int(f[3])
            marker_total = int(f[4])
            in_score = float(f[5]) if len(f) > 5 else float(in_count)
            per[cluster][clade] = (in_count, in_score)
            total[cluster] = marker_total
    return per, total


def load_guard(path):
    """From specificity_report.tsv (species rank): the cross-map guard footprint
    per cluster -- off_union (all off-target clades the cluster's copies hit),
    the set of mask-recovered clusters, and the set of clusters the guard actually
    evaluated (== emitted, length-passed). Empty/absent => guard not run."""
    off_union = defaultdict(set)
    recovered = set()
    queried = set()
    if not path or os.path.basename(path) == "NO_FILE" or not os.path.exists(path):
        return off_union, recovered, queried, False
    with open(path) as fh:
        head = fh.readline().rstrip("\n").split("\t")
        col = {n: i for i, n in enumerate(head)}
        if "cluster" not in col:
            return off_union, recovered, queried, False
        oc_i = col.get("offtarget_clades")
        rec_i = col.get("recovered")
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if col.get("rank") is not None and f[col["rank"]] != "species":
                continue
            cl = f[col["cluster"]]
            queried.add(cl)
            if oc_i is not None and oc_i < len(f) and f[oc_i]:
                off_union[cl] |= set(f[oc_i].split(";"))
            if rec_i is not None and rec_i < len(f) and f[rec_i] == "1":
                recovered.add(cl)
    return off_union, recovered, queried, True


def marker_set(clade_species, per, total, size, score, n_total, p, guard):
    """Final marker set for the merged clade, applying every pipeline filter:
    prevalence -> score+cap -> length+cross-map guard (off-targets inside the
    merged clade, or mask-recovered). Also returns the prevalence CORE set (used
    to pick the next neighbour)."""
    off_union, recovered, queried, guard_on = guard
    size_C = sum(size[s] for s in clade_species)
    score_C = sum(score[s] for s in clade_species)
    if size_C < p["min_clade_size"]:
        return set(), set()
    out_denom = n_total - size_C
    if out_denom <= 0:
        return set(), set()
    clade_set = set(clade_species)
    candidates = []       # (score_val, cluster) among prevalence survivors
    core = set()
    for cluster, spmap in per.items():
        inc = isc = 0.0
        hit = False
        for s in clade_species:
            v = spmap.get(s)
            if v:
                inc += v[0]
                isc += v[1]
                hit = True
        if not hit:
            continue
        in_prev = (size_C - (score_C - isc)) / size_C
        if in_prev > 1.0:
            in_prev = 1.0
        if in_prev >= p["min_in"]:
            core.add(cluster)
        if in_prev < p["min_in"]:
            continue
        out_prev = (total[cluster] - inc) / out_denom
        if out_prev > p["max_out"]:
            continue
        candidates.append((in_prev * (1.0 - out_prev) ** p["score_out_exp"], cluster))

    # SCORE's cap: keep the top max_per_clade by score BEFORE the guard (matching
    # the pipeline order), then apply length + cross-map guard.
    candidates.sort(key=lambda t: t[0], reverse=True)
    markers = set()
    for _, cl in candidates[: p["max_per_clade"]]:
        if guard_on:
            if cl not in queried:          # length-dropped or never emitted
                continue
            if not (cl in recovered or off_union.get(cl, set()) <= clade_set):
                continue                   # still cross-maps outside the merged clade
        markers.add(cl)
    return markers, core


def next_neighbour(clade_species, core, per, candidates):
    """Not-yet-merged candidate species sharing the most of the clade's CORE
    clusters (the genes whose out_prev it inflates). Greedy, one at a time."""
    best, best_ov = None, 0
    clade = set(clade_species)
    for b in candidates:
        if b in clade:
            continue
        ov = sum(1 for cl in core if b in per.get(cl, ()))
        if ov > best_ov:
            best, best_ov = b, ov
    return best, best_ov


def probe_species(sp, genus, genus_species, per, total, size, score, n_total, p,
                  guard, threshold, max_steps):
    """Greedy merge trajectory for one seed species. Returns a result dict."""
    clade = [sp]
    markers, core = marker_set(clade, per, total, size, score, n_total, p, guard)
    steps = [(sp, len(markers))]               # (seed-or-added species, markers)
    for _ in range(max_steps):
        if steps[-1][1] >= threshold:          # bar already cleared -- stop merging
            break
        nb, ov = next_neighbour(clade, core, per, genus_species)
        if not nb or ov == 0:                  # no overlapping species left
            break
        clade.append(nb)
        markers, core = marker_set(clade, per, total, size, score, n_total, p, guard)
        steps.append((nb, len(markers)))

    base = steps[0][1]
    # Prefer the SMALLEST clade that reaches the bar; else the peak-marker clade.
    reach_idx = next((i for i, (_, c) in enumerate(steps) if c >= threshold), None)
    if reach_idx is not None:
        rep, reached = reach_idx, "yes"
    else:
        rep = max(range(len(steps)), key=lambda i: steps[i][1])
        reached = "no"
    return {
        "species": sp,
        "genus": genus,
        "n_genomes": size.get(sp, 0),
        "baseline_markers": base,
        "merged_markers": steps[rep][1],
        "delta": steps[rep][1] - base,
        "reached_threshold": reached,
        "n_species_added": rep,
        "merged_clade": " + ".join(s for s, _ in steps[:rep + 1]),
        "trajectory": ">".join(str(c) for _, c in steps[:rep + 1]),
    }


def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--manifest", required=True)
    ap.add_argument("--counts", required=True, help="counts.tsv")
    ap.add_argument("--low_marker_species", required=True,
                    help="low_marker_species.tsv (flagged column)")
    ap.add_argument("--specificity", help="specificity_report.tsv (final-filter "
                    "guard; omit only if the cross-map guard did not run)")
    ap.add_argument("--threshold", type=int, default=50,
                    help="minimum requested markers; merging stops once reached")
    ap.add_argument("--min_in", type=float, default=0.80)
    ap.add_argument("--max_out", type=float, default=0.02)
    ap.add_argument("--min_clade_size", type=int, default=3)
    ap.add_argument("--max_per_clade", type=int, default=100)
    ap.add_argument("--score_out_exp", type=float, default=1.0)
    ap.add_argument("--max_steps", type=int, default=25)
    ap.add_argument("--outdir", default=".")
    args = ap.parse_args()

    p = dict(min_in=args.min_in, max_out=args.max_out,
             min_clade_size=args.min_clade_size, max_per_clade=args.max_per_clade,
             score_out_exp=args.score_out_exp)

    size, score, mani_genus, n_total = load_manifest(args.manifest)
    flagged, flag_genus = load_flagged(args.low_marker_species)
    guard = load_guard(args.specificity)
    if not guard[3]:
        sys.stderr.write("[merge_gain] WARNING: no specificity_report -- marker "
                         "counts are prevalence-only, NOT after the cross-map guard\n")

    cols = ["species", "genus", "n_genomes", "baseline_markers", "merged_markers",
            "delta", "reached_threshold", "n_species_added", "merged_clade",
            "trajectory"]
    out_path = os.path.join(args.outdir, "merge_gain.tsv")

    if not flagged:
        with open(out_path, "w") as fh:
            fh.write("\t".join(cols) + "\n")
        sys.stderr.write("[merge_gain] no flagged species -> empty merge_gain.tsv\n")
        return

    target_genera = {mani_genus.get(s, flag_genus.get(s)) for s in flagged}
    wanted = {s for s, g in mani_genus.items() if g in target_genera}
    sys.stderr.write(
        f"[merge_gain] {len(flagged)} flagged species in {len(target_genera)} "
        f"genera; loading counts for {len(wanted)} species...\n")
    per, total = load_counts(args.counts, wanted)

    rows = []
    for sp in flagged:
        genus = mani_genus.get(sp, flag_genus.get(sp))
        genus_species = {s for s, g in mani_genus.items() if g == genus}
        rows.append(probe_species(sp, genus, genus_species, per, total, size,
                                  score, n_total, p, guard, args.threshold,
                                  args.max_steps))

    rows.sort(key=lambda r: (r["reached_threshold"] != "yes", -r["delta"]))
    with open(out_path, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(str(r[c]) for c in cols) + "\n")

    n_reach = sum(1 for r in rows if r["reached_threshold"] == "yes")
    n_help = sum(1 for r in rows if r["delta"] > 0)
    sys.stderr.write(
        f"[merge_gain] wrote {out_path}: {len(rows)} species, {n_help} gain markers "
        f"by merging, {n_reach} reach the {args.threshold}-marker bar\n")


if __name__ == "__main__":
    main()

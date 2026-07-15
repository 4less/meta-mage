#!/usr/bin/env python3
"""
ANI species gap for low-marker species (aggregates a skani sparse table).

Given skani's all-vs-all sparse output (produced upstream, in the skani container)
and the manifest, this splits the pairwise ANI of each flagged species' genus into:
  * within  -- both genomes belong to the flagged species
  * between -- one genome is the flagged species, the other a sibling species
and answers the re-merge question: is the within-species ANI always tighter than
the between-species ANI (a clean gap), or do they overlap (the species may be one
taxon split in two)?

This step is pure stdlib-python (no skani here): the skani call lives in its own
process because its container has no python. Only same-genus pairs are kept.

Outputs (under --outdir):
    ani_pairs.tsv       -- focal_species, other_species, kind, ani, align_frac
                           (long form: every same-genus pair touching a flagged sp.)
    ani_gap_summary.tsv -- per flagged species: within min/median, nearest sibling
                           + its max between-ANI, the gap, overlap?, merge_candidate?
"""
import argparse
import os
import statistics
import sys
from collections import defaultdict


def load_manifest(path):
    """genome_id -> (genus, species); basename -> genome_id (for the sparse table)."""
    gid_info, name2gid = {}, {}
    with open(path) as fh:
        col = {n: i for i, n in enumerate(fh.readline().rstrip("\n").split("\t"))}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            genus, species = f[col["genus"]], f[col["species"]]
            if genus == "NA" or species == "NA":
                continue
            gid = f[col["genome_id"]]
            p = f[col["path"]]
            gid_info[gid] = (genus, species)
            name2gid[os.path.basename(p)] = gid
            name2gid[p] = gid
    return gid_info, name2gid


def load_flagged(path):
    """species set with flagged == yes, plus species -> genus."""
    flagged, sp_genus = set(), {}
    with open(path) as fh:
        col = {n: i for i, n in enumerate(fh.readline().rstrip("\n").split("\t"))}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            sp_genus[f[col["species"]]] = f[col["genus"]]
            if f[col["flagged"]] == "yes":
                flagged.add(f[col["species"]])
    return flagged, sp_genus


def lookup(name2gid, field):
    return (name2gid.get(field) or name2gid.get(os.path.basename(field))
            or name2gid.get(os.path.abspath(field)))


def parse_sparse(path, name2gid):
    """Yield (gid_a, gid_b, ani, af_max) from a skani -E sparse table.

    Columns: Ref_file  Query_file  ANI  Align_fraction_ref  Align_fraction_query
    (optional leading comment / header line). File fields are the genome paths we
    handed skani, mapped back to genome ids by basename.
    """
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith(("Ref_file", "#")):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 5:
                continue
            a, b = lookup(name2gid, f[0]), lookup(name2gid, f[1])
            if not a or not b:
                continue
            try:
                ani = float(f[2])
                af = max(float(f[3]), float(f[4])) / 100.0
            except ValueError:
                continue
            yield a, b, ani, af


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--manifest", required=True)
    ap.add_argument("--low_marker_species", required=True,
                    help="low_marker_species.tsv from low_marker_masking.py")
    ap.add_argument("--sparse", required=True, help="skani -E sparse table")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--ani", type=float, default=95.0,
                    help="between-ANI (%%) flagged as a merge candidate (default 95)")
    ap.add_argument("--min_af", type=float, default=0.5,
                    help="min align fraction for a pair to count (default 0.5)")
    args = ap.parse_args()

    gid_info, name2gid = load_manifest(args.manifest)
    flagged, sp_genus = load_flagged(args.low_marker_species)
    os.makedirs(args.outdir, exist_ok=True)

    pairs_path = os.path.join(args.outdir, "ani_pairs.tsv")
    gap_path = os.path.join(args.outdir, "ani_gap_summary.tsv")

    within = defaultdict(list)                       # focal sp -> [ani...]
    between = defaultdict(lambda: defaultdict(list))  # focal sp -> {other sp: [ani]}
    pair_rows = []
    for a, b, ani, af in parse_sparse(args.sparse, name2gid):
        if af < args.min_af:
            continue
        (ga, sa), (gb, sb) = gid_info[a], gid_info[b]
        if ga != gb:
            continue                                  # only within-genus pairs
        a_flag, b_flag = sa in flagged, sb in flagged
        if not (a_flag or b_flag):
            continue
        if sa == sb:
            if a_flag:
                within[sa].append(ani)
                pair_rows.append((sa, sb, "within", ani, af))
        else:
            if a_flag:
                between[sa][sb].append(ani)
                pair_rows.append((sa, sb, "between", ani, af))
            if b_flag:
                between[sb][sa].append(ani)
                pair_rows.append((sb, sa, "between", ani, af))

    with open(pairs_path, "w") as fh:
        fh.write("focal_species\tother_species\tkind\tani\talign_frac\n")
        for r in pair_rows:
            fh.write(f"{r[0]}\t{r[1]}\t{r[2]}\t{r[3]:.4f}\t{r[4]:.4f}\n")

    def fmt(x):
        return f"{x:.3f}" if x == x else "NA"        # NaN-safe (NaN != NaN)

    with open(gap_path, "w") as fh:
        fh.write("species\tgenus\tn_within\tmin_within\tmedian_within\tn_between\t"
                 "nearest_species\tmax_between\tgap\toverlap\tmerge_candidate\n")
        for sp in sorted(flagged):
            genus = sp_genus.get(sp, "NA")
            w = within.get(sp, [])
            btw = between.get(sp, {})
            n_within = len(w)
            min_w = min(w) if w else float("nan")
            med_w = statistics.median(w) if w else float("nan")
            n_between = sum(len(v) for v in btw.values())
            nearest, max_b = "", float("nan")
            if btw:
                nearest = max(btw, key=lambda s: max(btw[s]))
                max_b = max(btw[nearest])
            gap = (min_w - max_b) if (w and btw) else float("nan")
            overlap = ("yes" if (w and btw and max_b >= min_w) else
                       "no" if (w and btw) else "NA")
            merge = "yes" if (btw and max_b >= args.ani) else "no"
            fh.write(f"{sp}\t{genus}\t{n_within}\t{fmt(min_w)}\t{fmt(med_w)}\t"
                     f"{n_between}\t{nearest or 'NA'}\t{fmt(max_b)}\t{fmt(gap)}\t"
                     f"{overlap}\t{merge}\n")

    sys.stderr.write(f"[ani] {len(flagged)} flagged species, {len(pair_rows)} "
                     f"same-genus pairs. -> {gap_path}\n")


if __name__ == "__main__":
    main()

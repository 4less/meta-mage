#!/usr/bin/env python3
"""
ANI species gap for low-marker species (skani within- vs between-species).

For every flagged species (from low_marker_species.tsv) this computes, within the
species' genus, all pairwise genome ANI with skani and splits the pairs into:
  * within  -- both genomes belong to the flagged species
  * between -- one genome is the flagged species, the other a sibling species
and asks the re-merge question directly: is the within-species ANI always tighter
than the between-species ANI (a clean gap), or do they overlap (the species may be
one taxon split in two)?

skani triangle is run ONCE per genus that contains a flagged species (so a genus
with several flagged species is sketched once), reading gzipped FASTA directly.

Outputs (under --outdir):
    ani_pairs.tsv       -- focal_species, other_species, kind, ani, align_frac
                           (long form: every pair touching a flagged species)
    ani_gap_summary.tsv -- per flagged species: within min/median, nearest sibling
                           + its max between-ANI, the gap, overlap?, merge_candidate?

Requires: skani on PATH (or --skani /path/to/skani).
"""
import argparse
import os
import shutil
import statistics
import subprocess
import sys
from collections import defaultdict


def resolve_path(p, genome_dir):
    """Prefer a Nextflow-staged copy in genome_dir (by basename); else the
    manifest's original absolute path (works when the storage is bind-mounted)."""
    if genome_dir:
        local = os.path.join(genome_dir, os.path.basename(p))
        if os.path.exists(local):
            return local
    return p


def load_manifest(path, genome_dir):
    """genome_id -> (genus, species, fasta_path); also path -> genome_id."""
    gid_info, path2gid = {}, {}
    genus_members = defaultdict(list)          # genus -> [(gid, resolved_path)]
    with open(path) as fh:
        col = {n: i for i, n in enumerate(fh.readline().rstrip("\n").split("\t"))}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            genus, species = f[col["genus"]], f[col["species"]]
            if genus == "NA" or species == "NA":
                continue
            gid = f[col["genome_id"]]
            p = resolve_path(f[col["path"]], genome_dir)
            gid_info[gid] = (genus, species, p)
            path2gid[p] = gid
            path2gid[os.path.abspath(p)] = gid
            path2gid[os.path.basename(p)] = gid
            genus_members[genus].append((gid, p))
    return gid_info, path2gid, genus_members


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


def run_skani(skani, genome_list, out_sparse, threads, min_af, screen):
    subprocess.run(
        [skani, "triangle", "-l", genome_list, "-o", out_sparse, "-t", str(threads),
         "-E", "--min-af", str(min_af * 100.0), "-s", str(screen)],
        check=True, stdout=subprocess.DEVNULL)


def parse_sparse(path, path2gid):
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith(("Ref_file", "#")):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 5:
                continue
            a = (path2gid.get(f[0]) or path2gid.get(os.path.abspath(f[0]))
                 or path2gid.get(os.path.basename(f[0])))
            b = (path2gid.get(f[1]) or path2gid.get(os.path.abspath(f[1]))
                 or path2gid.get(os.path.basename(f[1])))
            if not a or not b:
                continue
            yield a, b, float(f[2]), max(float(f[3]), float(f[4])) / 100.0


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--manifest", required=True)
    ap.add_argument("--low_marker_species", required=True,
                    help="low_marker_species.tsv from low_marker_masking.py")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--genome_dir", default=".",
                    help="dir with (staged) genome files, matched by basename; "
                         "falls back to the manifest's absolute path")
    ap.add_argument("--ani", type=float, default=95.0,
                    help="between-ANI (%%) flagged as a merge candidate (default 95)")
    ap.add_argument("--min_af", type=float, default=0.5)
    ap.add_argument("--screen", type=float, default=80.0)
    ap.add_argument("--threads", type=int, default=os.cpu_count() or 4)
    ap.add_argument("--skani", default="skani")
    args = ap.parse_args()

    if not shutil.which(args.skani):
        sys.exit(f"error: skani not found ('{args.skani}'). Install it or pass --skani.")
    gid_info, path2gid, genus_members = load_manifest(args.manifest, args.genome_dir)
    flagged, sp_genus = load_flagged(args.low_marker_species)
    os.makedirs(args.outdir, exist_ok=True)

    # Genera that contain at least one flagged species (sketched once each).
    genera = {sp_genus[sp] for sp in flagged if sp_genus.get(sp, "NA") != "NA"}
    if not flagged:
        sys.stderr.write("[ani] no flagged species; nothing to do.\n")
        open(os.path.join(args.outdir, "ani_pairs.tsv"), "w").write(
            "focal_species\tother_species\tkind\tani\talign_frac\n")
        open(os.path.join(args.outdir, "ani_gap_summary.tsv"), "w").write(
            "species\tgenus\tn_within\tmin_within\tmedian_within\tn_between\t"
            "nearest_species\tmax_between\tgap\toverlap\tmerge_candidate\n")
        return

    # focal species -> within [ani...]; focal species -> {other_sp: [ani...]}
    within = defaultdict(list)
    between = defaultdict(lambda: defaultdict(list))
    pair_rows = []
    for genus in sorted(genera):
        members = genus_members.get(genus, [])
        if len(members) < 2:
            continue
        gdir = os.path.join(args.outdir, genus.replace(" ", "_"))
        os.makedirs(gdir, exist_ok=True)
        list_path = os.path.join(gdir, "genomes.txt")
        with open(list_path, "w") as fh:
            for _gid, p in members:
                fh.write(p + "\n")
        sparse = os.path.join(gdir, "skani.sparse.tsv")
        sys.stderr.write(f"[ani] {genus}: {len(members)} genomes\n")
        run_skani(args.skani, list_path, sparse, args.threads, args.min_af, args.screen)
        for a, b, ani, af in parse_sparse(sparse, path2gid):
            sa, sb = gid_info[a][1], gid_info[b][1]
            a_flag, b_flag = sa in flagged, sb in flagged
            if not (a_flag or b_flag):
                continue
            if sa == sb:
                within[sa].append(ani)
                pair_rows.append((sa, sb, "within", ani, af))
            else:
                if a_flag:
                    between[sa][sb].append(ani)
                    pair_rows.append((sa, sb, "between", ani, af))
                if b_flag:
                    between[sb][sa].append(ani)
                    pair_rows.append((sb, sa, "between", ani, af))

    with open(os.path.join(args.outdir, "ani_pairs.tsv"), "w") as fh:
        fh.write("focal_species\tother_species\tkind\tani\talign_frac\n")
        for r in pair_rows:
            fh.write(f"{r[0]}\t{r[1]}\t{r[2]}\t{r[3]:.4f}\t{r[4]:.4f}\n")

    def fmt(x):
        return f"{x:.3f}" if x == x else "NA"       # NaN-safe (NaN != NaN)

    with open(os.path.join(args.outdir, "ani_gap_summary.tsv"), "w") as fh:
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
            overlap = "yes" if (w and btw and max_b >= min_w) else (
                "no" if (w and btw) else "NA")
            merge = "yes" if (btw and max_b >= args.ani) else "no"
            fh.write(f"{sp}\t{genus}\t{n_within}\t{fmt(min_w)}\t{fmt(med_w)}\t"
                     f"{n_between}\t{nearest or 'NA'}\t{fmt(max_b)}\t{fmt(gap)}\t"
                     f"{overlap}\t{merge}\n")

    sys.stderr.write(f"[ani] {len(flagged)} flagged species across "
                     f"{len(genera)} genera. -> {args.outdir}/ani_gap_summary.tsv\n")


if __name__ == "__main__":
    main()

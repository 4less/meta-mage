#!/usr/bin/env python3
"""
Per-genus pairwise genome distances (skani) -- exploration for species re-merging.

Some species that survive the genome-count prefilter still yield zero markers
because every core gene is shared with a sibling in the same genus (out_prevalence
> max_out for the whole core). Whether two such species should really be one
taxon is a genome-distance question, not a gene-prevalence one, so this tool
computes all pairwise ANI within each genus with skani and rolls the genome-level
pairs up to a species x species view of merge candidates.

For every genus with >= 2 genomes it writes, under --outdir/<genus>/:
    genomes.txt          -- the genome FASTA paths fed to skani
    skani.sparse.tsv     -- raw skani all-vs-all edges (ANI + align fractions)
And across all genera:
    merge_candidates.tsv -- cross-species genome pairs at >= --ani / >= --min_af,
                            aggregated per (species_a, species_b): n pairs, mean
                            and max ANI, max align fraction. Sorted worst-first;
                            these are the pairs worth eyeballing for a merge.

skani is fast and reads gzipped FASTA directly. This does not modify the pipeline;
it is a standalone probe you point at an existing run's manifest.

Requires: skani on PATH (or --skani /path/to/skani).
"""
import argparse
import os
import shutil
import subprocess
import sys
from collections import defaultdict


def load_manifest(path):
    """Return [(genome_id, genus, species, fasta_path), ...] for placed genomes."""
    rows = []
    with open(path) as fh:
        col = {n: i for i, n in enumerate(fh.readline().rstrip("\n").split("\t"))}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            genus = f[col["genus"]]
            if genus == "NA":
                continue
            rows.append((f[col["genome_id"]], genus, f[col["species"]],
                         f[col["path"]]))
    return rows


def run_skani_triangle(skani, genome_list, out_sparse, threads, min_af, screen):
    """All-vs-all ANI for one genus -> sparse edge table (skani triangle -E)."""
    cmd = [skani, "triangle", "-l", genome_list, "-o", out_sparse,
           "-t", str(threads), "-E",              # -E: sparse (edge-list) output
           "--min-af", str(min_af * 100.0),       # skani wants a percentage
           "-s", str(screen)]                     # screen: skip pairs below this ANI
    subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL)


def parse_sparse(path, path2gid):
    """Yield (gid_a, gid_b, ani, af_max) from a skani -E sparse table.

    skani -E columns: Ref_file  Query_file  ANI  Align_fraction_ref  Align_fraction_query
    (a leading '#'-comment header line, then rows). File paths are what we passed
    in, so we map them straight back to genome ids.
    """
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("Ref_file") or line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 5:
                continue
            a = path2gid.get(f[0]) or path2gid.get(os.path.abspath(f[0]))
            b = path2gid.get(f[1]) or path2gid.get(os.path.abspath(f[1]))
            if not a or not b:
                continue
            ani = float(f[2])
            af = max(float(f[3]), float(f[4])) / 100.0
            yield a, b, ani, af


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--manifest", help="manifest.tsv (or use --base)")
    ap.add_argument("--base", help="nextflow outdir; manifest read from "
                    "<base>/manifest/manifest.tsv")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--genus", action="append",
                    help="limit to this genus (repeatable); default: all genera")
    ap.add_argument("--ani", type=float, default=95.0,
                    help="ANI (%%) at/above which a cross-species pair is a merge "
                         "candidate (default 95, the usual species boundary)")
    ap.add_argument("--min_af", type=float, default=0.5,
                    help="min alignment fraction for a pair to count (default 0.5)")
    ap.add_argument("--screen", type=float, default=80.0,
                    help="skani -s screen: skip pairs below this ANI (default 80)")
    ap.add_argument("--min_genomes", type=int, default=2,
                    help="skip genera with fewer genomes than this (default 2)")
    ap.add_argument("--threads", type=int, default=os.cpu_count() or 4)
    ap.add_argument("--skani", default="skani")
    args = ap.parse_args()

    manifest = args.manifest or (args.base and
                                 os.path.join(args.base, "manifest", "manifest.tsv"))
    if not manifest or not os.path.exists(manifest):
        sys.exit("error: give --manifest or --base pointing at a run with manifest.tsv")
    if not shutil.which(args.skani):
        sys.exit(f"error: skani not found ('{args.skani}'). Install it or pass --skani.")

    rows = load_manifest(manifest)
    by_genus = defaultdict(list)
    gid_species = {}
    path2gid = {}
    for gid, genus, species, path in rows:
        by_genus[genus].append((gid, path))
        gid_species[gid] = species
        path2gid[path] = gid
        path2gid[os.path.abspath(path)] = gid

    wanted = set(args.genus) if args.genus else None
    os.makedirs(args.outdir, exist_ok=True)

    # (species_a, species_b) -> [ani, ...] and running af max, over cross-species pairs.
    pair_anis = defaultdict(list)
    pair_afmax = defaultdict(float)
    done = 0
    for genus, members in sorted(by_genus.items()):
        if wanted and genus not in wanted:
            continue
        if len(members) < args.min_genomes:
            continue
        gdir = os.path.join(args.outdir, genus.replace(" ", "_"))
        os.makedirs(gdir, exist_ok=True)
        list_path = os.path.join(gdir, "genomes.txt")
        with open(list_path, "w") as fh:
            for _gid, path in members:
                fh.write(path + "\n")
        sparse = os.path.join(gdir, "skani.sparse.tsv")
        sys.stderr.write(f"[skani] {genus}: {len(members)} genomes\n")
        run_skani_triangle(args.skani, list_path, sparse, args.threads,
                           args.min_af, args.screen)
        for a, b, ani, af in parse_sparse(sparse, path2gid):
            sa, sb = gid_species.get(a, "NA"), gid_species.get(b, "NA")
            if sa == sb or "NA" in (sa, sb):
                continue           # within-species (or unplaced) pair: not a merge
            key = tuple(sorted((sa, sb)))
            pair_anis[key].append(ani)
            if af > pair_afmax[key]:
                pair_afmax[key] = af
        done += 1

    # Cross-species merge candidates.
    out = os.path.join(args.outdir, "merge_candidates.tsv")
    with open(out, "w") as fh:
        fh.write("species_a\tspecies_b\tn_genome_pairs\tmean_ani\tmax_ani\t"
                 "max_align_frac\tcandidate\n")
        cand = 0
        for key, anis in sorted(pair_anis.items(),
                                key=lambda kv: -max(kv[1])):
            mean_ani = sum(anis) / len(anis)
            max_ani = max(anis)
            afm = pair_afmax[key]
            is_cand = max_ani >= args.ani and afm >= args.min_af
            cand += is_cand
            fh.write(f"{key[0]}\t{key[1]}\t{len(anis)}\t{mean_ani:.2f}\t"
                     f"{max_ani:.2f}\t{afm:.3f}\t{'yes' if is_cand else 'no'}\n")

    sys.stderr.write(
        f"[skani] processed {done} genera; {len(pair_anis)} cross-species pairs, "
        f"{cand} flagged >= {args.ani}%% ANI. -> {out}\n")


if __name__ == "__main__":
    main()

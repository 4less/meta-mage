#!/usr/bin/env python3
"""
Pick the relaxed-tier candidates worth guarding for the adaptive rescue.

From the relaxed scored table (markers scored at the low in-prevalence floor),
keep only the rows that are:
  * for a FLAGGED species (still below the marker threshold after mask recovery,
    and whose merge probe does not reach the threshold), AND
  * NOT already guarded in the main run -- i.e. every candidate the guard has not
    yet had a verdict on.

The second test is against the emitted table (the markers that actually reached
CROSSMAP), NOT against in_prevalence. A high-prevalence candidate is only
"already covered" if it survived the per-clade score cap AND the min_gene_len
filter; for a species with more specific candidates than max_markers_per_clade,
most never reach the guard at all. Filtering on in_prev >= min_in_start assumed
the 0.8 run tested them and silently excluded the best untested candidates,
leaving relaxation to scrape accessory genes it was always going to fail on.

These are the extra candidates fed to EMIT_REPS -> CROSSMAP -> clean guard.
"""
import argparse


def load_flagged(low_marker, merge_gain):
    """Species that are flagged (below threshold) AND whose merge does not reach
    the threshold -- the ones the adaptive relaxation should try to rescue."""
    flagged = set()
    with open(low_marker) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        ci = {n: i for i, n in enumerate(header)}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if f[ci["flagged"]] == "yes":
                flagged.add(f[ci["species"]])
    # Drop species a merge would rescue (merge is preferred; don't also relax them).
    if merge_gain:
        try:
            fh = open(merge_gain)
        except OSError:
            fh = None
        if fh:
            with fh:
                header = fh.readline().rstrip("\n").split("\t")
                ci = {n: i for i, n in enumerate(header)}
                si, ri = ci.get("species"), ci.get("reached_threshold")
                if si is not None and ri is not None:
                    for line in fh:
                        f = line.rstrip("\n").split("\t")
                        if len(f) > max(si, ri) and f[ri] == "yes":
                            flagged.discard(f[si])
    return flagged


def load_guarded(path, rank):
    """(clade, cluster) pairs that already reached the guard in the main run."""
    guarded = set()
    if not path:
        return guarded
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        ci = {n: i for i, n in enumerate(header)}
        ri, ci_c, cl = ci.get("rank"), ci.get("clade"), ci.get("cluster")
        if None in (ri, ci_c, cl):
            return guarded
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) > max(ri, ci_c, cl) and f[ri] == rank:
                guarded.add((f[ci_c], f[cl]))
    return guarded


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--relaxed", required=True, help="markers_relaxed.tsv (scored at floor)")
    ap.add_argument("--low_marker", required=True, help="low_marker_species.tsv")
    ap.add_argument("--merge_gain", help="merge_gain.tsv (or NO_FILE)")
    ap.add_argument("--emitted", required=True,
                    help="markers.emitted.tsv -- candidates the main run actually "
                         "guarded; everything else is fair game for the rescue")
    ap.add_argument("--rank", default="species")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    mg = args.merge_gain
    if mg and mg.endswith("NO_FILE"):
        mg = None
    flagged = load_flagged(args.low_marker, mg)
    guarded = load_guarded(args.emitted, args.rank)

    kept = untested_hi = 0
    with open(args.relaxed) as fh, open(args.out, "w") as out:
        header = fh.readline()
        out.write(header)
        cols = header.rstrip("\n").split("\t")
        ip_i = cols.index("in_prevalence")
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if f[0] != args.rank or f[1] not in flagged:
                continue
            if (f[1], f[2]) in guarded:
                continue  # the guard already ruled on this one
            out.write(line if line.endswith("\n") else line + "\n")
            kept += 1
            try:
                if float(f[ip_i]) >= 0.8:
                    untested_hi += 1
            except (ValueError, IndexError):
                pass
    print(f"filter_relax: {len(flagged)} flagged species, {kept} extra candidates "
          f"({untested_hi} of them core-prevalence candidates the main run never "
          f"guarded)")


if __name__ == "__main__":
    main()

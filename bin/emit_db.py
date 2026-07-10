#!/usr/bin/env python3
"""
Finalise the classifier marker DB.

  * Annotate each sequence header with the marker's clade specificity and
    within-clade diversity, so the DB is self-describing.
  * Drop exact-duplicate CDS within a marker (identical rep copies add nothing).
    Non-identical copies are kept -- they cover within-clade nucleotide diversity,
    which is exactly what the diversity columns quantify.

Final header:
    >{cluster}|{rank}|{clade}|{genome_id} in={in_prev} out={out_prev} within_id={mean_within_id}
"""
import argparse
from collections import defaultdict


def load_table(path):
    """(rank,clade,cluster) -> dict of annotation fields."""
    info = {}
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        col = {n: i for i, n in enumerate(header)}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            key = (f[col["rank"]], f[col["clade"]], f[col["cluster"]])
            info[key] = {
                "in": f[col["in_prevalence"]],
                "out": f[col["out_prevalence"]],
                "within": f[col.get("mean_within_id", -1)] if "mean_within_id" in col else "NA",
            }
    return info


def iter_fasta(path):
    header = None
    buf = []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(buf)
                header = line[1:].strip()
                buf = []
            else:
                buf.append(line.strip())
        if header is not None:
            yield header, "".join(buf)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--table", required=True)
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--out_fasta", required=True)
    ap.add_argument("--out_table", required=True)
    args = ap.parse_args()

    info = load_table(args.table)
    seen = defaultdict(set)  # (rank,clade,cluster) -> set of seqs already written
    n_written = 0

    with open(args.out_fasta, "w") as out:
        for header, seq in iter_fasta(args.fasta):
            # header: rank|clade|cluster|genome_id
            parts = header.split("|")
            rank, clade, cluster, gid = (parts + ["", "", "", ""])[:4]
            key = (rank, clade, cluster)
            if seq in seen[key]:
                continue  # exact duplicate copy
            seen[key].add(seq)
            ann = info.get(key, {"in": "NA", "out": "NA", "within": "NA"})
            out.write(
                f">{cluster}|{rank}|{clade}|{gid} "
                f"in={ann['in']} out={ann['out']} within_id={ann['within']}\n"
            )
            for i in range(0, len(seq), 80):
                out.write(seq[i:i + 80] + "\n")
            n_written += 1

    # Marker-level summary table (one row per rank/clade/cluster kept).
    with open(args.out_table, "w") as out:
        out.write("cluster\trank\tclade\tn_seqs\tin_prevalence\tout_prevalence\tmean_within_id\n")
        for (rank, clade, cluster), seqs in seen.items():
            ann = info.get((rank, clade, cluster), {"in": "NA", "out": "NA", "within": "NA"})
            out.write(f"{cluster}\t{rank}\t{clade}\t{len(seqs)}\t"
                      f"{ann['in']}\t{ann['out']}\t{ann['within']}\n")

    print(f"emit_db: wrote {n_written} sequences across {len(seen)} markers")


if __name__ == "__main__":
    main()

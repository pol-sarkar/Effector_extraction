#!/usr/bin/env python3
"""Extract the best TBLASTN-supported effector sequence for each sample."""

from __future__ import annotations
import argparse
import csv
from collections import defaultdict
from pathlib import Path

def parse_args():
    p = argparse.ArgumentParser(description="Extract best TBLASTN-supported effector sequences per sample.")
    p.add_argument("--queries", required=True, help="Candidate effector protein FASTA.")
    p.add_argument("--genomes", required=True, help="Combined sample consensus nucleotide FASTA.")
    p.add_argument("--hits", required=True, help="TBLASTN outfmt 6 TSV.")
    p.add_argument("--metadata", required=True, help="TSV with sample, isolate, group columns.")
    p.add_argument("--outdir", required=True, help="Output directory.")
    return p.parse_args()

def read_fasta(path):
    seqs, header, parts = {}, None, []
    with open(path) as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    seqs[header] = "".join(parts)
                header = line[1:].split()[0]
                parts = []
            else:
                parts.append(line)
        if header is not None:
            seqs[header] = "".join(parts)
    return seqs

def read_metadata(path):
    metadata = {}
    with open(path, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"sample", "isolate", "group"}
        if not required.issubset(reader.fieldnames or []):
            raise ValueError("Metadata must contain: sample, isolate, group")
        for row in reader:
            metadata[row["sample"]] = {"isolate": row["isolate"], "group": row["group"]}
    return metadata

def reverse_complement(seq):
    table = str.maketrans("ACGTRYMKBDHVacgtrymkbdhvNn", "TGCAYRKMVHDBtgcayrkmvhdbNn")
    return seq.translate(table)[::-1]

def identify_sample(subject_id, sample_ids):
    for sample in sorted(sample_ids, key=len, reverse=True):
        if subject_id == sample or subject_id.startswith(sample + "_"):
            return sample
    return None

def main():
    args = parse_args()
    outdir = Path(args.outdir)
    nt_dir = outdir / "by_effector_nt"
    aa_dir = outdir / "by_effector_aa"
    nt_dir.mkdir(parents=True, exist_ok=True)
    aa_dir.mkdir(parents=True, exist_ok=True)

    queries = read_fasta(args.queries)
    genomes = read_fasta(args.genomes)
    metadata = read_metadata(args.metadata)
    sample_order = list(metadata)

    columns = ["qseqid","sseqid","pident","length","mismatch","qstart","qend","sstart","send","evalue","bitscore","sseq"]
    hits = []

    with open(args.hits) as handle:
        reader = csv.reader(handle, delimiter="\t")
        for row in reader:
            if len(row) < 12:
                continue
            r = dict(zip(columns, row[:12]))
            sample = identify_sample(r["sseqid"], metadata.keys())
            if sample is None or r["qseqid"] not in queries:
                continue

            r["sample"] = sample
            r["pident"] = float(r["pident"])
            r["length"] = int(r["length"])
            r["qstart"] = int(r["qstart"])
            r["qend"] = int(r["qend"])
            r["sstart"] = int(r["sstart"])
            r["send"] = int(r["send"])
            r["evalue"] = float(r["evalue"])
            r["bitscore"] = float(r["bitscore"])
            r["query_length_aa"] = len(queries[r["qseqid"]])
            r["qcov_pct"] = (abs(r["qend"] - r["qstart"]) + 1) / r["query_length_aa"] * 100.0
            hits.append(r)

    best = {}
    for r in hits:
        key = (r["qseqid"], r["sample"])
        rank = (r["bitscore"], r["qcov_pct"], r["pident"], r["length"])
        if key not in best or rank > best[key][0]:
            best[key] = (rank, r)

    records = []
    for (qid, sample), (_, r) in best.items():
        subject = r["sseqid"]
        if subject not in genomes:
            print(f"WARNING: subject not found: {subject}")
            continue

        start, end = min(r["sstart"], r["send"]), max(r["sstart"], r["send"])
        nt = genomes[subject][start-1:end]
        strand = "+"
        if r["sstart"] > r["send"]:
            nt = reverse_complement(nt)
            strand = "-"

        aa = r["sseq"].replace("-", "")
        records.append({
            "effector": qid,
            "sample": sample,
            "isolate": metadata[sample]["isolate"],
            "group": metadata[sample]["group"],
            "subject": subject,
            "pident": r["pident"],
            "qcov_pct": r["qcov_pct"],
            "alignment_length_aa": r["length"],
            "query_length_aa": r["query_length_aa"],
            "evalue": r["evalue"],
            "bitscore": r["bitscore"],
            "start": start,
            "end": end,
            "strand": strand,
            "nt_sequence": nt,
            "aa_hit_sequence": aa
        })

    sample_index = {s:i for i,s in enumerate(sample_order)}
    records.sort(key=lambda x: (x["effector"], sample_index.get(x["sample"], 10**9)))

    summary_fields = [
        "effector","sample","isolate","group","subject","pident","qcov_pct",
        "alignment_length_aa","query_length_aa","evalue","bitscore","start","end",
        "strand","nt_length","aa_hit_length"
    ]
    summary_path = outdir / "effector_sequences_summary.tsv"
    with open(summary_path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=summary_fields, delimiter="\t")
        writer.writeheader()
        for r in records:
            writer.writerow({
                "effector": r["effector"],
                "sample": r["sample"],
                "isolate": r["isolate"],
                "group": r["group"],
                "subject": r["subject"],
                "pident": f'{r["pident"]:.3f}',
                "qcov_pct": f'{r["qcov_pct"]:.2f}',
                "alignment_length_aa": r["alignment_length_aa"],
                "query_length_aa": r["query_length_aa"],
                "evalue": r["evalue"],
                "bitscore": r["bitscore"],
                "start": r["start"],
                "end": r["end"],
                "strand": r["strand"],
                "nt_length": len(r["nt_sequence"]),
                "aa_hit_length": len(r["aa_hit_sequence"])
            })

    with open(outdir / "all_effectors_nt.fasta", "w") as ntout, open(outdir / "all_effectors_aa_hits.fasta", "w") as aaout:
        for r in records:
            header = (
                f'{r["effector"]}|{r["sample"]}|{r["isolate"]}|{r["group"]}|'
                f'pident={r["pident"]:.2f}|qcov={r["qcov_pct"]:.1f}|'
                f'{r["subject"]}:{r["start"]}-{r["end"]}({r["strand"]})'
            )
            ntout.write(f">{header}\n{r['nt_sequence']}\n")
            aaout.write(f">{header}\n{r['aa_hit_sequence']}\n")

    grouped = defaultdict(list)
    for r in records:
        grouped[r["effector"]].append(r)

    for effector, subset in grouped.items():
        with open(nt_dir / f"{effector}_samples_nt.fasta", "w") as ntout, open(aa_dir / f"{effector}_samples_aa.fasta", "w") as aaout:
            aaout.write(f">REFERENCE|{effector}\n{queries[effector]}\n")
            for r in subset:
                header = f'{r["sample"]}|{r["isolate"]}|{r["group"]}|pident={r["pident"]:.2f}|qcov={r["qcov_pct"]:.1f}'
                ntout.write(f">{header}\n{r['nt_sequence']}\n")
                aaout.write(f">{header}\n{r['aa_hit_sequence']}\n")

    print(f"Selected effector x sample records: {len(records)}")
    print(f"Summary: {summary_path}")

if __name__ == "__main__":
    main()

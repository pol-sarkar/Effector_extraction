# CLso candidate effector recovery by TBLASTN

This repository documents a workflow for identifying and extracting candidate effector-associated sequences from sample-specific *Candidatus Liberibacter solanacearum* (CLso) consensus genomes.
Raw Illumina sequencing data used in this study are available through NCBI under BioProject PRJNA1399511.
## Overview

Candidate effector amino-acid sequences were searched against sample-specific consensus nucleotide sequences from 10 psyllid samples representing three geographic regions (Northwestern, Western, and Southwestern) using TBLASTN.

All computational analyses were performed on a Linux-based high-performance computing cluster (HPCC) using command-line tools. NCBI BLAST+ was used to construct the nucleotide BLAST database and perform TBLASTN searches, while custom Python scripts were used to process BLAST results, calculate query coverage, select the best-supported effector hit for each sample, and extract the corresponding nucleotide and translated amino-acid sequences.

```text
The workflow was:

Candidate effector proteins
        |
        v
Sample-specific CLso consensus genomes
        |
        v
Build local nucleotide BLAST database
        |
        v
TBLASTN protein-to-nucleotide search
        |
        v
Best hit for each effector × sample
        |
        +--> amino-acid percent identity (pident)
        +--> query coverage (qcov)
        +--> genomic coordinates and strand
        +--> nucleotide sequence extraction
        +--> translated amino-acid hit sequence
        |
        v
Effector-specific FASTA files
        |
        v
Multiple-sequence alignment and comparative analysis

TBLASTN was used because the query sequences are proteins and the target database contains nucleotide sequences. TBLASTN translates the nucleotide database in all six reading frames during the search.

## Repository structure

```text
clso-effector-tblastn/
├── README.md
├── extract_effector_sequences.py
├── sample_metadata.example.tsv
└── .gitignore
```

MEGA — downstream visualization and multiple-sequence alignment of selected effector sequences

## Requirements

- Python 3.8+
- NCBI BLAST+ (`makeblastdb`, `tblastn`)

The Python extraction script uses only the standard library.

## Inputs

### Candidate effector proteins


```text
effectors_query.fasta
```

### Sample-specific consensus genomes

Combine the consensus nucleotide sequences into one FASTA:

```text
all_consensus_genomes.fna
```

Subject identifiers should begin with the sample identifier in the metadata table, for example:

```text
>A_S1_NC_014774.1
```

### Sample metadata

File: sample_metadata.tsv

```text
File_name 	Sample_Name	Clso_Haplotype	Psyllid_Haplotype
A_S1	CSU9	CLsoG	Northwestern
B_S2	CSU14	CLsoG	Northwestern
G_S7	EOR13a	CLsoSumb2	Northwestern
H_S8	EOR18	CLsoSumb2	Northwestern
L_S12	DFS3	CLsoA	Western
M_S13	DFS10	CLsoA	Western
O_S15	JTL7	Southwestern	Southwestern
P_S16	JTL8	Southwestern	Southwestern
R_S18	ELL2	CLsoF	Northwestern
S_S19	ELL4	CLsoF	Northwestern

```

An example is provided as `sample_metadata.example.tsv`.

## 1. Build the nucleotide BLAST database

```bash
makeblastdb \
  -in all_consensus_genomes.fna \
  -dbtype nucl \
  -out all_consensus_db
```

## 2. Run TBLASTN

```bash
tblastn \
  -query effectors_query.fasta \
  -db all_consensus_db \
  -evalue 1e-5 \
  -outfmt "6 qseqid sseqid pident length mismatch qstart qend sstart send evalue bitscore sseq" \
  -out effector_tblastn_hits.tsv
```

The `pident` value is amino-acid identity because TBLASTN compares the protein query with translated nucleotide sequence.

## 3. Extract sample-specific sequences

```bash
python3 extract_effector_sequences.py \
  --queries effectors_query.fasta \
  --genomes all_consensus_genomes.fna \
  --hits effector_tblastn_hits.tsv \
  --metadata sample_metadata.tsv \
  --outdir results
```

For each effector × sample combination, the script retains the primary TBLASTN HSP ranked by bit score, then query coverage, amino-acid identity, and alignment length.

The corresponding nucleotide interval is extracted from the sample consensus genome. Reverse-strand hits are reverse complemented.

## Outputs

- `effector_sequences_summary.tsv`: BLAST metrics and genomic coordinates for each retained hit.
- `all_effectors_nt.fasta`: combined extracted nucleotide sequences.
- `all_effectors_aa_hits.fasta`: translated TBLASTN hit sequences.
- `by_effector_nt/`: one nucleotide FASTA per effector.
- `by_effector_aa/`: one amino-acid FASTA per effector, including the reference query.

The per-effector FASTA files can be imported individually into MEGA, MAFFT, MUSCLE, or other alignment software.

## Interpretation

Percent identity and query coverage should be interpreted together.

A hit with:

```text
pident = 100%, qcov = 100%
```

supports a full or near-full highly similar homolog.

A hit with:

```text
pident = 95%, qcov = 35%
```

represents only a partial match and should not be described as a complete effector solely from percent identity.

The amino-acid sequences written by this workflow are the translated TBLASTN HSPs. Partial hits may not represent complete predicted proteins. Full-length coding sequences should be confirmed independently before functional interpretation.

## Example 10-sample dataset

| Sample | Isolate | Group |
|---|---|---|
| A_S1 | CSU9 | CLsoG |
| B_S2 | CSU14 | CLsoG |
| G_S7 | EOR13a | CLsoSumb2 |
| H_S8 | EOR18 | CLsoSumb2 |
| L_S12 | DFS3 | CLsoA |
| M_S13 | DFS10 | CLsoA |
| O_S15 | JTL7 | Southwestern |
| P_S16 | JTL8 | Southwestern |
| R_S18 | ELL2 | CLsoF |
| S_S19 | ELL4 | CLsoF |

## Manuscript-ready methods summary

Candidate CLso effector protein sequences were searched against sample-specific CLso consensus nucleotide sequences using TBLASTN. Consensus sequences were formatted as a local nucleotide BLAST database using `makeblastdb`, and searches were performed with an E-value threshold of 1 × 10^-5. For each effector-sample combination, the highest-scoring homologous region was retained, and amino-acid percent identity and query coverage were recorded. The corresponding nucleotide sequence was extracted from the sample consensus genome using TBLASTN subject coordinates, with reverse complementation for negative-strand hits. Extracted sequences were organized by candidate effector for downstream multiple-sequence alignment and comparative analysis.

## Scope

This repository documents the targeted effector-recovery step. It does not recreate upstream trimming, host filtering, reference mapping, consensus generation, or assembly.

## Citation

Add the manuscript citation and/or DOI here when available.

## License

Before public release, add the license that best fits your project (for example MIT, BSD-3-Clause, or GPL-3.0).

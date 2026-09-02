Workflow
Candidate effector proteins
        ↓
TBLASTN against CLso-ZC1
(NC_014774.1)
        ↓
Identify approximate genomic locus
        ↓
Original trimmed Illumina reads
        ↓
Map each sample to the complete ZC1 genome
        ↓
Calculate nucleotide-by-nucleotide coverage
        ↓
Call sample-specific variants
        ↓
Mask insufficiently supported positions as N
        ↓
Extract read-supported effector nucleotide sequence
        ↓
Translate sufficiently complete sequences
        ↓
Protein alignment / sequence identity / phylogeny


# CLso candidate effector recovery by TBLASTN

This repository documents a workflow for identifying and extracting candidate effector-associated sequences from sample-specific *Candidatus Liberibacter solanacearum* (CLso) consensus genomes.
Raw sequencing data used in this study are available through NCBI under BioProject PRJNA1399511.
## Overview

Candidate effector amino-acid sequences were searched against sample-specific consensus nucleotide sequences from 10 psyllid samples representing three geographic regions (Northwestern, Western, and Southwestern) using TBLASTN. TBLASTN was used because the effector queries are amino acid sequences whereas the sample consensus genomes are nucleotide sequences. TBLASTN translates the nucleotide database in all six reading frames during the search and compares the resulting translated sequences with each protein query.
Following the initial TBLASTN search, downstream sequence extraction and comparative analyses were restricted to 10 selected psyllid samples representing three psyllid haplotypes/geographic groups: Northwestern, Western, and Southwestern. These samples also represented multiple CLso groups, including CLsoG, CLsoSumb2, CLsoA, Southwestern, and CLsoF.

All computational analyses were performed on a Linux-based high-performance computing cluster (HPCC) using command-line tools. NCBI BLAST+ was used to construct the nucleotide BLAST database and perform TBLASTN searches, while custom Python scripts were used to process BLAST results, calculate query coverage, select the best-supported effector hit for each sample, and extract the corresponding nucleotide and translated amino-acid sequences.

## Inputs
1. Candidate effector proteins: effectors_query.fasta (This FASTA file was used as the protein query input for TBLASTN.)

2. Sample-specific CLso consensus genomes: all_consensus_genomes.fna (These nucleotide sequences were used to construct the local BLAST database.)

The consensus sequences were not translated before the analysis. Translation of the nucleotide database was performed internally by TBLASTN in all six reading frames.

3. Sample metadata: sample_metadata.tsv

File_name   Sample_Name   CLso_Haplotype   Psyllid_Haplotype
A_S1        CSU9          CLsoG            Northwestern
B_S2        CSU14         CLsoG            Northwestern
G_S7        EOR13a        CLsoSumb2        Northwestern
H_S8        EOR18         CLsoSumb2        Northwestern
L_S12       DFS3          CLsoA            Western
M_S13       DFS10         CLsoA            Western
O_S15       JTL7          Southwestern     Southwestern
P_S16       JTL8          Southwestern     Southwestern
R_S18       ELL2          CLsoF            Northwestern
S_S19       ELL4          CLsoF            Northwestern

The metadata distinguish the CLso group detected in each sample from the psyllid haplotype/geographic group.

## Workflow
1. Build the nucleotide BLAST database

A local nucleotide database was generated from the combined sample consensus sequences using NCBI BLAST+:

makeblastdb \
  -in all_consensus_genomes.fna \
  -dbtype nucl \
  -out all_consensus_db

This generated the BLAST database files beginning with:
all_consensus_db.*
2. Run TBLASTN
Candidate effector proteins were searched against the nucleotide database using:

tblastn \
  -query effectors_query.fasta \
  -db all_consensus_db \
  -evalue 1e-5 \
  -outfmt "6 qseqid sseqid pident length mismatch qstart qend sstart send evalue bitscore sseq" \
  -out effector_tblastn_hits.tsv

The TBLASTN output contained:

Field	Description
qseqid	Candidate effector protein identifier
sseqid	Sample consensus sequence identifier
pident	Amino-acid percent identity
length	Length of the protein-level alignment
mismatch	Number of amino-acid mismatches
qstart, qend	Alignment coordinates in the effector query
sstart, send	Coordinates of the corresponding region in the nucleotide consensus sequence
evalue	BLAST E-value
bitscore	BLAST bit score
sseq	Translated subject sequence corresponding to the TBLASTN alignment

Because TBLASTN compares a protein query against translated nucleotide sequence, pident represents amino-acid identity, not nucleotide identity.

3. Select focal samples and extract effector sequences

The TBLASTN results were subsequently filtered to retain the 10 focal psyllid samples listed above.

For each candidate effector × sample combination, the best-supported TBLASTN hit was selected. Hits were ranked primarily by bit score, with query coverage, amino-acid identity, and alignment length used as additional criteria.

Query coverage was calculated as:

query coverage (%) = aligned query length / total query protein length × 100

The genomic coordinates reported by TBLASTN were then used to extract the corresponding nucleotide sequence from the appropriate sample consensus genome.

For negative-strand hits, the extracted nucleotide sequence was reverse complemented so that all sequences were reported in the coding orientation.

Python scripts were used to automate hit selection, calculation of query coverage, sequence extraction, strand correction, and generation of summary tables and FASTA files.

## Outputs:

The final 10-sample analysis is stored under:

Output/effectors/

Major outputs include:

all_effectors_10samples_aa_hits.fasta
all_effectors_10samples_nt.fasta
effector_best_hits_10samples.tsv
effector_hits_10samples_aa.fasta
effector_hits_10samples_nt.fasta
effector_sequences_10samples_summary.tsv
by_effector_aa/
by_effector_nt/
effector_sequences_10samples_summary.tsv

Summary information for the selected effector–sample hits, including BLAST similarity statistics, query coverage, genomic coordinates, and sequence information.

all_effectors_10samples_nt.fasta

Combined nucleotide sequences corresponding to the selected candidate effector regions across the 10 focal samples.

all_effectors_10samples_aa_hits.fasta

Combined translated amino-acid sequences corresponding to the selected TBLASTN hits.

by_effector_nt/

Contains a separate nucleotide FASTA file for each candidate effector across the focal samples.

by_effector_aa/

Contains a separate amino-acid FASTA file for each candidate effector. These files can be used individually for multiple-sequence alignment and comparative protein analysis.

## Interpretation

TBLASTN amino-acid identity (pident) and query coverage (qcov) was evaluated together. ( i.e. if pident = 100%, qcov   = 100%, the translated genomic region is identical to the reference effector protein across the complete query sequence. In contrast: pident = 95%, qcov   = 35% indicates high similarity across only a portion of the reference protein and should therefore be interpreted as a partial homologous match, rather than evidence for a complete effector protein.)

The amino-acid sequences recovered directly from TBLASTN represent the translated regions covered by individual BLAST alignments. Consequently, partial TBLASTN hits may not represent complete predicted proteins, and candidate genes of particular biological interest has to be further evaluated for complete coding sequence and open-reading-frame integrity.

## Output Data availability
https://drive.google.com/drive/folders/1FImwGDpTj1xvamab_ksXMS5gbJk7cWQh


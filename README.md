# clonal_scores
Perl code to find clonal sequences in B cell receptor sequencing data, defined as having the same IGHV usage, the same IGHJ usage, the same HCDR3 length, and at least X% similarity in HCDR3 sequence, where X is user defined.

## Usage
perl clonal_scores.pl sequence_data.csv X

where X is the minimum percent similarity in HCDR3 sequence for two sequences to be considered clonally related.

## Input file

The csv file should be formatted as follows:

Cell  Isotype Timepoint IGHV  IGHJ  HCDR3

atMBC IgG 0 IGHV1-58  IGHJ6 CAAGKYYYDSSGYSPYYYYYGMDVW

atMBC IgG 0 IGHV5-51  IGHJ6 CAGEQWKLRGVNYYFYGMDVW

## Output file

Running the script will generate an output file with extension "Results.txt" that lists the total cell counts in each pair of B cell subsets, the number of clonal connections, and the clonal score:

TP1	Ig1	Cell1	Connections1	Total count cells1	TP2	Ig2	Cell2	Connections2	Total count cells2	Clonal score

0	IgM	Naive_B_cell	169	5785	0	IgM	Naive_B_cell	169	5785	2.92

0	IgM	Naive_B_cell	3	5785	0	IgM	FcRL5-_classical_MBC	3	342	0.88



#!/usr/bin/env bash

# Modules to load
module load blast/2.13.0

# Start
echo date

# Working directory
cd /project/berglandlab/connor

# Blast sequence
blastn -query tcl2_pulex.fa \
-subject totalHiCwithallbestgapclosed.fa > \
tcl2_pulex.blast.out

# TCL2 gene
"Scaffold_9200_HRSCAF_10757:6676234-6686097"

# Extract translated fasta from genes
~/seqtk/seqtk subseq \
"Daphnia.proteins.aed.0.6.fasta" \
"toll/tcl2" > \
"toll/tcl2.aa.fa"

# Finish
echo "Finish"
echo date

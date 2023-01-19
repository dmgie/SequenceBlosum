#!/usr/bin/env bash

# Step 1. Obtain sequences via entrez script for protein, get fasta file for each temperature group
# ./entrez_get_fasta.sh


# Step 2. MSA using MEGA 
# Done in the GUI, saved as TempGroup_aligned.fas

# Step 3. Combine all fasta files into one
cat *.fasta > all.fasta

# Step 4. Get matrices from each MSA
for i in *.fas; do python3 script.py $i > matrices/$i.mat; done

# Step 5. Run clustalw2 on the combined fasta file using the each of the matrices 
# for i in matrices/*.mat; do clustalw2 -INFILE=all.fasta -MATRIX=$i -OUTFILE=$i.aln; done
# Or done using R script, with the "msa" library, and the command :
# msaClustalW("combined.fasta", type = "protein", substitutionMatrix = "matrices/<TempGroup>_Aligned.fas.matrix", verbose = TRUE)
# replacing <TempGroup> with the name of the temperature group

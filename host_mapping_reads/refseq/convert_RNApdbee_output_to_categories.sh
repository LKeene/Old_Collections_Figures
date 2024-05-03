#!/bin/bash

# this script converts the exported dot-bracket notation file 
# from RNApdbee to a file listing the category of each base in 
# the D. melanogaster rRNA sequence
#
# MDS 6/13/2023
#
#  Convert Structure into base-pairing / secondary structure info
#
# Structure: 
# PDB 4V6W
# Anger et al (2013) Nature
# PMID: 23636399
# high-resolution cryo-electron-microscopy structure

# main rRNA chains:
# B2: 18S SSU
# A5: 28S LSU
# A8: 5.8S 

# feed this structure into RNApdbee web interface:

# http://rnapdbee.cs.put.poznan.pl/
#
# On this website:
# 
# 1) Upload RNA 3D structure:
#    - from Protein Data Bank: 4V6W --> Get
#
# 2) Set parameters
#    - Include non-canonical ones [base pairs]: in text and visualization
# 
#  leave all other settings as defaults
#
#
# After run complete:
# - select all results
# - download selected results
#
# RNApdbee 2.0 
# T. Zok, M. Antczak, M. Zurkowski, M. Popenda, J. Blazewicz, R.W. Adamiak, M. Szachniuk. RNApdbee 2.0: multifunctional tool for RNA structure annotation. Nucleic Acids Research 46(W1), 2018, W30-W35, (doi:10.1093/nar/gky314)
# 
# 
# 

# separate out the output into b2.dot.txt and a5.dot.txt

# convert dot-bracket notation to interaction categories file using dot_bracket_to_bp script in this directory

./dot_bracket_to_bp -n B2  b2.dot.txt  > B2.category.txt
./dot_bracket_to_bp -n A5  a5.dot.txt  > A5.category.txt
cat A5.category.txt B2.category.txt  > interaction_categories.txt

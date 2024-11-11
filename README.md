# CandidateGeneLDWindow

## Clone this repository

git clone git@github.com:NikeeShrestha/CandidateGeneLDWindow.git

## Run this code by providing files (see example in the Data folder) to get candidate genes between the window provided for each SNP (see file: Data/LDwindowSignals.xlsx).

## Command for help

python CandidateGene.py or python CandidateGene.py -h

## Command example using files in the Data folder

python CandidateGene.py --gff3 Data/Sbicolor_730_v5.1.gene_exons.gff3 --LDWindow Data/LDwindowSignals.xlsx --annotationfile Data/Sbicolor_454_v3.1.1.annotation_info.txt --abundance Data/abundance.tsv

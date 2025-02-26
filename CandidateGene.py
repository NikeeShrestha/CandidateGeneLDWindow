from Parsefunctions import *

annotationfile=annotationparse('Data/Sbicolor_454_v3.1.1.annotation_info.txt')
##put k=0,don't change this value
LDfile('Data/LDwindowSignals.xlsx', chromcolumn=1, windowstartcolumn=2, windowendcolumn=3,
       snpcolumn=0, k=0, species='sorghum', pathtogff3='Data/Sbicolor_730_v5.1.gene_exons.gff3',annotationfile=annotationfile)
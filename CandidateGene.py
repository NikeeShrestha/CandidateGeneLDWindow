from Parsefunctions import *

annotationfile=annotationparse('/Users/nikeeshrestha/Desktop/CandidateGeneLDWindow/Data/Zmays_833_Zm-B73-REFERENCE-NAM-5.0.55.annotation_info.txt')
##put k=0,don't change this value
LDfile('/Users/nikeeshrestha/Desktop/candidategenes.xlsx', chromcolumn=3, windowstartcolumn=4, windowendcolumn=5,
       snpcolumn=4, k=0, species='maize', pathtogff3='/Users/nikeeshrestha/Desktop/CandidateGeneLDWindow/Data/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3',annotation=annotationfile)
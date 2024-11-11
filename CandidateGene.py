import importlib
import subprocess
import sys
import argparse

def import_or_install(package):
    try:
        importlib.import_module(package)
    except ImportError:
        print(f"{package} not found. Installing...")
        subprocess.check_call([sys.executable, "-m", "pip", "install", package])
    finally:
        globals()[package] = importlib.import_module(package)

# Example usage
# import_or_install('numpy') 
import_or_install('gffpandas')

import gffpandas.gffpandas as gffpd
import pandas as pd
import numpy as np

parser = argparse.ArgumentParser(description ='Candidate Gene Within Window')

parser.add_argument('--gff3', type = str, help ='full path to gff3 file')
parser.add_argument('--LDWindow', type = str, help='Excel file for SNP with LD window')
parser.add_argument("--annotationfile",type=str, default=None, help="Optional: Annotation file based on rice and arabidopsis; e.g. Sbicolor_454_v3.1.1.annotation_info.txt")
parser.add_argument('--abundance', type = str, default=None, help = 'Optional: TPM expression Data if available')

args = parser.parse_args()
if len(sys.argv) == 1:
    parser.print_help()
    sys.exit(1)
##GFF3 parse
GFF3path=args.gff3
annotation = gffpd.read_gff3(str(GFF3path))
annotation.df=annotation.df[annotation.df['type']=="gene"]
annotation.df["seq_id"]=annotation.df.seq_id.str.extract('(\d+)') ##change chr01 to 01 for comparison later
annotation.df.seq_id=pd.to_numeric(annotation.df.seq_id)
annotation.df.start=pd.to_numeric(annotation.df.start)
annotation.df.end=pd.to_numeric(annotation.df.end)
annotation.df['attributes'] = annotation.df['attributes'].str.extract(r'(?:ID=[^;]*;)?\s*Name=(.*?)(?:;|$)')
# annotation.df['attributes'] = annotation.df['attributes'].str.extract(r'(?:ID=[^;]*;)?\s*Name=(.*?)(?:;|$)')
# annotation.df

print(f'Loaded GFF3 file, Chormosomes (along with pseudo chromosomes) and Count: {annotation.df["seq_id"].value_counts().head(10)}')

##LD window parse
LDWindow=args.LDWindow
GWASLD=pd.read_excel(str(LDWindow))


# if GWASLD["Chr"].astype(str):
GWASLD["Chr"]=GWASLD.Chr.astype(str).str.extract('(\d+)')
GWASLD["Chr"]=pd.to_numeric(GWASLD["Chr"])
GWASLD.LDstart=pd.to_numeric(GWASLD.LDstart)
GWASLD.LDend=pd.to_numeric(GWASLD.LDend)
GWASLD.signal=pd.to_numeric(GWASLD.signal)
print(f'Loaded LD Window file, SNPs: {GWASLD["Signal"].value_counts()}')

###
if args.annotationfile:
   Info=pd.read_csv(str(args.annotationfile), sep="\t")
else:
   Info=None

if args.abundance:
   abundance=pd.read_csv(str(args.abundance), sep="\t")
   abundance['target_id'] = abundance['target_id'].str.replace(r'\.\d+$', '', regex=True)
else:
   abundance=None
   
###
k=0
for gwas, row1 in GWASLD.iterrows():
  chromosome1=row1["Chr"]
  LDstart_gwas=row1["LDstart"]
  LDend_gwas=row1["LDend"]
  signal_gwas=row1["signal"]
  print(chromosome1, signal_gwas,LDstart_gwas,LDend_gwas)

  Candidategene=[]

  for j, row2 in annotation.df.iterrows():
    chromosome2=row2["seq_id"]
    start=row2["start"]
    end=row2["end"]
    # if chromosome2!=chromosome1:
    #    continue
    # print(chromosome2, start,end)
    if (chromosome2==chromosome1 and ((LDstart_gwas) <= start <= (LDend_gwas) or (LDstart_gwas) <= end <= (LDend_gwas))):
      gene=row2["attributes"]
      startpositon=start
      endpsotion=end
      distancefromstart=start-signal_gwas
      distancefromend=end-signal_gwas
      totaldistance=(distancefromstart+distancefromend)/2
    #   print(gene, totaldistance)


      Candidategene.append(
            {'Candidate' : gene,
             'distance from signal': totaldistance
          }
      )
    Candidategenes=pd.DataFrame(Candidategene)
    # break
  k+=1

  print(k)

  if not Info.empty:

    Candidategenes = Candidategenes.merge(Info, how='left', left_on='Candidate', right_on='locusName')

# Dropping the extra 'locusName' column that was created during the merge
    Candidategenes = Candidategenes.drop(columns=['locusName'])
    # for i,row in Candidategenes.iterrows():
    #     gene=row["Candidate"]
    #     for j, row2 in Info.iterrows():
    #         if gene==row2["locusName"]:
    #             Candidategenes.at[i,"Best-hit-arabi-name"]=row2["Best-hit-arabi-name"]
    #             Candidategenes.at[i,"Arabdefine"]=row2["arabi-defline"]
    #             Candidategenes.at[i,"Best-hit-rice-name"]=row2["Best-hit-rice-name"]
    #             Candidategenes.at[i,"Ricedefine"]=row2["rice-defline"]
    #             break
  if not abundance.empty:
    Candidategenes = Candidategenes.merge(abundance, how='left', left_on='Candidate', right_on='target_id')
    Candidategenes = Candidategenes.drop(columns=['target_id'])
    # for i,row in Candidategenes.iterrows():
    #     gene=row["Candidate"]
    #     for j, row4 in abundance.iterrows():
    #         if gene==row4["target_id"]:
    #             Candidategenes.at[i,"tpm"]=row4["tpm"]
    #             break

  Candidategenes.to_excel("Chr"+str(chromosome1)+"_"+str(signal_gwas)+".xlsx", index=False)

  print(f'saved output for Chr{chromosome1}_{signal_gwas}')

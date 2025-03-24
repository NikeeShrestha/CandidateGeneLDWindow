import pandas as pd

def defparse(astr):
    dd = {}
    y = astr.strip().split(';')
    for a in y:
        b = a.split('=')
        dd[b[0]] = b[1]
    return dd

def annotationparse(afile):
    annotationfile={}
    with open(str(afile), 'r') as myfile:
        for line in myfile:
            myline=line.strip()
            if myline.startswith('#'): continue
            
            y=myline.split('\t')
            if len(y)<9: continue
            
            gene=y[1]
            
            # print(y)
            if len(y)>10 and y[10]:
                arabhit=y[10]
            else: arabhit=''
                
            if len(y)>12 and y[12]:
                arabdefine=y[12]
            else: arabdefine=''
            # if len(y)<2: continue
            
            if len(y)>13 and y[13]:
                ricehit=y[13]
            else: ricehit=''
            if  len(y)>15 and y[15]:
                ricedefine=y[15]
            else: ricedefine=''
                
            
            annotationfile[gene]={'Tairhit':arabhit, 'Tairdefine':arabdefine,
                                 'Ricehit':ricehit, 'Ricedefine':ricedefine}
    return annotationfile

def gff3parse(afile, species, intervalstart, intervalend, chrom, snpposition, annotation):
    
    chrom=chrom
    intervalstart=int(intervalstart)
    intervalend=int(intervalend)
    snpposiiton=int(snpposition)
    genes={}
    # count=0
    with open(str(afile), 'r') as mygff3:
        for line in mygff3:
            myline=line.strip()
            if myline.startswith('#'): continue
            
            y=myline.split('\t')
            
#             print(y)
#             count+=1
            
#             if count==50: break
            # print(y[2])
            if 'scaf' in y[0]: continue
            if y[2]!='gene': continue
            
            
            chromosome=int(y[0].lower().replace('chr',''))
            # print(chromosome)
            
            if chrom!=chromosome : continue
            
            defdict=defparse(y[-1])
            if species.lower()=='maize':
                gene=defdict['ID']
            elif species.lower()=='sorghum':
                gene=defdict['Name']
            
            genestart=int(y[3])
            geneend=int(y[4])
            midpoint=(genestart+geneend)/2
            
            
            if (intervalstart <= geneend and genestart <= intervalend):
                # print(gene)
                distancemin=int(abs(midpoint-snpposiiton))
                genes[gene]={'chrom': chromosome, 'start':genestart, 'end':geneend, 'distance_from_snp': distancemin}
                
                if annotation and gene in annotation:
                    # print('yes')
                    function=annotation[gene]
                    genes[gene].update(function)
            else:
                # print('no')
                continue
                
    return genes

def LDfile(afile, chromcolumn, windowstartcolumn, windowendcolumn, snpcolumn, k, species, pathtogff3,annotation):
    try:
        df=pd.read_excel(str(afile))
    except Exception:
        df=pd.read_csv(str(afile))

    # df=pd.read_excel(str(afile))
    # print(df.shape)
    
    if k < 0 or k > df.shape[0]-1: return
    
    row=df.iloc[k]
    # print(row[int(chromcolumn)])
    chrom=int(str(row[int(chromcolumn)]).lower().replace('chr',''))
    
    intervalstart=int(row[int(windowstartcolumn)])
    
    intervalend=int(row[int(windowendcolumn)])
    snp=row[int(snpcolumn)]
    
    # snpposition=int(str(row[int(snpcolumn)]).split('_')[1])
    snpposition=int(str(row[int(snpcolumn)]))
    
    genes=gff3parse(afile=pathtogff3, species=species, intervalstart=intervalstart,
                    intervalend=intervalend, chrom=chrom, snpposition=snpposition,
                    annotation=annotation)
    
    genesdf= pd.DataFrame.from_dict(genes, orient='index')
    # print(genesdf.head())

    genesdf.to_csv(f'Chrom_{chrom}_{snp}_candidategenes.csv', index_label='genes')
    
    # return genes,genesdf
    
    # print(k+1)
    LDfile(afile,chromcolumn, windowstartcolumn, windowendcolumn, snpcolumn, k+1, species, pathtogff3, annotation=annotation)
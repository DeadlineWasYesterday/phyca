import pandas as pd
import numpy as np
import os,sys
from Bio import SeqIO

file = sys.argv[1]
gf = 0.8
tf = 0.8

def afa_df(filepath):
    t = []
    for r in SeqIO.parse(filepath, 'fasta'):
        t.append([r.id]+list(r.seq))
    return pd.DataFrame(t)

t1 = afa_df(file)

#site gap filter
t2 = t1.replace('-',np.nan).dropna(thresh = gf*t1.shape[0], axis = 1)
#remove invariant sites
t2 = t2.loc[:,[0] + [i for i in t2.columns[1:] if t2[i].nunique() > 1]]
#taxon gap filter
t2 = t2.dropna(thresh = (tf*(t2.shape[1]-1)), axis = 0)

nt = t1.shape[0] #number of taxa
nf = t2.shape[1]-1 #no filter

t = [file, nt,nf]
#UB from 2 to 12
for x in range(2,15):
    t3 = t2.loc[:,[0] + [i2 for i2 in t2.columns[1:] if t2[i2].nunique() > x]] #UB 
    t.append(t3.shape[1]-1)
    t4 = t3.loc[:, [0] + [i2 for i2 in t3.columns[1:] if min(t3[i2].value_counts()) > 1]] #UB + NVC>1
    t.append(t4.shape[1]-1)
    t4 = t3.loc[:, [0] + [i2 for i2 in t3.columns[1:] if min(t3[i2].value_counts()) > 2]] #UB + NVC>2
    t.append(t4.shape[1]-1)
    t4 = t3.loc[:, [0] + [i2 for i2 in t3.columns[1:] if min(t3[i2].value_counts()) > 3]] #UB + NVC>3
    t.append(t4.shape[1]-1)
    
    t4 = t3.loc[:, [0] + [i2 for i2 in t3.columns[1:] if max(t3[i2].value_counts()) < t2.shape[0]*0.5]] #UB + XVC<50%
    t.append(t4.shape[1]-1)
    
    t4 = t3.loc[:, [0] + [i2 for i2 in t3.columns[1:] if max(t3[i2].value_counts()) < t2.shape[0]*0.2]] #UB + XVC<20%
    t.append(t4.shape[1]-1)
    
    t4 = t3.loc[:, [0] + [i2 for i2 in t3.columns[1:] if min(t3[i2].value_counts()) > 1]] #UB + NVC>1
    t5 = t4.loc[:, [0] + [i2 for i2 in t4.columns[1:] if max(t4[i2].value_counts()) < t2.shape[0]*0.5]] #UB + NVC>1 + XVC<50%
    t.append(t5.shape[1]-1)
    
    t4 = t3.loc[:, [0] + [i2 for i2 in t3.columns[1:] if min(t3[i2].value_counts()) > 1]] #UB + NVC>1
    t5 = t4.loc[:, [0] + [i2 for i2 in t4.columns[1:] if max(t4[i2].value_counts()) < t2.shape[0]*0.2]] #UB + NVC>1 + XVC<20%
    t.append(t5.shape[1]-1)
    
    t4 = t3.loc[:, [0] + [i2 for i2 in t3.columns[1:] if min(t3[i2].value_counts()) > 2]] #UB + NVC>2
    t5 = t4.loc[:, [0] + [i2 for i2 in t4.columns[1:] if max(t4[i2].value_counts()) < t2.shape[0]*0.5]] #UB + NVC>2 + XVC<50%
    t.append(t5.shape[1]-1)
    
    t4 = t3.loc[:, [0] + [i2 for i2 in t3.columns[1:] if min(t3[i2].value_counts()) > 2]] #UB + NVC>2
    t5 = t4.loc[:, [0] + [i2 for i2 in t4.columns[1:] if max(t4[i2].value_counts()) < t2.shape[0]*0.2]] #UB + NVC>2 + XVC<20%
    t.append(t5.shape[1]-1)

    t4 = t3.loc[:, [0] + [i2 for i2 in t3.columns[1:] if min(t3[i2].value_counts()) > 2]] #UB + NVC>3
    t5 = t4.loc[:, [0] + [i2 for i2 in t4.columns[1:] if max(t4[i2].value_counts()) < t2.shape[0]*0.15]] #UB + NVC>3 + XVC<15%
    t.append(t5.shape[1]-1)

print(' '.join([str(i) for i in t]))
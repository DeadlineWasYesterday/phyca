## nullify.py genome.fa full_table.tsv out.fa


import sys
import os
from Bio import SeqIO
import pandas as pd


t1 = pd.read_csv(sys.argv[2], sep = '\t', header = 0)
t1 = t1[(t1['Status'] == 'Single') | (t1['Status'] == 'Duplicated')]

      
t = []
for r in SeqIO.parse(sys.argv[1], 'fasta'):
    t2 = t1[t1['Sequence'] == r.id]
    ws = str(r.seq)
                
    before = ws.count('N')
                
    for s,e in t2[['Gene Start', 'Gene End']].astype(int).values:
        
        ws = ws[:min(s,e)] + 'N'* (max(s,e)-min(s,e)) + ws[max(s,e):]
                
    after = ws.count('N')

    #if after > before:
        #print(len(ws), before, after)
        #if (after - before) != ((t2['Gene End'] - t2['Gene Start']).astype(int)).sum():
            #print('Genes likely overlap')

    assert len(ws) == len(r.seq)

    t.append(('>'+r.id+'_null', ws))


with open(sys.argv[3], 'w') as f:
    for header, ws in t:
        f.write(header)
        f.write('\n')
        f.write(ws)
        f.write('\n')
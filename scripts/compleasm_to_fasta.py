#python mbfa.py full_table.tsv translated_protein.fa outdir prefix

import sys,os
import pandas as pd
from Bio import SeqIO

df = pd.read_csv(sys.argv[1], sep = '\t', header = 0)
df = df[(df['Status'] == 'Single') | (df['Status'] == 'Duplicated')]

df['GG'] = df['Best gene'] + '|' + df['Sequence'] + ':' + df['Gene Start'].astype(int).astype(str) + '-'+ df['Gene End'].astype(int).astype(str)

vals = df['GG'].values

t=[]
for r in SeqIO.parse(sys.argv[2], 'fasta'):
    if r.id in vals:
        t.append((r.id, str(r.seq)))

        
        
dfa = pd.DataFrame(t, columns = ['GG', 's']).drop_duplicates(['GG', 's'])
dfa = pd.merge(df, dfa[['GG', 's']], on = 'GG', how = 'left')
dfa = dfa.drop_duplicates(['GG'])
dfa = dfa.sort_values('Fraction',ascending = 0).drop_duplicates('Gene')
dfa['on'] = sys.argv[4]
dfa['h'] = '>' + dfa['on']

for i in dfa['Gene'].unique():
    assert dfa[dfa['Gene'] == i].shape[0] == 1
    dfa[dfa['Gene'] == i][['h','s']].dropna().to_csv('{0}/{1}.fa'.format(sys.argv[3],i), index = 0, header = None, sep = '\n')
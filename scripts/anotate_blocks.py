#annotate tree with BUSCO gene blocks

#usage
    #ant.py compiled_BUSCO_table compiled_null_table metadata treeleaves output

b = 2
    
import itertools, sys
import numpy as np
import pandas as pd

#read busco table
a1 = pd.read_csv(sys.argv[1]).drop_duplicates().reset_index(drop = 1)
a1 = a1[(a1['Status'] == 'Single') | (a1['Status'] == 'Duplicated')]
a1 = a1.drop('Unnamed: 0', axis = 1)
a1['Sequence'] = a1['Sequence']+a1['Assembly']
a1['Type'] = 't'
a1['pk1'] = a1['Sequence'] + '_' + a1['Gene Start'].astype(str)
a1['pk2'] = a1['Sequence'] + '_' + a1['Gene End'].astype(str)
a1 = a1.drop_duplicates('pk1').reset_index(drop = 1)
a1 = a1.drop_duplicates('pk2').reset_index(drop = 1)
a1 = a1.sort_values(['Assembly','Sequence','Gene Start']).reset_index(drop = 1)

#read nulls
x1 = pd.read_csv(sys.argv[2], sep = '\t').drop_duplicates().reset_index(drop = 1)
x1 = x1[(x1['Status'] == 'Single') | (x1['Status'] == 'Duplicated')]
x1['Sequence'] = x1['Sequence'].apply(lambda x: x.replace('_null',''))
x1['Sequence'] = x1['Sequence']+x1['Assembly']
x1['Type'] = 'n'
x1['pk1'] = x1['Sequence'] + '_' + x1['Gene Start'].astype(str)
x1['pk2'] = x1['Sequence'] + '_' + x1['Gene End'].astype(str)
x1 = x1.drop_duplicates('pk1').reset_index(drop = 1)
x1 = x1.drop_duplicates('pk2').reset_index(drop = 1)
x1 = x1.sort_values(['Assembly','Sequence','Gene Start']).reset_index(drop = 1)
#how many true counterparts the null gene has in the genome #exporting index in a1 was too slow
tmp2 = a1.groupby('Assembly')['Gene'].value_counts()
t1 = []
for a,g in x1[x1['Type'] == 'n'][['Assembly','Gene']].values:
    if (a,g) in tmp2.index:
        t1.append(tmp2.loc[a,g])
    else:
        t1.append(0)
x1['c'] = t1
x1 = pd.concat([a1,x1], axis = 0)
x1 = x1.sort_values(['Assembly','Sequence','Gene Start']).reset_index(drop = 1)

#compute blocks
def subs(block): #extract all consecutive subblocks
    t = []
    for i in range(len(block.split('_'))):
        for i2 in range(1+i,len(block.split('_'))+1):
            t.append('_'.join(block.split('_')[i:i2]))
    return t

def gc(block): #all possible gene combinations 2^n - 1
    block = block.split('_')
    t = []
    for L in range(len(block) + 1):
        for subset in itertools.combinations(block, L):
            if len(subset) > 0:
                t.append(list(subset))
    return t

d1 = a1.copy()
#append df in reverse orientation to make sure reverse blocks don't get lost
d1['Dir'] = 'f'
d1r = d1[::-1]
d1r['Dir'] = 'r'
d1 = pd.concat([d1,d1r],axis = 0).reset_index(drop = 1)

d1['gb'] = d1['Gene']
d1['ob'] = d1['Strand']
d1['tb'] = d1['Type']

#select most abundant block
for os in range(1,6): #block length
    d1['gb'] = d1['gb']+'_'+d1['Gene'].shift(-os)
    d2 = d1[d1['Sequence'] == d1['Sequence'].shift(-os)]    

qb = d2['gb'].value_counts().index[b] #highest incidence
print(qb)
print(d2['gb'].value_counts().iloc[b])

#use only contigs that have at least one gene from the block
y1 = x1.copy()
t = []
for i in qb.split('_'):
    t = t + list(set(x1[x1['Gene'] == g]['Sequence']))
d1 = d1[d1['Sequence'].isin(set(t))]

#append contigs containing null ##possible improvement: only append n genes left and right.
#computationally not possible for large blocks, ~>len(10)
c=0
for gs in gc(qb):
    n1 = y1.copy()
    for g in gs:
        n1 = n1[n1['Sequence'].isin(n1[n1['Gene'] == g]['Sequence'])]
        n1['Dir'] = 'f'     
    
    n1 = n1[(n1['Gene'].isin(gs)) | (n1['Type'] != 't')] #remove all null genes other than those in gs
    n1['Sequence'] = n1['Sequence'] + 'n' +str(c)
    n1r = n1[::-1]
    n1r['Dir'] = 'r'
    d1 = pd.concat([d1,n1,n1r],axis = 0).reset_index(drop = 1)
    c+=1
    
d1['c'] = d1['c'].fillna(0)

d1['gb'] = d1['Gene']
d1['ob'] = d1['Strand']
d1['tb'] = d1['Type']

for i in range(1,len(qb.split('_'))):
   
    d1['gb'] = d1['gb']+'_'+d1['Gene'].shift(-i)
    d1['ob'] = d1['ob']+d1['Strand'].shift(-i)
    d1['tb'] = d1['tb']+d1['Type'].shift(-i)


print(d1.shape)

d1 = d1[(d1['gb'] == qb) & (d1['Sequence'] == d1['Sequence'].shift(-i))]
    
    
print(d1.shape)
#metadata and tree
m1 = pd.read_csv(sys.argv[3], sep = '\t')
n1 = pd.read_csv(sys.argv[4], sep = '\t', header = None)

d1 = pd.merge(d1,m1, on = 'Assembly', how = 'left')

print(d1.shape)

gt = d1.groupby('genus')['tb'].value_counts().reset_index().drop_duplicates('genus')
gt = gt.rename(columns = {'count': 'type_count'})
go = d1.groupby('genus')['ob'].value_counts().reset_index().drop_duplicates('genus')
go = go.rename(columns = {'count': 'strand_count'})

print(gt.shape)
print(go.shape)

g = pd.merge(gt,go)[['genus','tb','ob', 'type_count', 'strand_count']]

n1 = n1[[1]]
n1.columns = ['genus']

pd.merge(n1,g,how = 'left').to_csv(sys.argv[5],index = 0, sep = '\t')

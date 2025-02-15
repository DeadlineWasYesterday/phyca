#merges ncbi assembly data json to taxonomy and exports as tsv
#usage python assemblydatajson_to_tsv.py assemblydatafile.jsonl taxtolin.csv output.tsv

import sys
import pandas as pd


z = pd.read_json(sys.argv[1], lines = True)

for i,r in z.iterrows():
    z.loc[i,'assemblyLevel'] = str(r['assemblyInfo'].get('assemblyLevel'))
    z.loc[i,'assemblyMethod'] = str(r['assemblyInfo'].get('assemblyMethod'))
    z.loc[i,'assemblyName'] = str(r['assemblyInfo'].get('assemblyName'))
    z.loc[i,'assemblyStatus'] = str(r['assemblyInfo'].get('assemblyStatus'))
    z.loc[i,'assemblyType'] = str(r['assemblyInfo'].get('assemblyType'))
    z.loc[i,'bioprojectAccession'] = str(r['assemblyInfo'].get('bioprojectAccession'))
    z.loc[i,'refseqCategory'] = str(r['assemblyInfo'].get('refseqCategory'))
    z.loc[i,'releaseDate'] = str(r['assemblyInfo'].get('releaseDate'))
    z.loc[i,'sequencingTech'] = str(r['assemblyInfo'].get('sequencingTech'))
    z.loc[i,'organismName'] = str(r['organism'].get('organismName'))
    z.loc[i,'taxId'] = str(r['organism'].get('taxId'))
    z.loc[i,'totalSequenceLength'] = str(r['assemblyStats'].get('totalSequenceLength'))
    z.loc[i,'totalNumberOfChromosomes'] = str(r['assemblyStats'].get('totalNumberOfChromosomes'))
    z.loc[i,'totalUngappedLength'] = str(r['assemblyStats'].get('totalUngappedLength'))
    z.loc[i,'contigL50'] = str(r['assemblyStats'].get('contigL50'))
    z.loc[i,'contigN50'] = str(r['assemblyStats'].get('contigN50'))
    z.loc[i,'gcCount'] = str(r['assemblyStats'].get('gcCount'))
    z.loc[i,'numberOfComponentSequences'] = str(r['assemblyStats'].get('numberOfComponentSequences'))
    z.loc[i,'numberOfContigs'] = str(r['assemblyStats'].get('numberOfContigs'))
    z.loc[i,'numberOfScaffolds'] = str(r['assemblyStats'].get('numberOfScaffolds'))
    z.loc[i,'scaffoldL50'] = str(r['assemblyStats'].get('scaffoldL50'))
    z.loc[i,'scaffoldN50'] = str(r['assemblyStats'].get('scaffoldN50'))
    z.loc[i,'numberOfScaffolds'] = str(r['assemblyStats'].get('numberOfScaffolds'))
    #z.loc[i,'stats'] = str(r['annotationInfo'].get('stats'))

y = pd.read_csv(sys.argv[2])

z['taxId'] = z['taxId'].astype(int)
z2 = z.merge(y[['tax_id', 'superkingdom', 'kingdom', 'subkingdom', 'superphylum', 'phylum', 'subphylum', 'superclass', 'class', 'subclass', 'superorder', 'order', 'suborder', 'superfamily', 'family', 'subfamily', 'genus', 'subgenus', 'species', 'subspecies']], left_on = 'taxId', right_on = 'tax_id')


print(z2['kingdom'].unique())
print(z2['kingdom'].value_counts())
print(z2['assemblyLevel'].value_counts())


z2['on'] = z2['genus'] + '_' + z2['species'] 


for i,r in z2.iterrows():
    if r['assemblyLevel'] == 'Complete Genome':
        z2.loc[i,'NumLev'] = 1
    elif r['assemblyLevel'] == 'Chromosome':
        z2.loc[i,'NumLev'] = 2
    elif r['assemblyLevel'] == 'Scaffold':
        z2.loc[i,'NumLev'] = 3
    elif r['assemblyLevel'] == 'Contig':
        z2.loc[i,'NumLev'] = 4


z2.to_csv(sys.argv[3], sep = '\t', index = 0)
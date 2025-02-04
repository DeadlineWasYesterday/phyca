#!/usr/bin/python3

import os, sys, shutil, argparse, subprocess, urllib, gzip, warnings
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import BioNick as bn
import matplotlib.ticker as ticker
from pkg_resources import resource_filename
from io import StringIO
from Bio import Phylo, SeqIO
from scipy.optimize import curve_fit
from Bio.Phylo.TreeConstruction import DistanceCalculator

from pandas.errors import SettingWithCopyWarning
warnings.simplefilter(action="ignore", category=SettingWithCopyWarning)

from scipy.optimize import OptimizeWarning
warnings.simplefilter(action="ignore", category=OptimizeWarning)



lineages = ['arthropoda','ascomycota','basidiomycota','chlorophyta','eudicots','fungi','liliopsida','metazoa','vertebrata','viridiplantae']


#remove when BioNick is updated
#####################
class node:
    def __init__(self, name, connections, branch_lengths):
        self.name = name
        self.connections = connections
        self.branch_lengths = branch_lengths
        self.num_con = len(self.connections)
        self.cbpairs = list(zip(connections,branch_lengths))
    def add_connection(self, connection, branch_length):      
        self.connections.append(connection)
        self.branch_lengths.append(branch_length)
        self.num_con = len(self.connections)
    def remove_connection(self,parent):
        self.cbpairs = [(x,y) for x,y in zip(self.connections,self.branch_lengths) if x!=parent]
        if len(self.connections) > 1:
            self.connections = list(list(zip(*self.cbpairs))[0])
            self.branch_lengths = list(list(zip(*self.cbpairs))[1])
        else:
            self.connections = []
            self.branch_length = []
    def set_parent(self,parent):
        self.parent = parent
    def expand(self):
        tmp1 = [a if isinstance(a,str) else 'intnode'+'%05.f'%a for a in self.connections]
        tmp2 = [str(a)+':'+'%g'%round(b,8) for a,b in zip(tmp1,self.branch_lengths)]
        return '('+','.join(tmp2)+')'
    def expand_rev(self):
        tmp1 = [a if isinstance(a,str) else 'intnode'+'%05.f'%a for a in self.connections]
        tmp2 = [str(a)+':'+'%g'%round(b,8) for a,b in zip(tmp1,self.branch_lengths)]
        return '('+','.join(tmp2[::-1])+')'

class newick:
    def __init__(self,text):
        self.text = text
        self.leaves = [x.split(':')[0].replace('(','') for x in text.split(',')]

class tree:
    def __init__(self, nodes):
        self.nodes = nodes

    def add_node(self,node):
        self.nodes.append(node)

    def remove_node(self,name):
        self.nodes = [x for x in self.nodes if x.name!=name]

    def get_node(self, name):
        for i in self.nodes:
            if i.name == name:
                return i
        return None

    def remove_biconnection(self,name1,name2):
        self.get_node(name1).remove_connection(name2)
        self.get_node(name2).remove_connection(name1)
        
    def num_tips(self):
        c=0
        for i in self.nodes:
            if len(i.connections) == 1:
                c+=1
        return c

    def num_internal_nodes(self):
        return self.num_nodes() - self.num_tips()
        
    def num_nodes(self):
        return (len(self.nodes))
    
    def listnodes(self):
        return [x.name for x in self.nodes]
    def listtips(self):
        return [i.name for i in self.nodes if len(i.connections) == 1]
    
    def unresolved_nodes(self):
        for i in self.nodes:
            if len(i.connections) != 1 and len(i.connections) != 3: 
                print(i.name)
                print(i.connections)
            
    def root_at_tip(self,tip):
        assert len(self.get_node(tip).connections) == 1
        self.root = self.get_node(tip).connections[0]
    def root_at_node(self,nodename):
        self.root = nodename
        pass

    def export_nw(self,nt,parent):
        if nt == '':
            #initiate from root
            nt = self.get_node(self.root).expand()
            parent = self.root
        #check for internal nodes in nt
        for i in newick(nt).leaves:
            if 'intnode' in i:
                node_label = int(i.replace('intnode',''))
                self.remove_biconnection(parent,node_label)
                nt = nt.replace(i,self.get_node(node_label).expand())
                for i2 in self.get_node(node_label).connections: #remove all connections to and from expanded node
                    if not isinstance(i2,str):
                        self.remove_biconnection(node_label,i2)
                #print(nt,i,node_label,parent)
                #return self.export_nw(nt,node_label)
        for i in newick(nt).leaves:
            if 'intnode' in i:
                node_label = int(i.replace('intnode',''))
                #print(nt,i,node_label,parent)
                return self.export_nw(nt,node_label)
        return nt
#####################



def format_stats(df,idf,total):
    #comprehensive
    odf = df['Gene'].value_counts().reset_index()['count'].value_counts().sort_index()
    odf.loc[0] = idf[~idf['Status'].isin(["Single","Duplicated"])].shape[0]
    odf.index.name = 'Copies'
    odf = odf.reset_index().sort_values('Copies')
    odf['Percentile'] = odf['count'].div(total)
    odf = odf.rename(columns = {'count':'Count'})

    #lame
    ldf = pd.DataFrame([('Complete', odf.loc[odf['Copies'] > 0,'Count'].sum(), odf.loc[odf['Copies'] > 0,'Percentile'].sum()),
                  ('Single', odf.loc[odf['Copies'] == 1,'Count'].sum(), odf.loc[odf['Copies'] == 1,'Percentile'].sum()),
                  ('Duplicated', odf.loc[odf['Copies'] > 1,'Count'].sum(), odf.loc[odf['Copies'] > 1,'Percentile'].sum()),
                  ('Fragmented', idf[idf['Status']=='Fragmented'].shape[0], idf[idf['Status']=='Fragmented'].shape[0]/total),
                  ('Interspaced',idf[idf['Status']=='Interspaced'].shape[0], idf[idf['Status']=='Interspaced'].shape[0]/total),
                  ('Missing',idf[idf['Status']=='Missing'].shape[0], idf[idf['Status']=='Missing'].shape[0]/total)],
                      columns = ['State','Count','Percentile'])

    return [odf, ldf]

def plot_bars(comp,lame,args):
    fig,ax = plt.subplots(1,2, figsize=(18, 8), gridspec_kw={'width_ratios': [4,12], 'wspace': 0.1})
    
    ax[0].bar(['Complete','Missing'],lame[lame['Type'] == 'CUSCO']['Percentile'].iloc[[1,5]]*100)
    ax[0].bar(['Complete'],lame[lame['Type'] == 'CUSCO']['Percentile'].iloc[[2]]*100, bottom = lame[lame['Type'] == 'CUSCO']['Percentile'].iloc[[1]]*100)
    ax[0].legend(['Single','Duplicated'])
    ax[0].set_xlabel('State',fontsize = 12)
    ax[0].set_ylabel('Percentile',fontsize = 12)
    ax[1].xaxis.set_major_locator(ticker.MaxNLocator(integer=True))
    ax[1].bar(comp[comp['Type'] == 'CUSCO']['Copies'],comp[comp['Type'] == 'CUSCO']['Percentile']*100)
    ax[1].set_xlabel('Copies',fontsize = 12)
    ax[1].set_ylabel('Percentile',fontsize = 12)
    plt.savefig('{0}/USCO_bars.pdf'.format(args.output), format = 'pdf')

def process_busco(bdf,cus,args):
    #remove and report bad hits
    idf = bdf.copy()
    btot,ctot = idf['Gene'].nunique(), len(set(idf['Gene'].unique()).intersection(set(cus['Gene'])))
    mtot = btot-ctot
    bdf = bdf[bdf['Status'].isin(["Single","Duplicated"])]
    bdf.loc[:,'fl'] = bdf.apply(lambda c: min(c['Gene Start'], c['Gene End']), axis = 1)
    bdf.loc[:,'fr'] = bdf.apply(lambda c: max(c['Gene Start'], c['Gene End']), axis = 1)
    bdf = bdf.sort_values(['Sequence', 'fl'])
    bdf['g-1'] = bdf['Gene'].shift(-1)
    bdf['fl-1'] = bdf['fl'].shift(-1)
    bdf['fso'] = bdf['fl-1'] - bdf['fr']
    bdf['s-1'] = bdf['Sequence'].shift(-1)
    df = bdf[~((bdf['Sequence'] == bdf['s-1']) & (bdf['Gene'] == bdf['g-1']) & (bdf['fso'] < 0))]

    if bdf[((bdf['Sequence'] == bdf['s-1']) & (bdf['fso'] < 0))].shape[0] != 0:
        tmp = bdf[((bdf['Sequence'] == bdf['s-1']) & (bdf['fso'] < 0))].shape[0]
        print("{0} BUSCO genes have overlapping sequences. Annotations of these genes may be questionable.".format(tmp*2))
    
    if bdf[((bdf['Sequence'] == bdf['s-1']) & (bdf['Gene'] == bdf['g-1']) & (bdf['fso'] < 0))].shape[0] != 0:
        tmp = bdf[((bdf['Sequence'] == bdf['s-1']) & (bdf['Gene'] == bdf['g-1']) & (bdf['fso'] < 0))].shape[0]
        print("Removed {0} identical BUSCO genes that were found in tandem with overlapping sequences.".format(tmp))
    
    #output stats
    ob,lb = format_stats(df,idf,btot)
    oc,lc = format_stats(df[df['Gene'].isin(cus['Gene'])], idf[idf['Gene'].isin(cus['Gene'])], ctot)
    om,lm = format_stats(df[~(df['Gene'].isin(cus['Gene']))], idf[~(idf['Gene'].isin(cus['Gene']))], mtot)
    ob['Type'],oc['Type'],om['Type'] = 'BUSCO','CUSCO','MUSCO'
    lb['Type'],lc['Type'],lm['Type'] = 'BUSCO','CUSCO','MUSCO'
    comp = pd.concat([ob,oc,om],axis=0)[['Type','Copies','Count','Percentile']]
    lame = pd.concat([lb,lc,lm],axis=0)[['Type','State','Count','Percentile']]
    comp.to_csv('{0}/USCO_copies.tsv'.format(args.output), sep = '\t', index = 0)
    lame.to_csv('{0}/USCO_stats.tsv'.format(args.output), sep = '\t', index = 0)
    
    #output graphs
    plot_bars(comp,lame,args)
    


#remove busco genes from assembly

def nullify(args):
    if args.compdir is None:
        print("Please specify compleasm output directory with -c.")
        sys.exit()
    
    t1 = pd.read_csv(args.compdir+'/'+args.lineage+'_odb10'+'/full_table.tsv', sep = '\t', header = 0)
    t1 = t1[(t1['Status'] == 'Single') | (t1['Status'] == 'Duplicated')]   

    t = []
    for r in SeqIO.parse(args.assembly, 'fasta'):
        t2 = t1[t1['Sequence'] == r.id]
        ws = str(r.seq)                    
        before = ws.count('N')
        for s,e in t2[['Gene Start', 'Gene End']].astype(int).values:
            
            ws = ws[:min(s,e)] + 'N'* (max(s,e)-min(s,e)) + ws[max(s,e):]
        after = ws.count('N')
        assert len(ws) == len(r.seq)
        t.append(('>'+r.id+'_null', ws))

    with open('{0}/{1}.null'.format(args.output,args.assembly), 'w') as f:
        for header, ws in t:
            f.write(header)
            f.write('\n')
            f.write(ws)
            f.write('\n')



#find the best reference
#export a synteny tree

#syntenic distance matrix from processed df
def syndm(jn):
    t = []
    ij = 0
    for ib,ic in jn.values:
        t2 = []
        ij2 = 0
        for ig,ih in jn.values:
            if ij2 < ij:
                t2.append(0)
                ij2+=1
                continue
    
            t2.append((len(ic.intersection(ih)) * 10000)/ len(ic.union(ih)))
            
            ij2+=1
        t.append([ib] + t2)
        ij+=1
    return t

#format distance matrix
def pdm(ba):
    ca = pd.DataFrame(ba)
    labs = ca[0].values
    ca = ca.iloc[:,1:]
    ca = pd.DataFrame(ca.values + ca.values.T - np.diag(np.diag(ca.values)))
    ca.index,ca.columns = labs,labs
    return ca

#distance matrix to newick
def nj(dm,sq):
    
    if dm.shape[0] == 2:
        sq.append((0,0,dm.columns[1],dm.iloc[0,1]))
        return sq
    n=dm.shape[0]
    #nj matrix
    dn = (dm*(n-2)).subtract(dm.sum().values,axis=0).subtract(dm.sum().values,axis=1)
    np.fill_diagonal(dn.values,0)
    #print(dn)
    #indices
    i,j = np.unravel_index(np.argmin(dn.values), dn.shape)
    #branches
    ib = 0.5*(dm.iloc[i,j]+abs(dm.sum(axis=0).iloc[i] - dm.sum(axis=0).iloc[j])/(dm.shape[0]-2))
    jb = 0.5*(dm.iloc[i,j]-abs(dm.sum(axis=0).iloc[i] - dm.sum(axis=0).iloc[j])/(dm.shape[0]-2))
    
    #update matrix
    t = []
    for k in range(n):
        t.append((dm.iloc[i,k]+dm.iloc[k,j]-dm.iloc[i,j])/2)

    dm.loc[len(sq)+10000,:] = t
    t.append(0)
    #print(dm)
    dm.loc[:,len(sq)+10000] = t
    #print(dm)

    #sequence
    sq.append((dm.index[i],ib,dm.columns[j],jb))
    
    dm = dm.drop(index = [dm.index[i],dm.index[j]], columns = [dm.columns[i],dm.columns[j]])
    np.fill_diagonal(dm.values,0)
    #print(dm)

    return nj(dm.copy(),sq)


def njtr(vt):
    tt = tree([])
    i=0
    for a,b,c,d in vt.values[:-1]:
        if not isinstance(a,str): #some kind of strange int vs numpyint instance issue on the second to last df entry
            a = a-10000
            nodet = tt.get_node(a)
            nodet.add_connection(i,b)
        else: #tip
            tt.add_node(node(a,[i],[b]))
    
        if not isinstance(c,str):
            c = c-10000
            nodet = tt.get_node(c)
            nodet.add_connection(i,d)
        else: #tip
            tt.add_node(node(c,[i],[d]))
        
        node1 = node(i,[],[])    
        node1.add_connection(a,b)
        node1.add_connection(c,d)
        tt.add_node(node1)
        i+=1
    
    #last connection
    a,b,c,d = vt.values[-1]
    if not isinstance(c,str): #internal #connect last two (?)
        tt.get_node(c-10000).add_connection(i-2,d)
        tt.get_node(i-2).add_connection(c-10000,d) #reciprocate
    else: #tip #connect to last
        tt.add_node(node(c,[i-1],[d]))
        tt.get_node(i-1).add_connection(c,d)
    return tt


#syntenice distance between two assemblies
def syndis(cmp1,cmp2,args):
    df1 = pd.read_csv('{0}/{1}_odb10/full_table.tsv'.format(cmp1, args.lineage), sep = '\t')
    df2 = pd.read_csv('{0}/{1}_odb10/full_table.tsv'.format(cmp2, args.lineage), sep = '\t')

    df1['Assembly'] = 'query'
    df2['Assembly'] = 'reference'

    df = pd.concat([df1,df2],axis = 0).reset_index(drop=True)

    j = syntree(df,args)
    
    return (len(j.values[0,1].intersection(j.values[1,1])) * 10000)/len(j.values[0,1].union(j.values[1,1]))



def syntree(df, args):
    df = df[df['Status'].isin(["Single","Duplicated"])]
    #left and right flank
    df.loc[:,'fl'] = df.apply(lambda c: min(c['Gene Start'], c['Gene End']), axis = 1)
    df.loc[:,'fr'] = df.apply(lambda c: max(c['Gene Start'], c['Gene End']), axis = 1)
    df = df.sort_values(['Assembly', 'Sequence', 'fl'])
    df['g-1'] = df['Gene'].shift(-1)
    df['fl-1'] = df['fl'].shift(-1)
    df['o-1'] = df['Strand'].shift(-1)
    df['s-1'] = df['Sequence'].shift(-1)
    df['fso'] = df['fl-1'] - df['fr']
    df = df[~((df['Sequence'] == df['s-1']) & (df['Gene'] == df['g-1']) & (df['fso'] < 0))]
    df['as'] = df['Assembly'] + '_' + df['Sequence']

    #do not remove singleton contigs if flag is passed to parser
    if not args.include_singleton_contigs:
        df = df[~df['as'].isin(df.groupby('as')['Gene'].count()[df.groupby('as')['Gene'].count() == 1].index)]

    #gene doublet
    df.loc[df['g-1'].notnull(), 'gd'] = df[df['g-1'].notnull()].apply(lambda c: str(sorted([c['Gene'], c['g-1']])), axis = 1)
    #doublet orientation
    df.loc[df['o-1'].notnull(), 'do'] = df[df['o-1'].notnull()].apply(lambda c: 's' if c['Strand'] == c['o-1'] else 'o', axis = 1)

    df = df[df['gd'].notnull()]
    df['agd'] = df['Assembly']+'_'+df['gd']
    #if ignoring orientation
    if args.ignore_orientation:
        #if including duplications
        if args.include_duplications:    
            df.loc[df['agd'].duplicated(keep=False), 'agd'] += df.groupby('agd').cumcount().add(1).astype(str)
        df['agdar'] = df['agd'].apply(lambda x: x.split('_')[-1])
        j = df.groupby('Assembly')['agdar'].agg(set).reset_index()
    else:
        df['agdo'] = df['agd']+df['do']
        #if including duplications
        if args.include_duplications:  
            df.loc[df['agdo'].duplicated(keep=False), 'agdo'] += df.groupby('agdo').cumcount().add(1).astype(str)
        df['agdoar'] = df['agdo'].apply(lambda x: x.split('_')[-1])
        j = df.groupby('Assembly')['agdoar'].agg(set).reset_index()

    return j
    
    
#use sepp to rapidly place on precomputed phylogeny  
#export synteny decay graph

def match_taxa(alnfile,treefile,outfilename):
    with open(treefile,'r') as f:
        tree = f.read()
    t=[]
    with gzip.open(alnfile, "rt") as handle:
        for entry in SeqIO.parse(handle,'fasta'):
            t.append((entry.id,str(entry.seq).replace('\n','')))
    aln = pd.DataFrame(t)
    leaves = bn.leaves(tree)
    headers = aln[0]
    shared = set(leaves).intersection(set(headers))
    treeout = bn.extract_subtree(bn.remove_node_labels(tree), shared)
    alnout = aln[aln[0].isin(shared)]

    with open(outfilename+'.tree','w') as f:
        f.write(treeout+';')
    
    alnout[0] = '>'+alnout[0]
    alnout.to_csv(outfilename+'.aln', sep ='\n', index = 0, header = None)

#extract busco gene sequence
def buscoseq(cmpdir,d2,gene,args):
    d2 = d2[(d2['Status'] == 'Single') | (d2['Status'] == 'Duplicated')].reset_index(drop = 1)
    #primary key to merge sequences
    d2['GG'] = d2['Best gene'] + '|' + d2['Sequence'] + ':' + d2['Gene Start'].astype(int).astype(str) + '-'+ d2['Gene End'].astype(int).astype(str)
    #remove alternate codon configurations
    d2 = d2.drop_duplicates('GG')
    tmp = d2.loc[d2['Gene'] == gene, 'GG'].iloc[0]
    t = []
    for r in SeqIO.parse(cmpdir+'/'+args.lineage+'_odb10/translated_protein.fasta', 'fasta'):
        if r.id == tmp:
            t.append((r.id, str(r.seq)))

    with open(args.output+'/temp.seq', 'w') as f:
        f.write('>'+'query_{0}'.format(gene))
        f.write('\n')
        f.write(t[0][1])



def phyca(args):

    if (args.assembly is None) and (args.compdir is None):
        print("Please specify either assembly with -a or Compleasm output directory with -c.")
        sys.exit()
        
    if not (args.lineage in lineages):
        print("Invalid BUSCO lineage. Supported lineages are (case-sensitive):")
        print(lineages)
        sys.exit()

    if not os.path.isdir(args.output):
        os.makedirs(args.output)

    # Nullify
    if args.nullify:
        nullify(args)
        print("BUSCO-depleted assembly saved at {0}/{1}.null".format(args.output,args.assembly))
        sys.exit()

    # Calculate distance
    if args.syndis:
        if (args.reference is None) and (args.rcompdir is None):
            print("Please either specify reference with -r or reference compleasm output directory with -m.")
            sys.exit()

        if args.compdir is None:
            if os.path.exists('cmout') and os.path.isdir('cmout'):
                print('Compleasm output from a previous run exists in ./cmout and will be overwritten.')
                shutil.rmtree('cmout')
            print('Compleasm run command: compleasm run -a {0} -l {1} -t {2} -o cmout'.format(args.assembly, args.lineage, str(args.threads)))
            cmplog = subprocess.run(["compleasm", "run", "-a", args.assembly, '-l', args.lineage, '-t', str(args.threads), '-o', 'cmout'],
                          capture_output = True, text = True)
            with open('{0}/compleasm.log'.format(args.output),'w') as f:
                f.write(cmplog.stdout)
            print('Compleasm run completed. Log saved in {0}/compleasm.log'.format(args.output))
            
            cmp1 = 'cmout'
        else:
            cmp1 = args.compdir

        if args.rcompdir is None:
            if os.path.exists('cmref') and os.path.isdir('cmref'):
                print('Compleasm output from a previous run exists in ./cmref and will be overwritten.')
                shutil.rmtree('cmref')
            print('Compleasm run command: compleasm run -a {0} -l {1} -t {2} -o cmref'.format(args.reference, args.lineage, str(args.threads)))
            cmplog = subprocess.run(["compleasm", "run", "-a", args.reference, '-l', args.lineage, '-t', str(args.threads), '-o', 'cmref'],
                          capture_output = True, text = True)
            with open('{0}/rcompleasm.log'.format(args.output),'w') as f:
                f.write(cmplog.stdout)
            print('Compleasm run completed. Log saved in {0}/rcompleasm.log'.format(args.output))
            cmp2 = 'cmref'
        else:
            cmp2 = args.rcompdir
        
        print("Syntenic identity (%): ", syndis(cmp1,cmp2,args)/100)
        sys.exit()      


    ### phyca
    #remove previous output
    if os.path.exists(args.output):
        print('phyca output from a previous run exists and will be overwritten. Use the -o option to set a new output directory.')
        shutil.rmtree(args.output)
    #create output directory   
    os.makedirs("{0}".format(args.output), exist_ok=True)
    if not args.reference is None: 
        pass

    cus = pd.read_csv(resource_filename(__name__, "CUS.tsv"), sep = '\t')
    if (not args.assembly is None) and (args.compdir is None): 

        if os.path.exists('cmout') and os.path.isdir('cmout'):
            print('Compleasm output from a previous run exists and will be overwritten.')
            shutil.rmtree('cmout')
        print('Compleasm run command: compleasm run -a {0} -l {1} -t {2} -o cmout'.format(args.assembly, args.lineage, str(args.threads)))
        cmplog = subprocess.run(["compleasm", "run", "-a", args.assembly, '-l', args.lineage, '-t', str(args.threads), '-o', 'cmout'],
                      capture_output = True, text = True)
        with open('{0}/compleasm.log'.format(args.output),'w') as f:
            f.write(cmplog.stdout)
        
        print('Compleasm run completed. Log saved in {0}/compleasm.log'.format(args.output))
        
        
        bdf = pd.read_csv('cmout/{0}_odb10/full_table.tsv'.format(args.lineage), sep = '\t')

        cmpdir = 'cmout'

    
    if not args.compdir is None:
        cmpdir = args.compdir
        bdf = pd.read_csv('{0}/{1}_odb10/full_table.tsv'.format(args.compdir, args.lineage), sep = '\t')

    process_busco(bdf,cus,args)

    
    #find the best reference
    #export a synteny tree
    fname = "{0}_ch.tsv.gz".format(args.lineage[:2])
    if not os.path.isfile(fname):
        print("Downloading chromosome-level BUSCO annotation table for {0}.".format(args.lineage))
        urllib.request.urlretrieve("https://ava.genome.arizona.edu/UniPhy/annotations/"+fname,fname)
    chdf = pd.read_csv(fname,sep = '\t', compression = 'gzip')
    
    fname = "{0}.tsv".format(args.lineage[:2])
    if not os.path.isfile(fname):
        print("Downloading genome metadata table for {0}.".format(args.lineage))
        urllib.request.urlretrieve("https://ava.genome.arizona.edu/UniPhy/metadata/"+fname,fname)
    mdf = pd.read_csv(fname,sep = '\t')

    mdf['label'] = mdf['on']+'_'+mdf['Assembly']
    chdf['Assembly'] = chdf['Assembly'].map(mdf.set_index('Assembly')['label'])
    
    bdf['Assembly'] = 'query'
    chdf = pd.concat([chdf,bdf],axis = 0).reset_index(drop=True)
    j = pdm(syndm(syntree(chdf,args)))

    #print best references
    print('Query most syntenic to the following references:')
    print(j[j['query']>0]['query'].sort_values()[::-1].iloc[1:6].reset_index().rename(columns = {'index':'genome','query':'syntenic_identity_score'}))
    j['query'].sort_values()[::-1].reset_index().rename(
        columns = {'index':'genome','query':'syntenic_identity_score'}).to_csv('{0}/SynIdentity.tsv'.format(args.output), sep = '\t', index = 0)
    
    xa = 10000 - j.copy() #similarity matrix to distance matrix
    tt = njtr(pd.DataFrame(nj(xa.copy(),[])))
    tt.root_at_tip('query')
    stre = tt.export_nw('','')

    #write synteny tree
    with open('{0}/synteny.tree'.format(args.output), 'w') as f:
        f.write(stre)

    #print synteny tree
    twn = stre
    from matplotlib.pyplot import figure
    figure(figsize=(max(5,len(bn.leaves(twn))/12), max(10,len(bn.leaves(twn))/5)), dpi=100)
    
    x = bn.draw_clad(bn.remove_node_labels(twn), dash = True, labels = True)
    plt.ylim(-1,len(bn.leaves(twn))+1)
    plt.gca().spines[['left','right', 'top']].set_visible(False)
    plt.gca().get_yaxis().set_visible(False)
    plt.xlabel('Substitutions/Site')
    plt.savefig('{0}/synteny_tree_dashed.pdf'.format(args.output), format = 'pdf', bbox_inches='tight')

    from matplotlib.pyplot import figure
    figure(figsize=(max(5,len(bn.leaves(twn))/12), max(10,len(bn.leaves(twn))/5)), dpi=100)
    
    x = bn.draw_clad(bn.remove_node_labels(twn), dash = False, labels = True)
    plt.ylim(-1,len(bn.leaves(twn))+1)
    plt.gca().spines[['left','right', 'top']].set_visible(False)
    plt.gca().get_yaxis().set_visible(False)
    plt.xlabel('Substitutions/Site')
    plt.savefig('{0}/synteny_tree.pdf'.format(args.output), format = 'pdf', bbox_inches='tight')

    
       
    #use sepp to rapidly place on precomputed phylogeny  
    #download alignment stats file
    if not os.path.isfile('aln_stats.tsv'):
        print("Downloading alignment stats table.")
        urllib.request.urlretrieve("https://ava.genome.arizona.edu/UniPhy/alignments/aln_stats.tsv","aln_stats.tsv")
    adf = pd.read_csv('aln_stats.tsv',sep = '\t')


    if args.lineage == 'metazoa':
        adf = adf[(adf['Lineage'] == args.lineage.capitalize()) & (adf['Type'] == 'CUSCO') & (adf['Mean Length'] > 350) ].sort_values('Coefficient of Variation')
    else:
        adf = adf[(adf['Lineage'] == args.lineage.capitalize()) & (adf['Type'] == 'CUSCO') & (adf['Mean Length'] > 700) ].sort_values('Coefficient of Variation')

    #make sure at least one placement gene is present in assembly  
    if len(set(bdf['Gene']).intersection(set(adf['Gene'].iloc[:20]))) == 0:
        print('Not enough BUSCO genes for rapid placement.')
        sys.exit()

    gene = adf[adf['Gene'].isin(bdf['Gene'])]['Gene'].iloc[0]
    #download stock gene alignment
    fname = '{0}.afa.gz'.format(gene)
    if not os.path.isfile(fname):
        print("Downloading stock alignment file.")
        urllib.request.urlretrieve("https://ava.genome.arizona.edu/UniPhy/alignments/{0}1/{1}".format(args.lineage[:2],fname),fname)

    #download tree
    fname = '{0}top40k.afa_1.treefile'.format(args.lineage[:2])
    if not os.path.isfile(fname):
        print("Downloading stock tree file.")
        urllib.request.urlretrieve("https://ava.genome.arizona.edu/UniPhy/trees/{0}".format(fname),fname)

    #download ml parameters file
    fname = 'RAxML_info.{0}'.format(gene)
    if not os.path.isfile(fname):
        print("Downloading MLinfo file.")
        urllib.request.urlretrieve("https://ava.genome.arizona.edu/UniPhy/mlinfo/{0}".format(fname),fname)

    #make sure leaves and headers are identical
    match_taxa('{0}.afa.gz'.format(gene), '{0}top40k.afa_1.treefile'.format(args.lineage[:2]), 'temp')

    buscoseq(cmpdir,bdf,gene,args)

    #run Sepp
    #run_sepp.py -t vi87347at33090.tree -r RAxML_info.87347at33090 -a vi87347at33090.aln -f testfraggene -m amino -o wtft
    # if os.path.isfile('{0}/seppout_placement.json'.format(args.output)):
    #     print("Removing previous seppout in {0}".format(args.output))

    print('Running SEPP. Cmd: '+' '.join(["run_sepp.py", "-t", "temp.tree", "-r", 'RAxML_info.{0}'.format(gene), 
                    '-a', "temp.aln", '-m', 'amino', '-f', args.output+"/temp.seq", '-d', args.output, '-o', 'seppout']))
    
    print('SEPP log will be saved at {0}/sepp.log'.format(args.output))
    
    sepplog = subprocess.run(["run_sepp.py", "-t", "temp.tree", "-r", 'RAxML_info.{0}'.format(gene), 
                    '-a', "temp.aln", '-m', 'amino', '-f', args.output+"/temp.seq", '-d', args.output, '-o', 'seppout'],
                  capture_output = True, text = True)

    with open('{0}/sepp.log'.format(args.output),'w') as f:
        f.write(sepplog.stdout)

    
    print('Running guppy. Cmd: guppy tog {0}/seppout_placement.json'.format(args.output))
    #run Guppy
    #guppy-64 tog wtfn_placement.json
    subprocess.run(["guppy", "tog", "{0}/seppout_placement.json".format(args.output)])
    shutil.move("seppout_placement.tog.tre", "{0}/placement.tree".format(args.output))


    print('Placement tree created.')
    
    #print placement tree

    with open("{0}/placement.tree".format(args.output), 'r') as f:
        twn = f.readlines()[0]
    from matplotlib.pyplot import figure
    figure(figsize=(max(5,len(bn.leaves(twn))/12), max(10,len(bn.leaves(twn))/5)), dpi=100)
    
    x = bn.draw_clad(bn.remove_node_labels(twn), dash = True, labels = True)
    plt.ylim(-1,len(bn.leaves(twn))+1)
    plt.gca().spines[['left','right', 'top']].set_visible(False)
    plt.gca().get_yaxis().set_visible(False)
    plt.xlabel('Substitutions/Site')
    plt.savefig('{0}/placement_tree_dashed.pdf'.format(args.output), format = 'pdf', bbox_inches='tight')

    
    from matplotlib.pyplot import figure
    figure(figsize=(max(5,len(bn.leaves(twn))/12), max(10,len(bn.leaves(twn))/5)), dpi=100)
    
    x = bn.draw_clad(bn.remove_node_labels(twn), dash = False, labels = True)
    plt.ylim(-1,len(bn.leaves(twn))+1)
    plt.gca().spines[['left','right', 'top']].set_visible(False)
    plt.gca().get_yaxis().set_visible(False)
    plt.xlabel('Substitutions/Site')
    plt.savefig('{0}/placement_tree.pdf'.format(args.output), format = 'pdf', bbox_inches='tight')


    print("Exporting synteny decay graph.")
    
    #export synteny decay graph
    figure(figsize=(6,6), dpi=100)
    p1 = Phylo.read('{0}/placement.tree'.format(args.output), 'newick')
    
    t = []
    for i in p1.get_terminals():
        t.append((i.name, p1.distance('query_{0}'.format(gene),i.name)))
    
    tmp = j.reset_index()[['index','query']].rename(columns = {'index':0})
    tmp[0] = tmp[0].apply(lambda x: x.split('_G')[0])
    rn = pd.merge(pd.DataFrame(t),tmp, how = 'inner')    
    xdata = rn[1].values
    ydata = rn['query'].values/10000
    plt.scatter(xdata, ydata, s = 5, color = 'blue')
    def func(x, a, b, c):
        return a * np.exp(-b * x) #- c*x
    popt, pcov = curve_fit(func, xdata, ydata)
    popt, pcov = curve_fit(func, xdata, ydata, bounds=(0, [10., 20., 1]), method = 'dogbox')
    plt.plot(#np.arange(min(xdata),max(xdata),(max(xdata)-min(xdata))/100), 
               #func(np.arange(min(xdata),max(xdata),(max(xdata)-min(xdata))/100), *popt), 
               np.arange(0,5,5/100), 
               func(np.arange(0,5,5/100), *popt), 
               '-', label='Fit: y=%5.3fe$^{%5.3fx}$' % tuple(popt)[:2], color = 'green')

    #ax[3].set_ylim(0,1)
    plt.ylabel('Syntenic BUSCO Connections (%)')
    plt.xlabel('Phylogenetic Distance (Substitutions/Site)')
    plt.yticks(np.arange(0,1.2,.2), labels = np.round(np.arange(0,120,20),2), rotation = 0)
    plt.legend()
    plt.savefig('{0}/SynDecay.pdf'.format(args.output), format = 'pdf', bbox_inches='tight')



def main():
    parser = argparse.ArgumentParser(description="phyca")
    parser.add_argument('-v', '--version', action='version', version='%(prog)s {version}'.format(version='__version__'))
    
    #group = parser.add_mutually_exclusive_group(required=True)
    parser.add_argument("-a", "--assembly", type=str, help="Assembly in FASTA format", metavar="assembly")
    parser.add_argument("-c", "--compdir", type=str, help="Compleasm output directory", metavar="compleasm_directory")

    parser.add_argument("-l",'--lineage', type=str, help="BUSCO lineage", metavar='lineage', required=True)
    
    parser.add_argument("-o",'--output', type=str, help="Output prefix", metavar='output', default='upout')
    parser.add_argument("-t",'--threads', type=int, help="Compleasm threads", metavar='threads', default=4)
    parser.add_argument("-r",'--reference', type=str, help="Reference assembly", metavar='reference')
    parser.add_argument("-m",'--rcompdir', type=str, help="Reference compleasm output directory", metavar='rcompleasm_directory')
    parser.add_argument("-n",'--nullify', action='store_true', help="Remove all BUSCO genes in assembly")
    parser.add_argument("-s",'--syndis', action='store_true', help="Compute syntenic distance from reference")
    parser.add_argument("-i","--ignore_orientation", action='store_true', 
                        help="Ignores orientation and only considers gene order for syntenic distances.")
    parser.add_argument("-d","--include_duplications", action='store_true', 
                        help="Duplicated gene pairs are considered distinct for syntenic distances.")
    parser.add_argument("-w","--include_singleton_contigs", action='store_true', 
                        help="Includes contigs with single genes for syntenic distances.")
    

    parser.set_defaults(func=phyca)

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit()

    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
# phyca: **phy**logeny and **c**ollinearity **a**ware assembly evaluation toolkit.

phyca is built around [Compleasm](https://github.com/huangnengCSU/compleasm/) utilizing the [NCBI Genome](https://www.ncbi.nlm.nih.gov/datasets/genome/) database. For a query assembly, phyca improves the precision of BUSCO/Compleasm annotations by up to 7%, makes syntenic comparisons to public reference genomes and rapidly places the assembly on a broad, precomputed phylogeny.

# Installation
```
pip install phyca
```

phyca is distributed through PyPI and github. A working installation of Compleasm (including SEPP and pplacer) is necessary to avail all functionality. I recommend creating a conda environment to install Compleasm first and installing phyca in that environment, e.g.,

```
# create environment
conda create -n phyca python=3.9.19
# install compleasm
conda install bioconda::compleasm=0.2.6
# install phyca
pip install phyca
```

Note that as of 02/03/2025, there is a known issue with pplacer and SEPP on Debian-based systems. A working solution is provided [here](https://github.com/smirarab/sepp/issues/140).

phyca has the following nonexhaustive dependency structure.
```
Python (tested with 3.9.19)
↓
│───numpy (tested with 2.0.2)
│───pandas (tested with 2.2.3)
|───matplotlib (tested with 3.9.4)
│───seaborn (tested with 0.13.2)
│───SciPy (tested with 1.13.1)
|───BioNick (tested with 0.0.7)
└───Compleasm (tested with 0.2.6)
        |─── hmmer (tested with 3.1b2)
        |─── miniprot (tested with 0.13-r248)
        |      └─── libgcc (tested with 14.2.0 under conda)
        └─── SEPP (tested with 4.4.0)
               └─── pplacer and guppy (v1.1.alpha19-0-g807f6f3) 
```

# Usage

phyca supports 10 BUSCO lineages: viridiplantae, liliopsida, eudicots, chlorophyta, fungi, ascomycota, basidiomycota, metazoa, arthropoda and vertebrata.

A simple run on a query assembly, would be:
```
phyca -a <assembly file> -l <lineage>
```
The Compleasm output folder can also be used as input if compleasm output was previously generated:
```
phyca -c <compleasm_direcoty> -l <lineage>
```

The above run will output BUSCO, CUSCO (Curated USCOs with higher precision) and MUSCO (remaining USCOs) statistics and graphs. It will compare the query to chromosome level genome assemblies from NCBI genome and output a table with a measure of synteny against each genome. It will output a Neighbor-Joining tree based on BUSCO synteny. Finally, it will place the assembly on a large precomputed phylogeny for the lineage and graph the observed decay in BUSCO synteny against inferred phylogenetic distance.


phyca can also be used to compute the syntenic distance between two assemblies with the -s flag. 
```
phyca -l <lineage> -s -a <assembly1> -r <assembly2>
```
The same comparison can be done by pointing to the compleasm output directoreis, if already available.
```
phyca -l <lineage> -s -c <assembly1_compdir> -m <assembly2_compdir>
```

# UniPhyDB
The bulk data used by phyca is hosted by [AGI](https://www.genome.arizona.edu/)'s [AVA cluster](https://www.genome.arizona.edu/services/instrumentation.html). Precomputed trees and more information is available on https://UniPhydb.github.io/ . 



# Example Output

USCO graph:

<img src="https://ava.genome.arizona.edu/UniPhy/web/USCO_bars.png" width="400">

Synteny decay plot:

<img src="https://ava.genome.arizona.edu/UniPhy/web/syndecay.png" width="400">


Placement tree snippet: 

<img src="https://ava.genome.arizona.edu/UniPhy/web/placement_snippet.png" width=400>
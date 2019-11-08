## Overview

Companion to "**[A genome-wide almanac of co-essential modules assigns function to uncharacterized genes](https://doi.org/10.1101/827071)**".

Contains code to generate co-essential gene pairs, co-essential modules, and modules with cancer type-specific dependencies. Coming soon: code to generate the two-dimensional layout (Fig. 1C). 

For the web tool associated with the paper, see **[coessentiality.net](http://coessentiality.net/)**.

## Code files

1. **[gene_pairs.py](https://github.com/kundajelab/coessentiality/blob/master/gene_pairs.py)**: generates co-essential gene pairs.
2. **[modules.py](https://github.com/kundajelab/coessentiality/blob/master/modules.py)**: generates co-essential modules using the gene pairs from #1.
3. **[cancer_type_dependencies.py](https://github.com/kundajelab/coessentiality/blob/master/cancer_type_dependencies.py)**: enumerates modules with cancer type-specific dependencies using the gene pairs and modules from #1 and #2.
4. **[load_screens.py](https://github.com/kundajelab/coessentiality/blob/master/load_screens.py)**: loads and bias-corrects CRISPR screens. Used by #1 and #3.

## Required external files

1. **[gene_effect.csv](https://ndownloader.figshare.com/files/12704099)**: CRISPR screens from the "DepMap Public 18Q3" release at **[https://depmap.org/portal/download/all/](https://depmap.org/portal/download/all/)**. Required for **gene_pairs.py** and **cancer_type_dependencies.py**.
2. **[sample_info.csv](https://ndownloader.figshare.com/files/12704612)**: metadata for the cell lines in gene_effect.csv. Required for **gene_pairs.py** and **cancer_type_dependencies.py**.
3. **[cluster_one-1.0.jar](https://www.paccanarolab.org/static_content/clusterone/cluster_one-1.0.jar)**: Java executable for ClusterONE. Required for **modules.py**.

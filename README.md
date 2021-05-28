## Overview

Companion to "**[A genome-wide almanac of co-essential modules assigns function to uncharacterized genes](https://doi.org/10.1101/827071)**".

Contains code to generate co-essential gene pairs, co-essential modules, and modules with cancer type-specific dependencies. Coming soon: code to generate the two-dimensional layout (Fig. 1C). 

For the web tool associated with the paper, see **[coessentiality.net](http://coessentiality.net/)**. If you would just like the final 17634 x 17634 matrix of p-values, you can download it [here](http://mitra.stanford.edu/bassik/coessentiality/GLS_p.npy). A corresponding 17634 x 17634 matrix with the sign of each correlation (1 = positive, -1 = negative) is downloadable [here](http://mitra.stanford.edu/bassik/coessentiality/GLS_sign.npy). The list of the 17634 genes that form the rows and columns of these matrices can be found [here](http://mitra.stanford.edu/bassik/coessentiality/genes.txt). These matrices are in NumPy's npy format and can be loaded with [`np.load`](https://numpy.org/doc/stable/reference/generated/numpy.load.html) in Python, the [RcppCNPy](https://dirk.eddelbuettel.com/code/rcpp.cnpy.html) library in R, and the [cnpy](https://github.com/rogersce/cnpy) library in C/C++. 

2D coordinates of the final gene layout can be downloaded [here](https://github.com/kundajelab/coessentiality/blob/master/vizdf.tsv). 

## Code files

1. **[gene_pairs.py](https://github.com/kundajelab/coessentiality/blob/master/gene_pairs.py)**: generates co-essential gene pairs.
2. **[modules.py](https://github.com/kundajelab/coessentiality/blob/master/modules.py)**: generates co-essential modules using the gene pairs from #1.
3. **[cancer_type_dependencies.py](https://github.com/kundajelab/coessentiality/blob/master/cancer_type_dependencies.py)**: enumerates modules with cancer type-specific dependencies using the gene pairs and modules from #1 and #2.
4. **[load_screens.py](https://github.com/kundajelab/coessentiality/blob/master/load_screens.py)**: loads and bias-corrects CRISPR screens. Used by #1 and #3.
5. **[generate_layout.tsv](https://github.com/kundajelab/coessentiality/blob/master/generate_layout.py)**: generates gene network for visualization, from the GLS p-value matrix and the gene modules from #2.

## Required external files

1. **[gene_effect.csv](https://ndownloader.figshare.com/files/12704099)**: CRISPR screens from the "DepMap Public 18Q3" release at **[https://depmap.org/portal/download/all/](https://depmap.org/portal/download/all/)**. Required for **gene_pairs.py** and **cancer_type_dependencies.py**.
2. **[sample_info.csv](https://ndownloader.figshare.com/files/12704612)**: metadata for the cell lines in gene_effect.csv. Required for **gene_pairs.py** and **cancer_type_dependencies.py**.
3. **[cluster_one-1.0.jar](https://www.paccanarolab.org/static_content/clusterone/cluster_one-1.0.jar)**: Java executable for ClusterONE. Required for **modules.py**.

import numpy as np, pandas as pd, subprocess
from statsmodels.stats.multitest import fdrcorrection

# Load GLS results

genes = pd.read_csv('genes.txt', header=None, squeeze=True)
GLS_p = pd.DataFrame(np.load('GLS_p.npy'), columns=genes, index=genes)

# Compute and save weights for ClusterONE

stacked_p = GLS_p.stack()
stacked_p = stacked_p[stacked_p.index.get_level_values(0) <
                      stacked_p.index.get_level_values(1)]
fdr = pd.Series(fdrcorrection(stacked_p)[1], index=stacked_p.index)
weights = 1 - fdr
weight_file = 'ClusterONE_weights.txt'
weights.to_csv(weight_file, sep=' ', header=None)

# Run ClusterONE at each d, and save results

for d in 0.2, 0.5, 0.9:
    subprocess.check_call(f'java -jar cluster_one-1.0.jar {weight_file} '
                          f'-d {d} -F csv > modules_d_{d}.csv', shell=True)
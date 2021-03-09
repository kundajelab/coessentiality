import numpy as np, pandas as pd, statsmodels.api as sm
from collections import defaultdict
from .load_screens import load_screens
from scipy.stats import cauchy
from statsmodels.stats.multitest import fdrcorrection

# Load screens and modules

screens = load_screens()
genes = screens.index
modules = np.array([
    module.split() for d in (0.2, 0.5, 0.9) for module in pd.read_csv(
        f'modules_d_{d}.csv', usecols=['Members'], squeeze=True)])

# Get cancer types; filter out cancer types that only appear once

cell_line_cancer_types = screens.columns.str.split('_').str[1:].str.join(
    ' ').str.capitalize().to_series()
cancer_type_frequencies = cell_line_cancer_types.value_counts()
singletons = cancer_type_frequencies.index[cancer_type_frequencies == 1]
non_singletons_mask = ~cell_line_cancer_types.isin(singletons).values
screens = screens.iloc[:, non_singletons_mask]
cell_line_cancer_types = cell_line_cancer_types[non_singletons_mask]

# One-hot encode

one_hot_cancer_types = pd.get_dummies(cell_line_cancer_types)
cancer_types = one_hot_cancer_types.columns

# Run multivariate GLS on each gene

cholsigmainv = np.linalg.cholesky(np.linalg.pinv(screens.cov())).T
inputs = cholsigmainv.dot(sm.add_constant(one_hot_cancer_types))
gene_ps = pd.DataFrame(index=genes, columns=cancer_types, dtype=float)
for gene in genes:
    output = cholsigmainv.dot(screens.loc[gene])
    model = sm.OLS(output, inputs).fit()
    gene_ps.loc[gene] = model.pvalues[1:]

# Combine GLS p-values with ACAT

ACAT = lambda ps: cauchy.sf(np.tan((0.5 - ps) * np.pi).mean())
results = defaultdict(list)
for gene_set in modules:
    gene_set_ps = gene_ps[genes.isin(gene_set)]
    for cancer_type in cancer_types:
        p_meta = ACAT(gene_set_ps[cancer_type].astype(np.float128))
        results['Module genes'].append(', '.join(gene_set))
        results['Cancer type'].append(cancer_type)
        results['p'].append(p_meta)
results = pd.DataFrame(results)

# FDR-correct by cancer type

results['FDR'] = results.p.groupby(results['Cancer type']).apply(
    lambda p: pd.Series(fdrcorrection(p)[1], index=p.index))

# Save significant module-cancer type dependencies

results[results.FDR < 0.1].sort_values('p').to_csv(
    'cancer_type_dependencies.tsv', sep='\t', index=False)
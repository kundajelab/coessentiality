import pandas as pd
from sklearn.decomposition import PCA

def load_screens():
    # Load screens
    screens = pd.read_csv('gene_effect.csv', index_col=0).T
    screens.index = screens.index.str.split(' ').str.get(0)
    # Map Broad ID to CCLE name
    cell_lines = pd.read_csv('sample_info.csv', index_col='Broad_ID',
                             usecols=['Broad_ID', 'CCLE_name'], squeeze=True)
    screens.columns = screens.columns.reindex(cell_lines)[0]
    # Bias-correct using "molecular function:olfactory receptor activity" genes
    olfactory_genes = pd.read_csv(
        'olfactory_genes.txt', header=None, squeeze=True)
    olfactory_data = screens.reindex(olfactory_genes).dropna()
    transformation = PCA(n_components=4)
    transformation.fit(olfactory_data)
    top_PC_effects = transformation.inverse_transform(
        transformation.transform(screens))
    screens -= top_PC_effects
    screens = screens.iloc[:, :-4]
    return screens
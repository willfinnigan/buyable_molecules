import pandas as pd
import numpy as np

def load_dfs(paths):
    print(f'loading {len(paths)} dataframes')
    list_dfs = []
    for path in paths:
        df = pd.read_csv(path, index_col=0)
        list_dfs.append(df)
    return list_dfs

def merge_dfs(list_dfs):
    print('merging dataframes')
    all_df = list_dfs.pop(0)
    for df in list_dfs:
        all_df = all_df.merge(df, on='SMILES', how='outer')
    return all_df

def fix_emol_id(id):
    if id is None:
        return None
    elif np.isnan(id) == True:
        return None
    else:
        return str(int(id))

if __name__ == "__main__":
    paths = ['mcule.csv',
             'molport.csv',
             'sigma.csv',
             'zincbb.csv',]

    dfs = load_dfs(paths)
    final_df = merge_dfs(dfs)
    print(f"{len(final_df.index)} rows")
    print(final_df.head())
    final_df.to_csv('final.csv')


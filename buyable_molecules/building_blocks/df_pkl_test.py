import pandas as pd
from pathlib import Path
from buyable_molecules.rdkit_pandas_functions import create_mol_col
import time
import dask.dataframe as dd
from dask.distributed import Client
import gc
from rdkit import Chem
import numpy as np

bb_folder = str(Path(__file__).parents[0]) + '/building_block_dfs'
final_building_block_df_path = f'{bb_folder}/final/final.csv'

def mol_col(row):
    try:
        return Chem.MolFromSmiles(row['SMILES'])
    except:
        print(f"error with smi {row['SMILES']}")
        return None

def test_func(smi):
    return smi

def test_loading_times(path):
    t0 = time.time()
    df = pd.read_csv(path)
    print(f"Time to load final csv = {round(time.time() - t0, 2)} seconds")
    print(df.head())


    t0 = time.time()
    df.to_pickle('final.pkl')
    print(f"Time to save final pkl = {round(time.time() - t0, 2)} seconds")

    t0 = time.time()
    df2 = pd.read_pickle('final.pkl')
    print(f"Time to load final pkl = {round(time.time() - t0, 2)} seconds")

def split_dfs(df, split=200, folder=f"{bb_folder}/final/final_split"):
    split_dfs = np.array_split(df, split)

    for i, s_df in enumerate(split_dfs):
        s_df.to_csv(f'{folder}/split_{i}.csv')

    split_dfs = None
    s_df = None
    gc.collect()

    return

def load_and_combine(split, folder=f"{bb_folder}/final/final_split/hdf"):
    df = pd.DataFrame()
    for i in range(split):
        print(i)
        df = df.append(pd.read_hdf(f"{folder}/split_{i}.h5", 'df'), ignore_index=True)

    return df

def create_mol_df_hdfs(num_splits, folder=f"{bb_folder}/final/final_split"):
    for i in range(num_splits):
        df = pd.read_csv(f"{folder}/split_{i}.csv", index_col=0, dtype={'SMILES': str,
                                                                          'molport_id': str,
                                                                          'mcule_id': str,
                                                                          'sigma_id': str,
                                                                          'zinc_id': str})
        df['mol'] = df['SMILES'].apply(create_mol_col)
        df.dropna(subset=['mol'], inplace=True)
        df.to_hdf(f"{folder}/hdf/split_{i}.h5", key='df', mode='w')


if __name__ == "__main__":
    split = 200

    """
    df = pd.read_csv(final_building_block_df_path, dtype={'SMILES': str, 'molport_id': str, 'mcule_id': str, 'sigma_id': str, 'zinc_id': str})
    split_dfs(df, split=split)

    del df
    gc.collect()
    """

    #create_mol_df_hdfs(split)

    #df = pd.read_hdf(f"{bb_folder}/final/final_split/hdf/split_0.h5", 'df')

    df = load_and_combine(split)
    print(df.head())

    #df.to_hdf(f'{bb_folder}/final/final_w_mols.hdf', key='df')

    """
    # test_loading_times(final_building_block_df_path)

    #df = pd.read_csv(final_building_block_df_path)
    #client = Client(memory_limit='4gb', n_workers=1)
    df = pd.read_csv(final_building_block_df_path, index_col=0, dtype={'SMILES': str,
                                                                      'molport_id': str,
                                                                      'mcule_id': str,
                                                                      'sigma_id': str,
                                                                      'zinc_id': str})
    dask_df = dd.from_pandas(df, chunksize=1000)
    df = None
    gc.collect()

    dask_df['mol'] = dask_df['SMILES'].map_partitions(test_func, meta=('mol', 'str'))
    #dask_df['mol'] = dask_df['SMILES'].apply(create_mol_col, meta=('mol', 'object'))
    #dask_df['mol'] = dask_df['SMILES'].map_partitions(create_mol_col, meta=('mol', 'object'))
    #res = dask_df.applymap(mol_col, meta=('mol', 'object'))
    dask_df = dask_df.dropna(subset=['mol'])
    dask_df.to_parquet('final.parquet', schema="infer")
    """




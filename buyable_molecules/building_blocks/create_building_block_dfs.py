from buyable_molecules.download_and_process import download_and_process
from pathlib import Path
import pandas as pd

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

def nan_to_none(df):
    df = df.where(pd.notnull(df), None)
    return df

def create_merged_df(paths):
    list_dfs = load_dfs(paths)
    merged_df = merge_dfs(list_dfs)
    merged_df = nan_to_none(merged_df)
    return merged_df


if __name__ == "__main__":

    save_folder = str(Path(__file__).parents[0]) + '/building_block_dfs'

    mcule_url = "https://mcule.s3.amazonaws.com/database/mcule_purchasable_building_blocks_210714.smi.gz"
    z_sigma_url = 'http://files.docking.org/catalogs/source/sialbb.src.txt'
    z_molport_url = 'http://files.docking.org/catalogs/source/molportbb.src.txt'
    z_bb_instock_url = 'http://files.docking.org/bb/current/In-Stock.txt'
    # emolecules_url = 'https://downloads.emolecules.com/free/2021-07-01/version.smi.gz'
    # emols links: https://reaxys.emolecules.com/cgi-bin/more?vid=477088

    #download_and_process(mcule_url, ["SMILES", "mcule_id"], None, '\t', f'{save_folder}/mcule.csv', download=True, unzip=True)
    download_and_process(z_sigma_url, ["SMILES", "sigma_id"], None, ' ', f'{save_folder}/sigma.csv', download=True, unzip=False)

    '''
    download_and_process(z_molport_url, ["SMILES", "molport_id"], None, '\t', f'{save_folder}/molport.csv', download=True, unzip=False)
    download_and_process(z_bb_instock_url, ["SMILES", "zinc_id", 'data_type', 'supplier', 'code'],
                         None, '\t', f'{save_folder}/zinc.csv', download=True, unzip=False)

    merged_df = create_merged_df([f'{save_folder}/mcule.csv',
                                  f'{save_folder}/sigma.csv',
                                  f'{save_folder}/molport.csv',
                                  f'{save_folder}/zinc.csv'])

    merged_df.to_csv(f"{save_folder}/final")
    '''

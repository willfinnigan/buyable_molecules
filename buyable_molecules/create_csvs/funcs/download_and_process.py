import pandas as pd
import numpy as np
import urllib.request
import gzip
import shutil
import os
from rdkit import Chem
from rdkit.Chem.SaltRemover import SaltRemover
from rdkit.Chem import PandasTools
remover = SaltRemover()
from rdkit import RDLogger
lg = RDLogger.logger()
lg.setLevel(RDLogger.ERROR)

from buyable_molecules.create_csvs.funcs import dataframe_functions

def load_data(path, cols, sep, header=None):
    print(f'Load path: {path}')
    data = pd.read_csv(path, sep=sep, header=header, names=cols)
    print(data.head())

    return data

def load_sdf(path):
    print(f'Load path: {path}')
    data = PandasTools.LoadSDF(path, smilesName='SMILES', molColName='Molecule',
                               includeFingerprints=False, removeHs=False, strictParsing=True)
    print(data.head())
    return data

def process_lifechem_sdf(df):
    def get_rdkit_smiles(mol):
        try:
            return Chem.MolToSmiles(mol)
        except:
            return None

    new_df = pd.DataFrame()
    df['Molecule'] = df['Molecule'].apply(dataframe_functions.remove_salt)
    new_df['SMILES'] = df['Molecule'].apply(get_rdkit_smiles)

    prices = []
    for index, row in df.iterrows():
        price = f"{row.get('Price_1', '')}$/{row.get('Pack_1', '')}"
        prices.append(price)
    new_df['lifechem_price'] = prices

    new_df['lifechem_id'] = df['IDNUMBER']

    new_df["SMILES"].replace('', np.nan, inplace=True)
    new_df = new_df.dropna(subset=['SMILES'])
    new_df = new_df.drop_duplicates()
    new_df = dataframe_functions.remove_duplicate_smiles(new_df)
    return new_df

def process_df_smiles(df_smiles):
    print('Data loaded - processing smiles..')
    print(f"Initial number of rows = {len(df_smiles.index)}")
    df_smiles["SMILES"] = df_smiles["SMILES"].apply(dataframe_functions.convert_to_rdkit)
    df_smiles["SMILES"].replace('', np.nan, inplace=True)
    df_smiles = df_smiles.dropna(subset=['SMILES'])
    df_smiles = df_smiles.drop_duplicates()
    df_smiles = dataframe_functions.remove_duplicate_smiles(df_smiles)
    print(f"Number of rows after processing = {len(df_smiles.index)}")

    return df_smiles



def download_file(url, filename):
    print(f"Downloading file from {url}")
    urllib.request.urlretrieve(url, filename)
    return filename

def download_and_unzip(url, filename='zip_file.gz', unzipped='unzipped.txt'):
    download_file(url, filename)
    print(f"Unzipping..")
    with gzip.open(filename, 'rb') as f_in:
        with open(unzipped, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
    return unzipped

def remove_file(path):
    if os.path.exists(path):
        os.remove(path)

def select_columns(df, columns):
    new_df = df[columns]
    return new_df

def download_and_process(url_or_path, columns, header, sep, save_to, download=True, unzip=True):
    if download == True:
        if unzip == True:
            filename = download_and_unzip(url_or_path, filename='zip_file.gz', unzipped='file_to_process.txt')
            remove_file('zip_file.gz')
        else:
            filename = download_file(url_or_path, 'file_to_process.txt')
    else:
        filename = url_or_path

    df = load_data(filename, columns, sep, header=header)

    if download == True:
        remove_file(filename)

    df = process_df_smiles(df)
    df.to_csv(save_to)

    return save_to

def download_and_process_lifechem(url, save_to):
    file = download_file(LIFECHEM, 'lifechem.sdf')
    df = load_sdf(file)
    remove_file(file)
    df = process_lifechem_sdf(df)
    df.to_csv(save_to)
    return save_to

if __name__ == "__main__":
    from buyable_molecules.urls import LIFECHEM
    file = download_file(LIFECHEM, 'lifechem.sdf')
    df = load_sdf(file)
    df = process_lifechem_sdf(df)



    print(df.head())

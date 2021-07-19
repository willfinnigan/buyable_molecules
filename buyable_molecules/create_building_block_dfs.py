import pandas as pd
import numpy as np
import urllib.request
import gzip
import shutil
import os
from rdkit import Chem
from rdkit.Chem.SaltRemover import SaltRemover
remover = SaltRemover()
from rdkit import RDLogger
lg = RDLogger.logger()
lg.setLevel(RDLogger.ERROR)

def convert_to_rdkit(smi):
    try:
        mol = Chem.MolFromSmiles(smi)
        mol = remove_salt(mol)
        new_smi = Chem.MolToSmiles(mol)
        return new_smi
    except:
        print(f"{smi} not accepted by rdkit")
        return None

def remove_salt(mol):
    try:
        mol = remover.StripMol(mol, dontRemoveEverything=True)
    except:
        pass
    return mol

def remove_duplicate_smiles(df):
    df = df.drop_duplicates(subset=['SMILES'])
    return df

def load_data(path, cols, sep, smi_col, header=None):
    print(f'Load path: {path}')
    data = pd.read_csv(path, sep=sep, header=header, names=cols)
    print(data.head())

    return data

def process_df_smiles(df_smiles):
    print('Data loaded - processing smiles..')
    print(f"Initial number of rows = {len(df_smiles.index)}")
    df_smiles["SMILES"] = df_smiles["SMILES"].apply(convert_to_rdkit)
    df_smiles["SMILES"].replace('', np.nan, inplace=True)
    df_smiles = df_smiles.dropna()
    df_smiles = df_smiles.drop_duplicates()
    df_smiles = remove_duplicate_smiles(df_smiles)
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

def get_mcule_file_and_process(url):
    download_and_unzip(url, filename='zip_file.gz', unzipped='mcule.smi')
    remove_file('zip_file.gz')
    df = load_data('mcule.smi', ["SMILES", "mcule_id"], '\t', "SMILES", header=None)
    df = process_df_smiles(df)
    df.to_csv('mcule.csv')
    remove_file('mcule.smi')

def get_zinc_sigma_file_and_process(url):
    filename = 'sialbb.src.txt'
    download_file(url, filename)
    df = load_data(filename, ["SMILES", "sigma_id"], ' ', "SMILES", header=None)
    df = process_df_smiles(df)
    df.to_csv('sigma.csv')
    remove_file(filename)

def get_zinc_molport_file_and_process(url):
    filename = 'molportbb.src.txt'
    download_file(url, filename)
    df = load_data(filename, ["SMILES", "molport_id"], '\t', "SMILES", header=None)
    df = process_df_smiles(df)
    df.to_csv('molport.csv')
    remove_file(filename)

def get_zinc_bb_file_and_process(url):
    filename = 'In-Stock.txt'
    download_file(url, filename)
    df = load_data(filename, ["SMILES", "zinc_id", 'data_type', 'supplier', 'code'], '\t', "SMILES", header=None)
    df = select_columns(df, ["SMILES", "zinc_id"])
    df = process_df_smiles(df)
    df.to_csv('zincbb.csv')
    remove_file(filename)

def get_emolecules_file_and_process(url):
    filename = 'emols.smi'
    download_and_unzip(url, filename='zip_file.gz', unzipped=filename)
    remove_file('zip_file.gz')
    df = load_data(filename, ['SMILES', 'version_id', 'emols_id'], ' ', "SMILES", header=0)
    df = select_columns(df, ["SMILES", "emols_id"])
    df = process_df_smiles(df)
    df.to_csv('emols.csv')
    remove_file(filename)


if __name__ == "__main__":
    mcule_url = "https://mcule.s3.amazonaws.com/database/mcule_purchasable_building_blocks_210714.smi.gz"
    z_sigma_url = 'http://files.docking.org/catalogs/source/sialbb.src.txt'
    z_molport_url = 'http://files.docking.org/catalogs/source/molportbb.src.txt'
    z_bb_instock_url = 'http://files.docking.org/bb/current/In-Stock.txt'
    emolecules_url = 'https://downloads.emolecules.com/free/2021-07-01/version.smi.gz'

    #get_mcule_file_and_process(mcule_url)
    #get_zinc_sigma_file_and_process(z_sigma_url)
    #get_zinc_molport_file_and_process(z_molport_url)
    #get_zinc_bb_file_and_process(z_bb_instock_url)
    get_emolecules_file_and_process(emolecules_url)

    # emols links: https://reaxys.emolecules.com/cgi-bin/more?vid=477088
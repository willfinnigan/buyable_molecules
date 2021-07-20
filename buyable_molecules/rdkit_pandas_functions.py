import pandas as pd
import numpy as np
import urllib.request
import gzip
import shutil
import os
from rdkit import Chem
from bson import Binary
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

def create_mol_col(smi):
    try:
        return Chem.MolFromSmiles(smi)
    except:
        print(f'error with smi {smi}')
        return None

def mol_binary(mol):
    try:
        return Binary(mol.ToBinary())
    except:
        return None

def get_inchi_key(mol):
    try:
        return Chem.MolToInchiKey(mol)
    except:
        return None

def load_mol_columns(df):
    df['mol'] = df['SMILES'].apply(create_mol_col)
    df = df.dropna(subset=['mol'])
    df['rdmol'] = df['mol'].apply(mol_binary)
    df = df.dropna(subset=['rdmol'])
    df['inchi_key'] = df['mol'].apply(get_inchi_key)
    df = df.dropna(subset=['inchi_key'])

    df = df.where(pd.notnull(df), None)

    return df
import mongoengine as db
from buyable_molecules.modals import BuildingBlock
import pandas as pd
import numpy as np
import gc

import urllib.request
from pathlib import Path
import gzip
import shutil
import os
import time
import pymongo
from pymongo import UpdateOne
from rdkit import Chem
from bson import Binary


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

    return df

def df_to_collection(df, collection):

    operations = []

    for index, row in df.iterrows():
        smi = row['SMILES']
        vendors = []
        if row['mcule_id'] != None:
            mcule_id = str(row['mcule_id'])
            vendors.append('mcule_bb')
        else:
            mcule_id = None

        if row['sigma_id'] != None:
            sigma_id = str(row['sigma_id'])
            vendors.append('sigma_bb')
        else:
            sigma_id = None

        if row['zinc_id'] != None:
            zinc_id = str(row['zinc_id'])
            vendors.append('zinc_bb_in_stock')
        else:
            zinc_id= None

        if row['molport_id'] != None:
            molport_id = str(row['molport_id'])
            vendors.append('molport_bb')
        else:
            molport_id = None

        try:
            inchi_key = row['inchi_key']
            mol_binary = row['inchi_key']

            op = UpdateOne({"index": inchi_key},
                           {"$set": {'smiles': smi,
                                     'index': inchi_key,
                                     'fingerprints': {},
                                     'rdmol': mol_binary,
                                     'vendors': vendors,
                                     'mcule_id': mcule_id,
                                     'sigma_id': sigma_id,
                                     'molport_id': molport_id,
                                     'zinc_id': zinc_id}}
                           )
            operations.append(op)
        except:
            print(f"Couldn't load smi {smi}")

    collection.bulk_write(operations)
    operations = None
    return


def data_to_mongo(df, chunk_size, collection):
    total_length = len(df.index)
    num_splits = int(total_length / chunk_size)

    split_dfs = np.array_split(df, num_splits)
    print(f"Split df into {num_splits} chunks")
    del df
    gc.collect()

    for i, s_df in enumerate(split_dfs):
        print(f"Processing data chunk {i} of {num_splits}")
        t0 = time.time()
        s_df = load_mol_columns(s_df)
        df_to_collection(s_df, collection)
        del s_df
        gc.collect()
        print(f"Time for chunk = {round(time.time() - t0, 1)} seconds")

if __name__ == '__main__':
    db.connect('molecules')
    BuildingBlock.drop_collection()

    client = pymongo.MongoClient('0.0.0.0', 27017)
    database = client['molecules']
    collection = database['building_block']

    df = pd.read_csv('final.csv', index_col=0)
    data_to_mongo(df, 10000, collection)

    test_mol = BuildingBlock.objects()[0].to_json()
    print(test_mol)


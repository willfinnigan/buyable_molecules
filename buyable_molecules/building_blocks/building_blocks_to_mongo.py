import mongoengine as db
import numpy as np
import gc
import time
from pathlib import Path
from buyable_molecules.modals import BuildingBlock
from buyable_molecules.rdkit_pandas_functions import load_mol_columns
import pandas as pd
import pymongo

final_building_block_df_path = str(Path(__file__).parents[0]) + '/building_block_dfs/final/final.csv'

def df_to_collection_via_mongoengine(df):

    to_insert = []

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

        inchi_key = row['inchi_key']
        mol_binary = row['rdmol']

        bb = BuildingBlock(smiles=smi,
                           index=inchi_key,
                           fingerprints={},
                           rdmol=mol_binary,
                           vendors=vendors,
                           mcule_id=mcule_id,
                           sigma_id=sigma_id,
                           molport_id=molport_id,
                           zinc_id=zinc_id)
        to_insert.append(bb)

        """op = UpdateOne({"index": inchi_key},
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
        operations.append(op)"""

    #collection.bulk_write(operations)
    #operations = None
    BuildingBlock.objects.insert(to_insert)
    to_insert = None
    bb = None
    gc.collect()
    return

def df_to_collection_via_pymongo(df, collection):

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

        inchi_key = row['inchi_key']
        mol_binary = row['rdmol']

        op = pymongo.UpdateOne({"index": inchi_key},
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

    collection.bulk_write(operations)
    operations = None
    op = None
    gc.collect()
    return

def data_to_mongo(df):

    t0 = time.time()
    df = load_mol_columns(df)
    df_to_collection_via_mongoengine(df)
    print(f"Time for df = {round(time.time() - t0, 1)} seconds")

def split_dfs(df, split=10):
    split_dfs = np.array_split(df, split)

    for i, s_df in enumerate(split_dfs):
        s_df.to_csv(f'final_split/split_{i}.csv')

    return

def data_to_mongo2(df, collection):
    df = load_mol_columns(df)
    df_to_collection_via_pymongo(df, collection)


if __name__ == '__main__':
    db.connect('molecules')
    '''
    BuildingBlock.drop_collection()

    df = pd.read_csv('final.csv', index_col=0)
    split=10
    split_dfs(df)
    del df
    gc.collect()

    for i in range(split):
        df = pd.read_csv(f"final_split/split_{i}.csv")
        data_to_mongo(df)
        del df
        gc.collect()
    '''
    t0 = time.time()
    bb = BuildingBlock.objects(db.Q(smiles='CO') & db.Q(vendors='mcule_bb'))
    #bb = BuildingBlock.objects(mcule_id='MCULE-1370061678')

    print(len(bb))

    print(bb[0].mcule_id)
    t1 = time.time()
    print(round(t1 - t0, 3))




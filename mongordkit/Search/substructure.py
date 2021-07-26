import pymongo
import rdkit
import math
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.rdmolops import PatternFingerprint
from mongordkit import Database
import time


def SubSearchNaive(pattern, mol_collection, chirality=False):
    """
    Search MOL_COLLECTION for molecules with PATTERN
    as a substructure.
    :param pattern: An rdmol object that represents a desired substructure pattern.
    :param mol_collection: A MongoDB collection.
    :param chilarity: Whether or not to include stereochemistry in search. Defaults to False.
    :return: A list of SMILES that have the desired substructure pattern.
    """
    results = []
    for molDoc in mol_collection.find():
        rdmol = Chem.Mol(molDoc['rdmol'])
        if rdmol.HasSubstructMatch(pattern, useChirality=chirality):
            results.append(molDoc['index'])
    return results


def AddPatternFingerprints(mol_collection):
    """
    Adds pattern fingerprints to each document in MOL_COLLECTION.
    :param mol_collection: A MongoDB collection.
    """

    print("Adding pattern fingerprints for substructure searching to molecules...")
    t0 = time.time()
    total_records = mol_collection.count_documents({})
    for i, moldoc in enumerate(mol_collection.find()):
        mol = Chem.Mol(moldoc['rdmol'])
        bit_vector = list(PatternFingerprint(mol).GetOnBits())
        count = len(bit_vector)
        mol_collection.update_one({'_id': moldoc['_id']}, {'$set': {'fingerprints.pattern_fp': {'bits': bit_vector, 'count': count}}})

        if i % 10000 == 0 and i != 0:
            elapsed = round(time.time()-t0, 1)
            estimated_total_time = (elapsed / i) * total_records
            time_remaining = estimated_total_time - elapsed
            print(f"{i} added in {elapsed} seconds - estimated time remaining = {round(time_remaining / 60, 2)} mins")

    return


def SubSearch(pattern, mol_collection, chirality=False):
    """
    Search MOL_COLLECTION for molecules with PATTERN
    as a substructure.
    :param pattern: An rdmol object that represents a desired substructure pattern.
    :param mol_collection: A MongoDB collection.
    :param chirality: Whether or not to include stereochemistry in search. Defaults to False.
    :return: A list of SMILES that have the desired substructure pattern.
    """
    results = []
    query_fp = list(PatternFingerprint(pattern).GetOnBits())
    qfp_len = len(query_fp)
    for molDoc in mol_collection.find({'fingerprints.pattern_fp.count': {'$gte': qfp_len},
                                     'fingerprints.pattern_fp.bits': {'$all': query_fp}
                                     }):
        rdmol = Chem.Mol(molDoc['rdmol'])
        if rdmol.HasSubstructMatch(pattern, useChirality=chirality):
            results.append(molDoc['index'])
    return results


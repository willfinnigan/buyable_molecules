from mongordkit.Search import similarity, substructure, utils
from mongordkit import Search
from mongordkit.Database import create, write
from rdkit import Chem
import rdkit
import pymongo
import time

if __name__ == '__main__':
    client = pymongo.MongoClient()
    database = client['molecules']
    bb_collection = database['building_block']
    mfp_collection = database['mfp_counts']

    #similarity.AddMorganFingerprints_new(bb_collection, mfp_collection)
    #similarity.create_counts(bb_collection, mfp_collection)
    #similarity.create_indexes(bb_collection)

    doc = bb_collection.find_one()['fingerprints']['morgan_fp']
    print(doc)

    q_mol = Chem.MolFromSmiles('Cc1ccccc1')

    # Perform a similarity search on TestDB for q_mol with a Tanimoto threshold of 0.4.
    t0 = time.time()
    results1 = similarity.SimSearch(q_mol, bb_collection, mfp_collection, 0.4)
    t1 = time.time()
    print(f"Time for SimSearch = {round(t1-t0,2)} seconds")

    # Do the same thing, but use the MongoDB Aggregation Pipeline.
    t0 = time.time()
    results2 = similarity.SimSearchAggregate(q_mol, bb_collection, mfp_collection, 0.4)
    t1 = time.time()
    print(f"Time for SimSearchAggregate = {round(t1 - t0, 2)} seconds")

    print('similaritySearch: {}'.format(results1))
    print('\n')
    print('similaritySearchAggregate: {}'.format(results2))

    substructure.AddPatternFingerprints(bb_collection)


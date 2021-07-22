from mongordkit.Search import similarity, substructure, utils
from mongordkit import Search
from mongordkit.Database import create, write
from rdkit import Chem
import rdkit
import pymongo

if __name__ == '__main__':
    client = pymongo.MongoClient()
    database = client['molecules']
    bb_collection = database['building_block']
    mfp_collection = database['mfp_counts']

    similarity.AddMorganFingerprints(bb_collection, mfp_collection)

    doc = bb_collection.find_one()['fingerprints']['morgan_fp']
    print(doc)

    q_mol = Chem.MolFromSmiles('Cc1ccccc1')

    # Perform a similarity search on TestDB for q_mol with a Tanimoto threshold of 0.4.
    results1 = similarity.SimSearch(q_mol, bb_collection, mfp_collection, 0.4)

    # Do the same thing, but use the MongoDB Aggregation Pipeline.
    results2 = similarity.SimSearchAggregate(q_mol, bb_collection, mfp_collection, 0.4)

    print('similaritySearch: {}'.format(results1))
    print('\n')
    print('similaritySearchAggregate: {}'.format(results2))


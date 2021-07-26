import pymongo
import mongoengine as db
import time
from buyable_molecules.modals import BuildingBlock, Emolecules

if __name__ == '__main__':
    db.connect('molecules')

    t0 = time.time()
    bb = BuildingBlock.objects(db.Q(smiles="CCCC=O"))
    # bb = BuildingBlock.objects(mcule_id='MCULE-1370061678')

    print(bb[0].to_json())
    t1 = time.time()
    print(round(t1 - t0, 3))


    client = pymongo.MongoClient()
    database = client['molecules']
    collection = database['building_block']

    collection.index_information()

    t0 = time.time()
    query = {"smiles": "O=CC(=O)Cc1ccccc1"}
    #query = {"sigma_id": "1084259"}
    docs = collection.find_one(query)

    print(docs)
    t1 = time.time()
    print(round(t1 - t0, 10))

    print(database.list_collection_names())

    collection2 = database['emolecules']


    info = collection.index_information()
    print(info)

    t0 = time.time()

    # query = {"sigma_id": "1084259"}
    docs = collection2.count_documents({})

    print(docs)
    t1 = time.time()
    print(round(t1 - t0, 10))

    t0 = time.time()
    query = {"smiles": "CCCCCC(=O)OC1CC[C@@]2(C)C(=CCC3C2CC[C@@]2(C)C3CC[C@@H]2[C@H](C)CCCC(C)C)C1"}
    docs = collection2.find_one(query)
    t1 = time.time()
    print(docs)
    print(round(t1 - t0, 10))

    t0 = time.time()
    docs = collection.find(batch_size=100)
    t1 = time.time()
    print(round(t1 - t0, 10))

    t0 = time.time()
    total_records = collection.count_documents({})
    t1 = time.time()
    print(round(t1 - t0, 10))
    print(total_records)
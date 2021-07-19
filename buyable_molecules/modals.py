import mongoengine as db

class BuildingBlock(db.DynamicDocument):
    smiles = db.StringField()
    rdmol = db.BinaryField()
    index = db.StringField()
    fingerprints = db.DictField()
    mcule_id = db.StringField()
    vendors = db.ListField(db.StringField())





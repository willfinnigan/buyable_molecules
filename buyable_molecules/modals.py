import mongoengine as db

class BuildingBlock(db.DynamicDocument):
    smiles = db.StringField()
    rdmol = db.BinaryField()
    index = db.StringField()
    fingerprints = db.DictField()
    mcule_id = db.StringField()
    sigma_id = db.StringField()
    molport_id = db.StringField()
    zinc_id = db.StringField()

    vendors = db.ListField(db.StringField())

    meta = {'indexes': ['smiles']}





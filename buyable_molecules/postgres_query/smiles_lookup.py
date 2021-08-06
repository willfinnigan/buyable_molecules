from buyable_molecules.create_postgres_db.funcs.db_setup import connect_to_db
from buyable_molecules.urls import VENDORS


def smiles_lookup(smiles, cursor, table='building_blocks', vendors=[], print_cmd=False):
    cmd = f"select * from {table} where smiles='{smiles}'"


    num_vendors = len(vendors)
    if num_vendors != 0:
        cmd += " and ("

    for i, vendor in enumerate(vendors):
        cmd += f"{vendor}_id is not NULL"

        if i != num_vendors-1:
            cmd += " or "
        else:
            cmd += ")"

    if print_cmd == True:
        print(cmd)

    cursor.execute(cmd)
    result = cursor.fetchone()
    return result

def get_vendor_smiles(cursor, vendor, table='building_blocks'):
    cmd = f"select * from {table} where {vendor}_id is not NULL"
    cursor.execute(cmd)
    result = cursor.fetchall()
    return result

def test_smi_fp(cursor, smi, table='building_blocks', fps='fps'):
    cmd = f"select {table}.smiles, bfp_to_binary_text({fps}.rdkit_fp) " \
          f"from {table} " \
          f"inner join {fps} on {table}.id={fps}.id " \
          f"where smiles='{smi}'"
    cursor.execute(cmd)
    result = cursor.fetchall()
    return result

def vendor_rows(cursor, vendor, table='building_blocks'):
    cmd = f"select count(*) from {table} where {vendor}_id is not NULL"
    cursor.execute(cmd)
    result = cursor.fetchone()

    return result

def get_vendor_counts(cursor, vendors=VENDORS, table='building_blocks'):
    for v in VENDORS:
        num = vendor_rows(cursor, v)
        print(f"{v} - {num[0]}")

if __name__ == '__main__':
    conn = connect_to_db()
    cursor = conn.cursor()
    result = smiles_lookup('CCCCO', cursor)
    print(result)

    get_vendor_counts(cursor)


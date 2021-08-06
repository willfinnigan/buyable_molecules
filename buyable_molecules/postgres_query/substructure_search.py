import time
from buyable_molecules.create_postgres_db.funcs.db_setup import connect_to_db
from rdkit import DataStructs

def substructure_search(smarts, vendors=[], cursor=None, table='building_blocks', mols='mols', fps='fps', print_cmd=False):
    if cursor is None:
        cursor = connect_to_db().cursor()

    cmd = f"select {table}.*, bfp_to_binary_text({fps}.rdkit_fp) from {table} " \
          f"inner join {fps} on {table}.id={fps}.id " \
          f"where {table}.id in " \
          f"(SELECT {mols}.id from {mols} " \
          f"where {mols}.m@>'{smarts}'::qmol) " \

    num_vendors = len(vendors)
    if num_vendors != 0:
        cmd += " and ("

    for i, vendor in enumerate(vendors):
        cmd += f"{vendor}_id is not NULL"

        if i != num_vendors - 1:
            cmd += " or "
        else:
            cmd += ")"

    if print_cmd == True:
        print(cmd)

    cursor.execute(cmd)
    results = cursor.fetchall()

    new_results = []
    for i, result in enumerate(results):
        fp = DataStructs.CreateFromBinaryText(bytes(result[-1]))
        new_results.append(list(result[0:-1]) + [fp])

    return new_results

def substructure_search_single_vendor(smarts, vendor, cursor=None, table='building_blocks', mols='mols', fps='fps', print_cmd=False):
    if cursor is None:
        cursor = connect_to_db().cursor()

    cmd = f"SELECT {table}.*, bfp_to_binary_text({fps}.rdkit_fp) FROM " \
                f"(select {mols}.* from {table} " \
                f"inner join {mols} on {table}.id={mols}.id " \
                f"where {table}.{vendor}_id is not NULL) sq_mol " \
          f"inner join {fps} on {fps}.id=sq_mol.id " \
          f"inner join {table} on {table}.id=sq_mol.id " \
          f"where sq_mol.m@>'{smarts}'::qmol"

    if print_cmd == True:
        print(cmd)

    cursor.execute(cmd)
    results = cursor.fetchall()

    new_results = []
    for i, result in enumerate(results):
        fp = DataStructs.CreateFromBinaryText(bytes(result[-1]))
        new_results.append(list(result[0:-1]) + [fp])

    return new_results


if __name__ == "__main__":
    t0 = time.time()
    result = sq_substructure_search('CC(=O)[OH]', vendor='sigma')
    print(result)
    print(f"{round(time.time()-t0,2)} seconds")



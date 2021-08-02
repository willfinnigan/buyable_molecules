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


if __name__ == "__main__":

    result = substructure_search('CC(=O)[OH]', vendors=['sigma'])
    print(result)



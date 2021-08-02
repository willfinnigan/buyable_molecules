import psycopg2
import time


def create_mol_table(conn, cursor, table='building_blocks', mols='mols', id='id'):
    cmd = f"select * into {mols} from (select id,mol_from_smiles(smiles::cstring) m from {table}) tmp where m is not null;"
    cursor.execute(cmd)
    conn.commit()

    cmd = f"create index molidx on {mols} using gist(m)"
    cursor.execute(cmd)
    conn.commit()

    cmd = f"CREATE INDEX {mols}_idx_id ON {mols} ({id})"
    cursor.execute(cmd)
    conn.commit()

def create_fingerprints(conn, cursor, mols='mols', fps='fps', id_col='id'):
    from rdkit import Chem
    from rdkit.Chem import rdFingerprintGenerator
    from rdkit import DataStructs

    fp_list = []
    fp_gen = rdFingerprintGenerator.GetRDKitFPGenerator()
    # grab the molecules, we're pulling them out in their pickled form:
    print('getting molecules..')
    cursor.execute(f'select id, mol_send(m) from {mols}')
    results = cursor.fetchall()
    num_results = len(results)
    print(f"{num_results} molecules retrieved")
    print(f"creating fingerprints..")
    t0 = time.time()
    for i, (id, pkl) in enumerate(results):
        if pkl is None: continue
        m = Chem.Mol(pkl.tobytes())
        fp = fp_gen.GetFingerprint(m)
        fp_list.append((id, fp))

        if i % 10000 == 0 and i != 0:
            elapsed = round(time.time() - t0, 1)
            estimated_total_time = (elapsed / i) * num_results
            time_remaining = estimated_total_time - elapsed
            print(f"{i} added in {elapsed} seconds - estimated time remaining = {round(time_remaining / 60, 2)} mins")

    print('done')

    print('writing fps..')
    cursor.execute(f'drop table if exists {fps}')
    cursor.execute(f'create table {fps} (id int, rdkit_fp bfp)')
    cursor.executemany(f'insert into {fps} values (%s,bfp_from_binary_text(%s))',
                      [(x, DataStructs.BitVectToBinaryText(y)) for x, y in fp_list])
    conn.commit()
    print('done')

    print('making indexes..')
    cmd = f"CREATE INDEX {fps}_idx_id ON {fps}({id_col});"
    cursor.execute(cmd)
    conn.commit()
    #cmd = f"create index {fps}_rdfp_idx on {fps} using gist(rdkit_fp);"
    #cursor.execute(cmd)
    #conn.commit()
    print('done')



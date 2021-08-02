from pathlib import Path
import pandas as pd
from buyable_molecules.create_postgres_db.funcs import db_setup, create_building_blocks_table, create_mol_and_fp_tables

SAVE_FOLDER = str(Path(__file__).parents[1]) + '/save_folder'

def create_db():
    print("creating db if it does not already exist")
    conn = db_setup.connect_to_db()
    cursor = conn.cursor()
    db_setup.install_rd_extension(conn, cursor)
    cursor.close()
    conn.close()

def create_building_blocks():
    print('Creating building_blocks table from df - this will take some time..')
    df = pd.read_csv(f"{SAVE_FOLDER}/final.csv", index_col=0, dtype={'SMILES': str, 'molport_id': str, 'mcule_id': str, 'sigma_id': str, 'zinc_id': str})
    df.rename(columns={"SMILES": "smiles"}, inplace=True)
    print(df.head())

    conn = db_setup.connect_to_db()
    cursor = conn.cursor()

    create_building_blocks_table.sqlalchemy_data_to_table(df)
    create_building_blocks_table.create_primary_key(conn, cursor)

    cursor.close()
    conn.close()


def create_mols():
    print('Creating mols and fps tables - this will take some time..')
    conn = db_setup.connect_to_db()
    cursor = conn.cursor()

    create_mol_and_fp_tables.create_mol_table(conn, cursor)
    create_mol_and_fp_tables.create_fingerprints(conn, cursor)

    cursor.close()
    conn.close()

def make_db():
    create_db()
    create_building_blocks()
    create_mols()

if __name__ == '__main__':
    make_db()


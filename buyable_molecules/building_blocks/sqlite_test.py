import sqlite3
import pandas as pd
from pathlib import Path
from buyable_molecules.rdkit_pandas_functions import create_mol_col

csv_folder = str(Path(__file__).parents[0]) + '/building_block_dfs'

def load_extension(connection):
    connection.enable_load_extension(True)
    connection.load_extension('chemicalite')
    connection.enable_load_extension(False)

def make_table(connection, table="BUYABLE"):
    c = connection.cursor()
    try:
        c.execute(f'''DROP TABLE {table}''')
        connection.commit
    except:
        pass

    c = connection.cursor()
    c.execute(f'''CREATE TABLE {table}
                 ([id] INTEGER PRIMARY KEY, [SMILES] text, [molecule] MOL, [sigma_id] text)''')
    connection.commit()

def df_to_sql2(connection, df, table='BUYABLE'):
    def iterate_df(df):
        for index, row in df.iterrows():
            yield row['SMILES'], row['sigma_id']

    with connection:
        connection.executemany(
            f"INSERT INTO {table}(smiles, sigma_id, molecule) "
            "VALUES(?1, ?2, mol_from_smiles(?1))", iterate_df(df))

def df_to_sql(connection, df, table='BUYABLE'):
    df.to_sql(table, connection, if_exists='replace', dtype={'SMILES': 'text', 'molecule': 'MOL', 'sigma_id': 'text'})

def create_index(db_path, name, table='BUYABLE'):
    conn = sqlite3.connect(db_path)
    c = conn.cursor()
    c.execute(f"CREATE INDEX index_{name} ON {table}({name})")
    conn.commit()



if __name__ == '__main__':
    connection = sqlite3.connect('testdb.sql')
    load_extension(connection)

    '''
    make_table(connection)
    df = pd.read_csv(f"{csv_folder}/sigma.csv", index_col=0)
    df['molecule'] = df['SMILES'].apply(create_mol_col)
    df = df.dropna()
    df = df[['SMILES', 'sigma_id']]
    print(df.head())

    df_to_sql2(connection, df, table='BUYABLE')
    '''

    #connection.execute("CREATE VIRTUAL TABLE str_idx_buyable_molecule " +
     #                  "USING rdtree(id, fp bits(2048))")


    c = connection.cursor()
    c.execute(
        "INSERT INTO str_idx_buyable_molecule(id, fp) " +
        "SELECT id, mol_pattern_bfp(molecule, 2048) FROM BUYABLE " +
        "WHERE molecule IS NOT NULL")
    connection.commit()



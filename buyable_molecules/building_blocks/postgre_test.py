import psycopg2
from psycopg2.extensions import ISOLATION_LEVEL_AUTOCOMMIT
from pathlib import Path
from io import StringIO
import pandas as pd
import sqlalchemy

bb_folder = str(Path(__file__).parents[0]) + '/building_block_dfs'
final_building_block_df_path = f'{bb_folder}/final/final.csv'


def create_db(conn_string, db_name):
    conn = psycopg2.connect(conn_string)
    conn.set_isolation_level(ISOLATION_LEVEL_AUTOCOMMIT)
    cursor = conn.cursor()
    cmd = f"create database {db_name}"
    cursor.execute(cmd)
    cursor.close()
    conn.close()
    print('database created')

def install_rd_extension(conn, cursor):
    cmd = "create extension if not exists rdkit;"
    cursor.execute(cmd)
    conn.commit()

def connect_to_db():
    try:
        conn = psycopg2.connect("user='postgres' dbname='molecules' host='0.0.0.0' password='rbc_postgres_1' port='5432'")
    except:
        create_db("user='postgres'host='0.0.0.0' password='rbc_postgres_1' port='5432'", 'molecules')
        conn = psycopg2.connect("user='postgres' dbname='molecules' host='0.0.0.0' password='rbc_postgres_1' port='5432'")

    return conn

def create_building_blocks_table(conn, cursor):
    # mcule_id,sigma_id,molport_id,zinc_id
    try:
        cmd = "CREATE TABLE building_blocks (id integer PRIMARY KEY, smiles VARCHAR(255), mcule_id VARCHAR(255), sigma_id VARCHAR(255), moleport_id VARCHAR(255), zinc_id VARCHAR(255));"
        cursor.execute(cmd)
        conn.commit()
        print('new table created')
    except:
        print('table not created')


def data_to_table(conn, df, table):
    buffer = StringIO()
    df.to_csv(buffer, index_label='id', header=False)
    buffer.seek(0)

    cursor = conn.cursor()
    try:
        cursor.copy_from(buffer, table, sep=',')
        conn.commit()
        print('data copy complete')
        cursor.close()
        return

    except (Exception, psycopg2.DatabaseError) as error:
        print(f"Error: {error}")
        conn.rollback()
        cursor.close()
        return

def sqlalchemy_data_to_table(df, table):
    engine = sqlalchemy.create_engine('postgresql://postgres:rbc_postgres_1@localhost:5432/molecules')
    df.to_sql(table, con=engine, index_label='id', if_exists='replace')

def create_smiles_table_primary_key(cursor, conn, table, id_col='id'):
    cmd = "alter table {table} add primary key ({id});"
    cursor.execute(cmd)
    conn.commit()

def create_mol_table(table):
    cmd = f"select * into mols from (select id,mol_from_smiles(SMILES::cstring) m from {table}) tmp where m is not null;"
    cursor.execute(cmd)
    conn.commit()

def create_mol_index():
    cmd = "create index molidx on mols using gist(m)"
    cursor.execute(cmd)
    conn.commit()

def table_columns(table, cursor):
    cursor.execute(f"Select * FROM {table} LIMIT 0")
    colnames = [desc[0] for desc in cursor.description]
    return colnames

if __name__ == '__main__':

    conn = connect_to_db()
    cursor = conn.cursor()
    install_rd_extension(conn, cursor)

    table = 'new_building_blocks'

    df = pd.read_csv(final_building_block_df_path, index_col=0, dtype={'SMILES': str, 'molport_id': str, 'mcule_id': str, 'sigma_id': str, 'zinc_id': str})
    df.rename(columns={"SMILES": "smiles"}, inplace=True)
    print(df.head())

    sqlalchemy_data_to_table(df, table)

    cols = table_columns(table, cursor)
    print(cols)

    create_mol_table(table)
    create_mol_index()

    cursor.close()
    conn.close()






import psycopg2
from psycopg2.extensions import ISOLATION_LEVEL_AUTOCOMMIT
from io import StringIO
import sqlalchemy

def sqlalchemy_data_to_table(df, table='building_blocks'):
    print("data to postgress..")
    engine = sqlalchemy.create_engine('postgresql://postgres:rbc_postgres_1@localhost:5432/molecules')
    df.to_sql(table, con=engine, index_label='id', if_exists='replace')
    print('done')

def create_primary_key(conn, cursor, table='building_blocks', id_col='id'):
    print('Create primary key and index..')
    cmd = f"alter table {table} add primary key ({id_col});"
    cursor.execute(cmd)
    conn.commit()
    cmd = f"CREATE INDEX {table}_idx_id ON {table} ({id_col});"
    cursor.execute(cmd)
    conn.commit()
    print('done')

def create_smiles_index(conn, cursor, table='building_blocks', smi_col='smiles'):
    cmd = f"CREATE INDEX {table}_idx_smiles ON {table} ({smi_col});"
    cursor.execute(cmd)
    conn.commit()
    print('done')

if __name__ == '__main__':
    from buyable_molecules.create_postgres_db.funcs import db_setup
    conn = db_setup.connect_to_db()
    cursor = conn.cursor()

    create_smiles_index(conn, cursor)


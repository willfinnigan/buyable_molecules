import psycopg2
from psycopg2.extensions import ISOLATION_LEVEL_AUTOCOMMIT
from buyable_molecules.config import POSTGRES_PASSWORD, USER, HOST, PORT, DB_NAME

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
        conn = psycopg2.connect(
            f"user='{USER}' dbname='{DB_NAME}' host='{HOST}' password='{POSTGRES_PASSWORD}' port='{PORT}'")
    except:
        create_db(f"user='{USER}' host='{HOST}' password='{POSTGRES_PASSWORD}' port='{PORT}'", DB_NAME)
        conn = psycopg2.connect(
            f"user='{USER}' dbname='{DB_NAME}' host='{HOST}' password='{POSTGRES_PASSWORD}' port='{PORT}'")

    cursor = conn.cursor()
    install_rd_extension(conn, cursor)
    cursor.close()

    return conn

def table_columns(table, cursor):
    cursor.execute(f"Select * FROM {table} LIMIT 0")
    colnames = [desc[0] for desc in cursor.description]
    return colnames
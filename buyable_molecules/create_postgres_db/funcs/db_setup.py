import psycopg2
from psycopg2.extensions import ISOLATION_LEVEL_AUTOCOMMIT
from io import StringIO
import sqlalchemy

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
            "user='postgres' dbname='molecules' host='0.0.0.0' password='rbc_postgres_1' port='5432'")
    except:
        create_db("user='postgres'host='0.0.0.0' password='rbc_postgres_1' port='5432'", 'molecules')
        conn = psycopg2.connect(
            "user='postgres' dbname='molecules' host='0.0.0.0' password='rbc_postgres_1' port='5432'")

    cursor = conn.cursor()
    install_rd_extension(conn, cursor)
    cursor.close()

    return conn

def table_columns(table, cursor):
    cursor.execute(f"Select * FROM {table} LIMIT 0")
    colnames = [desc[0] for desc in cursor.description]
    return colnames
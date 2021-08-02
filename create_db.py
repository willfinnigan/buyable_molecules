from buyable_molecules.create_csvs.create_buyable_csv import create_buyable_csv
from buyable_molecules.create_postgres_db.create_postgres_db import make_db

if __name__ == "__main__":
    create_buyable_csv()
    make_db()



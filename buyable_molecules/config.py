import os

USER = os.environ.get('PS_USER') or "postgres"
POSTGRES_PASSWORD = os.environ.get('PS_PASSWORD') or "secret_password"
HOST = os.environ.get('PS_HOST') or "0.0.0.0"
PORT = os.environ.get('PS_PORT') or "5432"
DB_NAME = os.environ.get('PS_DB_NAME') or "molecules"
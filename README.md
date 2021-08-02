This is a quick module for creating and then querying a database of buyable building block molecules from a selection of vendors.

The mcule building blocks, zinc "in-stock", sigma building blocks and molport building blocks are included.  
All catalogs except for mcule are downloaded from Zinc.

This module requires a postgres database running on the default port, with the rdkit extension available.

One option is to use a docker container (eg from https://github.com/mcs07/docker-postgres-rdkit)
`docker run -p 5432:5432 -e POSTGRES_PASSWORD=secret_password -d mcs07/postgres-rdkit`

Settings for the postgres connection can be changed in config.py, or by using environment variables.  Otherwise default settings are used.

edit urls.py to change to files that are downloaded (or if they go out of date)

Run create_db.py to create the database (will take a long to time to run)  

Run buyable_molecules.substructure_search(smarts, vendors=[]) to perform a substructure search. Optionally filter by vendor.  This function takes about 5 minutes on my 2015 macbookpro.

Rdkit fingerprints are returned with the substructure search.

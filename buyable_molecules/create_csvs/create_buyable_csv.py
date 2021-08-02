from buyable_molecules.create_csvs.funcs.download_and_process import download_and_process, remove_file
from buyable_molecules.create_csvs.funcs import dataframe_functions
from pathlib import Path
from buyable_molecules import urls

SAVE_FOLDER = str(Path(__file__).parents[1]) + '/save_folder'
COLUMNS = ['SMILES', 'mcule_id', 'sigma_id', 'molport_id', 'zinc_id', 'apollo_id', 'flurochem_id']

def create_buyable_csv():
    mcule = download_and_process(urls.MCULE_URL, ["SMILES", "mcule_id"], None, '\t', f'{SAVE_FOLDER}/mcule.csv', download=True, unzip=True)
    sigma = download_and_process(urls.SIGMA_URL, ["SMILES", "sigma_id"], None, ' ', f'{SAVE_FOLDER}/sigma.csv', download=True, unzip=False)
    molport = download_and_process(urls.MOLPORT_URL, ["SMILES", "molport_id"], None, '\t', f'{SAVE_FOLDER}/molport.csv', download=True, unzip=False)
    zinc = download_and_process(urls.ZINC_IN_STOCK_URL, ["SMILES", "zinc_id", 'data_type', 'supplier', 'code'],
                         None, '\t', f'{SAVE_FOLDER }/zinc.csv', download=True, unzip=False)
    apollo = download_and_process(urls.APOLLO_BB_URL, ["SMILES", "apollo_id"], None, ' ', f'{SAVE_FOLDER}/apollo.csv',
                                 download=True, unzip=False)
    fluro = download_and_process(urls.FLUROCHEM_URL, ["SMILES", "flurochem_id"], None, ' ', f'{SAVE_FOLDER}/fluro.csv',
                                  download=True, unzip=False)

    merged_df = dataframe_functions.create_merged_df([mcule, sigma, molport, zinc, apollo, fluro])

    merged_df = merged_df[COLUMNS]
    merged_df.to_csv(f"{SAVE_FOLDER}/final.csv")

    remove_file(mcule)
    remove_file(sigma)
    remove_file(molport)
    remove_file(zinc)
    remove_file(apollo)
    remove_file(fluro)

if __name__ == "__main__":
    print('creating buyable csv..')
    create_buyable_csv()
    print('..done')



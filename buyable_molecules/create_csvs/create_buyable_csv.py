from buyable_molecules.create_csvs.funcs.download_and_process import download_and_process
from buyable_molecules.create_csvs.funcs import dataframe_functions
from pathlib import Path
from buyable_molecules import urls

SAVE_FOLDER = str(Path(__file__).parents[1]) + '/save_folder'

def create_buyable_csv():
    mcule = download_and_process(urls.MCULE_URL, ["SMILES", "mcule_id"], None, '\t', f'{SAVE_FOLDER}/mcule.csv', download=True, unzip=True)
    sigma = download_and_process(urls.SIGMA_URL, ["SMILES", "sigma_id"], None, ' ', f'{SAVE_FOLDER}/sigma.csv', download=True, unzip=False)
    molport = download_and_process(urls.MOLPORT_URL, ["SMILES", "molport_id"], None, '\t', f'{SAVE_FOLDER}/molport.csv', download=True, unzip=False)
    zinc = download_and_process(urls.ZINC_IN_STOCK_URL, ["SMILES", "zinc_id", 'data_type', 'supplier', 'code'],
                         None, '\t', f'{SAVE_FOLDER }/zinc.csv', download=True, unzip=False)

    merged_df = dataframe_functions.create_merged_df([mcule, sigma, molport, zinc])

    merged_df = merged_df[['SMILES', 'mcule_id', 'sigma_id', 'molport_id', 'zinc_id']]
    merged_df.to_csv(f"{SAVE_FOLDER}/final.csv")

if __name__ == "__main__":
    print('creating buyable csv..')
    create_buyable_csv()
    print('..done')



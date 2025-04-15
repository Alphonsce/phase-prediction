import os
import time
import random
import pandas as pd
from tqdm import tqdm
import requests

import numpy as np
import os
import rootutils

rootutils.setup_root(os.path.abspath('./'), indicator=".project-root", pythonpath=True, dotenv=True, cwd=True)

from src.cif_utils import get_smiles_from_id_on_web, get_smiles_from_chemical_name

USER_AGENTS = [
    "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/117.0.0.0 Safari/537.36",
    "Mozilla/5.0 (Macintosh; Intel Mac OS X 13_2) AppleWebKit/537.36 (KHTML, like Gecko) Firefox/116.0",
    "Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/119.0.0.0 Safari/537.36",
    # Add more user agents if needed
]

if __name__ == "__main__":
    ## Loading all ids:
    try:
        cifs_path = "/Users/aleksandr.varlamov/cif/all_cifs"
        cifs = sorted(os.listdir(cifs_path))
        cifs_id = list(map(lambda x: x.split(".")[0], cifs))
    except:
        cifs_id = np.load("cifs_id.npy").tolist()
    
    csv_path = "cod_parsed.csv"

    if os.path.exists(csv_path):
        results_df = pd.read_csv(csv_path)
    else:
        results_df = pd.DataFrame(columns=["id", "smiles", "common_name", "chemical_name"])

    for i, id in tqdm(enumerate(cifs_id), total=len(cifs_id)):
        if int(id) in results_df['id'].values:
            continue

        try:
            smiles, common_name, chemical_name = get_smiles_from_id_on_web(id)
            if smiles is None and common_name is not None:
                try:
                    smiles = get_smiles_from_chemical_name(common_name)
                    print(f"Successfully got smiles from name for {id}")
                except Exception as e:
                    if chemical_name is not None:
                        try:
                            smiles = get_smiles_from_chemical_name(chemical_name)
                            print(f"Successfully got smiles from name for {id}")
                        except Exception as e:
                            print(f"Failed to get smiles from name for {id}: {e}")
                            continue
                
            new_row = pd.DataFrame({"id": [id], "smiles": [smiles], "common_name": [common_name], "chemical_name": [chemical_name]})
            results_df = pd.concat([results_df, new_row], ignore_index=True)
            results_df.to_csv(csv_path, index=False)
            time.sleep(random.uniform(0.2, 0.5))  # Random sleep between 0.5 and 2.5 seconds

        except Exception as e:
            print(f"Failed for {id}: {e}")
            time.sleep(random.uniform(2, 5))  # Sleep longer after a failure
            continue



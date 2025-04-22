import os
import time
import random
import pandas as pd
from tqdm import tqdm
import concurrent.futures
import numpy as np
import rootutils
import argparse

# Set up the root environment as before
rootutils.setup_root(os.path.abspath('./'),
                       indicator=".project-root",
                       pythonpath=True,
                       dotenv=True,
                       cwd=True)

from src.cif_utils import get_smiles_from_id_on_web, get_smiles_from_chemical_name

USER_AGENTS = [
    "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/117.0.0.0 Safari/537.36",
    "Mozilla/5.0 (Macintosh; Intel Mac OS X 13_2) AppleWebKit/537.36 (KHTML, like Gecko) Firefox/116.0",
    "Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/119.0.0.0 Safari/537.36",
    "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Edge/120.0.0.0 Safari/537.36",
    "Mozilla/5.0 (Macintosh; Intel Mac OS X 14_0) AppleWebKit/605.1.15 (KHTML, like Gecko) Version/17.0 Safari/605.1.15",
    "Mozilla/5.0 (iPad; CPU OS 17_0 like Mac OS X) AppleWebKit/605.1.15 (KHTML, like Gecko) Version/17.0 Mobile/15E148 Safari/604.1",
    "Mozilla/5.0 (iPhone; CPU iPhone OS 17_0 like Mac OS X) AppleWebKit/605.1.15 (KHTML, like Gecko) Version/17.0 Mobile/15E148 Safari/604.1",
    "Mozilla/5.0 (Linux; Android 13; SM-S908B) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/112.0.0.0 Mobile Safari/537.36",
    "Mozilla/5.0 (Windows NT 10.0; Win64; x64; rv:109.0) Gecko/20100101 Firefox/119.0",
    "Mozilla/5.0 (X11; Ubuntu; Linux x86_64; rv:109.0) Gecko/20100101 Firefox/119.0",
    "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/118.0.0.0 Safari/537.36 OPR/104.0.0.0",
    "Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/118.0.0.0 Safari/537.36 OPR/104.0.0.0",
    "Mozilla/5.0 (Windows NT 10.0; WOW64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/118.0.0.0 Safari/537.36 Vivaldi/6.2.3105.58"
]

# Define a function to process a single CIF id
def process_id(cif_id):
    try:
        # Get smiles, common_name, and chemical_name using the first method
        smiles, common_name, chemical_name = get_smiles_from_id_on_web(cif_id)
        
        # If the first method didn't return smiles, try using chemical names
        if smiles is None and common_name is not None:
            try:
                smiles = get_smiles_from_chemical_name(common_name)
                print(f"Successfully got smiles from name for {cif_id}")
            except Exception as e:
                # If common_name lookup fails and chemical_name is available, try that
                if chemical_name is not None:
                    try:
                        smiles = get_smiles_from_chemical_name(chemical_name)
                        print(f"Successfully got smiles from chemical_name for {cif_id}")
                    except Exception as e:
                        print(f"Failed to get smiles from name for {cif_id}: {e}")
                        return {"id": cif_id, "smiles": smiles, "common_name": common_name, "chemical_name": chemical_name}
                else:
                    return {"id": cif_id, "smiles": smiles, "common_name": common_name, "chemical_name": chemical_name}

        # Mimic a random sleep to respect rate limits
        time.sleep(random.uniform(0.7, 1.0))
        
        return {"id": cif_id, "smiles": smiles, "common_name": common_name, "chemical_name": chemical_name}
    
    except Exception as e:
        print(f"Failed for {cif_id}: {e}")
        time.sleep(random.uniform(2, 5))  # Sleep longer after failure
        return None

def get_args():
    parser = argparse.ArgumentParser(description='Parse COD HTML pages in parallel to extract chemical information')
    parser.add_argument('--max-workers', type=int, default=16, help='Maximum number of worker threads (default: 4)')
    parser.add_argument('--csv-path', type=str, default="cod_parsed.csv", help='Path to the output CSV file (default: cod_parsed.csv)')
    parser.add_argument('--write-every', type=int, default=20, help='Write to CSV every N processed entries (default: 100)')
    args = parser.parse_args()
    return args

if __name__ == "__main__":
    args = get_args()
    max_workers = args.max_workers
    csv_path = args.csv_path
    write_every_steps = args.write_every
    
    # Load the list of cif ids
    try:
        cifs_path = "/Users/aleksandr.varlamov/cif/all_cifs"
        cifs = sorted(os.listdir(cifs_path))
        cifs_id = [x.split(".")[0] for x in cifs]
    except Exception as e:
        cifs_id = np.load("cifs_id.npy").tolist()
    
    if os.path.exists(csv_path):
        results_df = pd.read_csv(csv_path)
    else:
        results_df = pd.DataFrame(columns=["id", "smiles", "common_name", "chemical_name"])

    # Filter IDs that have not been processed already
    processed_ids = set(results_df['id'].astype(str).values)
    ids_to_process = [cif_id for cif_id in cifs_id if cif_id not in processed_ids]
    
    # Counter for periodic CSV writing
    count = 0
    with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as executor:
        future_to_id = {executor.submit(process_id, cif_id): cif_id for cif_id in ids_to_process}
        
        for future in tqdm(concurrent.futures.as_completed(future_to_id), total=len(future_to_id), desc="Processing IDs"):
            cif_id = future_to_id[future]
            try:
                result = future.result()
                if result is not None:
                    new_row = pd.DataFrame([result])
                    results_df = pd.concat([results_df, new_row], ignore_index=True)
                    count += 1
                    
                    # Write to CSV every 100 iterations.
                    if count % write_every_steps == 0:
                        results_df.to_csv(csv_path, index=False)
                        print(f"Flushed to CSV after {count} updates.")
            except Exception as e:
                print(f"Error processing {cif_id}: {e}")

    # Final write to CSV to make sure all results are saved.
    results_df.to_csv(csv_path, index=False)
    print("Processing complete.")

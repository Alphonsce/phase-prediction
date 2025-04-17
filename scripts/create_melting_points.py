import os
import pandas as pd
from tqdm import tqdm
from pymatgen.io.cif import CifParser
import rootutils

rootutils.setup_root(os.path.abspath('./'), indicator=".project-root", pythonpath=True, dotenv=True, cwd=True)

from src.cif_utils import cif_from_file

if __name__ == "__main__":
    cifs_path = "/Users/aleksandr.varlamov/cif/all_cifs"
    cifs = sorted(os.listdir(cifs_path))
    cifs_id = list(map(lambda x: x.split(".")[0], cifs))

    # Check if melting_points.csv already exists and load it
    if os.path.exists('melting_points.csv'):
        existing_df = pd.read_csv('melting_points.csv')
        existing_ids = set(existing_df['id'].astype(str))
        print(f"Found existing melting_points.csv with {len(existing_ids)} entries")
    else:
        existing_df = pd.DataFrame(columns=['id', 'T'])
        existing_ids = set()
        print("No existing melting_points.csv found, creating new file")
        existing_df.to_csv('melting_points.csv', index=False)

    # Process CIFs and write data after each iteration
    batch_size = 100
    current_batch = []
    processed_count = 0

    for i, id in enumerate(tqdm(cifs_id, total=len(cifs_id))):
        # Skip if ID already exists in the CSV
        if id in existing_ids:
            continue
            
        path = os.path.join(cifs_path, f"{id}.cif")
        try:
            cif = cif_from_file(path)
            # Check if melting point exists in the cif dictionary
            melting_point = cif.get('_chemical_melting_point', None)
            # Replace '?' with None
            if melting_point == '?':
                melting_point = None
            current_batch.append({'id': id, 'T': melting_point})
        except Exception as e:
            print(f"Error reading {path}: {e}")
            current_batch.append({'id': id, 'T': None})
            continue
        
        # Write data after each batch
        if len(current_batch) >= batch_size:
            new_df = pd.DataFrame(current_batch)
            # Append to existing CSV file
            new_df.to_csv('melting_points.csv', mode='a', header=False, index=False)
            processed_count += len(current_batch)
            print(f"Added {len(current_batch)} entries to melting_points.csv (Total: {processed_count})")
            # Update existing IDs to avoid duplicates
            existing_ids.update(new_df['id'].astype(str))
            current_batch = []
    
    # Write any remaining entries
    if current_batch:
        new_df = pd.DataFrame(current_batch)
        new_df.to_csv('melting_points.csv', mode='a', header=False, index=False)
        processed_count += len(current_batch)
        print(f"Added final {len(current_batch)} entries to melting_points.csv (Total: {processed_count})")
    
    if processed_count > 0:
        print(f"Successfully added {processed_count} new entries to melting_points.csv")
    else:
        print("No new entries to add")
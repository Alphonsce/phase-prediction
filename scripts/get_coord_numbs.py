import os
import argparse
import pandas as pd
from tqdm.auto import tqdm
from concurrent.futures import ProcessPoolExecutor, as_completed

import rootutils
rootutils.setup_root(
    os.path.abspath('./'),
    indicator='.project-root',
    pythonpath=True,
    dotenv=True,
    cwd=True
)
from src.cif_utils import get_fixed_length_descriptor
# from src.__cutoffs import cutoffs_gpt as cutoffs
from src.__vdw_cutoffs import vdw_cutoffs as cutoffs

# ------------ Worker Function ------------
def process_entry(entry, cutoffs, struct_repeat):
    """
    Worker function to read a CIF, compute the descriptor,
    and package it up as a dict. On failure, returns a dict
    with an 'error' key.
    """
    cif_path = entry['cif_path']
    try:
        coord_numbs = get_fixed_length_descriptor(
            cif_path,
            cutoffs,
            use_interior=True,
            struct_repeat=struct_repeat,
        )
        row = {
            'id': entry['id'],
            'smiles': entry['can_smiles'],
            **coord_numbs.to_dict()
        }
        return {'result': row}
    except Exception as e:
        return {'error': cif_path, 'exception': str(e)}
    
def create_arguments():
    meta_csv = "data_cod/cod_bradley_merged.csv"
    cifs_dir = "cifs"
    
    parser = argparse.ArgumentParser(
        description='Parallel CIF descriptor computation.'
    )
    parser.add_argument('--meta-csv', type=str, default=meta_csv,
                        help='Path to metadata CSV file')
    parser.add_argument('--cifs-dir', type=str, default=cifs_dir,
                        help='Directory containing .cif files')
    parser.add_argument('--result-csv', type=str, default='coord_numbs.csv',
                        help='Output CSV path')
    parser.add_argument('--save-every', type=int, default=5,
                        help='Save results after this many entries processed')
    parser.add_argument('--struct-repeat', nargs=3, type=int, default=[2, 2, 2],
                        metavar=('X', 'Y', 'Z'),
                        help='Structure repetition factors (3 ints)')
    parser.add_argument('--num-workers', type=int, default=8,
                        help='Number of parallel worker processes')
    return parser.parse_args()

def main():
    args = create_arguments()

    # Read metadata
    meta_df = (pd.read_csv(args.meta_csv)
               .drop_duplicates()
               .sort_values('id'))
    meta_df['cif_path'] = meta_df['id'].astype(str).apply(
        lambda x: os.path.join(args.cifs_dir, f"{x}.cif")
    )
    entries = meta_df.to_dict('records')

    # Prepare or load result DataFrame
    if os.path.exists(args.result_csv):
        cn_df = pd.read_csv(args.result_csv)
        processed = len(cn_df)
        # Skip entries already processed
        existing_ids = set(cn_df['id'])
        entries = [e for e in entries if e['id'] not in existing_ids]
    else:
        cn_df = pd.DataFrame()
        processed = 0

    error_cifs = []
    saved_count = processed

    # Set up ProcessPoolExecutor
    with ProcessPoolExecutor(max_workers=args.num_workers) as executor:
        futures = {
            executor.submit(
                process_entry, entry, cutoffs, tuple(args.struct_repeat)
            ): idx for idx, entry in enumerate(entries)
        }

        for future in tqdm(as_completed(futures), total=len(entries),
                           desc='Processing CIFs'):
            out = future.result()
            if 'result' in out:
                cn_df = pd.concat([cn_df, pd.DataFrame([out['result']])],
                                  ignore_index=True)
                saved_count += 1
            else:
                error_cifs.append(out['error'])
                print(f"Failed: {out['error']} — {out['exception']}")

            # Periodic save
            if saved_count and saved_count % args.save_every == 0:
                cn_df.to_csv(args.result_csv, index=False,
                            float_format='%.3f')
                print(f"Saved {saved_count} records so far.")

    # Final save and summary
    cn_df.to_csv(args.result_csv, index=False, float_format='%.3f')
    print(f"Done. Total processed: {saved_count}/{processed + len(entries)}.")
    if error_cifs:
        print(f"Encountered errors for {len(error_cifs)} files:")
        for p in error_cifs:
            print('  -', p)

if __name__ == '__main__':
    main()

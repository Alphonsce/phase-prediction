from pymatgen.io.cif import CifParser, CifFile, CifBlock
import json
import pubchempy as pcp
from rdkit import Chem

import requests
from bs4 import BeautifulSoup
import pandas as pd
import random

USER_AGENTS = [
    "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/117.0.0.0 Safari/537.36",
    "Mozilla/5.0 (Macintosh; Intel Mac OS X 13_2) AppleWebKit/605.1.15 (KHTML, like Gecko) Version/16.3 Safari/605.1.15",
    "Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/119.0.0.0 Safari/537.36",
]

def get_canonical_smiles(row: pd.Series) -> pd.Series:
    """
    Get the canonical SMILES string for a row in a DataFrame
    """
    try:
        mol = Chem.MolFromSmiles(row["smiles"])
        canonical_smiles = Chem.MolToSmiles(mol, isomericSmiles=True, canonical=True)
        row["can_smiles"] = canonical_smiles
        return row
    except Exception as e:
        row["can_smiles"] = None
        return row

def get_smiles_from_chemical_name(chemical_name: str) -> str:
    """
    Get the SMILES string for a given chemical name using PubChem.
    """
    compound_name = chemical_name
    compounds = pcp.get_compounds(compound_name, 'name')

    if compounds:
        smiles = compounds[0].isomeric_smiles
        return smiles
    else:
        raise ValueError(f"Compound not found: {chemical_name}")

def get_smiles_from_id_on_web(cod_id: str) -> str:
    """
    Get the SMILES string and common name for a given COD ID from the COD web page.
    """
    url = f"https://qiserver.ugr.es/cod/{cod_id}.html"
    headers = {
        "User-Agent": random.choice(USER_AGENTS),
        "Accept-Language": "en-US,en;q=0.9",
        "Referer": "https://qiserver.ugr.es/",
        "DNT": "1",
        "Connection": "keep-alive"
    }

    try:
        response = requests.get(url, headers=headers, timeout=10)
        response.raise_for_status()
    except Exception as e:
        raise Exception(f"Failed to fetch page: {url} with error: {e}")

    soup = BeautifulSoup(response.text, 'html.parser')
    tables = soup.find_all('table')
    
    smiles = None
    common_name = None

    for table in tables:
        for row in table.find_all('tr'):
            header = row.find('th')
            data = row.find('td')
            if not header or not data:
                continue
            if 'SMILES' in header.text:
                smiles = data.text.strip()
            elif 'Common name' in header.text:
                common_name = data.text.strip()
    
    return smiles, common_name

def cif_from_file(file_path: str):
    """
    Read a .cif file and return a dict with all fields
    """
    parser = CifParser(file_path)
    data_dict = list(parser.as_dict().values())[0]  # only data, no id
    
    return data_dict

def dict_to_cif_block(dict_cif: dict, id: str) -> CifBlock:
    block = CifBlock(dict_cif, loops=[], header=id)
    return block

# ----- Coordination number functions -----
from pymatgen.core import Structure
import numpy as np
import pandas as pd
from collections import defaultdict

def is_interior(site, frac_margin):
    """
    Check if a site is fully inside the cell by a margin in fractional coords.
    """
    f = site.frac_coords
    return np.all(f > frac_margin) and np.all(f < 1 - frac_margin)


def compute_avg_cn_df(
        cif_path: str,
        raw_cutoffs: dict,
        use_interior: bool = False,
        struct_repeat: tuple = (1, 1, 1)
    ) -> pd.DataFrame:
    """
    Compute average coordination numbers for each central->neighbor pair.

    Parameters:
    - cif_path: path to CIF file
    - raw_cutoffs: dict of (element1, element2) -> cutoff distance (Å)
    - use_interior: if True, filter out atoms near cell boundaries for each cutoff
    - struct_repeat: how many times to repeat the unit cell in (a,b,c)

    Returns:
    - DataFrame with columns: central, neighbor, avg_cn, n_sites
    """
    # Normalize cutoffs to unordered keys for symmetry
    cutoffs = {frozenset((e1, e2)): cutoff
               for (e1, e2), cutoff in raw_cutoffs.items()}

    struct = Structure.from_file(cif_path)
    struct = struct * struct_repeat
    elems = sorted({site.specie.symbol for site in struct})

    # Precompute per-cutoff fractional margins if interior filtering
    frac_margins = {}
    if use_interior:
        # Use minimal lattice parameter to convert Å to fractional
        min_axis = min(struct.lattice.abc)
        for pair, cutoff in cutoffs.items():
            frac_margins[pair] = cutoff / min_axis

    # Collect counts
    coord_numbers = defaultdict(lambda: defaultdict(list))
    for site in struct:
        central = site.specie.symbol
        # Loop over each defined pair and its cutoff
        for pair, cutoff in cutoffs.items():
            if central not in pair:
                continue
            # if using interior, skip sites too close to boundary for this cutoff
            if use_interior:
                if not is_interior(site, frac_margins[pair]):
                    continue
            # Determine partner (handles homo- and heteronuclear)
            if len(pair) == 1:
                partner = central
            else:
                partner = next(iter(pair - {central}))
            # Count neighbors of this partner type within cutoff
            neighs = struct.get_neighbors(site, cutoff)
            cnt = sum(1 for n in neighs if n.specie.symbol == partner)
            coord_numbers[central][partner].append(cnt)

    # Build full DataFrame (including homo- and heteronuclear)
    rows = []
    for central in elems:
        for partner in elems:
            counts = coord_numbers[central].get(partner, [])
            avg_cn = float(np.mean(counts)) if counts else 0.0
            n_sites = len(counts)
            rows.append({
                "central": central,
                "neighbor": partner,
                "avg_cn": avg_cn,
                "n_sites": n_sites
            })
    return pd.DataFrame(rows)


def get_fixed_length_descriptor(
        cif_path: str,
        raw_cutoffs: dict,
        use_interior: bool = False,
        struct_repeat: tuple = (1, 1, 1)
    ) -> pd.Series:
    """
    Generate a fixed-length avg_cn descriptor vector from a CIF.

    Features are ordered as [E1->E2, E2->E1] for each (E1, E2) in sorted raw_cutoffs.
    If a pair is not present, avg_cn = 0.

    Returns:
    - pandas Series indexed by feature names "E1–E2" with avg_cn values
    """
    # Compute per-pair stats
    df = compute_avg_cn_df(
        cif_path,
        raw_cutoffs,
        use_interior=use_interior,
        struct_repeat=struct_repeat
    )

    # Build ordered feature list (both directions)
    features = []
    for (e1, e2) in sorted(raw_cutoffs.keys()):
        features.append(f"{e1}–{e2}")
        features.append(f"{e2}–{e1}")

    # Map observed avg_cn
    avg_cn_map = {f"{row.central}–{row.neighbor}": row.avg_cn
                  for _, row in df.iterrows()}
    descriptor = {feat: avg_cn_map.get(feat, 0.0) for feat in features}

    return pd.Series(descriptor, name="avg_cn_descriptor")
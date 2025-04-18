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
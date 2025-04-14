from CifFile import ReadCif
import json
import pubchempy as pcp
from rdkit import Chem

import requests
from bs4 import BeautifulSoup
import pandas as pd

def get_canonical_smiles(row: pd.Series) -> pd.Series:
    """
    Get the canonical SMILES string for a row in a DataFrame
    """
    mol = Chem.MolFromSmiles(row["smiles"])
    canonical_smiles = Chem.MolToSmiles(mol, isomericSmiles=True, canonical=True)
    row["can_smiles"] = canonical_smiles
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
    
def read_cif(file_path: str) -> dict:
    """
    Read a .cif file and return a dict with all fields
    """
    cf = ReadCif(file_path)
    block = cf.first_block()
    return {key: block[key] for key in block.keys()}

def get_smiles_from_id_on_web(cod_id: str) -> str:
    """
    Get the SMILES string for a given COD ID from the COD web page.
    """
    url = f"https://www.crystallography.net/cod/{cod_id}.html"
    response = requests.get(url)
    if response.status_code != 200:
        raise Exception(f"Failed to fetch page: {url}")
    
    soup = BeautifulSoup(response.text, 'html.parser')
    
    # Find the table containing the SMILES string
    tables = soup.find_all('table')
    for table in tables:
        for row in table.find_all('tr'):
            header = row.find('th')
            if header and 'SMILES' in header.text:
                data = row.find('td')
                if data:
                    return data.text.strip()
    
    raise Exception(f"SMILES string not found for {cod_id}")
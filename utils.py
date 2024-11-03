import pandas as pd
import numpy as np
import rdkit
import matplotlib.pyplot as plt

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors, Descriptors, Lipinski
from rdkit.Chem import AllChem
from rdkit.Chem import rdFingerprintGenerator

from sklearn.preprocessing import normalize
from sklearn.decomposition import PCA

from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score, roc_auc_score, r2_score, balanced_accuracy_score
from sklearn.metrics import mean_absolute_error, mean_squared_error, root_mean_squared_error

def create_df(data_path, columns) -> pd.DataFrame:
    df = pd.read_csv(data_path, sep=' ', header=None, names=columns)
    df["label"] = df['label'].replace({1: 0, 2: 1})
    return df

def compute_fingerprints(smiles, radius, nbits):
    mol = Chem.MolFromSmiles(smiles)
    mfpgen = rdFingerprintGenerator.GetMorganGenerator(radius=radius, fpSize=nbits)
    return mfpgen.GetFingerprint(mol)

def compute_descriptors(smiles, descriptors):
    mol = Chem.MolFromSmiles(smiles)

    avail_descriptors = {
        'MolWt': lambda mol: Descriptors.MolWt(mol),
        'LogP': lambda mol: Descriptors.MolLogP(mol),
        'TPSA': lambda mol: Descriptors.TPSA(mol),
        'NumRotatableBonds': lambda mol: Descriptors.NumRotatableBonds(mol),
        'NumHDonors': lambda mol: Descriptors.NumHDonors(mol),
        'NumHAcceptors': lambda mol: Descriptors.NumHAcceptors(mol),
        'FractionCSP3': lambda mol: Lipinski.FractionCSP3(mol),
        'NumAromaticRings': lambda mol: Descriptors.NumAromaticRings(mol),
        'FractionRotatableBonds': lambda mol: rdMolDescriptors.CalcFractionCSP3(mol),

        # ----------
        'NumHBD': lambda mol: rdMolDescriptors.CalcNumHBD(mol),
        "NumHeavyAtoms": lambda mol: rdMolDescriptors.CalcNumHeavyAtoms(mol),
        # ----------
        
        'NumHBA': lambda mol: rdMolDescriptors.CalcNumHBA(mol),
        'NumRings': lambda mol: rdMolDescriptors.CalcNumRings(mol),
        'NumHeteroatoms': lambda mol: rdMolDescriptors.CalcNumHeteroatoms(mol),
        'Chi0v': lambda mol: rdMolDescriptors.CalcChi0v(mol),
        'Chi1v': lambda mol: rdMolDescriptors.CalcChi1v(mol),
        'Chi2v': lambda mol: rdMolDescriptors.CalcChi2v(mol),
    }

    X = []
    for desc_name in descriptors:
        X.append(avail_descriptors[desc_name](mol))
    
    return pd.Series(X)
    
def create_data(df, descriptors: list, apply_norm=False, radius=2, nbits=2048, fingerprints_pca=True, pca_dim=32) -> tuple:
    '''
    Arguments:
    -------
    radius: radius for Morgan fingerprints
    nbits: number of bits for Morgan fingerprints

    fingerprints_pca: use PCA on fingerprints features or not
    pca_dim: obvious

    Returns: (X_fps, X_at, y)
    ------
    X_fps: fingerprint features
    X_at: non-fingerprint atomic features
    y: label
    '''
    smiles = df['smiles']
    y = df['label'].values

    df['fingerprints'] = smiles.apply(compute_fingerprints, args=(radius, nbits,))
    df = df.dropna(subset=['fingerprints'])
    X_fps = np.array([np.array(fp) for fp in df['fingerprints']])
    if fingerprints_pca:
        pca = PCA(n_components=pca_dim)
        X_fps = pca.fit_transform(X_fps)

    X_at = smiles.apply(compute_descriptors, args=(descriptors,))
    X_at = np.array(X_at)

    if apply_norm:
        X_fps = normalize(X_fps)
        X_at = normalize(X_at)

    return X_fps, X_at, y

def eval_metrics(y_true, y_pred, type="classification"):
    if type == "classification":
        return {
            "ACC": accuracy_score(y_true, y_pred),
            "BAL-ACC": balanced_accuracy_score(y_true, y_pred),
            "ROC-AUC": roc_auc_score(y_true, y_pred),
            "R2": r2_score(y_true, y_pred)
        }
    return {
        "MSE": mean_squared_error(y_true, y_pred),
        "RMSE": root_mean_squared_error(y_true, y_pred),
        "MAE": mean_absolute_error(y_true, y_pred),
        "R2": r2_score(y_true, y_pred)
    }
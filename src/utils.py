import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import rdkit
import rootutils
from rdkit import Chem
from rdkit.Chem import (AllChem, Descriptors, Fragments, Lipinski,
                        rdFingerprintGenerator, rdMolDescriptors)
from sklearn.decomposition import PCA
from sklearn.metrics import (accuracy_score, balanced_accuracy_score,
                             mean_absolute_error, mean_squared_error, r2_score,
                             roc_auc_score)
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import normalize

from src.avail_descriptors import AVAIL_DESCRIPTORS

from tqdm.auto import tqdm

rootutils.setup_root(os.path.abspath('./'), indicator=".project-root", pythonpath=True, dotenv=True, cwd=True)
tqdm.pandas()

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
    mol = Chem.AddHs(mol)          # Add hydrogens
    AllChem.EmbedMolecule(mol)      # Generate a 3D conformer

    X = []
    for desc_name in descriptors:
        try:
            X.append(AVAIL_DESCRIPTORS[desc_name](mol))
        except:
            X.append(None)
    
    return pd.Series(X)
    
def create_data(df, descriptors: list, create_fingerprints=True, apply_norm=False, radius=2, nbits=2048, fingerprints_pca=True, pca_dim=32) -> tuple:
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

    if create_fingerprints:
        df['fingerprints'] = smiles.apply(compute_fingerprints, args=(radius, nbits,))
        df = df.dropna(subset=['fingerprints'])
        X_fps = np.array([np.array(fp) for fp in df['fingerprints']])
        if fingerprints_pca:
            pca = PCA(n_components=pca_dim)
            X_fps = pca.fit_transform(X_fps)
    else:
        X_fps = None

    X_at = smiles.progress_apply(compute_descriptors, args=(descriptors,))
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
        "MAE": mean_absolute_error(y_true, y_pred),
        "R2": r2_score(y_true, y_pred)
    }

def plot_pred_true(pred, true):
    plt.figure(figsize=(6, 4))

    plt.scatter(
        pred, true
    )

    plt.xlabel("T pred", fontsize=16)
    plt.ylabel("T true", fontsize=16)

    plt.legend(fontsize=16)

def plot_importance(importances, labels):
    plt.figure(figsize=(6, 4))

    plt.bar(range(len(importances)), importances, tick_label=labels)

    plt.grid(alpha=0.4)
    plt.xticks(rotation=90, fontsize=12)
    plt.ylabel("Feature Importance", fontsize=16)

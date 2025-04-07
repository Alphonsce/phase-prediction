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
                             roc_auc_score, f1_score)
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
    try:
        mol = Chem.MolFromSmiles(smiles)
        mfpgen = rdFingerprintGenerator.GetMorganGenerator(radius=radius, fpSize=nbits)
        return mfpgen.GetFingerprint(mol)
    except Exception as e:
        with open('error_smiles.txt', 'a') as f:
            f.write(f"Error in fingerprints: {smiles}\n{str(e)}\n")
        return None

def compute_descriptors(smiles, descriptors):
    try:
        mol = Chem.MolFromSmiles(smiles)
        mol = Chem.AddHs(mol)          # Add hydrogens
        # AllChem.EmbedMolecule(mol)      # Generate a 3D conformer, takes too long, sometimes failes for no reason

        X = []
        for desc_name in descriptors:
            X.append(AVAIL_DESCRIPTORS[desc_name](mol))
        
        return pd.Series(X)
    except Exception as e:
        with open('error_smiles.txt', 'a') as f:
            f.write(f"Error in descriptors: {smiles}\n{str(e)}\n")
        return pd.Series([None] * len(descriptors))
    
def create_data(df, descriptors: list, create_fingerprints=True, apply_norm=False, radius=2, nbits=2048, fingerprints_pca=True, pca_dim=32, temp_column=False) -> tuple:
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
 
    if temp_column:
        temp = df['T'].values

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
        if create_fingerprints:
            X_fps = normalize(X_fps)
        X_at = normalize(X_at)

    if temp_column:
        return X_fps, X_at, y, temp
    return X_fps, X_at, y

def create_or_load_data(df, data_args, save_or_load_dir, load_data=True):
    if not os.path.exists(save_or_load_dir):
        os.makedirs(save_or_load_dir)

    x_path = os.path.join(save_or_load_dir, 'X.npy')
    labels_path = os.path.join(save_or_load_dir, 'labels.npy')
    temp_path = os.path.join(save_or_load_dir, 'temp.npy')

    if load_data and all(os.path.exists(path) for path in [x_path, labels_path, temp_path]):
        print(f"Loading data from {save_or_load_dir}")
        X = np.load(x_path)
        labels = np.load(labels_path)
        temp = np.load(temp_path)

    elif load_data:
        raise ValueError(f"Data not found in {save_or_load_dir}")
    
    else:
        print(f"Creating data in {save_or_load_dir}")
        _, X, labels, temp = create_data(df, **data_args)
        np.save(x_path, X)
        np.save(labels_path, labels)
        np.save(temp_path, temp)
    
    return X, labels, temp

def eval_metrics(y_true, y_pred, type="classification"):
    if type == "classification":
        return {
            "ACC": accuracy_score(y_true, y_pred),
            "F1": f1_score(y_true, y_pred, average="binary"),
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
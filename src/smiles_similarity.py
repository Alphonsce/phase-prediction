from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import DataStructs





# Example reference set (smiles_list1) and candidate set (smiles_list2)



smiles_list1 = [
    "CCO",      # ethanol
    "CCN",      # ethylamine
    "CCOCC",    # diethyl ether
    "CC(C)O"    # isopropanol
]

f = open("merged_pub_bradley",'r')
smiles_list1=[]

for line in f:
    s=line.split()
    if s[2]=='2' and s[1].find(".")==-1 and s[1].find("+")==-1:
        smiles_list1.append(s[1])
f.close()


smiles_list2 = [
    "CCC",      # propane
    "CCOC",     # ethyl methyl ether
    "CCOC(C)=O",# ethyl acetate
    "CCCN",     # propylamine
    "COC"       # dimethyl ether
]

f = open("small_sim_brad",'r')
smiles_list2=[]

for line in f:
    s=line.split()
    if s[1].find(".")==-1 and s[1].find("+")==-1:
        smiles_list2.append(s[1])
f.close()

# Parameters:
# r2: similarity threshold (molecules are considered similar if Tanimoto similarity >= r2)
# k1: minimum number of similar molecules (from smiles_list1) that must be found for a candidate to pass.
r2 = 0.4
k1 = 3

# Fingerprint parameters
radius = 2
nBits = 2048

def generate_fps(smiles_list):
    """Converts a list of SMILES to RDKit molecules and returns a list of (smiles, mol, fingerprint) tuples."""
    results = []
    for smi in smiles_list:
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            print(f"Warning: Could not parse SMILES: {smi}")
            continue
        fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius=radius, nBits=nBits)
        results.append((smi, mol, fp))
    return results

# Generate fingerprints for both lists
ref_fps = generate_fps(smiles_list1)
cand_fps = generate_fps(smiles_list2)

# For each candidate molecule, count how many molecules in the reference set are similar (similarity >= r2)
selected_candidates = []

for cand_smi, cand_mol, cand_fp in cand_fps:
    similar_count = 0
    similar_refs = []
    for ref_smi, ref_mol, ref_fp in ref_fps:
        sim = DataStructs.TanimotoSimilarity(cand_fp, ref_fp)
        if sim >= r2:
            similar_count += 1
            similar_refs.append((ref_smi, sim))
    if similar_count >= k1:
        selected_candidates.append({
            "Candidate_SMILES": cand_smi,
            "Similar_Count": similar_count,
            "Similar_References": similar_refs
        })

# Print the results
print(f"Candidates from smiles_list2 that have at least {k1} molecules from smiles_list1 with similarity >= {r2}:")
for entry in selected_candidates:
    print("\nCandidate:", entry["Candidate_SMILES"])
    print("Number of similar reference molecules:", entry["Similar_Count"])
    print("Similar reference molecules (SMILES, similarity):")
    for ref in entry["Similar_References"]:
        print(f"  {ref[0]}  ({ref[1]:.3f})")


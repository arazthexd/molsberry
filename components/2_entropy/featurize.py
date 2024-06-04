import argparse
from pandarallel import pandarallel
import pandas as pd
from rdkit import Chem
from rdkit.Chem import PandasTools
from utils import (
    CalcRingFlexibility, 
    CalcHBondFoldability, 
    CalcPiPiStackFoldability, 
    CalcNumQueryMatches
)

SGS_SMARTS = [
    "[NX3]!@;-[CX3]=[OX1]",
    "[O]!@;-[CX3]=[OX1]",
    "[NX3]!@;-[CX3]=[SX1]"
]
SGS_QUERY = [Chem.MolFromSmarts(smarts) for smarts in SGS_SMARTS]

METHYL_QUERY = Chem.MolFromSmarts("[CH3]")
ROTOR_QUERY = Chem.MolFromSmarts("[!$(*#*)&!D1]-!@[!$(*#*)&!D1]")

# Function to calculate features
def calculate_features(row):

    mol = row["Mol"]
    return pd.Series({
        "RingFlexibility": CalcRingFlexibility(mol),
        "HBondFoldability": CalcHBondFoldability(mol),
        "PiPiStackFoldability": CalcPiPiStackFoldability(mol),
        "NumMethyls": CalcNumQueryMatches(mol, METHYL_QUERY),
        "NumRotors": CalcNumQueryMatches(mol, ROTOR_QUERY),
        "NumSpecifiedGroups": sum(CalcNumQueryMatches(mol, q) for q in SGS_QUERY)
    })

def smiles_to_mol(smi):
    return Chem.MolFromSmiles(smi)

if __name__ == "__main__":

    # Initialize parallel processing
    pandarallel.initialize(progress_bar=True)

    # Set up the argument parser
    parser = argparse.ArgumentParser(description="Featurize SMILES data.")
    parser.add_argument("--input", type=str, help="Path to the input file (csv)")
    parser.add_argument("--output", type=str, default="features/features.csv", help="Path to save the output")
    parser.add_argument("--include_conf_entropy", action="store_true", help="Whether to include the conf entropy in features file")
    parser.add_argument("--include_smiles", action="store_true", help="Whether to include smiles representation in features file")
    parser.add_argument("--feature_set", type=str, default="default", help="Feature set to calculate")
    parser.add_argument("--specific_groups", type=str, default="", help="Path to a txt or smi file with SMARTS patterns")
    args = parser.parse_args()

    # Read the CSV files
    df = pd.read_csv(f"{args.input}")

    # Add RDKit molecules
    print("\nConstructing Molecules From Smiles...")
    df["Mol"] = df["SMILES"].parallel_apply(smiles_to_mol)

    # Calculate features in parallel
    print("\nFeaturizing Molecules...")
    train_features = df.parallel_apply(calculate_features, axis=1)

    if args.include_smiles:
        train_features = pd.concat([df["SMILES"], train_features], axis=1)

    if args.include_conf_entropy:
        train_features = pd.concat([train_features, df["ConfEntropy"]], axis=1)

    # Save the features
    print("\nSaving Features...")
    train_features.to_csv(f"{args.output}", index=False)

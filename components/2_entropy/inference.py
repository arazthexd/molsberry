import argparse
import joblib

import numpy as np
import pandas as pd
from pandarallel import pandarallel
pandarallel.initialize(progress_bar=True)

from rdkit import Chem, RDLogger
from rdkit.Chem import PandasTools
RDLogger.DisableLog('rdApp.*')   

from featurize import calculate_features

parser = argparse.ArgumentParser(description="Predict Entropy For Ligands")
parser.add_argument("--input", type=str, help="Path to the input file (.sdf)")
parser.add_argument("--model", type=str, help="Path to the model file (.pkl)")
parser.add_argument("--output", type=str, help="Path to the output file with entropies (.csv)")
args = parser.parse_args()

print("\nLoading Model and SDF...")
model = joblib.load(args.model)

df: pd.DataFrame = PandasTools.LoadSDF(args.input, removeHs=False, molColName="Mol")
print(df.head())

print("\nCalculating Features...")
features = df.parallel_apply(calculate_features, axis=1)
df = pd.concat([df, features], axis=1)

print("\n\nPredicting Entropies...")
entropies = np.maximum(model.predict(np.log(features.values+1)), 0)
df["PredEntropy"] = pd.Series(entropies)

print("\nSaving Entropies in a CSV File...\n")
df["PredEntropy"].to_csv(args.output, index=False)
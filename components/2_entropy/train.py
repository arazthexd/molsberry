import argparse
import joblib

import pandas as pd
import numpy as np

from sklearn.linear_model import LinearRegression
from sklearn.neural_network import MLPRegressor
from sklearn.metrics import mean_absolute_error

parser = argparse.ArgumentParser(description="Train and Validate Entropy Model.")
parser.add_argument("--train", type=str, help="Path to the train input file (.csv)")
parser.add_argument("--holdout", type=str, help="Path to the holdout input file (.csv)")
parser.add_argument("--save", type=str, help="Path for saving the model (.pkl)")
args = parser.parse_args()

feature_names = ["RingFlexibility", "HBondFoldability", "PiPiStackFoldability", "NumMethyls", "NumRotors", "NumSpecifiedGroups"]

print("\nLoading Data...")
df_train = pd.read_csv(args.train)
x_train = np.log(df_train[feature_names].values + 1)
y_train = df_train["ConfEntropy"].values

df_holdout = pd.read_csv(args.holdout)
x_holdout= np.log(df_holdout[feature_names].values + 1)
y_holdout = df_holdout["ConfEntropy"].values

print("\nTraining Model...")
model = MLPRegressor((8), activation="tanh", max_iter=500)
model.fit(x_train, y_train)

print("train_mae: ", mean_absolute_error(y_train, model.predict(x_train)))
print("holdout_mae: ", mean_absolute_error(y_holdout, model.predict(x_holdout)))
print("train_score: ", model.score(x_train, y_train))
print("holdout_score: ", model.score(x_holdout, y_holdout))

print("\nSaving Model...")
joblib.dump(model, args.save)
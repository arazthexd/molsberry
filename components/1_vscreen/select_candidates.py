import argparse

import pandas as pd

from rdkit import Chem, RDLogger
from rdkit.Chem import PandasTools
RDLogger.DisableLog('rdApp.*') 

if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        prog="Docked Ligand Selector",
    )
    parser.add_argument("filename", help="path to docked ligands file (.sdf)")
    parser.add_argument("outfile", help="path to selected docked ligands file (.sdf)")
    args = parser.parse_args()

    df = PandasTools.LoadSDF(args.filename)
    n = len(df)

    df["ConfidentAffinity"] = pd.to_numeric(df["CNNaffinity"]) - pd.to_numeric(df["CNNaffinity_variance"]) * 3
    df["minimizedAffinity"] = pd.to_numeric(df["minimizedAffinity"])
    df_top10 = df.sort_values("ConfidentAffinity", ascending=False).iloc[:n//10]
    df["UIDNumTop10CNN"] = df["UniqueID"].apply(lambda uid: (df_top10["UniqueID"] == uid).sum())
    df["IDNumTop10CNN"] = df["ID"].apply(lambda uid: (df_top10["ID"] == uid).sum())
    df_top1 = df_top10.sort_values("CNNaffinity", ascending=False).iloc[:n//100]
    df["UIDNumTop1CNN"] = df["UniqueID"].apply(lambda uid: (df_top1["UniqueID"] == uid).sum())
    df["IDNumTop1CNN"] = df["ID"].apply(lambda uid: (df_top1["ID"] == uid).sum())
    def selection(row):
        s1 = row["IDNumTop1CNN"] > 1 or row["IDNumTop10CNN"] > 5 and row["minimizedAffinity"] < -4
        sorted_df = df[df["ID"] == row["ID"]].sort_values("CNNaffinity", ascending=False)
        return s1 and row.name == sorted_df.iloc[0].name
    df["Select"] = df.apply(selection, axis=1)
    df_write = df[df["Select"]]

    PandasTools.WriteSDF(df_write, args.outfile, idName="ID", properties=[
        "UniqueID", "SMILES", "minimizedAffinity", "ConfidentAffinity", "CNNaffinity", 
        "UIDNumTop10CNN", "IDNumTop10CNN", "UIDNumTop1CNN", "IDNumTop1CNN"])

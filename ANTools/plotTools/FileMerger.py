
import os
import json
import pickle as pk
import uproot3
import pandas as pd


def save_dfs(dfs, file_name):
    with open(file_name, 'wb') as f:
        pk.dump(dfs, f)

def load_dfs(file_name):
    with open(file_name, 'rb') as f:
        dfs = pk.load(f)
    return dfs

def read_single_rootfile_topd(filename, treepath,):
    tree = uproot3.open(filename)[treepath]
    df = tree.pandas.df(branches=tree.keys())
    print(f"{filename} has been opend!")
    return df

def read_samples(file_dir, treepath):
    
    df = pd.concat([read_single_rootfile_topd(file_dir + filename, treepath) for filename in os.listdir(file_dir) if filename[-5:] == ".root"])
    df = df.reset_index()
    del df["entry"]
    
    return df
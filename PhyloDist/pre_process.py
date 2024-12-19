import pandas as pd
import numpy as np
from sklearn import preprocessing
from Bio import Phylo


def input(x):
    dfx = pd.read_csv(x, sep=",", index_col=0)
    return dfx


def is_binary(df):
    binary = df.map(lambda x: x in [0, 1]).all().all()
    return binary


def to_binary(dfx, treshold=60):
    binary = is_binary(dfx)
    if binary is False:
        for i in range(0, len(dfx)):
            for j in range(0, len(dfx.columns)):
                if dfx.iloc[i][j] > treshold:
                    dfx.iloc[i][j] = 1
                else:
                    dfx.iloc[i][j] = 0
    return dfx
    
def normalize(dfx):
    binary = is_binary(dfx)
    if binary is False:
        for i in dfx.index:
            dfx.loc[i] = preprocessing.normalize(dfx.loc[i])
        return dfx
    else:
        return dfx


def transition_vector(dfx, binary):
    if binary:
        tv = np.zeros((len(dfx.index),(len(dfx.columns))))
        for j in range(0, len(dfx)):
            vec = dfx.iloc[j]
            pos_ori = np.asarray(vec[:-1])
            pos_new = np.asarray(vec[1:])
            trans_vec = pos_new-pos_ori
            tv[j] = np.pad(trans_vec, (1,0))
        tv = pd.DataFrame(tv, index=dfx.index, columns=dfx.columns)
        return tv
    else:
        return "Need binary profiles, you can use to_binary()"


def order_by_tree(dfx, tree):
    phylo = Phylo.read(tree, "newick")
    print(tree)
    leaf = phylo.get_terminals()
    print(leaf)
    ordered_df = pd.DataFrame()
    for l in leaf:  # noqa: E741
        ordered_df[str(l.name)] = dfx[str(l.name)]
    return ordered_df
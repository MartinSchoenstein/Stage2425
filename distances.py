import pandas as pd
from sklearn.metrics import jaccard_score
from scipy.spatial import distance
from scipy import stats
from sklearn.metrics import mutual_info_score
from sklearn import preprocessing
from newick import read
from Bio import Phylo


def distance_profiles(x, y, method):
    if isinstance(x, str):
        dfx, dfy = input(x, y)
    else:
        dfx, dfy = x, y
    if method == "Jaccard":
        return jaccard(dfx, dfy)
    if method == "Hamming":
        return hamming(dfx, dfy)
    if method == "Pearson":
        return pearson(dfx, dfy)
    if method == "MI":
        return mi(dfx, dfy)


def input(x, y):
    dfx = pd.read_csv(x, sep="\t", index_col=0)
    dfy = pd.read_csv(y, sep="\t", index_col=0)
    global binary
    binary = is_binary(dfx)
    return (dfx, dfy)


def jaccard(dfx, dfy):
    if binary == False:
        return "Binary profiles only ; use to_Binary() fonction"
    else:
        jaccard_distance = pd.DataFrame(index=list(dfx.index), columns=list(dfy.index))
        for i in dfx.index:
            for j in dfy.index:
                jaccard_distance.loc[i, j] = 1 - jaccard_score(dfx.loc[i], dfy.loc[j])
    print("Jaccard Distance :")
    return jaccard_distance


def hamming(dfx, dfy):
    if binary == False:
        return "Binary profiles only ; use to_Binary() fonction"
    else:
        hamming_distance = pd.DataFrame(index=dfx.index, columns=dfy.index)
        for i in dfx.index:
            for j in dfy.index:
                hamming_distance.loc[i, j] = distance.hamming(dfx.loc[i], dfy.loc[j])
    print("Hamming Distance :")
    return hamming_distance


def pearson(dfx, dfy):
    if binary == False:
        print(
            "You use continuous profiles, think about normalize your datas with normalize() function"
        )
    pearson_results = pd.DataFrame(index=dfx.index, columns=dfy.index)
    for i in dfx.index:
        for j in dfy.index:
            pearson_correlation = stats.pearsonr(dfx.loc[i], dfy.loc[j])
            pearson_results.loc[i, j] = [pearson_correlation[0], pearson_correlation[1]]
    print("Pearson Correlation and associated p-values :")
    return pearson_results


def mi(dfx, dfy):
    if binary == False:
        print(
            "You use continuous profiles, think about normalize your datas with normalize() function"
        )
    mi_distance = pd.DataFrame(index=dfx.index, columns=dfy.index)
    for i in dfx.index:
        for j in dfy.index:
            mi_distance.loc[i, j] = mutual_info_score(dfx.loc[i], dfy.loc[j])
    print("Mutual Information :")
    return mi_distance


def to_binary(path, treshold=60):
    df = pd.read_csv(path, sep="\t", index_col=0)
    global binary
    binary = is_binary(df)
    if binary is False:
        for i in range(0, len(df)):
            for j in range(0, len(df.columns)):
                if df.iloc[i][j] > treshold:
                    df.iloc[i][j] = 1
                else:
                    df.iloc[i][j] = 0
        return df
    else:
        return "Already binary profiles"


# def normalize(df):
#   for i in range(0, len(df)):
#      for j in range(0, len(df.columns)):
#         df.iloc[i][j] = df.iloc[i][j] / max(df.iloc[i])
# return df


def normalize(path):
    df = pd.read_csv(path, sep="\t", index_col=0)
    global binary
    binary = is_binary(df)
    if binary is False:
        for i in df.index:
            df.loc[i] = preprocessing.normalize(df.loc[i])
        return df
    else:
        return "Need continous profiles"


def transition_vector(path):
    pp = pd.read_csv(path, sep="\t", index_col=0)
    global binary
    binary = is_binary(pp)
    if binary is True:
        tv = [0]
        for i in range(1, len(pp.columns)):
            if pp.iloc[0][i] == pp.iloc[0][i - 1]:
                tv.append(0)
            if pp.iloc[0][i] > pp.iloc[0][i - 1]:
                tv.append(1)
            if pp.iloc[0][i] < pp.iloc[0][i - 1]:
                tv.append(-1)
        return tv
    else:
        return "Need binary profiles, you can use to_binary()"


def is_binary(df):
    for x in df.icol[0]:
        if x not in [0, 1]:
            global binary
            binary = False
            return binary
        else:
            binary = True
            return binary


def order_by_tree(x, tree):
    if isinstance(x, str):
        dfx = pd.read_csv(x, sep="\t", index_col=0)
    else:
        dfx = x
    phylo = Phylo.read(tree, "newick")
    leaf = phylo.get_terminals()
    ordered_df = pd.DataFrame()
    for l in leaf:
        ordered_df[str(l)] = dfx[str(l)]
    return ordered_df

import pandas as pd
import numpy as np
from sklearn.metrics import jaccard_score
from scipy.spatial import distance
from scipy import linalg, stats
from sklearn.metrics import mutual_info_score
from sklearn import preprocessing
from newick import read
from Bio import Phylo


def distance_profiles(x, y, method, parameter = False):
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
    if method == "cotransition":
        return cotransition(dfx, dfy, parameter)
    if method == "SVDphy":
        return SVD_phy(A, p)


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


def cotransition(tvx, tvy, parameter = False):
    print("Take care to use transition vectors and not classic profiles")
    if binary == False:
        return "Binary profiles only ; use to_Binary() fonction"
    cotransition_scores = pd.DataFrame(index=tvx.index, columns=tvy.index)
    for i in tvx.index:
        for j in tvy.index:
            t1 = 0
            t2 = 0
            c = 0
            d = 0
            k = 0
            last_transition_x = ""
            last_transition_y = ""
            for x in range(0, len(tvx.columns)):
                if parameter == False or last_transition_x != x -1:
                    if tvx.loc[i][x] != 0:
                        t1 = t1 + 1
                if parameter == False or last_transition_y != x -1:
                    if tvy.loc[j][x] !=0:
                        t2 = t2 + 1
                if (last_transition_x != x - 1 and last_transition_y != x -1) or parameter == False:
                    if tvx.loc[i][x] != 0 and tvy.loc[j][x] != 0:
                        if tvx.loc[i][x] == tvy.loc[j][x]:
                            c = c + 1
                        else:
                            d = d + 1
            k = c -d
            cotransition_scores.loc[i, j] = k / (t1 + t2 - abs(k))
    print("Cotransition score :")
    return cotransition_scores


def SVD_phy(a, p):

    U, S, V = np.linalg.svd(,False) #SVD de la matrice A
    k = int(p * len(U)) #nb de colonnes à garder dans U
    U_truncated = U[:, :k] #ajout des colonnes U

    # normalisation
    for i in range(0, len(U_truncated)):
        for j in range(0, k):
            U_truncated[i][j] = U_truncated[i][j] / np.max(U_truncated[i])

    # Calculate the distance euclidienne
    SVDphy_distance = np.zeros((len(U_truncated), len(U_truncated)))
    for i in range(0, len(U_truncated)):
        for j in range(0, len(U_truncated)):
            SVDphy_distance[i][j] = np.sqrt(np.sum((U_truncated[i] - U_truncated[j])**2))
    
    row_labels = a.index
    col_labels = a.index
    SVDphy_distance = pd.DataFrame(SVDphy_distance, index=row_labels, columns=col_labels)

    return SVDphy_distance


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
        tv = pd.DataFrame(index=pp.index, columns=pp.columns)
        for j in range(0, len(pp)):
            tv.iloc[j][0] = 0
            for i in range(1, len(pp.columns)):
                if pp.iloc[j][i] == pp.iloc[j][i - 1]:
                    tv.iloc[j][i] = 0
                if pp.iloc[j][i] > pp.iloc[j][i - 1]:
                    tv.iloc[j][i] = 1
                if pp.iloc[j][i] < pp.iloc[j][i - 1]:
                    tv.iloc[j][i] = -1
        return tv
    else:
        return "Need binary profiles, you can use to_binary()"


def is_binary(df):
    for x in df.iloc[0]:
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


import pandas as pd
from scipy import stats
from scipy.spatial import distance
from sklearn.metrics import jaccard_score, mutual_info_score


def distance_profiles(x, y, method):
    dfx, dfy = input(x, y)
    if method == "Jaccard":
        return Jaccard(dfx, dfy)
    if method == "Hamming":
        return Hamming(dfx, dfy)
    if method == "Pearson":
        return pearson(dfx, dfy)
    if method == "MI":
        return calc_mi(dfx, dfy)


def input(x, y):
    dfx = pd.read_csv(x, sep="\t", index_col=0)
    dfy = pd.read_csv(y, sep="\t", index_col=0)
    for x in dfx.icol[0]:
        if x not in [0, 1]:
            global binary
            binary = False
            break
        else:
            binary = True
    return (dfx, dfy)


def Jaccard(dfx, dfy):  # noqa: N802
    if binary is False:
        return "Binary profiles only ; use to_Binary() fonction"
    else:
        jaccard_distance = pd.DataFrame(index=list(dfx.index), columns=list(dfy.index))
        for i in dfx.index:
            for j in dfy.index:
                jaccard_distance.loc[i, j] = 1 - jaccard_score(dfx.loc[i], dfy.loc[j])
    print("Jaccard Distance :")
    return jaccard_distance


def Hamming(dfx, dfy):  # noqa: N802
    if binary is False:
        return "Binary profiles only ; use to_Binary() fonction"
    else:
        hamming_distance = pd.DataFrame(index=dfx.index, columns=dfy.index)
        for i in dfx.index:
            for j in dfy.index:
                hamming_distance.loc[i, j] = distance.hamming(dfx.loc[i], dfy.loc[j])
    print("Hamming Distance :")
    return hamming_distance


def pearson(dfx, dfy):
    if binary is False:
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


def calc_mi(dfx, dfy):
    if binary is False:
        print(
            "You use continuous profiles, think about normalize your datas with normalize() function"
        )
    mi_distance = pd.DataFrame(index=dfx.index, columns=dfy.index)
    for i in dfx.index:
        for j in dfy.index:
            mi_distance.loc[i, j] = mutual_info_score(dfx.loc[i], dfy.loc[j])
    print("Mutual Information :")
    return mi_distance


def to_binary(df, treshold=60):
    for i in range(0, len(df)):
        for j in range(0, len(df.columns)):
            if df.iloc[i][j] > treshold:
                df.iloc[i][j] = 1
            else:
                df.iloc[i][j] = 0
    return df


def normalize(df):
    for i in range(0, len(df)):
        for j in range(0, len(df.columns)):
            df.iloc[i][j] = df.iloc[i][j] / max(df.iloc[i])
    return df

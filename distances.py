import pandas as pd
import numpy as np
from sklearn.metrics import jaccard_score
from scipy.spatial import distance
from scipy import stats
from sklearn.metrics import mutual_info_score
from sklearn import preprocessing
from newick import read
from Bio import Phylo


def distance_profiles(
    method,                         
    x,                              #profils phylogénétiques
    y=None,                         #2eme matrice de profils, facultative, si on veut comparer spécifiquement un aux autres par exemple
    successive_transitions=True,    #option pour le cotransition score, True = normal, False = ne prend pas en compte les transitions successives (dans un profil:  0010 ne compte que pour 1 evenement de transition par exemple)
    confidence=1.5,                 #option pour le pcs, poids positif pour les transitions doubles (0011 dans un profil par exemple)
    penalty=0.6,                    #option pour le pcs, penalité pour les transitions non communes
    truncation=0.6,                 #option pour svdphy, proportion de la matrice u conservée après la troncature
):
    if isinstance(x, str):
        dfx = input(x)
    else:
        dfx = x
    if isinstance(y, str):
        dfy = input(y)
    else:
        dfy = y
    if method == "Jaccard" or method == "jaccard":
        return jaccard(dfx, dfy)
    if method == "Hamming" or method == "hamming":
        return hamming(dfx, dfy)
    if method == "Pearson" or method == "pearson":
        return pearson(dfx, dfy)
    if method == "MI" or method == "mi":
        return mi(dfx, dfy)
    if method == "cotransition" or method == "Cotransition":
        return cotransition(dfx, dfy, successive_transitions)
    if method == "pcs" or method == "PCS":
        return pcs(dfx, dfy, confidence, penalty)
    if method == "svd_phy" or method == "SVD_phy":
        return svd_phy(dfx, truncation)


def input(x):
    dfx = pd.read_csv(x, sep="\t", index_col=0)
    global binary
    binary = is_binary(dfx)
    return dfx


def jaccard(dfx, dfy):
    if binary == False:
        return "Binary profiles only ; use to_Binary() fonction"
    else:
        if dfy is not None:             #si 2 matrices en entrée
            jaccard_distance = pd.DataFrame(
                index=list(dfx.index), columns=list(dfy.index)
            )
            for i in dfx.index:
                for j in dfy.index:
                    jaccard_distance.loc[i, j] = 1 - jaccard_score(
                        dfx.loc[i], dfy.loc[j]
                    )
        else:                           #si une seule matrice en entrée
            jaccard_distance = pd.DataFrame(
                index=list(dfx.index), columns=list(dfx.index)
            )
            for i in dfx.index:
                for j in dfx.index:
                    jaccard_distance.loc[i, j] = 1 - jaccard_score(
                        dfx.loc[i], dfx.loc[j]
                    )
    print("Jaccard Distance :")
    return jaccard_distance


def hamming(dfx, dfy):
    if binary == False:
        return "Binary profiles only ; use to_Binary() fonction"
    else:
        if dfy is not None:             #si 2 matrices en entrée
            hamming_distance = pd.DataFrame(index=dfx.index, columns=dfy.index)
            for i in dfx.index:
                for j in dfy.index:
                    hamming_distance.loc[i, j] = distance.hamming(
                        dfx.loc[i], dfy.loc[j]
                    )
        else:                           #si une seule matrice en entrée
            hamming_distance = pd.DataFrame(index=dfx.index, columns=dfx.index)
            for i in dfx.index:
                for j in dfx.index:
                    hamming_distance.loc[i, j] = distance.hamming(
                        dfx.loc[i], dfx.loc[j]
                    )
    print("Hamming Distance :")
    return hamming_distance


def pearson(dfx, dfy):
    if binary == False:
        print(
            "You use continuous profiles, think about normalize your datas with normalize() function"
        )
    if dfy is not None:             #si 2 matrices en entrée
        pearson_results = pd.DataFrame(index=dfx.index, columns=dfy.index)
        for i in dfx.index:
            for j in dfy.index:
                pearson_correlation = stats.pearsonr(dfx.loc[i], dfy.loc[j])
                pearson_results.loc[i, j] = [
                    pearson_correlation[0],         #valeur de la correlation
                    pearson_correlation[1],         #p-valeur associée
                ]
    else:                          #si une seule matrice en entrée
        pearson_results = pd.DataFrame(index=dfx.index, columns=dfx.index)
        for i in dfx.index:
            for j in dfx.index:
                pearson_correlation = stats.pearsonr(dfx.loc[i], dfx.loc[j])
                pearson_results.loc[i, j] = [
                    pearson_correlation[0],         #valeur de la correlation
                    pearson_correlation[1],         #p-valeur associée
                ]
    print("Pearson Correlation and associated p-values :")
    return pearson_results


def mi(dfx, dfy):
    if binary == False:
        print(
            "You use continuous profiles, think about normalize your datas with normalize() function"
        )
    if dfy is not None:         #si 2 matrices en entrée
        mi_distance = pd.DataFrame(index=dfx.index, columns=dfy.index)
        for i in dfx.index:
            for j in dfy.index:
                mi_distance.loc[i, j] = mutual_info_score(dfx.loc[i], dfy.loc[j])
    else:                       #si une seule matrice en entrée
        mi_distance = pd.DataFrame(index=dfx.index, columns=dfx.index)
        for i in dfx.index:
            for j in dfx.index:
                mi_distance.loc[i, j] = mutual_info_score(dfx.loc[i], dfx.loc[j])
    print("Mutual Information :")
    return mi_distance


def cotransition(tvx, tvy, successive_transitions=True):
    print("Take care to use transition vectors and not classic profiles, ordered by a tree")
    if binary == False:
        return "Binary profiles only ; use to_Binary() fonction"
    if tvy is not None:                 #si 2 matrices en entrée
        cotransition_scores = pd.DataFrame(index=tvx.index, columns=tvy.index)
        for i in tvx.index:
            for j in tvy.index:
                t1 = 0              #nombre de transitions dans un profil
                t2 = 0              #nombre de transitions dans l'autre
                c = 0               #nombre de transitions communes allant dans le même sens 
                d = 0               #nombre de transitions communes n'allant pas dans le même sens
                k = 0               # k = c -d
                last_transition_x = ""   #index de la dernière transition détectée dans un profil
                last_transition_y = ""   #index de la dernière transition détectée dans l'autre
                for x in range(0, len(tvx.columns)):
                    if successive_transitions == True or last_transition_x != x - 1:           #successive_transitions = False alors on regarde si on a pas déjà compté une transition juste avant
                        if tvx.loc[i][x] != 0:
                            last_transition_x = x
                            t1 = t1 + 1
                    if successive_transitions == True or last_transition_y != x - 1:           #successive_transitions = False alors on regarde si on a pas déjà compté une transition juste avant
                        if tvy.loc[j][x] != 0:
                            last_transition_y = x
                            t2 = t2 + 1
                    if (
                        last_transition_x != x - 1 and last_transition_y != x - 1              
                    ) or successive_transitions == True:
                        if tvx.loc[i][x] != 0 and tvy.loc[j][x] != 0:       #transitions communes (au niveau de la même espèce dans le profil)
                            if tvx.loc[i][x] == tvy.loc[j][x]:              #transitions dans le même sens
                                c = c + 1
                            else:                                           #transitions dans deux sens différents
                                d = d + 1
                k = c - d
                if t1 == 0 and t2 == 0:                                     #si on a aucune transitions dans tout le profil, on pourra pas faire la divison ensuite, on retourne None pour la distance
                    cotransition_scores.loc[i, j] = None
                else:
                    cotransition_scores.loc[i, j] = k / (t1 + t2 - abs(k))
    else:                       #si une seule matrice en entrée
        cotransition_scores = pd.DataFrame(index=tvx.index, columns=tvx.index)
        for i in tvx.index:
            for j in tvx.index:
                t1 = 0
                t2 = 0
                c = 0
                d = 0
                k = 0
                last_transition_x = ""
                last_transition_y = ""
                for x in range(0, len(tvx.columns)):
                    if successive_transitions == True or last_transition_x != x - 1:
                        if tvx.loc[i][x] != 0:
                            last_transition_x = x
                            t1 = t1 + 1
                    if successive_transitions == True or last_transition_y != x - 1:
                        if tvx.loc[j][x] != 0:
                            last_transition_y = x
                            t2 = t2 + 1
                    if (
                        last_transition_x != x - 1 and last_transition_y != x - 1
                    ) or successive_transitions == True:
                        if tvx.loc[i][x] != 0 and tvx.loc[j][x] != 0:
                            if tvx.loc[i][x] == tvx.loc[j][x]:
                                c = c + 1
                            else:
                                d = d + 1
                k = c - d
                print(i, j, x, t1, t2, c, d, k)
                if t1 == 0 and t2 == 0:
                    cotransition_scores.loc[i, j] = None
                else:
                    cotransition_scores.loc[i, j] = k / (t1 + t2 - abs(k))
    print("Cotransition score :")
    return cotransition_scores


def pcs(tvx, tvy, confidence=1.5, penalty=0.6):
    print("Take care to use transition vectors and not classic profiles, ordered by a tree")
    if binary == False:
        return "Binary profiles only ; use to_Binary() fonction"
    if tvy is not None:                 #si 2 matrices en entrée
        pcs_scores = pd.DataFrame(index=tvx.index, columns=tvy.index)
        for i in tvx.index:
            for j in tvy.index:
                match_1 = 0             #nombre de transitions simples (01 ou 10) identiques 
                mismatch_1 = 0          #nombre de transitions simples (01 ou 10) non partagées
                match_2 = 0             #nombre de transitions doubles (0011 ou 1100) identiques 
                mismatch_2 = 0          #nombre de transitions doubles (0011 ou 1100) non partagées
                for x in range(1, len(tvx.columns)):
                    if tvx.loc[i][x] != 0 or tvy.loc[j][x] != 0:        #si on a une transition sur au moins un des deux profils
                        if tvx.loc[i][x] == tvy.loc[j][x]:              #si les transitions sont égales et partagées
                            if len(tvx.columns) - 1 - x > 0:            #si on est pas au dernier élément du profil (car à ce moment la on ne pourrait plus faire index+1)
                                if (                                    #vérifie si on est dans une situation de transitions doubles pour les deux profils
                                    tvx.loc[i][x - 1]
                                    == tvx.loc[i][x + 1]
                                    == tvy.loc[j][x - 1]
                                    == tvy.loc[j][x + 1]
                                    == 0
                                ):
                                    match_2 = match_2 + 1
                                else:                                   #sinon c'est juste un match simple
                                    match_1 = match_1 + 1
                            else:
                                match_1 = match_1 + 1
                        elif len(tvx.columns) - 1 - x > 0:              #si les transitions ne sont pas partagées ou égales et qu'on est pas au dernier élément
                            if (                                        #si on avait une transition double sur un des deux profils
                                tvx.loc[i][x] != 0
                                and tvx.loc[i][x - 1] == tvx.loc[i][x + 1] == 0
                            ) or (
                                tvy.loc[j][x] != 0
                                and tvy.loc[j][x - 1] == tvy.loc[j][x + 1] == 0
                            ):
                                mismatch_2 = mismatch_2 + 1
                            else:                                       #sinon c'est juste un mismatch simple
                                mismatch_1 = mismatch_1 + 1
                        else:
                            mismatch_1 = mismatch_1 + 1
                pcs_scores.loc[i, j] = (
                    (match_1)
                    + (match_2 * confidence)
                    - penalty * (mismatch_1 + mismatch_2 * confidence)
                )
    else:                                       #si une seule matrice en entrée
        pcs_scores = pd.DataFrame(index=tvx.index, columns=tvx.index)
        for i in tvx.index:
            for j in tvx.index:
                match_1 = 0
                mismatch_1 = 0
                match_2 = 0
                mismatch_2 = 0
                for x in range(1, len(tvx.columns)):
                    if tvx.loc[i][x] != 0 or tvx.loc[j][x] != 0:
                        if tvx.loc[i][x] == tvx.loc[j][x]:
                            if len(tvx.columns) - 1 - x > 0:
                                if (
                                    tvx.loc[i][x - 1]
                                    == tvx.loc[i][x + 1]
                                    == tvx.loc[j][x - 1]
                                    == tvx.loc[j][x + 1]
                                    == 0
                                ):
                                    match_2 = match_2 + 1
                                else:
                                    match_1 = match_1 + 1
                            else:
                                match_1 = match_1 + 1
                        elif len(tvx.columns) - 1 - x > 0:
                            if (
                                tvx.loc[i][x] != 0
                                and tvx.loc[i][x - 1] == tvx.loc[i][x + 1] == 0
                            ) or (
                                tvx.loc[j][x] != 0
                                and tvx.loc[j][x - 1] == tvx.loc[j][x + 1] == 0
                            ):
                                mismatch_2 = mismatch_2 + 1
                            else:
                                mismatch_1 = mismatch_1 + 1
                        else:
                            mismatch_1 = mismatch_1 + 1
                pcs_scores.loc[i, j] = (
                    (match_1)
                    + (match_2 * confidence)
                    - penalty * (mismatch_1 + mismatch_2 * confidence)
                )
    print("PCS :")
    return pcs_scores


def svd_phy(a, p):
    u, s, v = np.linalg.svd(a, False)  # SVD de la matrice A
    k = int(p * len(u))  # nb de colonnes à garder dans U
    u_truncated = u[:, :k]  # ajout des colonnes U

    # normalisation
    for i in range(0, len(u_truncated)):
        for j in range(0, k):
            u_truncated[i][j] = u_truncated[i][j] / np.max(u_truncated[i])

    # Calculate the distance euclidienne
    svdphy_distance = np.zeros((len(u_truncated), len(u_truncated)))
    for i in range(0, len(u_truncated)):
        for j in range(0, len(u_truncated)):
            svdphy_distance[i][j] = np.sqrt(
                np.sum((u_truncated[i] - u_truncated[j]) ** 2)
            )

    row_labels = a.index
    col_labels = a.index
    svdphy_distance = pd.DataFrame(
        svdphy_distance, index=row_labels, columns=col_labels
    )

    return svdphy_distance


def to_binary(path, treshold=60):                       #transforme les profils continues en des profils binaires
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


def normalize(path):                                        #normalise des profils continues
    df = pd.read_csv(path, sep="\t", index_col=0)
    global binary
    binary = is_binary(df)
    if binary is False:
        for i in df.index:
            df.loc[i] = preprocessing.normalize(df.loc[i])
        return df
    else:
        return "Need continous profiles"


def transition_vector(path):                                #génére des vecteurs de transitions à partir de profils binaires (001101 --> 0010-11)
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


def order_by_tree(x, tree):                                #ordone un profile en fonction des feuilles d'un arbre sous format newick (les id pour les espèces dans le profils doivent être les même que dans l'arbre)
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

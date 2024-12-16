import pandas as pd
import numpy as np
from sklearn.metrics import jaccard_score
from scipy.spatial import distance
from scipy import stats
from sklearn.metrics import mutual_info_score
from sklearn import preprocessing
from Bio import Phylo
from sqlite3 import Binary
from scipy import stats
from scipy.stats import fisher_exact
import time


def distance_profiles(
    method,
    x,
    y=None,
    subset_prot=None,
    successive_transitions=True,
    confidence=1.5,
    penalty=0.6,
    truncation=None,
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
        return SVD_phy(dfx, truncation, subset_prot)


def input(x):
    dfx = pd.read_csv(x, sep=",", index_col=0)
    global binary
    binary = is_binary(dfx)
    return dfx


def jaccard(dfx, dfy):
    if binary == False:
        return "Binary profiles only ; use to_Binary() fonction"
    else:
        symetry = False
        if dfy is None:
            dfy = dfx
            symetry=True
        jaccard_distance = np.zeros((len(dfx.index), len(dfy.index)))
        for a,i in enumerate(dfx.index): # Éviter les doublons si symétrique
            query = dfx.loc[i]
            for b,j in enumerate(dfy.index):
                if symetry and b<a:
                    continue
                #score_temp = 1 - jaccard_score(query, dfy.loc[j])
                score_temp = distance.jaccard(query, dfy.loc[j])
                jaccard_distance[a, b] = score_temp
                if symetry:
                    jaccard_distance[b, a] = score_temp
    jaccard_distance = pd.DataFrame(jaccard_distance, index=dfx.index, columns=dfy.index)
    jaccard_distance.to_csv('jaccard_distance.csv', index=True)
    return jaccard_distance


def hamming(dfx, dfy):
    if True == False:
        return "Binary profiles only ; use to_Binary() fonction"
    else:
        symetry = False
        if dfy is None:
            dfy = dfx
            symetry=True
        hamming_distance = np.zeros((len(dfx.index), len(dfy.index)))
        for a,i in enumerate(dfx.index): # Éviter les doublons si symétrique
            query = dfx.loc[i]
            for b,j in enumerate(dfy.index):
                if symetry and b<a:
                    continue
                score_temp = distance.hamming(query, dfy.loc[j])
                hamming_distance[a,b] = score_temp
                if symetry: 
                    hamming_distance[b,a] = score_temp
    hamming_distance = pd.DataFrame(hamming_distance, index=dfx.index, columns=dfy.index)
    hamming_distance.to_csv('hamming_distance.csv', index=True)
    return hamming_distance


def pearson(dfx, dfy):
    if binary == False:
        print(
            "You use continuous profiles, think about normalize your datas with normalize() function"
        )
    symetry = False
    if dfy is None:
        dfy = dfx
        symetry=True
    pearson_results = np.zeros((len(dfx.index), len(dfy.index)))
    for a,i in enumerate(dfx.index): # Éviter les doublons si symétrique
        query = dfx.loc[i]
        for b,j in enumerate(dfy.index):
            if symetry and b<a:
                continue
            pearson_correlation = stats.pearsonr(query, dfy.loc[j])
            if not np.isnan(pearson_correlation[0]):
                score_temp = 1 - pearson_correlation[0]
                pearson_results[a, b] = score_temp
                if symetry:
                    pearson_results.loc[b, a] = score_temp
            else:
                pearson_results[a, b] = np.nan
                if symetry:
                    pearson_results.loc[b, a] = np.nan
    pearson_results = pd.DataFrame(pearson_results, index=dfx.index, columns=dfy.index)
    pearson_results.to_csv('pearson_correlation.csv', index=True)
    return pearson_results


def mi(dfx, dfy):
    if binary == False:
        print(
            "You use continuous profiles, think about normalize your datas with normalize() function"
        )
    symetry = False
    if dfy is None:
        dfy = dfx
        symetry=True
    mi_distance = np.zeros((len(dfx.index), len(dfy.index)))
    for a,i in enumerate(dfx.index):
        #query = dfx.loc[i]
        for b,j in enumerate(dfy.index):
            if symetry and b<a:
                continue
            score_temp = 1 - mutual_info_score(dfx.loc[i], dfy.loc[j])
            mi_distance[a, b] = score_temp
            if symetry:
                mi_distance[b, a] = score_temp
    mi_distance = pd.DataFrame(mi_distance, index=dfx.index, columns=dfy.index)
    mi_distance.to_csv('mi_distance.csv', index=True)
    return mi_distance


def cotransition(dfx, dfy, successive_transitions=True):
    #if binary == False:
    #    return "Binary profiles only ; use to_Binary() fonction"
    tvx = transition_vector(dfx)
    print(len(tvx))
    symetry = False
    if dfy is None:
        tvy = tvx
        dfy = dfx
        symetry=True
    else:
        tvy = transition_vector(dfy)
    cotransition_scores = np.zeros((len(dfx.index), len(dfy.index)))
    p_values = np.zeros((len(dfx.index), len(dfy.index)))
    print(len(p_values))
    for a,i in enumerate(tvx.index):
        #start_row = time.time()
        for b,j in enumerate(tvy.index):
            if symetry and b<a:
                continue
            t1 = 0
            t2 = 0
            c = 0
            d = 0
            k = 0
            last_transition_x = ""
            last_transition_y = ""
            v1 = tvx.loc[i]
            v2 = tvy.loc[j]
            t1 = np.count_nonzero(v1)
            t2 = np.count_nonzero(v2)
            nonz = (v1 != 0) & (v2 != 0)
            c = np.sum(v1[nonz]==v2[nonz])
            d = np.count_nonzero(v1[nonz]-v2[nonz])
            #for x in range(0, len(tvx.columns)):
                #if successive_transitions == True or last_transition_x != x - 1:
                #    if tvx.loc[i][x] != 0:
                #            t1 = t1 + 1
                #if successive_transitions == True or last_transition_y != x - 1:
                #    if tvy.loc[j][x] != 0:
                #            t2 = t2 + 1
                #if (
                #    last_transition_x != x - 1 and last_transition_y != x - 1
                #) or successive_transitions == True:
                #    if tvx.loc[i][x] != 0 and tvy.loc[j][x] != 0:
                #        if tvx.loc[i][x] == tvy.loc[j][x]:
                #            c = c + 1
                #        else:
                #            d = d + 1
            k = c - d
            #print(i, j, len(tvx.columns), t1, t2, c, d, k)
            if t1 == 0 and t2 == 0:
                cotransition_scores[a, b] = None
                if symetry:
                    cotransition_scores[b, a] = None
            else:
                score_temp = k / (t1 + t2 - abs(k))
                cotransition_scores[a, b] = score_temp
                if symetry:
                    cotransition_scores[b, a] = score_temp
            #print(time.time() - start_row)

                #tableau de contingence:
            contingency_table = [[abs(k),t1-abs(k)], [t2-abs(k),(len(tvx.columns))-t1-t2+abs(k)]]
            score = fisher_exact(contingency_table, alternative="greater")
            p_values[a, b] = score.pvalue
            if symetry:
                p_values[b, a] = score.pvalue
    print("Cotransition score :")
    cotransition_scores = pd.DataFrame(cotransition_scores, index=dfx.index, columns=dfy.index)
    p_values = pd.DataFrame(p_values, index=dfx.index, columns=dfy.index)
    cotransition_scores.to_csv('cotransition_scores.csv', index=True)
    p_values.to_csv('p_values.csv', index=True)
    return cotransition_scores, p_values



def pcs(dfx, dfy, confidence=1.5, penalty=0.6):
    print("Take care to use transition vectors and not classic profiles")
    if binary == False:
        return "Binary profiles only ; use to_Binary() fonction"
    tvx = transition_vector(dfx)
    symetry = False
    if dfy is None:
        tvy = tvx
        symetry=True
    else:
        tvy = transition_vector(dfy)
    pcs_scores = np.zeros((len(dfx.index), len(dfy.index)))
    for a,i in enumerate(tvx.index):
        start = time.time()
        for b,j in enumerate(tvy.index):
            if symetry and b<a:
                continue
            match_1 = 0
            mismatch_1 = 0
            match_2 = 0
            mismatch_2 = 0
            for x in range(1, len(tvx.columns)):
                if tvx.loc[i][x] != 0 or tvy.loc[j][x] != 0:
                    if tvx.loc[i][x] == tvy.loc[j][x]:
                        if len(tvx.columns) - 1 - x > 0:
                            if (
                                tvx.loc[i][x - 1]
                                == tvx.loc[i][x + 1]
                                == tvy.loc[j][x - 1]
                                == tvy.loc[j][x + 1]
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
                            tvy.loc[j][x] != 0
                            and tvy.loc[j][x - 1] == tvy.loc[j][x + 1] == 0
                        ):
                            mismatch_2 = mismatch_2 + 1
                        else:
                            mismatch_1 = mismatch_1 + 1
                    else:
                        mismatch_1 = mismatch_1 + 1
            score_temp = (
                (match_1) + (match_2 * confidence) - penalty * (mismatch_1 + mismatch_2 * confidence)
            )
            pcs_scores[a, b] = score_temp
            if symetry:
                pcs_scores[b, a] = score_temp
            if j%100==0:
                print(time.time()-start)
    pcs_scores = pd.DataFrame(pcs_scores, index=dfx.index, columns=dfy.index)
    pcs_scores.to_csv('pcs_scores.csv', index=True)
    return pcs_scores


def SVD_phy(A,  p, subset_prot=None): 
    start = time.time()
    U, S, V = np.linalg.svd(A,False) #SVD de la matrice A
    k = int(p * len(U)) #nb de colonnes à garder dans U
    U_truncated = U[:, :k] #ajout des colonnes U
    
    print(time.time()-start)
    print(U.shape)
    print(U_truncated.shape)


    # normalisation
    #for i in range(0, len(U_truncated)):
        #for j in range(0, k):
        #U_truncated[i][j] = U_truncated[i][j] / np.max(U_truncated[i])
    # Calculate the distance euclidienne
    if subset_prot!=None:
        index_U = A.index
        U = pd.DataFrame(U_truncated, index=index_U)
        subset_U = U.loc[subset_prot]
        row_labels  = subset_U.index
        subset_U = subset_U.to_numpy()
        print(subset_U.shape)
    else:
        subset_U = U_truncated
        row_labels = A.index
    col_labels = A.index
    
    SVDphy_distance = np.zeros((len(subset_U), len(U_truncated)))

    for i in range(0, len(subset_U)):
        start = time.time()
        for j in range(0, len(U_truncated)):
            SVDphy_distance[i][j] = np.sqrt(np.sum((subset_U[i] - U_truncated[j])**2))
        if i%100==0:
            print(time.time()-start)

    SVDphy_distance = pd.DataFrame(SVDphy_distance, index=row_labels, columns=col_labels)
    SVDphy_distance.to_csv('SVDphy_distance.csv', index=True)
    return SVDphy_distance

   


def to_binary(dfx, treshold=60):
    global binary
    binary = is_binary(dfx)
    if binary is False:
        for i in range(0, len(dfx)):
            for j in range(0, len(dfx.columns)):
                if dfx.iloc[i][j] > treshold:
                    dfx.iloc[i][j] = 1
                else:
                    dfx.iloc[i][j] = 0
        return dfx
    else:
        return print("Already binary profiles")


def normalize(dfx):
    global binary
    binary = is_binary(dfx)
    if binary is False:
        for i in dfx.index:
            dfx.loc[i] = preprocessing.normalize(dfx.loc[i])
        return dfx
    else:
        return dfx


def transition_vector(dfx):
    global binary
    binary = is_binary(dfx)
    if binary is True:
        tv = np.zeros((len(dfx.index),(len(dfx.columns))))
        start = time.time()
        for j in range(0, len(dfx)):
            vec = dfx.iloc[j]
            pos_ori = np.asarray(vec[:-1])
            pos_new = np.asarray(vec[1:])
            trans_vec = pos_new-pos_ori
            tv[j] = np.pad(trans_vec, (1,0))
            #tv[j][0] = 0
            #for i in range(1, len(dfx.columns)):
            #    curr = dfx.iloc[j][i]
            #    prev = dfx.iloc[j][i-1]
            #    if curr == prev:
            #        tv[j][i] = 0
            #    elif curr > prev:
            #        tv[j][i] = 1
            #    elif curr < prev:
            #        tv[j][i] = -1
            if j%100==0:
                print(time.time()-start)
        tv = pd.DataFrame(tv, index=dfx.index, columns=dfx.columns)
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


def order_by_tree(dfx, tree):
    phylo = Phylo.read(tree, "newick")
    print(tree)
    leaf = phylo.get_terminals()
    print(leaf)
    ordered_df = pd.DataFrame()
    for l in leaf:
        ordered_df[str(l.name)] = dfx[str(l.name)]
    return ordered_df

if __name__ == "__main__":
    dfx = input("strigamia-acuminata_profiles.tsv.short")
    Jresuts = distance_profiles("jaccard", dfx)
    print(Jresuts)

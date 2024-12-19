import pandas as pd
import numpy as np
from scipy.spatial import distance
from scipy import stats
from sklearn.metrics import mutual_info_score
from scipy.stats import fisher_exact
from PhyloDist.pre_process import is_binary, transition_vector


def jaccard(dfx, dfy):
    binary_x = is_binary(dfx)
    symetry = False
    if dfy is None:
        dfy = dfx
        symetry=True
        binary_y = binary_x
    else:
        binary_y = is_binary(dfy)
    if not binary_x or not binary_y:
        print("Binary profiles only ; use to_Binary() fonction")
    
    jaccard_distance = np.zeros((len(dfx.index), len(dfy.index)))
    for a,i in enumerate(dfx.index): # Éviter les doublons si symétrique
        query = dfx.loc[i]
        for b,j in enumerate(dfy.index):
            if symetry and b<a:
                continue
            score_temp = distance.jaccard(query, dfy.loc[j])
            jaccard_distance[a, b] = score_temp
            if symetry:
                jaccard_distance[b, a] = score_temp
    index=dfx.index 
    columns=dfy.index
    jaccard_distance = pd.DataFrame(jaccard_distance, index=index, columns=columns)
    return jaccard_distance


def hamming(dfx, dfy):
    binary_x = is_binary(dfx)
    symetry = False
    if dfy is None:
        dfy = dfx
        symetry=True
        binary_y = binary_x
    else:
        binary_y = is_binary(dfy)
    if not binary_x or not binary_y:
        print("Binary profiles only ; use to_Binary() fonction")
    
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
    index=dfx.index 
    columns=dfy.index
    hamming_distance = pd.DataFrame(hamming_distance, index=index, columns=columns)
    return hamming_distance


def pearson(dfx, dfy):
    binary_x = is_binary(dfx)
    symetry = False
    if dfy is None:
        dfy = dfx
        symetry=True
        binary_y = binary_x
    else:
        binary_y = is_binary(dfy)
    if not binary_x or not binary_y:
        print("You use continuous profiles, think about normalize your datas with normalize() function")
    
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
                    pearson_results[b, a] = score_temp
            else:
                pearson_results[a, b] = np.nan
                if symetry:
                    pearson_results[b, a] = np.nan
    index=dfx.index 
    columns=dfy.index
    pearson_results = pd.DataFrame(pearson_results, index=index, columns=columns)
    return pearson_results


def mi(dfx, dfy):
    binary_x = is_binary(dfx)
    symetry = False
    if dfy is None:
        dfy = dfx
        symetry=True
        binary_y = binary_x
    else:
        binary_y = is_binary(dfy)
    if not binary_x or not binary_y:
        print("You use continuous profiles, think about normalize your datas with normalize() function")

    mi_distance = np.zeros((len(dfx.index), len(dfy.index)))
    for a,i in enumerate(dfx.index):
        query = dfx.loc[i]
        for b,j in enumerate(dfy.index):
            if symetry and b<a:
                continue
            score_temp = 1 - mutual_info_score(query, dfy.loc[j])
            mi_distance[a, b] = score_temp
            if symetry:
                mi_distance[b, a] = score_temp
    index=dfx.index 
    columns=dfy.index
    mi_distance = pd.DataFrame(mi_distance, index=index, columns=columns)
    return mi_distance


def cotransition(dfx, dfy = None, successive_transitions=True):
    binary_x = is_binary(dfx)
    tvx = transition_vector(dfx, binary_x)
    symetry = False
    if dfy is None:
        symetry=True
        dfy = dfx
        binary_y = binary_x
        tvy = tvx
    else:
        binary_y = is_binary(dfy)
        tvy = transition_vector(dfy, binary_y)
    if not binary_x or not binary_y:
        print("Binary profiles only ; use to_Binary() fonction")
   
    cotransition_scores = np.zeros((len(dfx.index), len(dfy.index)))
    p_values = np.zeros((len(dfx.index), len(dfy.index)))
    for a,i in enumerate(tvx.index):
        for b,j in enumerate(tvy.index):
            if symetry and b<a:
                continue
            t1 = 0
            t2 = 0
            c = 0
            d = 0
            k = 0
            v1 = tvx.loc[i]
            v2 = tvy.loc[j]
            t1 = np.count_nonzero(v1)
            t2 = np.count_nonzero(v2)
            nonz = (v1 != 0) & (v2 != 0)
            c = np.sum(v1[nonz]==v2[nonz])
            d = np.count_nonzero(v1[nonz]-v2[nonz])
            k = c - d
            if t1 == 0 and t2 == 0:
                cotransition_scores[a, b] = None
                if symetry:
                    cotransition_scores[b, a] = None
            else:
                score_temp = k / (t1 + t2 - abs(k))
                cotransition_scores[a, b] = score_temp
                if symetry:
                    cotransition_scores[b, a] = score_temp

                #tableau de contingence:
            contingency_table = [[abs(k),t1-abs(k)], [t2-abs(k),(len(tvx.columns))-t1-t2+abs(k)]]
            score = fisher_exact(contingency_table, alternative="greater")
            p_values[a, b] = score.pvalue
            if symetry:
                p_values[b, a] = score.pvalue
    cotransition_scores = pd.DataFrame(cotransition_scores, index=dfx.index, columns=dfy.index)
    p_values = pd.DataFrame(p_values, index=dfx.index, columns=dfy.index)
    return cotransition_scores, p_values



def pcs(dfx, dfy, confidence=1.5, penalty=0.6):
    binary_x = is_binary(dfx)
    tvx = transition_vector(dfx, binary_x)
    symetry = False
    if dfy is None:
        symetry=True
        dfy = dfx
        binary_y = binary_x
        tvy = tvx
    
    else:
        binary_y = is_binary(dfy)
        tvy = transition_vector(dfy, binary_y)
    if not binary_x or not binary_y:
        print("Binary profiles only ; use to_Binary() fonction")
   
    pcs_scores = np.zeros((len(dfx.index), len(dfy.index)))
    for a,i in enumerate(tvx.index):
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
    pcs_scores = pd.DataFrame(pcs_scores, index=dfx.index, columns=dfy.index)
    return pcs_scores


def SVD_phy(dfx,  p, subset_prot=None): 
    U, S, V = np.linalg.svd(dfx,False) #SVD de la matrice dfx
    k = int(p * len(U)) #nb de colonnes à garder dans U
    U_truncated = U[:, :k] #ajout des colonnes U
    
    if subset_prot is not None:
        index_U = dfx.index
        U = pd.DataFrame(U_truncated, index=index_U)
        subset_U = U.loc[subset_prot]
        row_labels  = subset_U.index
        subset_U = subset_U.to_numpy()
        print(subset_U.shape)
    else:
        subset_U = U_truncated
        row_labels = dfx.index
    col_labels = dfx.index
    
    SVDphy_distance = np.zeros((len(subset_U), len(U_truncated)))
    print(SVDphy_distance.shape)

    for i in range(0, len(subset_U)):
        for j in range(0, len(U_truncated)):
            SVDphy_distance[i][j] = np.sqrt(np.sum((subset_U[i] - U_truncated[j])**2))

    SVDphy_distance = pd.DataFrame(SVDphy_distance, index=row_labels, columns=col_labels)
    return SVDphy_distance


from PhyloDist.distances import SVD_phy, cotransition, jaccard, hamming, mi, pcs, pearson
from PhyloDist.output import to_csv
import PhyloDist.pre_process as pre_process

def distance_profiles(
    method,
    x,
    y=None,
    subset_prot=None,
    successive_transitions=True,
    confidence=1.5,
    penalty=0.6,
    truncation=None,
    csv=False
):
    
    dfx = pre_process.input(x)
    if y is None:
        dfy = None
    else:
        dfy = pre_process.input(y)
    p_value = None
    if method == "Jaccard" or method == "jaccard":
        result = jaccard(dfx, dfy)
    if method == "Hamming" or method == "hamming":
        result = hamming(dfx, dfy) 
    if method == "Pearson" or method == "pearson":
        result = pearson(dfx, dfy)
    if method == "MI" or method == "mi":
        result = mi(dfx, dfy)
    if method == "cotransition" or method == "Cotransition":
        result, p_value = cotransition(dfx, dfy, successive_transitions)
    if method == "pcs" or method == "PCS":
        result = pcs(dfx, dfy, confidence, penalty)
    if method == "svd_phy" or method == "SVD_phy":
        result = SVD_phy(dfx, truncation, subset_prot)
    if csv is True:
        if p_value != None:
            return to_csv(result, name=f"{method}.csv"), to_csv(p_value, name=f"{method}_pvalue.csv")
        else:
            return to_csv(result, name=f"{method}.csv")
    else:
        return result
from PhyloDist.distances import jaccard, hamming, mi, pearson
from PhyloDist.output import to_csv
from PhyloDist.pre_process import input, is_binary, binary

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
        print(binary)
    else:
        dfx = x
        binary = True
    if isinstance(y, str):
        dfy = input(dfy)
    else:
        dfy = y
        binary = True
    if method == "Jaccard" or method == "jaccard":
        Jaccard, index, columns = jaccard(dfx, dfy)
        return to_csv(Jaccard, index, columns)
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
import pandas as pd
import numpy as np
from scipy import linalg

dfx = pd.read_csv("strigamia-acuminata_profiles.tsv.short", sep="\t", index_col=0)
#dfy = pd.read_csv("strigamia-acuminata_seq2.tsv.short", sep="\t", index_col=0)

"""
def SVD_phy(A, p):

    U, S, Vh = linalg.svd(A) #SVD de la matrice A

    k = int(p * len(U)) #nb de colonnes à garder dans U

    dfx_truncated = U[:, :k] #ajout des colonnes U

    return dfx_truncated


result = SVD_phy(dfx, 0.4)

print("SVP-phy", result)
"""
U, S, Vh = linalg.svd(dfx) #SVD de la matrice A
print("Dimensions de U :", U.shape)
k = int(0.4 * len(U)) #nb de colonnes à garder dans U
print(k)
dfx_truncated = U[:, :k] #ajout des colonnes U
print("Dimensions de dfx_truncated :", dfx_truncated.shape)

np.savetxt('U.tsv', U, delimiter='\t')
np.savetxt('dfx_truncated.tsv', dfx_truncated, delimiter='\t')
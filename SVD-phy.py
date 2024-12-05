import pandas as pd
import numpy as np
from scipy import linalg

dfx = pd.read_csv("strigamia-acuminata_profiles.tsv.short", sep="\t", index_col=0)


def SVD_phy(A, p):

    #SVD de la matrice A
    U, S, V = linalg.svd(A,False)

    #nb de colonnes à garder dans U
    k = int(p * len(U))

    #ajout des colonnes U
    U_truncated = U[:, :k]

    # normalisation
    for i in range(0, len(U_truncated)):
        for j in range(0, k):
            U_truncated[i][j] = U_truncated[i][j] / np.max(U_truncated[i])

    # Calculate the distance euclidienne
    SVDphy_distance = np.zeros((len(U_truncated), len(U_truncated)))
    for i in range(0, len(U_truncated)):
        for j in range(0, len(U_truncated)):
            SVDphy_distance[i][j] = np.sqrt(np.sum((U_truncated[i] - U_truncated[j])**2))
    
    row_labels = A.index
    col_labels = A.index
    SVDphy_distance = pd.DataFrame(SVDphy_distance, index=row_labels, columns=col_labels)

    return SVDphy_distance

result = SVD_phy(dfx, 0.6)
print(result.shape)
print(result)


# afficher les 10 meilleurs scores de la première ligne
first_row = result.iloc[0]
top_10_indices = np.argsort(first_row)[:10]
top_10_scores = first_row[top_10_indices]

print("Top 10 scores de la première ligne:")
for i, score in zip(top_10_indices, top_10_scores):
    print(f"Index: {i}, Score: {score}")

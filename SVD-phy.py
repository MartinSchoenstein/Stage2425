import pandas as pd
from scipy import linalg

dfx = pd.read_csv("strigamia-acuminata_profiles.tsv.short", sep="\t", index_col=0)

U, S, Vh = linalg.svd(dfx)

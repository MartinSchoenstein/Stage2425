import pandas as pd
import random

df = pd.read_csv("strigamia-acuminata_profiles.tsv.short", sep="\t", index_col=0)

for i in range(0, len(df)):
    for j in range(0, len(df.columns)):
        if df.iloc[i][j] == 1:
            df.iloc[i][j] = random.randint(60, 10000)
        if df.iloc[i][j] == 0:
            df.iloc[i][j] = random.randint(1, 60)
print(df)

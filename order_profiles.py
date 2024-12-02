import pandas as pd
from newick import read
from Bio import Phylo

phylo_file = "/gstock/MetaInvert/martin/tree_ortho_eukaryota_myriapoda.nwk"
prof_file = "/home/schoenstein/stage2425/strigamia-acuminata_profiles.tsv.short"
df = pd.read_csv(prof_file, sep="\t", index_col=0)
phylo = Phylo.read(phylo_file, "newick")
leaf = phylo.get_terminals()
new_df = pd.DataFrame()

for x in leaf:
    new_df[str(x)] = df[str(x)]

new_df.to_csv(
    "/home/schoenstein/stage2425/ordered_strigamia-acuminata_profiles.tsv.short",
    sep="\t",
)

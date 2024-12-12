import pandas as pd
from distances import distance_profiles, input, transition_vector, order_by_tree
import time

#input
dfx = input("profiles_HUMAN.txt")
lprotein = ['IFT43_HUMAN','IF122_HUMAN','TT21B_HUMAN','IF140_HUMAN','WDR19_HUMAN','WDR35_HUMAN',
            'IFT20_HUMAN','IFT22_HUMAN','IFT25_HUMAN','IFT27_HUMAN','IFT46_HUMAN','IFT52_HUMAN',
            'IFT57_HUMAN','IFT74_HUMAN','IFT80_HUMAN','IFT81_HUMAN','IFT88_HUMAN','IF172_HUMAN',
            'IFT56_HUMAN','TRAF3_HUMAN','CLUA1_HUMAN','BBS1_HUMAN','BBS2_HUMAN',
            'ARL6_HUMAN','BBS4_HUMAN','BBS5_HUMAN','BBS7_HUMAN','TTC8_HUMAN','PTHB1_HUMAN','BBIP1_HUMAN',
            'VPS51_HUMAN','VPS52_HUMAN','VPS53_HUMAN','VPS54_HUMAN']
dfy = dfx.loc[lprotein]
print(dfx['Komagataella phaffii (strain GS115 / ATCC 20864)'])
#newik
ordered_dfx = order_by_tree(dfx, "newick_quote.txt")
#ordered_dfy = order_by_tree(dfy, "newick_quote.txt")
#print("order by tree done")

#Jaccard = distance_profiles("jaccard", dfy, dfx)
#Hamming = distance_profiles("hamming", dfy, dfx)
#Pearson = distance_profiles("pearson", dfy, dfx)
#MI = distance_profiles("mi", dfy, dfx)
#cotransition = distance_profiles("cotransition", ordered_dfy, ordered_dfx)
#PCS = distance_profiles("pcs", dfy, dfx)


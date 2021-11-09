from sys import argv
script,folder=argv

import pandas as pd
import os
filename = os.listdir(folder)
print(filename)
#filename.remove(".gmt.gmt")
all_gene = []
for item in filename:
    if "grey" not in item:
        print(item)
        genes = pd.read_table(folder+item)
        genes = genes.iloc[:,0].str.split(".",expand=True)[0]
        #print(item)
        all_gene.append(item+"\t"+"NA"+"\t"+"\t".join(list(genes))+"\n")

with open(folder+"gmt.gmt","w") as f:
    for item in all_gene:
        f.write(item)


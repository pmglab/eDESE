# -*- coding: utf-8 -*-
"""
Created on Wed Aug 18 13:00:30 2021

@author: DEII
"""

import os
import pandas as pd
all_brain = os.listdir("E:\\000research\\000eqtl_transcript_gwas\\MCGA_result\\Revised_eLife\\WGCNA\\Brain_tmm_normalized_geneExpression\\")

brain = pd.read_table("E:\\000research\\000eqtl_transcript_gwas\\MCGA_result\\Revised_eLife\\WGCNA\\Brain_tmm_normalized_geneExpression\\"+all_brain[0],sep="\t")
common_gene = list(brain.Name)
for item in all_brain[1:]:
    brain = pd.read_table("E:\\000research\\000eqtl_transcript_gwas\\MCGA_result\\Revised_eLife\\WGCNA\\Brain_tmm_normalized_geneExpression\\"+item,sep="\t")
    gene = list(brain.Name)
    common_gene = [item for item in gene if item in common_gene]

for item in all_brain:
    brain = pd.read_table("E:\\000research\\000eqtl_transcript_gwas\\MCGA_result\\Revised_eLife\\WGCNA\\Brain_tmm_normalized_geneExpression\\"+item,sep="\t")
    brain1 = brain[brain.Name.isin(common_gene)]
    brain1.to_csv("E:\\000research\\000eqtl_transcript_gwas\\MCGA_result\\Revised_eLife\\WGCNA\\Brain_tmm_normalized_geneExpression\\"+item+"concensus_gene.txt",sep="\t",index=False)    
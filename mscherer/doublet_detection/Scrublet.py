#!/usr/bin/env python
# coding: utf-8

# In[4]:


import scanpy as sc
import scrublet as scr
import pandas as pd


# In[5]:


dat = pd.read_table('/users/lvelten/project/Methylome/analysis/missionbio/tapestri/Sample4_80_percent/tsv/Sample4_80_percent.barcode.cell.distribution.tsv',sep='\t')

adata = sc.AnnData(dat)


# In[6]:


adata


# In[ ]:


res = scr.Scrublet(adata)


# In[ ]:


doublet_scores, predicted_doublets = res.scrub_doublets()


# In[ ]:


doublet_scores.to_csv('/users/lvelten/project/Methylome/analysis/missionbio/tapestri/Sample4_80_percent/tsv/doublet_scores.csv')
predicted_doublets.to_csv('/users/lvelten/project/Methylome/analysis/missionbio/tapestri/Sample4_80_percent/tsv/predicted_doublets.csv')

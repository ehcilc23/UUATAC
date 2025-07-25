import numpy as np
import scanpy as sc
import pandas as pd
import datatable
import os
import gc
# import dask.dataframe
import pandas.util.testing as tm

### Sample cells 
df_id = pd.read_csv('cellspass_QC_UUATAC.csv')
df_id['cellsPass'] = df_id['cellsPass'].str.replace(r'.*#','')
df_id.head()
df_id.shape

df_id_use = df_id.sample(frac = 1) 
df_id_use = df_id_use.iloc[0:3000,]
df_id_use.head()
df_id_use.shape

df_id_use.to_csv('cellspass_QC_UUATAC_sample3k_celluse.csv', index=False)

### get select cells
df = datatable.fread("UUATAC_mouse_brain.sort.bed.gz", sep = '\t').to_pandas()
df.columns = ['chr', 'start', 'end', 'barcode', 'Frag']

df = df.loc[df['barcode'].isin(df_id_use['cellsPass'])]
df.shape 
df.to_csv('cellspass_QC_UUATAC_sample3k_fragmentuse.tsv', sep='\t', header=False, index=False)

# sort -k1,1V -k2,2n -k3,3n cellspass_QC_UUATAC_sample3k_fragmentuse.tsv | bgzip > cellspass_QC_UUATAC_sample3k_fragmentuse.tsv.gz

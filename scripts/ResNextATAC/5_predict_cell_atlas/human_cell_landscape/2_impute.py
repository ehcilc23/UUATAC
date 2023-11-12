import h5py, os, argparse, logging, time
import numpy as np
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument("--region",dest="region")
#parser.add_argument("--anno", dest="anno")
args = parser.parse_args()
table = pd.read_table(args.region, index_col=0, header=None)

def one_hot(seq):
    seq_len = len(seq.item(0))
    seqindex = {'A':0, 'C':1, 'G':2, 'T':3, 'a':0, 'c':1, 'g':2, 't':3}
    seq_vec = np.zeros((len(seq),seq_len,4), dtype='bool')
    for i in range(len(seq)):
        thisseq = seq.item(i)
        for j in range(seq_len):
            try:
                seq_vec[i,j,seqindex[thisseq[j]]] = 1
            except:
                pass
    return seq_vec
seq = one_hot(table.values)

X=seq
X = X.swapaxes(-1,1).astype(np.float32)
X.shape

import sys
sys.path.append('/file/path/prefix/arch/')
from ResNeXt_modified import *

model_new = resnext34()

import h5py, os, argparse, logging, time



import torch
from torch import nn
from torch.optim import Adam
from torch.utils.data import DataLoader

load_params = torch.load("/file/path/prefix/best_model.pth")
model_params = model_new.state_dict()


model_params.keys()

# extract the corresponding params
co_params = {k: v for k, v in load_params.items() if k in model_params.keys()}

model_params.update(co_params)
model_new.load_state_dict(model_params)

cell_embed=load_params['fc1.dense.weight']
cell_embed.shape

device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
model_new.eval()

model_new.to(device)

w=cell_embed.cpu().numpy()
import scipy.sparse
from scipy import sparse
import anndata
n_samples = X.shape[0]
from tqdm import tqdm
from numpy import count_nonzero
from scipy.sparse import csr_matrix 
bs = 100
y_pred = []
for i in tqdm(range(0, n_samples, bs)):#
    X_batch = X[i:i+bs, ...]
    X_batch = torch.from_numpy(X_batch).to(device)
    y_pred_batch = model_new.forward(X_batch).squeeze().cpu().data.numpy()
   # y_pred_batch= np.where(y_pred_batch > t, 1, 0).astype(np.float32)
    y_pred_batch=y_pred_batch.astype(np.float32)
    y_pred_batch=np.dot(y_pred_batch,w.T)
    y_pred_batch=np.divide(1,1+np.exp(-y_pred_batch))
    tmp= np.where(y_pred_batch>0.99, 1, 0).astype(np.float32)
    tmp1=tmp[np.sum(tmp,axis=1)>1500]
    #print(tmp.shape)
    #s=1-count_nonzero(A)/A.size
    #print(s)
    ad1=anndata.AnnData(tmp1)
    ad1.X=sparse.csr_matrix(ad1.X)
    y_pred.append(ad1)
    
ad=anndata.concat(y_pred,join='outer')
print(ad.X.shape)
tmp=ad.X

tmp=sparse.csr_matrix(tmp)

from scipy.io import mmwrite
mmwrite('./pmat_pred_imputed.mtx',tmp)

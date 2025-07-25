import sys
sys.path.append("/software/path/prefix/NvTK/")
import h5py, os, argparse, logging, time
import numpy as np
import pandas as pd
import scipy.sparse
from scipy import sparse
import anndata
from tqdm import tqdm
import matplotlib.pyplot as plt

import torch
from torch import nn
from torch.optim import Adam
from torch.utils.data import DataLoader

import NvTK
from NvTK import Trainer
from NvTK.Evaluator import calculate_correlation, calculate_pr, calculate_roc
from NvTK.Explainer import get_activate_W, meme_generate, save_activate_seqlets, calc_frequency_W

from NvTK.Explainer import seq_logo, plot_seq_logo
from NvTK.Modules import BasicPredictor

# set_all_random_seed
NvTK.set_random_seed()
NvTK.set_torch_seed()
NvTK.set_torch_benchmark()
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")


## define model
sys.path.append("/software/path/prefix/arch/")
from ResNeXt_conv1_128_btnk_2dense import *
n_tasks = 50040
model = resnext34(num_classes = n_tasks)

# define criterion
criterion = nn.BCELoss().to(device)

# define optimizer
optimizer = Adam(model.parameters(), lr=1e-3, betas=(0.9, 0.999), eps=1e-08, weight_decay=0,)

# define trainer
trainer = Trainer(model, criterion, optimizer, device, 
                    patience=10, tasktype='binary_classification', metric_sample=100,
                    use_tensorbord=True)


## reload best model
model = trainer.load_best_model('/file/path/prefix/best_model.pth')
model.eval()

# load input region
parser = argparse.ArgumentParser()
parser.add_argument("--region",dest="region")
parser.add_argument("--anno", dest="anno")
args = parser.parse_args()
table = pd.read_table(args.region, index_col=0, header=None)
blank = 'N'*500

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
X = seq
X = X.swapaxes(-1,1).astype(np.float32)
del seq
n_samples = X.shape[0]
print(n_samples)

# load cell annotation
anno = pd.read_table(args.anno,sep='\t')

## model prediction for batches
bs = 2000
y_pred = []
os.makedirs("Predict_track", exist_ok=True)
for i in tqdm(range(0, n_samples, bs)):#
    X_batch = X[i:i+bs, ...]
    X_batch = torch.from_numpy(X_batch).to(device)
    y_pred_batch = model.forward(X_batch).cpu().data.numpy()
    pmat=y_pred_batch.astype(np.float32)
    blank_idx = (table[[1]][i:i+bs] == blank).values.flatten()
    pmat[blank_idx,:] = 0
    pmat=pd.DataFrame(pmat.T)
    pmat.index=anno.id
    pmat=pmat.groupby(anno.id.values).mean()
    ad=anndata.AnnData(sparse.csr_matrix(pmat.T))
    y_pred.append(ad)

ad=anndata.concat(y_pred,join='outer')
anno_id=pmat.index
anno_id
del pmat

for j in range(ad.X.shape[-1]):
        predictions_test=ad.X[:,j].toarray().flatten()
        df = pd.DataFrame(np.column_stack((table.index, predictions_test)))
        df.to_csv("./Predict_track/merged_"+anno_id[j]+".txt", index=False, header=False, sep='\t')

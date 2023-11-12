import sys
sys.path.append("/software/path/prefix/NvTK/")
import h5py, os, argparse, logging, time

import numpy as np
import pandas as pd

import torch
from torch import nn
from torch.optim import Adam
from torch.utils.data import DataLoader

import NvTK
from NvTK import Trainer
from NvTK.Evaluator import calculate_correlation, calculate_pr, calculate_roc
from NvTK.Explainer import get_activate_W, meme_generate, save_activate_seqlets, calc_frequency_W

import matplotlib.pyplot as plt
from NvTK.Explainer import seq_logo, plot_seq_logo

#from NvTK import resnet18
from NvTK.Modules import BasicPredictor
# set_all_random_seed
NvTK.set_random_seed()
NvTK.set_torch_seed()
NvTK.set_torch_benchmark()
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

## define model
sys.path.append("/file/path/prefix/arch/")
from ResNeXt_conv1_128_btnk_2dense import *
n_tasks=50040
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

#load

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

import pandas as pd
anno = pd.read_table(args.anno,sep='\t')
anno

X=seq
X = X.swapaxes(-1,1).astype(np.float32)
X.shape
#X = torch.tensor(X).unsqueeze(2)
#X.shape
#X=X.numpy()
#X.shape
import scipy.sparse
from scipy import sparse

import anndata
n_samples = X.shape[0]
print(n_samples)
from tqdm import tqdm
bs = 2000
y_pred = []
os.makedirs("Predict_track", exist_ok=True)
t=0
for i in tqdm(range(0, n_samples, bs)):#
    X_batch = X[i:i+bs, ...]
    X_batch = torch.from_numpy(X_batch).to(device)
    y_pred_batch = model.forward(X_batch).cpu().data.numpy()
    ad=anndata.AnnData(sparse.csr_matrix(y_pred_batch.astype(np.float32)))
    blank_idx = (table[[1]][i:i+bs] == blank).values.flatten()
    pmat=ad.X
    pmat[blank_idx,:] = 0
    pmat1=pd.DataFrame(pmat.toarray().T)
#anno.anno=anno.subcluster 
    pmat1.index=anno.id
    pmat1=pmat1.groupby(anno.id.values).mean()
    predictions_test=pmat1.T.values
    for j in range(predictions_test.shape[-1]):
        p_test = predictions_test[:, j].flatten()
        df = pd.DataFrame(np.column_stack((table.index[i:i+bs], p_test)))
        df.to_csv("./Predict_track/test_predict_"+pmat1.index[j]+'_'+str(t)+".txt", index=False, header=False, sep='\t')
    t+=1



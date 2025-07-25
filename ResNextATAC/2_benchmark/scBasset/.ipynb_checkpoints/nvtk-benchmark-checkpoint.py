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
import seaborn as sns

from NvTK.Explainer import seq_logo, plot_seq_logo



# set_all_random_seed
NvTK.set_random_seed()
NvTK.set_torch_seed()
NvTK.set_torch_benchmark()

# help-functions
def choose_activation(activation1):
    activation_args={}
    if activation1 == "ReLU":
        activation = nn.ReLU
    elif activation1 == "Sigmoid":
        activation = nn.Sigmoid
    elif activation1 == "Tanh":
        activation = nn.Tanh
    elif activation1 == "LeakyReLU":
        activation = nn.LeakyReLU
        activation_args={"negative_slope":0.2}
    elif activation1 == "Exp":
        from NvTK.Modules import EXP
        activation = EXP
    return activation, activation_args

# logging
os.makedirs("Log", exist_ok=True)
logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s %(filename)s[line:%(lineno)d] %(levelname)s %(message)s',
                    datefmt='%a, %d %b %Y %H:%M:%S',
                    filename=time.strftime('./Log/log_Nvwa.%m%d.%H:%M:%S.txt'),
                    filemode='w')

# args
parser = argparse.ArgumentParser()
parser.add_argument("data")
parser.add_argument("--tasktype", dest="tasktype", default='regression', type=str)

# parser.add_argument("--mode", dest="mode", default="train")
parser.add_argument("--gpu-device", dest="device_id", default="0")
parser.add_argument("--use_tensorbord", dest="use_tensorbord", default=True, type=bool)

parser.add_argument("--subset_task_by", dest="subset_task_by", default=None, type=str)
parser.add_argument("--subset_task", dest="subset_task", default=None, type=str)
parser.add_argument("--sample_size", dest="sample_size", default=None, type=int)
parser.add_argument("--pos", dest="pos", default=None, type=int)
parser.add_argument("--left_pos", dest="left_pos", default=None, type=int)
parser.add_argument("--right_pos", dest="right_pos", default=None, type=int)

parser.add_argument("--numFiltersConv1", dest="numFiltersConv1", default=128, type=int)
parser.add_argument("--filterLenConv1", dest="filterLenConv1", default=15, type=int)
parser.add_argument("--pooltype", dest="pooltype", default='avg', type=str)
parser.add_argument("--Pool1", dest="Pool1", default=15, type=int)
parser.add_argument("--activation", dest="activation1", default='ReLU', type=str)
# parser.add_argument("--globalpoolsize", dest="globalpoolsize", default=None, type=int)
parser.add_argument("--use_BN", dest="use_BN", default=False, type=bool)

parser.add_argument("--use_CBAM", dest="use_CBAM", default=False, type=bool)
parser.add_argument("--use_Transformer", dest="use_transformer", action="store_true", default=False)
parser.add_argument("--use_ResNet", dest="use_ResNet", default=None, type=int)
parser.add_argument("--use_DeepCNN", dest="use_DeepCNN", default=1, type=int)
# parser.add_argument("--use_CharCNN", dest="use_CharCNN", action="store_true", default=False)

parser.add_argument("--tower_by", dest="tower_by", default="Celltype")
parser.add_argument("--tower_hidden", dest="tower_hidden", default=None, type=int)

parser.add_argument("--l1_weight", dest="l1_weight", default=0, type=float)
parser.add_argument("--use_focal_loss", dest="use_focal_loss", action="store_true", default=False)
parser.add_argument("--use_scMTLoss", dest="use_scMTLoss", action="store_true", default=False)
parser.add_argument("--scMTLoss_level", dest="scMTLoss_level", default="Celltype", type=str)

parser.add_argument("--use_data_rc_augment", dest="use_data_rc_augment", action="store_true", default=False)
parser.add_argument("--use_data_shift_augment", dest="use_data_shift_augment", action="store_true", default=False)
parser.add_argument("--use_hard_mining_augment", dest="use_hard_mining_augment", action="store_true", default=False)

parser.add_argument("--lr", dest="lr", default=1e-5, type=float)
parser.add_argument("--EPOCH", dest="EPOCH", default=1000, type=int)
parser.add_argument("--patience", dest="patience", default=10, type=int)
parser.add_argument("--batch_size", dest="batch_size", default=8, type=int)
parser.add_argument("--metric_sample", dest="metric_sample", default=100, type=int)

args = parser.parse_args()
logging.info(args)

## change device
os.environ["CUDA_VISIBLE_DEVICES"] = args.device_id
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")


## unpack dataset
h5file = h5py.File(args.data, 'r')
X = h5file["pmat"]["X"][:].swapaxes(-1,1).astype(np.float32)
peak_idx = h5file['pmat']['pmat_sc']['i'][:]
cell_idx = h5file['pmat']['pmat_sc']['j'][:]
x = h5file['pmat']['pmat_sc']['x'][:]
dim = h5file['pmat']['pmat_sc']['dim'][:]
y = np.zeros((dim[0], dim[1]), dtype = np.float32)
y[peak_idx, cell_idx] = x
features = h5file["pmat"]["peak"][:]
h5file.close()

# unpack anno
n_tasks = y.shape[-1]
mask = features[:,-1].astype(str)
train_idx, val_idx, test_idx = mask=='train', mask=='val', mask=='test'
x_train, x_val, x_test = X[train_idx], X[val_idx], X[test_idx]
y_train, y_val, y_test = y[train_idx], y[val_idx], y[test_idx]
# np.sum(train_idx), np.sum(val_idx), np.sum(test_idx)


# define data loader
batch_size = args.batch_size
train_loader = DataLoader(list(zip(x_train, y_train)), batch_size=batch_size,
                            shuffle=True, num_workers=0, drop_last=False, pin_memory=True)
validate_loader = DataLoader(list(zip(x_val, y_val)), batch_size=batch_size, 
                            shuffle=False, num_workers=0, drop_last=False, pin_memory=True)
test_loader = DataLoader(list(zip(x_test, y_test)), batch_size=batch_size, 
                            shuffle=False, num_workers=0, drop_last=False, pin_memory=True)

# # define model
# pooltype = nn.MaxPool1d if args.pooltype == 'max' else nn.AvgPool1d
# GAP = 64 if args.tasktype == 'regression' else 8
# activation, activation_args = choose_activation(args.activation1)
# model = BaselineCNN(output_size=n_tasks, 
#                     emb_planes=args.numFiltersConv1, kernel_size=args.filterLenConv1, bn=args.use_BN,
#                     activation=activation, activation_args=activation_args,
#                     pool=pooltype, pool_args={"kernel_size":args.Pool1},
#                     n_deep_convlayers=args.use_DeepCNN, use_CBAM=args.use_CBAM, GAP=GAP,
#                     tasktype=args.tasktype)

# if args.use_ResNet:
#     from BaselineModel import get_resnet
#     model = get_resnet(output_size=n_tasks, layers=args.use_ResNet, tasktype=args.tasktype)

# if args.use_transformer:
#     from BaselineModel import get_transformer
#     model = get_transformer(x_test[:2].shape, output_size=n_tasks, tasktype=args.tasktype)




## define model
from Model import *
model = scBasset(n_cells= n_tasks)


logging.info(model.__str__())

optimizer = Adam(model.parameters(), lr=args.lr, betas=(0.9, 0.999), eps=1e-08, weight_decay=0,)

if args.tasktype == "regression":
    criterion = nn.MSELoss().to(device)
elif args.tasktype == "binary_classification":
    criterion = nn.BCELoss().to(device)
elif args.tasktype == "classification":
    criterion = nn.CrossEntropyLoss().to(device)

trainer = Trainer(model, criterion, optimizer, device, 
                    patience=args.patience, tasktype=args.tasktype, metric_sample=args.metric_sample,
                    use_tensorbord=args.use_tensorbord)

# train
trainer.train_until_converge(train_loader, validate_loader, test_loader, EPOCH=args.EPOCH)

# reload best model
model = trainer.load_best_model()
model.eval()

os.makedirs("./Test", exist_ok=True)
# predict test-set
_, _, test_predictions, test_targets = trainer.predict(test_loader)
pd.DataFrame(test_targets).to_csv("./Test/test_targets.csv")
pd.DataFrame(test_predictions).to_csv("./Test/test_predictions.csv")
# pd.DataFrame(test_targets, columns=task_name, index=test_gene).to_csv("./Test/test_targets.csv")
# pd.DataFrame(test_predictions, columns=task_name, index=test_gene).to_csv("./Test/test_predictions.csv")

# metric test-set
if args.tasktype == 'regression':
    correlation, pvalue = calculate_correlation(test_targets, test_predictions, method="pearson")
    pcc = [correlation[k] for k in correlation.keys() if k not in ["macro", "micro"]] # dict keys ordered by default in py3.7+

    correlation, pvalue = calculate_correlation(test_targets, test_predictions, method="spearman")
    scc = [correlation[k] for k in correlation.keys() if k not in ["macro", "micro"]] # dict keys ordered by default in py3.7+
    
    pd.DataFrame({"pcc":pcc, "scc":scc}).T.to_csv("./Test/Metric.csv")

else:
    fpr, tpr, roc_auc = calculate_roc(test_targets, test_predictions)
    auroc = [roc_auc[k] for k in roc_auc.keys() if k not in ["macro", "micro"]] # dict keys ordered by default in py3.7+

    p, r, average_precision = calculate_pr(test_targets, test_predictions)
    aupr = [average_precision[k] for k in average_precision.keys() if k not in ["macro", "micro"]] # dict keys ordered by default in py3.7+
    
    pd.DataFrame({"auroc":auroc, "aupr":aupr}).T.to_csv("./Test/Metric.csv")
    
# plot 
metric = pd.read_csv('./Test/Metric.csv',index_col=0)
metric = metric.T
metric

plt.figure(dpi = 100)
sns.boxplot(y = metric['auroc'],
            orient = 'v',
            showfliers=True,# don't show outliers
            fliersize=3, # outlier size
            boxprops = {'color':'black',# box line color
        'facecolor':'#9983BD'# box color
        },
        )
plt.xlabel('median=%0.2f' % np.median(metric['auroc']))
plt.savefig('./Figures/auroc.pdf')
# plt.show()

plt.figure(dpi = 100)
sns.boxplot(y = metric['aupr'],
            orient = 'v',
            showfliers=True,# don't show outliers
            fliersize=3, # outlier size
            boxprops = {'color':'black',# box line color
        'facecolor':'#F49B7C'# box color
        },
        )
plt.xlabel('median=%0.2f' % np.median(metric['aupr']))
plt.savefig('./Figures/aupr.pdf')
# plt.show()


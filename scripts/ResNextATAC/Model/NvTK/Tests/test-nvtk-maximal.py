import h5py, os, argparse, logging, time

import numpy as np
import pandas as pd

import torch
from torch import nn
from torch.optim import Adam
from torch.utils.data import DataLoader

from NvTK import Trainer
from NvTK.Model.ConvModel import CNN
from NvTK.Evaluator import calculate_roc, calculate_pr
from NvTK.Explainer import get_activate_W, meme_generate, save_activate_seqlets

# logging
os.makedirs("./Log", exist_ok=True)
logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s %(filename)s[line:%(lineno)d] %(levelname)s %(message)s',
                    datefmt='%a, %d %b %Y %H:%M:%S',
                    filename=time.strftime('./Log/log_nvtk_minimal.%m%d.%H:%M:%S.txt'),
                    filemode='w')

# args
parser = argparse.ArgumentParser()
parser.add_argument("data")
parser.add_argument("--gpu-device", dest="device_id", default="0")
args = parser.parse_args()
logging.info(args)

## change device
os.environ["CUDA_VISIBLE_DEVICES"] = args.device_id
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

# unpack datasets
h5file = h5py.File(args.data, 'r')
anno = h5file["annotation"][:]
x_train = h5file["train_data"][:].astype(np.float32)
y_train = h5file["train_label"][:].astype(np.float32)
x_val = h5file["val_data"][:].astype(np.float32)
y_val = h5file["val_label"][:].astype(np.float32)
x_test = h5file["test_data"][:].astype(np.float32)
y_test = h5file["test_label"][:].astype(np.float32)
h5file.close()

# unpack anno
n_tasks = anno.shape[0]
task_name = anno[:,0]

# define data loader
batch_size = 16
train_loader = DataLoader(list(zip(x_train, y_train)), batch_size=batch_size,
                            shuffle=True, num_workers=2, drop_last=False, pin_memory=True)
validate_loader = DataLoader(list(zip(x_val, y_val)), batch_size=batch_size, 
                            shuffle=False, num_workers=2, drop_last=False, pin_memory=True)
test_loader = DataLoader(list(zip(x_test, y_test)), batch_size=batch_size, 
                            shuffle=False, num_workers=2, drop_last=False, pin_memory=True)

# define model
model = CNN(output_size=n_tasks, tasktype='binary_classification')
optimizer = Adam(model.parameters(), lr=1e-4, betas=(0.9, 0.999), eps=1e-08, weight_decay=0,)
criterion = nn.BCELoss().to(device)
trainer = Trainer(model, criterion, optimizer, device, tasktype='binary_classification')

# train
trainer.train_until_converge(train_loader, validate_loader, test_loader, EPOCH=500)
# reload best model
model = trainer.load_best_model()
# predict test-set
_, _, test_predictions, test_targets = trainer.predict(test_loader)
# metric test-set
fpr, tpr, roc_auc = calculate_roc(test_targets, test_predictions)
auroc = [roc_auc[k] for k in roc_auc.keys() if k not in ["macro", "micro"]] # dict keys ordered by default in py3.7+

p, r, average_precision = calculate_pr(test_targets, test_predictions)
aupr = [average_precision[k] for k in average_precision.keys() if k not in ["macro", "micro"]] # dict keys ordered by default in py3.7+

pd.DataFrame({"auroc":auroc, "aupr":aupr}, index=task_name).T.to_csv("Metric.csv")

# explain based on feature-map
W = get_activate_W(model, model.Embedding.conv, test_loader, motif_width=3)
meme_generate(W, output_file='meme.txt', prefix='Filter_')

save_activate_seqlets(model, model.Embedding.conv, test_loader, threshold=0.999,
                        out_fname='seqlets.fasta', motif_width=3)
os.system('homer seqlets.fasta')

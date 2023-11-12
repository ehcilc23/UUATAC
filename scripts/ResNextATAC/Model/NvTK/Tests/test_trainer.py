import torch
from torch import nn
import torch.nn.functional as F
import torch.utils.data as Data

import matplotlib.pyplot as plt
from sklearn.model_selection import train_test_split
import numpy as np

from NvTK.Trainer import Trainer

np.random.seed(666)
X = np.linspace(-1, 1, 1000)
y = np.power(X, 2) + 0.1 * np.random.normal(0, 1, X.size)
print(X.shape)
print(y.shape)
plt.scatter(X, y)

X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.3, random_state=1024)

X_train = torch.from_numpy(X_train).type(torch.FloatTensor)
X_train = torch.unsqueeze(X_train, dim=1)  #转换成二维
y_train = torch.from_numpy(y_train).type(torch.FloatTensor)
y_train = torch.unsqueeze(y_train, dim=1)

X_test = torch.from_numpy(X_test).type(torch.FloatTensor)
X_test = torch.unsqueeze(X_test, dim=1)  #转换成二维

#将数据装载dataloader中, 对数据进行分批训练
torch_data  = Data.TensorDataset(X_train, y_train)
loader = Data.DataLoader(dataset=torch_data, batch_size=50, shuffle=True)
next(iter(loader))[0].shape, next(iter(loader))[-1].shape

#创建自己的nn
class Net(nn.Module):
    def __init__(self):
        super(Net, self).__init__()
        self.hidden = nn.Linear(1, 20)
        self.predict = nn.Linear(20, 1)
 
    def forward(self, x):
        x = F.relu(self.hidden(x))
        x = self.predict(x)
        return x
    
    def save(self, fname):
        torch.save(self.parameters, fname)

net = Net()
device = torch.device("cpu")
optimizer = torch.optim.Adam(net.parameters(), lr=1e-3, betas=(0.9, 0.999), eps=1e-08, weight_decay=0,)
criterion = torch.nn.MSELoss().to(device)
trainer = Trainer(net, criterion, optimizer, device, pred_prob=False)

trainer.train_until_converge(loader, loader, loader, 100, resume=False, verbose_step=5)

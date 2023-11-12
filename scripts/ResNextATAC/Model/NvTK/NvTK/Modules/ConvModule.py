import logging

import torch
from torch import nn

from .BasicModule import BasicModule


class ShallowWideConvModule(nn.Module):
    def __init__(self, numFiltersConv1=40, filterLenConv1=5,
                        numFiltersConv2=44, filterLenConv2=15,
                        numFiltersConv3=44, filterLenConv3=25):
        
        super().__init__()
        self.conv1 = nn.Conv1d(in_channels=4, out_channels=numFiltersConv1, 
                                kernel_size=filterLenConv1, padding=(filterLenConv1 - 1) // 2)
        self.conv2 = nn.Conv1d(in_channels=4, out_channels=numFiltersConv2, 
                                kernel_size=filterLenConv2, padding=(filterLenConv2 - 1) // 2)
        self.conv3 = nn.Conv1d(in_channels=4, out_channels=numFiltersConv3, 
                                kernel_size=filterLenConv3, padding=(filterLenConv3 - 1) // 2)

    def forward(self, x):
        logging.debug(x.shape)

        out1 = self.conv1(x)
        logging.debug(out1.shape)
        
        out2 = self.conv2(x)
        logging.debug(out2.shape)

        out3 = self.conv3(x)
        logging.debug(out3.shape)
        
        out = torch.cat([out1, out2, out3], dim=1)
        logging.debug(out.shape)

        return out


class DeepConvModule(nn.Module):
    def __init__(self, in_planes, out_planes, Layers=4, hidden_channels=256, kernel_size=3,
                    bn=True, activation=nn.ReLU, activation_args={}, 
                    pool=nn.AvgPool1d, pool_args={'kernel_size': 3}):
        super().__init__()
        layers = []
        layers.append(nn.Conv1d(in_planes, hidden_channels, kernel_size=kernel_size))
        
        for _ in range(Layers):
            layers.append(nn.Conv1d(hidden_channels, hidden_channels, kernel_size=kernel_size))
            if bn:
                layers.append(nn.BatchNorm1d(hidden_channels, eps=1e-5, momentum=0.01, affine=True))
            if activation:
                layers.append(activation(**activation_args))
            if pool:
                layers.append(pool(**pool_args))

        layers.append(nn.Conv1d(hidden_channels, out_planes, kernel_size=kernel_size))

        self.deepconv = nn.Sequential(*layers)

    def forward(self, x):
        x = self.deepconv(x)
        logging.debug(x.shape)

        return x



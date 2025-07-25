""" Reimplement of Publicated Model Archiectures applicable in NvTK.
This module provides 

1.  `NINCNN` - Network in Network CNN model Archiecture

2. `DeepSEA` - DeepSEA architecture (Zhou & Troyanskaya, 2015).

3. `Beluga` - DeepSEA architecture used in Expecto (Zhou & Troyanskaya, 2019).

4. `DanQ` - DanQ architecture (Quang & Xie, 2016).

5. `Basset` - Basset architecture (Kelley, 2016).

and supporting methods.
"""

import logging
import numpy as np

import torch
import torch.nn as nn

from ..Modules import BasicModule

__all__ = ['NINCNN', 'DeepSEA', 'Beluga', 'DanQ', 'Basset']

# TODO update Nvwa model
# class Nvwa(BasicModule):
#     def __init__(self, sequence_length, n_genomic_features):
#         super().__init__()

class NINCNN(BasicModule):
    """
        @misc{https://doi.org/10.48550/arxiv.1312.4400,
            doi = {10.48550/ARXIV.1312.4400},
            url = {https://arxiv.org/abs/1312.4400},
            author = {Lin, Min and Chen, Qiang and Yan, Shuicheng},
            keywords = {Neural and Evolutionary Computing (cs.NE), Computer Vision and Pattern Recognition (cs.CV), Machine Learning (cs.LG), FOS: Computer and information sciences, FOS: Computer and information sciences},
            title = {Network In Network},
            publisher = {arXiv},
            year = {2013},
            copyright = {arXiv.org perpetual, non-exclusive license}
        }
    """
    def __init__(self, sequence_length, n_genomic_features):
        super().__init__()
        self.conv1 = nn.Sequential(
            nn.Conv1d(4, 512, 3, 1),
            nn.ReLU(),
            nn.Conv1d(512, 512, 1, 1),
            nn.ReLU(),
            nn.Conv1d(512, 512, 1, 1),
            nn.ReLU(),
            nn.AvgPool1d(3),
        )
        self.conv2 = nn.Sequential(
            nn.Conv1d(512, 768, 3, 1),
            nn.ReLU(),
            nn.Conv1d(768, 768, 1, 1),
            nn.ReLU(),
            nn.Conv1d(768, 768, 1, 1),
            nn.ReLU(),
            nn.AvgPool1d(3),
        )
        self.conv3 = nn.Sequential(
            nn.Conv1d(768, 1024, 3, 1),
            nn.ReLU(),
            nn.Conv1d(1024, 1024, 1, 1),
            nn.ReLU(),
            nn.Conv1d(1024, 1024, 1, 1),
            nn.ReLU(),
            nn.AvgPool1d(3),
        )
        self.classifier = nn.Conv1d(1024, n_genomic_features, 1, 1)
        self.GAP = nn.AdaptiveAvgPool1d(1)
        
    def forward(self, x):
        logging.debug(x.shape)
        x = self.conv1(x)
        logging.debug(x.shape)
        x = self.conv2(x)
        logging.debug(x.shape)
        x = self.conv3(x)
        logging.debug(x.shape)
        x = self.classifier(x)
        logging.debug(x.shape)
        x = self.GAP(x)
        logging.debug(x.shape)
        x = x.view(x.size(0), -1)
        logging.debug(x.shape)
        return x


class DeepSEA(BasicModule):
    """
    DeepSEA architecture (Zhou & Troyanskaya, 2015).
    """
    def __init__(self, sequence_length, n_genomic_features):
        """
        Parameters
        ----------
        sequence_length : int
        n_genomic_features : int
        """
        super(DeepSEA, self).__init__()
        conv_kernel_size = 8
        pool_kernel_size = 4

        self.conv_net = nn.Sequential(
            nn.Conv1d(4, 320, kernel_size=conv_kernel_size),
            nn.ReLU(inplace=True),
            nn.MaxPool1d(
                kernel_size=pool_kernel_size, stride=pool_kernel_size),
            nn.Dropout(p=0.2),

            nn.Conv1d(320, 480, kernel_size=conv_kernel_size),
            nn.ReLU(inplace=True),
            nn.MaxPool1d(
                kernel_size=pool_kernel_size, stride=pool_kernel_size),
            nn.Dropout(p=0.2),

            nn.Conv1d(480, 960, kernel_size=conv_kernel_size),
            nn.ReLU(inplace=True),
            nn.Dropout(p=0.5))

        reduce_by = conv_kernel_size - 1
        pool_kernel_size = float(pool_kernel_size)
        self.n_channels = int(
            np.floor(
                (np.floor(
                    (sequence_length - reduce_by) / pool_kernel_size)
                 - reduce_by) / pool_kernel_size)
            - reduce_by)
        self.classifier = nn.Sequential(
            nn.Linear(960 * self.n_channels, n_genomic_features),
            nn.ReLU(inplace=True),
            nn.Linear(n_genomic_features, n_genomic_features),
            nn.Sigmoid())

    def forward(self, x):
        """Forward propagation of a batch.
        """
        out = self.conv_net(x)
        reshape_out = out.view(out.size(0), 960 * self.n_channels)
        predict = self.classifier(reshape_out)
        return predict

    def criterion(self):
        """
        The criterion the model aims to minimize.
        """
        return nn.BCELoss()

    def get_optimizer(self, lr):
        """
        The optimizer and the parameters with which to initialize the optimizer.
        At a later time, we initialize the optimizer by also passing in the model
        parameters (`model.parameters()`). We cannot initialize the optimizer
        until the model has been initialized.
        """
        return (torch.optim.SGD,
                {"lr": lr, "weight_decay": 1e-6, "momentum": 0.9})


class LambdaBase(nn.Sequential):
    def __init__(self, fn, *args):
        super(LambdaBase, self).__init__(*args)
        self.lambda_func = fn

    def forward_prepare(self, input):
        output = []
        for module in self._modules.values():
            output.append(module(input))
        return output if output else input

class Lambda(LambdaBase):
    def forward(self, input):
        return self.lambda_func(self.forward_prepare(input))

class Beluga(BasicModule):
    """
    DeepSEA architecture used in Expecto (Zhou & Troyanskaya, 2019).
    """
    def __init__(self, sequence_length, n_genomic_features):
        super(Beluga, self).__init__()
        conv_kernel_size = 8
        pool_kernel_size = 8
        n_hiddens = 32

        reduce_by = (conv_kernel_size - 1) * 2 # conv twice
        self.n_channels = int(
            np.floor(
                (np.floor(
                    (sequence_length - reduce_by) / pool_kernel_size)
                 - reduce_by) / pool_kernel_size)
            - reduce_by)

        self.model = nn.Sequential(
            nn.Sequential(
                nn.Conv2d(4,320,(1, conv_kernel_size)),
                nn.ReLU(),
                nn.Conv2d(320,320,(1, conv_kernel_size)),
                nn.ReLU(),
                nn.Dropout(0.2),
                nn.MaxPool2d((1, pool_kernel_size),(1, pool_kernel_size)),
                nn.Conv2d(320,480,(1, conv_kernel_size)),
                nn.ReLU(),
                nn.Conv2d(480,480,(1, conv_kernel_size)),
                nn.ReLU(),
                nn.Dropout(0.2),
                nn.MaxPool2d((1, pool_kernel_size),(1, pool_kernel_size)),
                nn.Conv2d(480,640,(1, conv_kernel_size)),
                nn.ReLU(),
                nn.Conv2d(640,640,(1, conv_kernel_size)),
                nn.ReLU(),
            ),
            nn.Sequential(
                nn.Dropout(0.5),
                Lambda(lambda x: x.view(x.size(0),-1)),
                nn.Sequential(Lambda(lambda x: x.view(1,-1) if 1==len(x.size()) else x ),nn.Linear(640 * self.n_channels, n_hiddens)),
                nn.ReLU(),
                nn.Sequential(Lambda(lambda x: x.view(1,-1) if 1==len(x.size()) else x ),nn.Linear(n_hiddens, n_genomic_features)),
            ),
            nn.Sigmoid(),
        )

    def forward(self, x):
        x = x.unsqueeze(2) # update 2D sequences
        return self.model(x)


class DanQ(nn.Module):
    """
    DanQ architecture (Quang & Xie, 2016).
    """
    def __init__(self, sequence_length, n_genomic_features):
        """
        Parameters
        ----------
        sequence_length : int
            Input sequence length
        n_genomic_features : int
            Total number of features to predict
        """
        super(DanQ, self).__init__()
        self.nnet = nn.Sequential(
            nn.Conv1d(4, 320, kernel_size=26),
            nn.ReLU(inplace=True),
            nn.MaxPool1d(
                kernel_size=13, stride=13),
            nn.Dropout(0.2))

        self.bdlstm = nn.Sequential(
            nn.LSTM(
                320, 320, num_layers=1, batch_first=True, bidirectional=True))

        self._n_channels = np.floor(
            (sequence_length - 25) / 13)
        self.classifier = nn.Sequential(
            nn.Dropout(0.5),
            nn.Linear(self._n_channels * 640, 925),
            nn.ReLU(inplace=True),
            nn.Linear(925, n_genomic_features),
            nn.Sigmoid())

    def forward(self, x):
        """Forward propagation of a batch.
        """
        out = self.nnet(x)
        reshape_out = out.transpose(0, 1).transpose(0, 2)
        out, _ = self.bdlstm(reshape_out)
        out = out.transpose(0, 1)
        reshape_out = out.contiguous().view(
            out.size(0), 640 * self._n_channels)
        predict = self.classifier(reshape_out)
        return predict

    def criterion(self):
        return nn.BCELoss()

    def get_optimizer(self, lr):
        return (torch.optim.RMSprop, {"lr": lr})


class Basset(BasicModule):
    '''Basset architecture (Kelley, 2016).
    Deep convolutional neural networks for DNA sequence analysis.
    The architecture and optimization parameters for the DNaseI-seq compendium analyzed in the paper.
    '''
    def __init__(self, sequence_length, n_genomic_features):
        super(Basset, self).__init__()

        self.model = nn.Sequential(
            nn.Sequential(
                nn.Conv2d(4,300,(1, 19)),
                nn.BatchNorm2d(300),
                nn.ReLU(),
                nn.MaxPool2d((1, 8),(1, 8)),
                nn.Conv2d(300,200,(1, 11)),
                nn.BatchNorm2d(200),
                nn.ReLU(),
                nn.MaxPool2d((1, 8),(1, 8)),
                nn.Conv2d(200,200,(1, 7)),
                nn.BatchNorm2d(200),
                nn.ReLU(),
                nn.MaxPool2d((1, 8),(1, 8)),
            ),
            nn.Sequential(
                Lambda(lambda x: x.view(x.size(0),-1)),
                nn.Sequential(Lambda(lambda x: x.view(1,-1) if 1==len(x.size()) else x ),nn.Linear(4800, 1000)),
                nn.BatchNorm1d(1000),
                nn.Dropout(0.3),
                nn.ReLU(),
                nn.Linear(1000, 32),
                nn.BatchNorm1d(32),
                nn.Dropout(0.3),
                nn.ReLU(),
                nn.Linear(32, n_genomic_features),
            ),
            nn.Sigmoid(),
        )

    def forward(self, x):
        x = x.unsqueeze(2) # update 2D sequences
        return self.model(x)
    
    def architecture(self):
        d = {'conv_filters1':300,
            'conv_filters2':200,
            'conv_filters3':200,
            'conv_filter_sizes1':19,
            'conv_filter_sizes2':11,
            'conv_filter_sizes3':7,
            'pool_width1':3,
            'pool_width2':4,
            'pool_width3':4,
            'hidden_units1':1000,
            'hidden_units2':32,
            'hidden_dropouts1':0.3,
            'hidden_dropouts2':0.3,
            'learning_rate':0.002,
            'weight_norm':7,
            'momentum':0.98}
        return d


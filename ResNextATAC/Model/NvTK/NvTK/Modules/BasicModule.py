"""Basic module in NvTK.
This module provides 

1.  `BasicModule` class - the general abstract class

2.  `BasicConv1d` class - Basic Convolutional Module (1d)

3.  `BasicRNNModule` class - Basic RNN(LSTM) Module in batch-first style

4.  `BasicLinearModule`

5.  `BasicPredictor` Module

6.  `BasicLoss` Module

and supporting methods.
"""

import random, time, logging
import numpy as np

import torch
from torch import nn

__all__ = ["BasicModule", "BasicConv1d", "BasicRNNModule", "BasicLinearModule", "BasicPredictor", "BasicLoss", "Flatten"]


class BasicModule(nn.Module):
    """Basic module class in NvTK."""

    def __init__(self):
        super(BasicModule,self).__init__()

    def initialize_weights(self):
        """initialize module parameters.

        Conv module weight will be initialize in xavier_normal_,
        bias will be initialize in zero_

        Linear module weight will be initialize in xavier_normal_,
        bias will be initialize in zero_

        BatchNorm module weight will be initialize in constant = 1,
        bias will be initialize in constant = 0

        LSTM module weight will be initialize in orthogonal_
        """
        for m in self.modules():
            if isinstance(m, (nn.Conv1d, nn.Conv2d)):
                nn.init.xavier_normal_(m.weight) 
                if m.bias is not None:
                    m.bias.data.zero_()
                logging.debug("init Conv param...")

            elif isinstance(m, nn.Linear):
                nn.init.xavier_normal_(m.weight) 
                if m.bias is not None:
                    m.bias.data.zero_()
                logging.debug("init Linear param...")

            elif isinstance(m, (nn.BatchNorm1d, nn.BatchNorm2d)):
                nn.init.constant_(m.weight, 1)
                nn.init.constant_(m.bias, 0)
                logging.debug("init BatchNorm param...")

            elif isinstance(m, nn.LSTM):
                nn.init.orthogonal_(m.all_weights[0][0])
                nn.init.orthogonal_(m.all_weights[0][1])
                nn.init.orthogonal_(m.all_weights[1][0])
                nn.init.orthogonal_(m.all_weights[1][1])
                logging.debug("init LSTM param...")

    def initialize_weights_from_pretrained(self, pretrained_net_fname):
        """initialize module weights from pretrained model

        Parameters
        ----------
        pretrained_net_fname : str
            the pretrained model file path (e.g. `checkpoint.pth`).
        """
        pretrained_dict = torch.load(pretrained_net_fname)
        net_state_dict = self.state_dict()
        pretrained_dict = {k: v for k, v in pretrained_dict.items() if k in net_state_dict}
        net_state_dict.update(pretrained_dict)
        self.load_state_dict(net_state_dict)
        logging.info("params loaded from: %s" % pretrained_net_fname)

    def load(self, path):
        """load module weights from saved model 

        Parameters
        ----------
        path : str
            the saved model file path (e.g. `checkpoint.pth`).
        """
        self.load_state_dict(torch.load(path))

    def save(self, fname=None):
        """save module weights to file

        Parameters
        ----------
        fname : str, optional
            Specify the saved model file path.
            Default is "None". Saved file will be formatted as "model.time.pth".
        """
        if fname is None:
            fname = time.strftime("model" + '%m%d_%H:%M:%S.pth')
        torch.save(self.state_dict(), fname)
        return fname

    # TODO (jiaqili, 20220417): test the tensor flow in module forward
    def test(self, input_size):
        device = list(self.parameters)[0].device
        x = torch.zeros(input_size).to(device)
        self.forward(x)
        logging.info("Test: all the tensor flow shape reported")


class Flatten(nn.Module):
    """Flatten Module: flatten the tensor as (batch_size, -1)."""
    def forward(self, x):
        return x.view(x.size(0), -1)


class EXP(nn.Module):
    """Exp Module: calculate the exp of tensor as `x.exp()`"""
    def forward(self, x):
        return x.exp()


class BasicConv1d(BasicModule):
    """Basic Convolutional Module (1d) in NvTK.

    Parameters
    ----------
    in_planes : int
        Number of input channels
    out_planes : int
        Number of output channels produced by the convolution
    kernel_size : int, optional
        Size of the convolving kernel
    conv_args : dict, optional
        Other convolutional args, Default is dict().
        Will be pass to `torch.nn.Conv1d(**conv_args)`
        (e.g. `conv_args={'dilation':1}`)
    bn : bool, optional
        Whether to use BatchNorm1d, Default is True.
    activation : nn.Module, optional
        Activation Module, Default is nn.ReLU.
    activation_args : dict, optional
        Other activation args, Default is dict().
        Will be pass to `activation(**activation_args)`
        (e.g. `activation=nn.LeakyReLU, activation_args={'p':0.2}`)
    dropout : bool, optional
        Whether to use Dropout, Default is True.
    dropout_args : dict, optional
        Dropout args, Default is {'p':0.5}.
        Will be pass to `nn.Dropout(**dropout_args)` if dropout
        (e.g. `dropout=True, dropout_args={'p':0.5}`)
    pool : nn.Module, optional
        Pool Module (1d), Default is nn.AvgPool1d.
    pool_args : dict, optional
        Other pool args, Default is {'kernel_size': 3}.
        Will be pass to `pool(**pool_args)`
        (e.g. `pool=nn.AvgPool1d, pool_args={'kernel_size': 3}`)

    Attributes
    ----------
    in_channels : int

    out_channels : int

    conv : nn.Conv1d
        The convolutional neural network component of the model.
    bn : nn.BatchNorm1d
        The Batch Normalization 
    activation : nn.Module
        The activation Module
    dropout : nn.Dropout
        The Dropout Module
    pool : nn.Module
        The pool Module

    Tensor flows
    ----------
    -> conv(x)

    -> bn(x) if bn
    
    -> activation(x) if activation
    
    -> dropout(x) if dropout
    
    -> pool(x) if pool

    """
    def __init__(self, in_planes, out_planes, kernel_size=3, conv_args={}, 
                    bn=True, 
                    activation=nn.ReLU, activation_args={}, 
                    dropout=True, dropout_args={'p':0.5},
                    pool=nn.AvgPool1d, pool_args={'kernel_size': 3}):
        super().__init__()
        self.in_channels = in_planes
        self.out_channels = out_planes
        self.conv = nn.Conv1d(in_planes, out_planes, kernel_size=kernel_size, **conv_args)
        self.bn = nn.BatchNorm1d(out_planes, eps=1e-5, momentum=0.01, affine=True) if bn else None
        self.activation = activation(**activation_args) if activation is not None else None
        self.dropout = nn.Dropout(**dropout_args) if dropout else None
        self.pool = pool(**pool_args) if pool else None

    def forward(self, x):
        x = self.conv(x)
        if self.bn is not None:
            x = self.bn(x)
        if self.activation is not None:
            x = self.activation(x)
        if self.dropout is not None:
            x = self.dropout(x)
        if self.pool is not None:
            x = self.pool(x)
        logging.debug(x.shape)
        return x


class BasicRNNModule(BasicModule):
    """
    Basic RNN(LSTM) Module in batch-first style
    """
    def __init__(self, LSTM_input_size=512, LSTM_hidden_size=512, LSTM_hidden_layes=2):
        super().__init__()
        self.rnn_hidden_state = None
        self.rnn = nn.LSTM(
            input_size=LSTM_input_size,
            hidden_size=LSTM_hidden_size,
            num_layers=LSTM_hidden_layes,
            batch_first=True, # batch, seq, feature
            bidirectional=True,
        )
    def forward(self, input):
        output, self.rnn_hidden_state = self.rnn(input, None)
        logging.debug(output.shape)
        return output


class BasicLinearModule(BasicModule):
    """
    Basic Linear Module in NvTK.

    Parameters
    ----------
    input_size : int
        Number of input size
    output_size : int
        Number of output size produced by the Linear
    bias : bool, optional
        Bias of the Linear, Default is True.
        It could be False when use BatchNorm.
    bn : bool, optional
        Whether to use BatchNorm1d, Default is True.
    activation : nn.Module, optional
        Activation Module, Default is nn.ReLU.
    activation_args : dict, optional
        Other activation args, Default is dict().
        Will be pass to `activation(**activation_args)`
        (e.g. `activation=nn.LeakyReLU, activation_args={'p':0.2}`)
    dropout : bool, optional
        Whether to use Dropout, Default is True.
    dropout_args : dict, optional
        Dropout args, Default is {'p':0.5}.
        Will be pass to `nn.Dropout(**dropout_args)` if dropout
        (e.g. `dropout=True, dropout_args={'p':0.5}`)

    Attributes
    ----------
    input_size : int

    output_size : int

    linear : nn.Linear
        The Linear neural network component
    bn : nn.BatchNorm1d
        The Batch Normalization 
    activation : nn.Module
        The activation Module
    dropout : nn.Dropout
        The Dropout Module

    Tensor flows
    ----------
    -> linear(x)

    -> bn(x) if bn
    
    -> activation(x) if activation
    
    -> dropout(x) if dropout

    """
    def __init__(self, input_size, output_size, bias=True, bn=True, 
                    activation=nn.ReLU, activation_args={}, 
                    dropout=True, dropout_args={'p':0.5}):
        super().__init__()
        self.input_size = input_size
        self.output_size = output_size
        self.linear = nn.Linear(input_size, output_size, bias=bias)
        self.bn = nn.BatchNorm1d(output_size, eps=1e-5, momentum=0.01, affine=True) if bn else None
        self.activation = activation(**activation_args) if activation is not None else None
        self.dropout = nn.Dropout(**dropout_args) if dropout else None

    def forward(self, x):
        x = self.linear(x)
        if self.bn is not None:
            x = self.bn(x)
        if self.activation is not None:
            x = self.activation(x)
        if self.dropout is not None:
            x = self.dropout(x)
        logging.debug(x.shape)
        return x


class BasicPredictor(BasicModule):
    """BasicPredictor Module in NvTK.
    
    BasicPredictor support task types of 'none', 'binary_classification', 'classification', 'regression';
    1. 'none' : nullify the whole BaiscPredictor with identity
    2. 'binary_classification' : activate with Sigmoid
    3. 'classification' : activate with Softmax(dim=1)
    4. 'regression' : Identity

    Parameters
    ----------
    input_size : int
        Number of input size
    output_size : int
        Number of output size (task numbers)
    tasktype : str, optional
        Specify the task type, Default is "binary_classification".
        (e.g. `tasktype="regression"`)

    Attributes
    ----------
    supported_tasks : currently supported task types
    tasktype : task type of Predictor
    input_size : int
    Map : nn.Linear
        The Linear Module Mapping input to output.
    Pred: nn.Module
        The Activation Module in specified task type.

    Tensor flow
    ----------
    -> Map(x)

    -> Pred(x)

    """
    def __init__(self, input_size, output_size, tasktype='binary_classification'):
        super().__init__()
        self.supported_tasks = ['none', 'binary_classification', 'classification', 'regression']

        self.input_size = input_size
        self.tasktype = tasktype
        
        self.Map = nn.Linear(input_size, output_size, bias=True)
        self.switch_task(tasktype) # init self.Map

    def forward(self, x):
        return self.Pred(self.Map(x))

    def switch_task(self, tasktype):
        """switch to specified task type

        Parameters
        ----------
        tasktype : str
            Specify the task type (e.g. `tasktype="regression"`)
        """

        msg = 'tasktype: %s not supported, check the document' % tasktype
        assert tasktype in self.supported_tasks, msg

        if tasktype == 'none':
            self.Map = nn.Sequential()
            self.Pred = nn.Sequential()
        elif tasktype == 'classification':
            self.Pred = nn.Softmax(dim=1)
        elif tasktype == 'binary_classification':
            self.Pred = nn.Sigmoid()
        elif tasktype == 'regression':
            self.Pred = nn.Sequential()

        self.tasktype = tasktype
    
    # TODO
    # def finetune(self, output_size, tasktype='binary_classification'):
    #     self.Map = nn.Linear(input_size, output_size, bias=True)
    #     self.Pred = nn.Sequential()

    #     self.tasktype = tasktype
    #     self.switch_task(tasktype)

    def current_task(self):
        """return current task type"""
        return self.tasktype

    def remove(self):
        """Predictor.remove: replace predictor with null Sequential,
        same as switch_task('none')
        """
        self.switch_task('none')


class BasicLoss(nn.Module):
    """BasicLoss Module in NvTK.
    
    BasicLoss support task types of 'binary_classification', 'classification', 'regression';
    1. 'binary_classification' : BCELoss function
    2. 'classification' : CrossEntropyLoss function
    3. 'regression' : MSELoss function

    Parameters
    ----------
    tasktype : str, optional
        Specify the task type, Default is "binary_classification".
        (e.g. `tasktype="regression"`)
    reduction : str, optional
        Specifies the reduction to apply to the output: `'none'` | `'mean'` | `'sum'`.

    Attributes
    ----------
    supported_tasks : currently supported task types
    tasktype : task type of Predictor
    loss : loss function

    """
    def __init__(self, tasktype='binary_classification', reduction='mean'):
        super().__init__()
        self.supported_tasks = ['binary_classification', 'classification', 'regression']

        self.tasktype = tasktype
        self.reduction = reduction

        self.switch_task(tasktype) # init self.loss

    def forward(self, pred, target):
        return self.loss(pred, target)

    def switch_task(self, tasktype):
        """switch to specified task type

        Parameters
        ----------
        tasktype : str
            Specify the task type (e.g. `tasktype="regression"`)
        """

        msg = 'tasktype: %s not supported, check the document' % tasktype
        assert tasktype in self.supported_tasks, msg

        if tasktype == 'classification':
            self.loss = nn.CrossEntropyLoss(reduction=self.reduction)
        elif tasktype == 'binary_classification':
            self.loss = nn.BCELoss(reduction=self.reduction)
        elif tasktype == 'regression':
            self.loss = nn.MSELoss(reduction=self.reduction)

        self.tasktype = tasktype

"""Sequence Embedding module in NvTK.
This module provides 

1.  `BasicConvEmbed` class - Basic Convolutional Embedding Module (1d)

2.  `RevCompConvEmbed` class - Convolutional Embedding Module considering Reverse-Complement Sequence 

3.  `CharConvModule` class - Wide and shallow Charactor-level Convolution Module

and supporting methods.
"""

# Code:   jiaqili@zju.edu

import logging

import torch
import torch.nn as nn

__all__ = ["BasicConvEmbed", "RevComp", "RevCompConvEmbed", "CharConvModule"]

class BasicConvEmbed(nn.Module):
    ''' Basic Convolutional Embedding Module in NvTK.
    Embed Sequence using Convolution Layer.

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

    '''
    def __init__(self, out_planes, kernel_size=3, in_planes=4, 
                    conv_args={'stride':1, 'padding':0, 'dilation':1, 'groups':1}, 
                    bn=False, activation=nn.ReLU, activation_args={}, 
                    pool=nn.AvgPool1d, pool_args={'kernel_size': 3}):
        super().__init__()
        self.out_channels = out_planes
        self.conv = nn.Conv1d(in_planes, out_planes, kernel_size=kernel_size, 
                        **conv_args)
        self.bn = nn.BatchNorm1d(out_planes, eps=1e-5, momentum=0.01, affine=True) if bn else None
        self.activation = activation(**activation_args) if activation else None
        self.pool = pool(**pool_args) if pool else None

    def forward(self, x):
        x = self.conv(x)
        if self.bn is not None:
            x = self.bn(x)
        if self.activation is not None:
            x = self.activation(x)
        if self.pool is not None:
            x = self.pool(x)
        logging.debug(x.shape)
        return x


class RevComp(nn.Module):
    """Reverse Complement of onehot Sequence"""
    def forward(self, x):
        if len(x.shape) == 3:
            return x.flip([1, 2]) # (batchsize, 4, seqlen)
        elif len(x.shape) == 3:
            return x.flip([1, -1]) # (batchsize, 4, 1, seqlen)


class RevCompConvEmbed(nn.Module):
    ''' Reverse Complement Convolutional Embedding Module in NvTK.
    Embed Sequence and Reverse Complement Sequence using Convolution Layer.

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
    RevCompConvEmbed : BasicConvEmbed
        Basic Convolutional Embedding Module in NvTK
    RevComp : nn.Module
        Reverse Complement of onehot Sequence

    Tensor flows
    ----------
    -> x1 = RevComp(x) 

    -> x1 = RevCompConvEmbed(x1)
    
    -> x2 = RevCompConvEmbed(x)
    
    -> x1 + x2
        
    '''
    def __init__(self, out_planes, kernel_size=3, in_planes=4, 
                    conv_args={'stride':1, 'padding':0, 'dilation':1, 'groups':1}, 
                    bn=False, activation=nn.ReLU, activation_args={}, 
                    pool=nn.AvgPool1d, pool_args={'kernel_size': 3}):
        super().__init__()
        self.RevCompConvEmbed = BasicConvEmbed(out_planes, kernel_size=kernel_size, in_planes=in_planes, 
                    conv_args=conv_args, 
                    bn=bn, activation=activation, activation_args=activation_args, 
                    pool=pool, pool_args=pool_args)
        self.RevComp = RevComp()

    def forward(self, x):
        fmap1 = self.RevCompConvEmbed(x)
        fmap2 = self.RevCompConvEmbed(self.RevComp(x))
        return fmap1 + fmap2


class CharConvModule(nn.Module):
    '''
    Embed Sequence using wide and shallow CharConvolution Layer.
    '''
    def __init__(self, numFiltersConv1=40, filterLenConv1=5,
                        numFiltersConv2=44, filterLenConv2=15,
                        numFiltersConv3=44, filterLenConv3=25,
                        bn=False, activation=nn.ReLU, activation_args={}, 
                        pool=nn.AvgPool1d, pool_args={'kernel_size': 3}):
        
        super().__init__()
        self.conv1 = nn.Conv1d(in_channels=4, out_channels=numFiltersConv1, 
                                kernel_size=filterLenConv1, padding=(filterLenConv1 - 1) // 2)
        self.conv2 = nn.Conv1d(in_channels=4, out_channels=numFiltersConv2, 
                                kernel_size=filterLenConv2, padding=(filterLenConv2 - 1) // 2)
        self.conv3 = nn.Conv1d(in_channels=4, out_channels=numFiltersConv3, 
                                kernel_size=filterLenConv3, padding=(filterLenConv3 - 1) // 2)

        out_planes = numFiltersConv1 + numFiltersConv2 + numFiltersConv3
        self.bn = nn.BatchNorm1d(out_planes, eps=1e-5, momentum=0.01, affine=True) if bn else None
        self.activation = activation(**activation_args) if activation is not None else None
        self.pool = pool(**pool_args) if pool else None

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

        if self.bn is not None:
            out = self.bn(out)
        if self.activation is not None:
            out = self.activation(out)
        if self.pool is not None:
            out = self.pool(out)
        logging.debug(out.shape)

        return out


"""Convolutional Models in NvTK.
This module provides 

1.  `CNN` class - Convolutional Model in NvTK

2.  `CAN` class - Convolutional Attention Model using CBAM

"""

# Code:   jiaqili@zju.edu

from collections import OrderedDict

import torch
from torch import nn

from .BasicModel import BasicModel
from ..Modules import BasicConvEmbed, RevCompConvEmbed, CharConvModule
from ..Modules import BasicConv1d, Flatten, BasicLinearModule, BasicPredictor
from ..Modules import CBAM


class CNN(BasicModel):
    '''Convolutional Model in NvTK.

    contain Embedding, Encoder, Decoder, Predictor
    '''
    def __init__(self, output_size, 
                    out_planes=128, kernel_size=3, in_planes=4, 
                    conv_args={'stride':1, 'padding':0, 'dilation':1, 'groups':1, 'bias':True}, 
                    bn=False, activation=nn.ReLU, activation_args={}, 
                    pool=nn.AvgPool1d, pool_args={'kernel_size': 3},
                    tasktype='regression'):
        super().__init__()
        self.Embedding = BasicConvEmbed(out_planes=out_planes, 
                    kernel_size=kernel_size, in_planes=in_planes, conv_args=conv_args, 
                    bn=bn, activation=activation, activation_args=activation_args, 
                    pool=pool, pool_args=pool_args)
        self.Encoder = nn.Sequential(OrderedDict([
                        ('Conv', BasicConv1d(in_planes=out_planes, out_planes=256)),
                        ('GAP', nn.AdaptiveAvgPool1d(8)),
                        ('Flatten', Flatten())
                        ]))
        self.Decoder = BasicLinearModule(256 * 8, 256)
        self.Predictor = BasicPredictor(256, output_size, tasktype=tasktype)


class CAN(BasicModel):
    '''
    Covolution Attention Model in NvTK
    contain Embedding, Encoder(CBAM), Decoder, Predictor
    '''
    def __init__(self, output_size, 
                    out_planes=128, kernel_size=3, in_planes=4, 
                    conv_args={'stride':1, 'padding':0, 'dilation':1, 'groups':1, 'bias':False}, 
                    bn=True, activation=nn.ReLU, activation_args={}, 
                    pool=nn.AvgPool1d, pool_args={'kernel_size': 3},
                    tasktype='regression'):
        super().__init__()
        self.Embedding = BasicConvEmbed(out_planes=out_planes, 
                    kernel_size=kernel_size, in_planes=in_planes, conv_args=conv_args, 
                    bn=bn, activation=activation, activation_args=activation_args, 
                    pool=pool, pool_args=pool_args)
        self.Encoder = nn.Sequential(OrderedDict([
                        ('Conv', BasicConv1d(in_planes=128, out_planes=256)),
                        ('Attention', CBAM(256)),
                        ('GAP', nn.AdaptiveAvgPool1d(8)),
                        ('Flatten', Flatten())
                        ]))
        self.Decoder = BasicLinearModule(256 * 8, 256)
        self.Predictor = BasicPredictor(256, output_size, tasktype=tasktype)



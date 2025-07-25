"""Basic Model in NvTK.
This module provides 

1.  `BasicModel` class - the general abstract class

and supporting methods.
"""

# Code: jiaqili@zju.edu

import logging

import torch
from torch import nn

from ..Modules import BasicModule

# TODO maybe not suitable for probability models
class BasicModel(BasicModule):
    """Basic Model class in NvTK.
    Prototype for a sequence-based deep-learning model. 
    BasicModel contains Embedding, Encoder, Decoder, Predictor layers.
    
    Embedding : embed the sequence into vectors

    Encoder : encode the inputs into feature-maps

    Decoder : decode the encoded inputs (Flattened) into higher feature-maps

    Predictor : mapp the decoded feature-maps into task-specific space 
    and make prediction

    Tensor flows
    ------------
    -> Embedding(x)

    -> Encoder(x)

    -> Flatten(x)

    -> Decoder(x)

    -> Predictor(x)

    """
    def __init__(self):
        super().__init__()
        self.Embedding = nn.Sequential()
        self.Encoder = nn.Sequential()
        self.Decoder = nn.Sequential()
        self.Predictor = nn.Sequential()

    def forward(self, x):
        embed = self.Embedding(x)
        logging.debug(embed.shape)

        fmap = self.Encoder(embed)
        logging.debug(fmap.shape)

        if len(fmap.shape) > 2:
            fmap = fmap.reshape((fmap.size(0), -1))
            logging.warning("fmap after Encoder reshaped as (batchsize, -1), \n \
                            Add Flatten module in Encoder to deprecate this warning")

        fmap = self.Decoder(fmap)
        logging.debug(fmap.shape)

        out = self.Predictor(fmap)
        logging.debug(out.shape)

        return out


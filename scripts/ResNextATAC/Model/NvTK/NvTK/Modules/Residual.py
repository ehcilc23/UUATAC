'''
    @inproceedings{li2019selective,
        title={Selective Kernel Networks},
        author={Li, Xiang and Wang, Wenhai and Hu, Xiaolin and Yang, Jian},
        journal={IEEE Conference on Computer Vision and Pattern Recognition},
        year={2019}
    }

    @inproceedings{li2019spatial,
        title={Spatial Group-wise Enhance: Enhancing Semantic Feature Learning in Convolutional Networks},
        author={Li, Xiang and Hu, Xiaolin and Xia, Yan and Yang, Jian},
        journal={arXiv preprint arXiv:1905.09646},
        year={2019}
    }

    @inproceedings{li2019understanding,
        title={Understanding the Disharmony between Weight Normalization Family and Weight Decay: e-shifted L2 Regularizer},
        author={Li, Xiang and Chen, Shuo and Yang, Jian},
        journal={arXiv preprint arXiv:},
        year={2019}
    }

    @inproceedings{li2019generalization,
        title={Generalization Bound Regularizer: A Unified Framework for Understanding Weight Decay},
        author={Li, Xiang and Chen, Shuo and Gong, Chen and Xia, Yan and Yang, Jian},
        journal={arXiv preprint arXiv:},
        year={2019}
    }

'''

# Code:   https://github.com/implus/PytorchInsight/blob/master/classification/models/imagenet/resnet_cbam.py
# Note:   modified for onehot sequence input
#         add attention module using CBAM

__all__  = ["BasicBlock", "Bottleneck", "ResNet", "ResidualNet"]

import torch
import torch.nn as nn
from torch.nn import init
import torch.nn.functional as F

from .BasicModule import BasicModule

class BasicConv(nn.Module):
    def __init__(self, in_planes, out_planes, kernel_size, stride=1, padding=0, 
                        dilation=1, groups=1, relu=True, bn=True, bias=False):
        super(BasicConv, self).__init__()
        self.out_channels = out_planes
        self.conv = nn.Conv2d(in_planes, out_planes, kernel_size=kernel_size, 
                                stride=stride, padding=padding, dilation=dilation,
                                groups=groups, bias=bias)
        self.bn = nn.BatchNorm2d(out_planes,eps=1e-5, momentum=0.01, affine=True) if bn else None
        self.relu = nn.ReLU() if relu else None

    def forward(self, x):
        x = self.conv(x)
        if self.bn is not None:
            x = self.bn(x)
        if self.relu is not None:
            x = self.relu(x)
        return x

class Flatten(nn.Module):
    def forward(self, x):
        return x.view(x.size(0), -1)

class ChannelGate(nn.Module):
    def __init__(self, gate_channels, reduction_ratio=16, pool_types=['avg', 'max']):
        super(ChannelGate, self).__init__()
        self.gate_channels = gate_channels
        self.mlp = nn.Sequential(
            Flatten(),
            nn.Linear(gate_channels, gate_channels // reduction_ratio),
            nn.ReLU(),
            nn.Linear(gate_channels // reduction_ratio, gate_channels)
            )
        self.pool_types = pool_types
        self.avgpool = nn.AdaptiveAvgPool2d(1)
        self.maxpool = nn.AdaptiveMaxPool2d(1)
        self.sigmoid = nn.Sigmoid()

    def forward(self, x):
        channel_att_sum = None
        for pool_type in self.pool_types:
            if pool_type=='avg':
                avg_pool = self.avgpool(x)
                channel_att_raw = self.mlp( avg_pool )
            elif pool_type=='max':
                max_pool = self.maxpool(x)
                channel_att_raw = self.mlp( max_pool )

            if channel_att_sum is None:
                channel_att_sum = channel_att_raw
            else:
                channel_att_sum = channel_att_sum + channel_att_raw

        scale = self.sigmoid(channel_att_sum).unsqueeze(2).unsqueeze(3).expand_as(x)
        return x * scale


class ChannelPool(nn.Module):
    def forward(self, x):
        return torch.cat( (torch.max(x,1)[0].unsqueeze(1), torch.mean(x,1).unsqueeze(1)), dim=1 )

class SpatialGate(nn.Module):
    def __init__(self):
        super(SpatialGate, self).__init__()
        kernel_size = 7
        self.compress = ChannelPool()
        self.spatial = BasicConv(2, 1, (1, kernel_size), stride=1, padding=(0, (kernel_size-1) // 2), relu=False)
        self.sigmoid = nn.Sigmoid()
    def forward(self, x):
        x_compress = self.compress(x)
        x_out = self.spatial(x_compress)
        scale = self.sigmoid(x_out) # broadcasting
        return x * scale

class CBAM(nn.Module):
    def __init__(self, gate_channels, reduction_ratio=16, pool_types=['avg', 'max'], no_spatial=False):
        super(CBAM, self).__init__()
        self.ChannelGate = ChannelGate(gate_channels, reduction_ratio, pool_types)
        self.no_spatial=no_spatial
        if not no_spatial:
            self.SpatialGate = SpatialGate()
    def forward(self, x):
        if len(x.shape) == 3 and x.shape[1] == 4: # (B, 4, L)
            x = x.unsqueeze(2) # (B, 4, 1, L) update 2D sequences
        x_out = self.ChannelGate(x)
        if not self.no_spatial:
            x_out = self.SpatialGate(x_out)
        return x_out

def conv1x3(in_planes, out_planes, stride=1):
    "1x3 convolution with padding"
    return nn.Conv2d(in_planes, out_planes, kernel_size=(1, 3), stride=stride,
                     padding=(0,1), bias=False)

class BasicBlock(nn.Module):
    '''Basic residual Block with conv1x3.

    Parameters
    ----------
    inplanes : int
        Number of input channels
    planes : int
        Number of output channels produced by the convolution
    stride : int, optional
        Number of stride, Default is 1.
    use_cbam : bool, optional
        Whether to use CBAM, Default is False.

    Attributes
    ----------
    conv1 : conv1x3
        The convolutional neural network component of the model.
    bn1 : nn.BatchNorm2d
        The Batch Normalization 
    conv2 : conv1x3
        The convolutional neural network component of the model.
    bn2 : nn.BatchNorm2d
        The Batch Normalization 
    relu : nn.ReLU
        The relu activation Module
    pool : nn.Module
        The pool Module
    cbam : CBAM
        Convolutional Block Attention Module

    Tensor flows
    ----------
    -> residual = x
    
    -> relu(bn1(conv1(x)))
    
    -> bn2(conv2(x))

    -> cbam(x) if use_cbam

    -> relu(x + residual)

    '''
    expansion = 1

    def __init__(self, inplanes, planes, stride=1, downsample=None, use_cbam=False):
        super(BasicBlock, self).__init__()
        self.conv1 = conv1x3(inplanes, planes, stride)
        self.bn1 = nn.BatchNorm2d(planes)
        self.relu = nn.ReLU(inplace=True)
        self.conv2 = conv1x3(planes, planes)
        self.bn2 = nn.BatchNorm2d(planes)
        self.downsample = downsample
        self.stride = stride

        if use_cbam:
            self.cbam = CBAM( planes, 16 )
        else:
            self.cbam = None

    def forward(self, x):
        residual = x

        out = self.conv1(x)
        out = self.bn1(out)+1
        out = self.relu(out)

        out = self.conv2(out)
        out = self.bn2(out)+1

        if self.downsample is not None:
            residual = self.downsample(x)

        if not self.cbam is None:
            out = self.cbam(out)

        out += residual
        out = self.relu(out)

        return out


class Bottleneck(nn.Module):
    '''Basic residual Bottleneck with conv2d.

    Parameters
    ----------
    inplanes : int
        Number of input channels
    planes : int
        Number of output channels produced by the convolution
    stride : int, optional
        Number of stride, Default is 1.
    use_cbam : bool, optional
        Whether to use CBAM, Default is False.

    Attributes
    ----------
    conv1 : nn.Conv2d
        The convolutional neural network component of the model.
    bn1 : nn.BatchNorm2d
        The Batch Normalization 
    conv2 : nn.Conv2d
        The convolutional neural network component of the model.
    bn2 : nn.BatchNorm2d
        The Batch Normalization 
    conv3 : nn.Conv2d
        The convolutional neural network component of the model.
    bn3 : nn.BatchNorm2d
        The Batch Normalization 
    relu : nn.ReLU
        The relu activation Module
    pool : nn.Module
        The pool Module
    cbam : CBAM
        Convolutional Block Attention Module
        
    Tensor flows
    ----------
    -> residual = x
    
    -> relu(bn1(conv1(x)))
    
    -> relu(bn2(conv2(x)))

    -> bn3(conv3(x))

    -> cbam(x) if use_cbam

    -> relu(x + residual)

    '''

    expansion = 4

    def __init__(self, inplanes, planes, stride=1, downsample=None, use_cbam=False):
        super(Bottleneck, self).__init__()
        self.conv1 = nn.Conv2d(inplanes, planes, kernel_size=1, bias=False)
        self.bn1 = nn.BatchNorm2d(planes)
        self.conv2 = nn.Conv2d(planes, planes, kernel_size=(1,3), stride=stride,
                               padding=(0,1), bias=False)
        self.bn2 = nn.BatchNorm2d(planes)
        self.conv3 = nn.Conv2d(planes, planes * 4, kernel_size=1, bias=False)
        self.bn3 = nn.BatchNorm2d(planes * 4)
        self.relu = nn.ReLU(inplace=True)
        self.downsample = downsample
        self.stride = stride

        if use_cbam:
            self.cbam = CBAM( planes * 4, 16 )
        else:
            self.cbam = None

    def forward(self, x):
        residual = x

        out = self.conv1(x)
        out = self.bn1(out)+1
        out = self.relu(out)

        out = self.conv2(out)
        out = self.bn2(out)+1
        out = self.relu(out)

        out = self.conv3(out)
        out = self.bn3(out)

        if self.downsample is not None:
            residual = self.downsample(x)

        if not self.cbam is None:
            out = self.cbam(out)

        out += residual
        out = self.relu(out)

        return out


class ResNet(BasicModule):
    '''ResNet module in NvTK.

    Parameters
    ----------
    block : BasicBlock, Bottleneck
        One of the Residual Block Module
    layers : list of int
        List of Number of output channels in ResNet.
        (e.g. `layers=[2, 2, 2, 2]`)
    network_type : str
        Type of ResNet, Default is None.
        Should be one of `'ImageNet', 'CIFAR10', 'CIFAR100', None`
        (e.g. `network_type='ImageNet'`)
    num_classes : int, optional
        Number of output classes
    att_type : str, optional
        Whether to use CBAM_ResNet, Default is None.

    Attributes
    ----------
    inplanes : int

    network_type : str

    conv1 : nn.Conv2d
        The convolutional neural network component of the model.
    bn1 : nn.BatchNorm2d
        The Batch Normalization 
    conv2 : nn.Conv2d
        The convolutional neural network component of the model.
    layer1 : nn.Sequential
        The Residual convolutional layer
    bam1 : CBAM
        The Convolutional Block Attention Module
    layer2 : nn.Sequential
        The Residual convolutional layer
    bam2 : CBAM
        The Convolutional Block Attention Module
    layer3 : nn.Sequential
        The Residual convolutional layer
    bam3 : CBAM
        The Convolutional Block Attention Module
    layer4 : nn.Sequential
        The Residual convolutional layer
    fc : nn.Linear
        The Full connectted layer
        
    Tensor flows
    ----------
    -> x if network_type is None

    -> relu(bn1(conv1(x))) if network_type is CIFAR
    
    -> maxpool(relu(bn1(conv1(x)))) if network_type is ImageNet
    
    -> layer1(x)
    
    -> bam1(x) if bam1

    -> layer2(x)
    
    -> bam2(x) if bam2

    -> layer3(x)
    
    -> bam3(x) if bam3

    -> layer4(x)
    
    -> Faltten(x)

    -> fc(x) if fc

    '''
    
    def __init__(self, block, layers, network_type=None, num_classes=None, att_type=None):
        self.inplanes = 64
        super(ResNet, self).__init__()
        self.network_type = network_type
        # different model config between ImageNet and CIFAR 
        if network_type is None:
            self.conv1 = None
        elif network_type == "ImageNet":
            self.conv1 = nn.Conv2d(4, 64, kernel_size=(1, 7), stride=1, padding=(0,3), bias=False)
            self.maxpool = nn.MaxPool2d(kernel_size=(1, 3), stride=1, padding=(0,1))
            self.avgpool = nn.AdaptiveAvgPool2d(1)
        elif network_type == "CIFAR10" or network_type == "CIFAR100" :
            self.conv1 = nn.Conv2d(4, 64, kernel_size=(1, 3), stride=1, padding=(0,1), bias=False)

        self.bn1 = nn.BatchNorm2d(64)
        self.relu = nn.ReLU(inplace=True)

        if att_type=='CBAM':
            self.bam1 = CBAM(64*block.expansion)
            self.bam2 = CBAM(128*block.expansion)
            self.bam3 = CBAM(256*block.expansion)
        else:
            self.bam1, self.bam2, self.bam3 = None, None, None

        self.layer1 = self._make_layer(block, 64,  layers[0], att_type=att_type)
        self.layer2 = self._make_layer(block, 128, layers[1], stride=1, att_type=att_type)
        self.layer3 = self._make_layer(block, 256, layers[2], stride=1, att_type=att_type)
        self.layer4 = self._make_layer(block, 512, layers[3], stride=1, att_type=att_type)
        
        if num_classes is not None:
            self.fc = nn.Linear(512 * block.expansion, num_classes)
        else:
            self.fc = None

        init.kaiming_normal(self.fc.weight)
        for key in self.state_dict():
            if key.split('.')[-1]=="weight":
                if "conv" in key:
                    init.kaiming_normal(self.state_dict()[key], mode='fan_out')
                if "bn" in key:
                    if "SpatialGate" in key:
                        self.state_dict()[key][...] = 0
                    else:
                        self.state_dict()[key][...] = 1
            elif key.split(".")[-1]=='bias':
                self.state_dict()[key][...] = 0

    def _make_layer(self, block, planes, blocks, stride=1, att_type=None):
        downsample = None
        if stride != 1 or self.inplanes != planes * block.expansion:
            downsample = nn.Sequential(
                nn.Conv2d(self.inplanes, planes * block.expansion,
                          kernel_size=1, stride=stride, bias=False),
                nn.BatchNorm2d(planes * block.expansion),
            )

        layers = []
        layers.append(block(self.inplanes, planes, stride, downsample, use_cbam=att_type=='CBAM'))
        self.inplanes = planes * block.expansion
        for i in range(1, blocks):
            layers.append(block(self.inplanes, planes, use_cbam=att_type=='CBAM'))

        return nn.Sequential(*layers)

    def forward(self, x):
        if self.network_type is not None: # which should embed sequence directly
            if len(x.shape) == 3: # (B, 4, L)
                x = x.unsqueeze(2) # (B, 4, 1, L) update 2D sequences
            
            x = self.conv1(x)
            x = self.bn1(x)
            x = self.relu(x)

        if self.network_type == "ImageNet":
            x = self.maxpool(x)

        x = self.layer1(x)
        if not self.bam1 is None:
            x = self.bam1(x)

        x = self.layer2(x)
        if not self.bam2 is None:
            x = self.bam2(x)

        x = self.layer3(x)
        if not self.bam3 is None:
            x = self.bam3(x)

        x = self.layer4(x)

        if self.network_type == "ImageNet":
            x = self.avgpool(x)
        else:
            x = F.avg_pool2d(x, 4)
        x = x.view(x.size(0), -1)

        if not self.fc is None:
            x = self.fc(x)

        return x
        

def ResidualNet(network_type, depth, num_classes, att_type):
    assert network_type in ["ImageNet", "CIFAR10", "CIFAR100"], "network type should be ImageNet or CIFAR10 / CIFAR100"
    assert depth in [18, 34, 50, 101], 'network depth should be 18, 34, 50 or 101'

    if depth == 18:
        model = ResNet(BasicBlock, [2, 2, 2, 2], network_type, num_classes, att_type)

    elif depth == 34:
        model = ResNet(BasicBlock, [3, 4, 6, 3], network_type, num_classes, att_type)

    elif depth == 50:
        model = ResNet(Bottleneck, [3, 4, 6, 3], network_type, num_classes, att_type)

    elif depth == 101:
        model = ResNet(Bottleneck, [3, 4, 23, 3], network_type, num_classes, att_type)

    return model


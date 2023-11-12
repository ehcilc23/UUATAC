import torch
from torch import nn

from NvTK.Model import BasicModel
from NvTK.Modules import BasicConv1d, BasicConvEmbed, Flatten, BasicLinearModule, BasicPredictor
from NvTK.Modules import CBAM, CharConvModule

from collections import OrderedDict


class BaselineCNN(BasicModel):
    '''
    Baseline CNN (Default two-layer CNN; 
        Embedding: BasicConvEmbed(N128_L15_P15_Pa_ReLU_BN)
        Encoder: BasicConv1d(N256|CBAM) * n_deep_convlayers; GAP{64|Expr, 1|Epig}; Flatten
        Decoder: MTLnone
        Predictor: regression/binary_classification)
    BasicConvEmbedding could be tuned afterward.
    '''
    def __init__(self, output_size, 
                    emb_planes=128, kernel_size=15, in_planes=4, 
                    conv_args={'stride':1, 'padding':0, 'dilation':1, 'groups':1}, 
                    bn=False, activation=nn.ReLU, activation_args={}, 
                    pool=nn.AvgPool1d, pool_args={'kernel_size': 15},
                    n_deep_convlayers=1, use_CBAM=False, GAP=64, tasktype='regression'):
        super().__init__()
        self.Embedding = BasicConvEmbed(out_planes=emb_planes, 
                    kernel_size=kernel_size, in_planes=in_planes, conv_args=conv_args, 
                    bn=bn, activation=activation, activation_args=activation_args, 
                    pool=pool, pool_args=pool_args)

        encoder_layers = OrderedDict([('Conv', BasicConv1d(in_planes=emb_planes, out_planes=256))])
        if use_CBAM:
            encoder_layers['CBAM'] = CBAM(256)
        if n_deep_convlayers > 1:
            for n in range(1, n_deep_convlayers):
                encoder_layers['Conv_'+str(n)] = BasicConv1d(256, 256, pool=False)
                if use_CBAM:
                    encoder_layers['CBAM_'+str(n)] = CBAM(256)
        encoder_layers['GAP'] = nn.AdaptiveAvgPool1d(GAP) # (batch_size, 256, GAP)
        encoder_layers['Flatten'] = Flatten() # (batch_size, 256*GAP)
        self.Encoder = nn.Sequential(encoder_layers)

        self.Decoder = BasicLinearModule(256 * GAP, 256)
        self.Predictor = BasicPredictor(256, output_size, tasktype=tasktype)


def switch_to_charCNN(model):
    model.Embedding = CharConvModule(numFiltersConv1=44, filterLenConv1=5,
                                    numFiltersConv2=40, filterLenConv2=15,
                                    numFiltersConv3=44, filterLenConv3=25)
    return model


def get_resnet(output_size, layers=18, tasktype='regression'):
    assert layers in [18, 50, 101]
    from NvTK import resnet18, resnet50, resnet101
    if layers == 18:
        model = resnet18(output_size)
    elif layers == 50:
        model = resnet50(output_size)
    elif layers == 101:
        model = resnet101(output_size)
    model.fc = BasicPredictor(512, output_size, tasktype=tasktype)
    return model


def get_transformer(input_size, output_size, tasktype='regression'):
    from NvTK import TransformerEncoder
    model = BasicModel()
    model.Embedding = BasicConvEmbed(out_planes=128, kernel_size=15, conv_args={'padding':7}, 
                                    activation=nn.LeakyReLU, activation_args={"negative_slope":0.2}, 
                                    pool=nn.AvgPool1d, pool_args={'kernel_size': 2})
    data = torch.zeros(input_size)
    _, d_model, seq_len = model.Embedding.forward(data).shape
    model.Encoder = TransformerEncoder(seq_len, d_model=d_model, n_layers=2, n_heads=2, d_ff=32, embedding=nn.Sequential())
    # model.Decoder = BasicLinearModule(d_model, 128)
    model.Predictor = BasicPredictor(d_model, output_size, tasktype=tasktype)
    return model


# test 
# for channel in (8, 32, 128, 512):
#     model = BaselineExprCNN(output_size=100, out_planes=channel)
#     print(model)


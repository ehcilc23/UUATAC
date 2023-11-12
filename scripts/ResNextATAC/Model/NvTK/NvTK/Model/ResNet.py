from ..Modules import ResidualNet

def cbam_resnet18(n_features=1000, **kwargs):
    model = ResidualNet('ImageNet', 18, n_features, 'CBAM')
    return model

def cbam_resnet34(n_features=1000, **kwargs):
    model = ResidualNet('ImageNet', 34, n_features, 'CBAM')
    return model

def cbam_resnet50(n_features=1000, **kwargs):
    model = ResidualNet('ImageNet', 50, n_features, 'CBAM')
    return model

def cbam_resnet101(n_features=1000, **kwargs):
    model = ResidualNet('ImageNet', 101, n_features, 'CBAM')
    return model

def resnet18(n_features=1000, **kwargs):
    model = ResidualNet('ImageNet', 18, n_features, None)
    return model

def resnet34(n_features=1000, **kwargs):
    model = ResidualNet('ImageNet', 34, n_features, None)
    return model

def resnet50(n_features=1000, **kwargs):
    model = ResidualNet('ImageNet', 50, n_features, None)
    return model

def resnet101(n_features=1000, **kwargs):
    model = ResidualNet('ImageNet', 101, n_features, None)
    return model
    
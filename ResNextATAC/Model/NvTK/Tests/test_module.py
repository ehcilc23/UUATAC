import torch

from NvTK.Modules import *

x = torch.rand((2, 32, 100)) # (batchsize, channel, seqlen)
module = BasicConv1d(32, 128, 3)
o = module.forward(x)
o.shape

x = torch.rand((2, 100, 32)) # batch, seq, feature
module = BasicRNNModule(32, 128, 2)
o = module.forward(x)
assert o.shape == (2, 100, 128*2) # batch, seq, feature*2

x = torch.rand((2, 100)) # (batchsize, -1)
module = BasicLinearModule(100, 5)
o = module.forward(x)
assert o.shape == (2, 5)

x = torch.rand((2, 5)) # (batchsize, 10)
module = BaiscPredictor()
o = module.forward(x)
assert o.shape == (2, 5)

module.switch_task('classification')
o = module.forward(x)
assert o.shape == (2, 5)

module.switch_task('regression')
o = module.forward(x)
assert o.shape == (2, 5)

module.switch_task('none')
o = module.forward(x)
assert o.shape == (2, 5)

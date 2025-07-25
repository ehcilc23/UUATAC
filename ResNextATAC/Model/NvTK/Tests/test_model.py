import torch

from NvTK.Model import *

x = torch.rand((2, 4, 1000)) # (batchsize, 4, seqlen)

model = CNN(10)
o = model(x)
assert o.shape == (2, 10)

model = CAN(10)
o = model(x)
assert o.shape == (2, 10)

model = resnet18(10)
o = model(x)
assert o.shape == (2, 10)

model = NINCNN(1000, 10)
o = model(x)
assert o.shape == (2, 10)

model = DeepSEA(1000, 10)
o = model(x)
assert o.shape == (2, 10)

model = Beluga(1000, 10)
o = model(x)
assert o.shape == (2, 10)

model = DanQ(1000, 10)
o = model(x)
assert o.shape == (2, 10)

model = Basset(1000, 10)
o = model(x)
assert o.shape == (2, 10)

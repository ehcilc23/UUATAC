# NvTK

Source code used for ```Systematic evaluation of deep learning for mapping sequence to single-cell data using NvTK```

NvTK (NvwaToolKit), is a systemmatic and easy-using deep learning software in genomics. NvTK support modern deep learning achitectures in genomics, such as Residual Module, ResNet, Attention Module, CBAM, Transformer and so on. 


## Requirements
- Python packages
```
python>=3.7
numpy
pandas>=0.21
matplotlib==3.0.*
# h5py > 2.10 may returns b'strings' reading h5file
h5py=2.10.0
tqdm
scikit-learn>=0.21.2
# torch >=1.10.1 support tensorboard, and ModifyOutputHook
torch>=1.10.1
tensorboard=2.7.0
captum=0.5.0
networkx
# higher version of scipy do not support resize
pillows
```

- external softwares (optional)
```
# meme could align the deep learning motif with known TFBS database
meme-5.4.1
# homer2 could search the motif in activated seqlets
homer2
```
<!-- biopython-1.79 -->


"""Feature map based model interpretation methods in NvTK.
"""

import torch
import logging
import numpy as np

from .Motif import normalize_pwm, threshold_nrom

__all__ = ["get_activate_W", "get_fmap", "get_activate_W_from_fmap",
    "get_activate_sequence_from_fmap", "save_activate_seqlets"]

def _get_W_from_conv(model, motif_width=5, normalize=True, device=torch.device("cuda")):
    '''
    Experimental function!
    get motif directly from convolution parameters, 
    PWM were extracted from `model.Embedding.conv`
    '''
    x_tensor = torch.zeros((4, 4, motif_width)).to(device)
    x_tensor[0,0,:] = 1
    x_tensor[1,1,:] = 1
    x_tensor[2,2,:] = 1
    x_tensor[3,3,:] = 1

    try:
        fmap = model.Embedding.conv(x_tensor).data.cpu().numpy()
    except AttributeError:
        logging.error("Check if you model have model.Embedding.conv attr?")
        raise AttributeError

    W = fmap.swapaxes(0, 1).clip(0)
    if normalize:
        W = np.array([normalize_pwm(pwm) for pwm in W])
    return W


# hook
class ActivateFeaturesHook():
    def __init__(self, module):
        self.hook = module.register_forward_hook(self.hook_fn)
    def hook_fn(self, module, input, output):
        self.features = output.cpu().data.numpy()#.mean(-1)
    def get_features(self):
        return self.features
    def close(self):
        self.hook.remove()


def get_fmap(model, hook_module, data_loader, device=torch.device("cuda")):
    """Get feature map of input data at model.hook_module

    Parameters
    ----------
    model : 
        model
    hook_module : int
        hook_module
    data_loader : torch.Data.Dataloader
        input data
    device : torch.device, optional
        torch.device, Default is `torch.device("cuda")`.

    Returns
    ----------
    fmap : np.ndarr
        feature map of input data at model.hook_module
    X : np.ndarr
        input data
    """
    fmap, X = [], []
    model.eval()
    with torch.no_grad():
        activations = ActivateFeaturesHook(hook_module)
        for x_tensor, _ in data_loader:
            x_tensor = x_tensor.to(device)
            _ = model(x_tensor)
            X.append(x_tensor.cpu().numpy())
            fmap.append(activations.get_features())
        fmap = np.vstack(fmap)
        X = np.vstack(X)
        activations.close()
    return fmap, X


def get_activate_W_from_fmap(fmap, X, pool=1, threshold=0.99, motif_width=10, pad=0, axis=1):
    """Get activated motif pwm from feature map

    Parameters
    ----------
    fmap : np.ndarr
        feature map of input data at model.hook_module
    X : np.ndarr
        input data
    pool : int
        input data
    threshold : floor
        threshold determine the activated sites in feature map
    motif_width : int
        width of motif, the width region sequence of activated sites 
        will be normalized as counts

    Returns
    ----------
    W : np.ndarr
        array of activated motif pwm, 
        shape of W (n_filters, 4, motif_width)
    """

    motif_nb = fmap.shape[1]
    X_dim, seq_len = X.shape[1], X.shape[-1]

    W=[]
    for filter_index in range(motif_nb):
        # find regions above threshold
        data_index, pos_index = np.where(fmap[:,filter_index,:] > np.max(fmap[:,filter_index,:], axis=axis, keepdims=True)*(threshold*threshold_nrom))

        seq_align = []; count_matrix = []
        for i in range(len(pos_index)):
            # pad 1-nt
            start = pos_index[i] - 1 - pad
            end = start + motif_width + 2
            # handle boundary conditions
            if end > seq_len:
                end = seq_len
                start = end - motif_width - 2 
            if start < 0:
                start = 0 
                end = start + motif_width + 2 

            seq = X[data_index[i], :, start*pool:end*pool]
            seq_align.append(seq)
            count_matrix.append(np.sum(seq, axis=0, keepdims=True))

        seq_align = np.array(seq_align)
        count_matrix = np.array(count_matrix)

        # normalize counts
        seq_align = (np.sum(seq_align, axis=0)/np.sum(count_matrix, axis=0))*np.ones((X_dim, (motif_width+2)*pool))
        seq_align[np.isnan(seq_align)] = 0
        W.append(seq_align)

    W = np.array(W)
    return W


def get_activate_W(model, hook_module, data, pool=1, pad=0, threshold=0.99, motif_width=20, axis=1):
    """Get activated motif pwm of input data at model.hook_module

    Parameters
    ----------
    model : 
        model
    hook_module : int
        hook_module
    data_loader : torch.Data.Dataloader
        input data
    device : torch.device, optional
        torch.device, Default is `torch.device("cuda")`.
    pool : int
        input data
    threshold : floor
        threshold determine the activated sites in feature map
    motif_width : int
        width of motif, the width region sequence of activated sites 
        will be normalized as counts

    Returns
    ----------
    W : np.ndarr
        array of activated motif pwm, 
        shape of W (n_filters, 4, motif_width)
    """
    fmap, X = get_fmap(model, hook_module, data)
    W = get_activate_W_from_fmap(fmap, X, pool, threshold, motif_width, pad=pad, axis=axis)
    return W


def onehot2seq(gene_seq, gene_name, out_fname):
    d = {0:'A', 1:'C', 2:'G', 3:'T'}
    s = ''
    for i, fas in zip(gene_name, map(lambda y: ''.join(map(lambda x:d[x], np.where(y.T==1)[-1])), gene_seq)):
        s += '>'+str(i)+'\n'
        s += fas+'\n'
    with open(out_fname, 'w') as fh:
        fh.write(s)


def get_activate_sequence_from_fmap(fmap, X, pool=1, threshold=0.99, motif_width=40):
    """Get activated sequence from feature map.
    Seqlets could be further analyzed by bioinformatic softwares, 
    such as Homer2.

    Parameters
    ----------
    fmap : np.ndarr
        feature map of input data at model.hook_module
    X : np.ndarr
        input data
    pool : int
        input data
    threshold : floor
        threshold determine the activated sites in feature map
    motif_width : int
        width of motif, the width region sequence of activated sites 
        will be normalized as counts

    Returns
    ----------
    W : list
        list of activated motif seqlets, 
        shape of W (n_filters, 4, motif_width)
    M : list
        Seqlet Names, defined as "Motif_Act"
    """

    motif_nb = fmap.shape[1]
    seq_len = X.shape[-1]

    W, M = [], []
    for filter_index in range(motif_nb):
        # find regions above threshold
        data_index, pos_index = np.where(fmap[:,filter_index,:] > np.max(fmap[:,filter_index,:], axis=1, keepdims=True)*threshold)

        for i in range(len(pos_index)):
            # handle boundary conditions
            start = pos_index[i] - 1
            end = pos_index[i] + motif_width + 2
            if end > seq_len:
                end = seq_len
                start= end - motif_width - 2 
            if start < 0:
                start = 0 
                end = start + motif_width + 2

            seq = X[data_index[i], :, start*pool:end*pool]
            W.append(seq)
            M.append('_'.join(("Motif", str(filter_index), "Act", str(i))))

    return W, M


def save_activate_seqlets(model, hook_module, data, out_fname, pool=1, threshold=0.99, motif_width=40):
    """Save activated Seqlets pwm from feature map
    Seqlets could be further analyzed by bioinformatic softwares, 
    such as Homer2.

    Parameters
    ----------
    model : 
        model
    hook_module : int
        hook_module
    data_loader : torch.Data.Dataloader
        input data
    out_fname : str
        output file name
    device : torch.device, optional
        torch.device, Default is `torch.device("cuda")`.
    pool : int
        input data
    threshold : floor
        threshold determine the activated sites in feature map
    motif_width : int
        width of motif, the width region sequence of activated sites 
        will be normalized as counts
    """

    fmap, X = get_fmap(model, hook_module, data)
    gene_seq, gene_name = get_activate_sequence_from_fmap(fmap, X, pool=pool, threshold=threshold, motif_width=motif_width)
    onehot2seq(gene_seq, gene_name, out_fname)


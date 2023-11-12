import numpy as np

def seq_rc_augment(X, y):
    X_rc = np.array([seq[::-1,::-1] for seq in X]) # reverse_comp
    y_rc = y # keep same label
    return np.vstack((X_rc, X)), np.vstack((y_rc, y))

def shift_sequence(seq, shift, pad_value=0.25):
    """Shift a sequence left or right by shift_amount.

    Args:
    seq: [batch_size, seq_depth, seq_length] sequence
    shift: signed shift value (int)
    pad_value: value to fill the padding (primitive or scalar)
    """
    if len(seq.shape) != 3:
        raise ValueError('input sequence should be rank 3')
    input_shape = seq.shape

    pad = pad_value * np.ones_like(seq[:, :, 0:np.abs(shift)])

    def _shift_right(_seq):
        # shift is positive
        sliced_seq = _seq[:, :, :-shift:]
        return np.concatenate([pad, sliced_seq], axis=-1)

    def _shift_left(_seq):
        # shift is negative
        sliced_seq = _seq[:, :, -shift:]
        return np.concatenate([sliced_seq, pad], axis=-1)

    if shift > 0:
        sseq = _shift_right(seq)
    else:
        sseq = _shift_left(seq)

    return sseq

def seq_shift_augment(X, y, shift=None):
    if shift is None:
        shift = X.shape[-1] // 4
        
    X_rc_left = shift_sequence(X, shift) # reverse_comp
    X_rc_right = shift_sequence(X, -shift) # reverse_comp

    return np.vstack((X_rc_left, X_rc_right, X)), np.vstack([y]*3)


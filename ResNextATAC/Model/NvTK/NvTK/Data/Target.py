

def onehot_encode(label):
    from sklearn.preprocessing import label_binarize
    return label_binarize(label, classes=range(np.max(label)+1))

def map_prob2label(y_pred_prob, map_fn=np.argmax):
    assert isinstance(y_pred_prob, np.ndarray)
    return np.array(list(map(map_fn, y_pred_prob)))


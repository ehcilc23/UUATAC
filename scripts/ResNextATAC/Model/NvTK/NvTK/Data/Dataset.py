import h5py
import numpy as np

from torch.utils.data import DataLoader

def unpack_datasets(fname):
    # unpack datasets
    h5file = h5py.File(fname, 'r')
    anno = h5file["annotation"][:]
    x_train = h5file["train_data"][:].astype(np.float32)
    y_train = h5file["train_label"][:].astype(np.float32)
    x_val = h5file["val_data"][:].astype(np.float32)
    y_val = h5file["val_label"][:].astype(np.float32)
    x_test = h5file["test_data"][:].astype(np.float32)
    y_test = h5file["test_label"][:].astype(np.float32)
    h5file.close()

    return anno, x_train, y_train, x_val, y_val, x_test, y_test


def generate_dataloader(x_train, y_train, x_val, y_val, x_test, y_test, batch_size = 16):
    # define data loader
    train_loader = DataLoader(list(zip(x_train, y_train)), batch_size=batch_size,
                                shuffle=True, num_workers=2, drop_last=False, pin_memory=True)
    validate_loader = DataLoader(list(zip(x_val, y_val)), batch_size=batch_size, 
                                shuffle=False, num_workers=2, drop_last=False, pin_memory=True)
    test_loader = DataLoader(list(zip(x_test, y_test)), batch_size=batch_size, 
                                shuffle=False, num_workers=2, drop_last=False, pin_memory=True)

    return train_loader, validate_loader, test_loader


def generate_dataloader_from_datasets(fname, batch_size = 16):
    anno, x_train, y_train, x_val, y_val, x_test, y_test = unpack_datasets(fname)
    train_loader, validate_loader, test_loader = generate_dataloader(x_train, y_train, x_val, y_val, x_test, y_test, batch_size=batch_size)
    return train_loader, validate_loader, test_loader

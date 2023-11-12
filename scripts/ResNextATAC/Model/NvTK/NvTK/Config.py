# TODO whether to define pipeline module like sklearn
import json
import logging
from copy import deepcopy

import NvTK
from NvTK.Data.Dataset import generate_dataloader_from_datasets

from torch import nn
from torch import optim

def load_config_from_json(fname):
    with open(fname) as f:
        config = json.load(f)
    return config


def dump_config_to_json(config, fname='config_dump.json'):
    with open(fname, 'w') as f:
        f.write(json.dumps(config, indent=4))


def generate_dataloader_from_config(config):
    data_type = config['data'].get('type', 'h5Dataset')
    fname = config['data']['fpath']
    train_loader, validate_loader, test_loader = generate_dataloader_from_datasets(fname, batch_size = 16)
    return train_loader, validate_loader, test_loader


def parse_model_args(args):
    '''help-function of get_model_from_config'''
    args = deepcopy(args)
    for k, v in args.items():
        if k in ["pool", "activation"]:
            if hasattr(nn, v):
                args[k] = getattr(nn, v)
            elif hasattr(NvTK, v):
                args[k] = getattr(NvTK, v)
            else:
                logging.error("model args[%s]=%s not valid!"%(k,v))
                raise ValueError
    return args

def get_model_from_config(config):
    model_type = config['model'].get('type', 'CNN')
    model_args = config['model'].get('args', None)
    model_args = parse_model_args(model_args) if model_args else dict()

    if hasattr(NvTK, model_type):
        model = getattr(NvTK, model_type)(**model_args)
    return model


def get_optimizer_from_config(config, model):
    if 'optimizer' in config:
        optimizer_type = config['optimizer']['type']
        args = config['optimizer'].get('args', {"lr":1e-4})
    else:
        optimizer_type = 'Adam'
        args = {"lr":1e-4}

    if hasattr(optim, optimizer_type):
        optimizer = getattr(optim, optimizer_type)(model.parameters(), **args)
    return optimizer


def get_criterion_from_config(config):
    if 'criterion' in config:
        criterion_type = config['criterion']['type']
        args = config['criterion'].get('args', {})
    if 'tasktype' in config:
        tasktype = config['tasktype']

    if tasktype == 'regression':
        criterion_type = 'MSELoss'
        args = {}
    else:
        criterion_type = 'BCELoss'
        args = {}

    if hasattr(nn, criterion_type):
        criterion = getattr(nn, criterion_type)(**args)
    elif hasattr(NvTK.Modules.Loss, criterion_type):
        criterion = getattr(NvTK.Modules.Loss, criterion_type)(**args)
    return criterion


def parse_trainer_args(config):
    trainer_args = {}
    if "trainer" in config:
        if "args" in config["trainer"]:
            trainer_args = config["trainer"]["args"]
    return trainer_args


def parse_modes_from_config(config):
    if 'modes' in config:
        return [k for k in config['modes']]
    else:
        return ["hpo", "train", "evaluate", "explain"]

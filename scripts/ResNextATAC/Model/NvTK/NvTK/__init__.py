from .Modules import *

from .Model import *

from .Trainer import Trainer
from .Evaluator import *
from .Explainer import *

# from .Config import *
from .nvtk import *

__author__ = "Jiaqili@zju.edu.cn"

__email__ = "Jiaqili@zju.edu.cn"

__version__ = "0.0"


def set_random_seed(random_seed = 12):
    '''set random_seed'''
    import random
    import numpy as np
    random.seed(random_seed)
    np.random.seed(random_seed)


def set_torch_seed(random_seed = 12):
    '''set torch random_seed'''
    import torch
    torch.random.seed = random_seed
    torch.manual_seed(random_seed)
    torch.cuda.manual_seed(random_seed) 
    torch.cuda.manual_seed_all(random_seed)


def set_torch_benchmark():
    '''set torch benchmark'''
    torch.backends.cudnn.benchmark = True
    torch.backends.cudnn.deterministic = True


def initialize_logger(output_path, verbosity=1):
    """
    Initializes the logger for NvTK.
    This function can only be called successfully once.
    If the logger has already been initialized with handlers,
    the function exits. Otherwise, it proceeds to set the
    logger configurations.

    Parameters
    ----------
    output_path : str
        The path to the output file where logs will be written.

    verbosity : int, {2, 1, 0}
        Default is 1. The level of logging verbosity to use.

            * 0 - Only warnings will be logged.
            * 1 - Information and warnings will be logged.
            * 2 - Debug messages, information, and warnings will all be logged.

    """
    logger = logging.getLogger("nvtk")
    # check if logger has already been initialized
    if len(logger.handlers):
        return

    if verbosity == 0:
        logger.setLevel(logging.WARN)
    elif verbosity == 1:
        logger.setLevel(logging.INFO)
    elif verbosity == 2:
        logger.setLevel(logging.DEBUG)

    file_formatter = logging.Formatter(
        "%(asctime)s - %(levelname)s - %(message)s")

    file_handle = logging.FileHandler(output_path)
    file_handle.setFormatter(file_formatter)
    logger.addHandler(file_handle)

    stdout_formatter = logging.Formatter(
        "%(asctime)s - %(message)s")

    stdout_handle = logging.StreamHandler(sys.stdout)
    stdout_handle.setFormatter(stdout_formatter)
    stdout_handle.setLevel(logging.INFO)
    logger.addHandler(stdout_handle)
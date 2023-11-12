"""
This module provides 

1.  `NvTKCLI` class - the client class

2.  `execute_nvtk_hpo` function - execute Hyper Parameter Optimization in NvTK

3.  `execute_nvtk_train` function - execute model training in NvTK

4.  `execute_nvtk_evaluate` function - execute model evaluation in NvTK

5.  `execute_nvtk_explain` function - execute model explaination in NvTK

and supporting methods.
"""

import sys
import torch

from .Config import *
from .Trainer import Trainer

class NvTKCLI(object):
    def __init__(self, config):
        super().__init__()

        if isinstance(config, str):
            self.config = load_config_from_json(config)
        else:
            self.config = config
        self.modes = parse_modes_from_config(self.config)
        self.device = torch.device("cuda") # TODO

        self.train_loader, self.validate_loader, self.test_loader = generate_dataloader_from_config(self.config)

        self.model = get_model_from_config(self.config)
        self.optimizer = get_optimizer_from_config(self.config, self.model)
        self.criterion = get_criterion_from_config(self.config)

        trainer_args = parse_trainer_args(self.config)
        self.trainer = Trainer(self.model, 
                            self.criterion, 
                            self.optimizer, 
                            self.device, **trainer_args)

    def execute(self):
        '''execute nvtk'''
        for mode in self.modes:
            assert mode in ["hpo", "train", "evaluate", "explain"]

        for mode in self.modes:
            if mode == "hpo":
                execute_nvtk_hpo()
            if mode == "train":
                self.model, self.trainer = execute_nvtk_train(self.trainer, 
                                            self.train_loader, self.validate_loader, self.test_loader)
            elif mode == "evaluate":
                execute_nvtk_evaluate(self.trainer, self.test_loader)
            elif mode == "explain":
                motif_width =30
                execute_nvtk_explain(self.model, self.test_loader, motif_width)


def execute_nvtk_hpo(search_space=None):
    from .Architecture import hyperparameter_tune
    hyperparameter_tune(search_space=search_space, num_samples=10, max_num_epochs=10, gpus_per_trial=1)


def execute_nvtk_train(trainer, train_loader, validate_loader, test_loader):
    # train
    trainer.train_until_converge(train_loader, validate_loader, test_loader, EPOCH=500)
    # reload best model
    model = trainer.load_best_model()
    return model, trainer


def execute_nvtk_evaluate(trainer, test_loader, task_name=None):
    import pandas as pd
    from .Evaluator import calculate_roc, calculate_pr

    # reload best model
    model = trainer.load_best_model()

    # predict test-set
    _, _, test_predictions, test_targets = trainer.predict(test_loader)
    # metric test-set
    fpr, tpr, roc_auc = calculate_roc(test_targets, test_predictions)
    auroc = [roc_auc[k] for k in roc_auc.keys() if k not in ["macro", "micro"]] # dict keys ordered by default in py3.7+

    p, r, average_precision = calculate_pr(test_targets, test_predictions)
    aupr = [average_precision[k] for k in average_precision.keys() if k not in ["macro", "micro"]] # dict keys ordered by default in py3.7+
    return pd.DataFrame({"auroc":auroc, "aupr":aupr}, index=task_name).T.to_csv("Metric.csv")


def execute_nvtk_explain(model, test_loader, motif_width):
    from .Explainer import get_activate_W, meme_generate, save_activate_seqlets
    # explain based on feature-map
    W = get_activate_W(model, model.Embedding.conv, test_loader, motif_width=motif_width)
    meme_generate(W, output_file='meme.txt', prefix='Filter_')

    save_activate_seqlets(model, model.Embedding.conv, test_loader, threshold=0.999,
                            out_fname='seqlets.fasta', motif_width=motif_width)



if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("data")
    parser.add_argument("--config", dest="config", default=None, type=str)
    args = parser.parse_args()

    if args.config is not None:
        nvtk = NvTKCLI(args.config)
        nvtk.execute()

'''
@article{liaw2018tune,
    title={Tune: A Research Platform for Distributed Model Selection and Training},
    author={Liaw, Richard and Liang, Eric and Nishihara, Robert
            and Moritz, Philipp and Gonzalez, Joseph E and Stoica, Ion},
    journal={arXiv preprint arXiv:1807.05118},
    year={2018}
}
'''
__all__ = ["hyperparameter_tune"]

import torch
from torch import nn
from torch.optim import Adam

from ray import tune
from ray.tune.schedulers import ASHAScheduler

from .Model import CNN
from .Trainer import Trainer

# wrap data loading and training in functions,
# make some network parameters configurable,
# add checkpointing (optional),
# and define the search space for the model tuning

def train_model(config, train_loader=None, validate_loader=None):

    # TODO support more NvTK architectures
    # Hyperparameters
    model = CNN(**config, output_size=n_tasks, tasktype="binary_classification").to(device)
    
    optimizer = Adam(model.parameters(), lr=1e-3)
    trainer = Trainer(model, criterion, optimizer, device)
    
    # Train for 10 EPOCH
    for epoch in range(10):
        # train epoch
        trainer.train_per_epoch(train_loader, epoch)
        
        # validation metrics
        _, val_loss, val_pred_prob, val_target_prob = trainer.predict(validate_loader)
        val_metric = trainer.evaluate(val_pred_prob, val_target_prob)
                
        # Send the current training result back to Tune
        tune.report(acc=val_metric, loss=val_loss)


def hyperparameter_tune(search_space=None, num_samples=10, max_num_epochs=10, gpus_per_trial=1, **kwargs):
    """Hyper Parameter Tune in NvTK.

    Currently, it only support NvTK.CNN architectures.
    The search_space define the Ray Tune’s search space, 
    where Tune will now randomly sample a combination of parameters.
    It will then train a number of models in parallel 
    and find the best performing one among these.
    We also use the ASHAScheduler which will terminate bad performing trials early.

    Parameters
    ----------
    search_space : dict, optional
        Ray Tune’s search space, Default is None.
        Here is an example:
        `search_space={
            "out_planes": tune.grid_search([32, 128, 512]),
            "kernel_size": tune.grid_search([5, 15, 25]),
            "bn": tune.choice([True, False])
        }`
    num_samples : int
        Number of sampled searching trials, Default is 10.
    max_num_epochs : int
        max number of epochs in ASHAScheduler, Default is 10.
    gpus_per_trial : int, floor
        specify the number of GPUs, Default is 1.

    Returns
    ----------
    best_trial
    """
    if search_space is None:
        search_space={
            "out_planes": tune.grid_search([32, 128, 512]),
            "kernel_size": tune.grid_search([5, 15, 25])
        }
    scheduler = ASHAScheduler(metric="loss", mode="min", max_t=max_num_epochs,
            grace_period=1, reduction_factor=2)
    result = tune.run(
        tune.with_parameters(train_model, 
                            train_loader=train_loader, 
                            validate_loader=validate_loader),
        config=search_space,
        metric="loss",
        mode="min",
        num_samples=num_samples,
        scheduler=scheduler
    )

    best_trial = result.get_best_trial("loss", "min", "last")
    print("Best trial config: {}".format(best_trial.config))
    print("Best trial final validation loss: {}".format(
        best_trial.last_result["loss"]))
    print("Best trial final validation accuracy: {}".format(
        best_trial.last_result["accuracy"]))
    print(best_trial)


if __name__ == "main":
    hyperparameter_tune(num_samples=2, max_num_epochs=2, gpus_per_trial=0)

import os, time, shutil, logging, pickle, random
import numpy as np
import torch

from .Evaluator import calculate_correlation, calculate_roc

class Trainer(object):
    """Model Trainer in NvTK.
    
    Trainer class could train and validate a NvTK or pytorch-based model.
    Trainer saved the model parameters to `Log` dirs, overwriting it 
    if the latest validation performance is better than 
    the previous best-performing model.

    Trainer maintain a Log dict to record the loss and metric 
    during training, Trainer will save the Log to `Figures`.

    Parameters
    ----------
    model : NvTK.model, torch.nn.module
        The NvTK model architecture, 
        or torch module container of parameters and methods.
    criterion : torch.nn.module
        The criterion of loss function,
        criterion was optimized during training.
    optimizer : torch.optim
        The optimizer to update the weights of Parameter,
        for minimizing the criterion during training.
    device : torch.device
        Context-manager that changes the selected device,
        representing the device on which the model is or will be allocated,
        ('cuda' or 'cpu').
    tasktype : str, optional
        Specify the task type, Default is "regression".
        (e.g. `tasktype="binary_classification"`)
    metric_sample : int, optional
        The number of sampled tasks to metric. Default is 100.
        For multi-task learning with more than `metric_sample` task, 
        we will sample tasks when calculating the metric values.
    item_sample : int, optional
        The number of sampled items to metric. Default is 50000.
        For dataset with more than `item_sample` items, 
        we will sample items when calculating the metric values.
    patience : int, optional
        The number of patience in early stopping methods, Default is 10.
        Early stopping could not work with a very large patience,
        set `patience=np.inf` to omit early stopping.
    resume : bool, optional
        Whether to resume model training, Default is False.
    resume_dirs : str, optional
        The directory in which to resume model training, Default is False.
    use_tensorboard : bool, optional
        Whether to use tensorboard. Default is False.
    tensorboard_args : dict, optional
        tensorboard arguments.


    Attributes
    ----------
    model : 
        The container of parameters and methods.
    device : str
        An object representing the device on which the model is or will be allocated.
        ('cuda' or 'cpu').
    criterion : 
        The criterion the model aims to minimize.
    optimizer : 
        The algorithm that updates the weights of Parameter during the backward step.
    patience : int
        The number of epochs to be trained after activating early stopping. 
    pred_prob : bool
        To judge the task type.
    metric_sample : int
        The number of metric sample.
    item_sample : int
        The number of item sample.
    logs: dict
        logs maintains loss and metric information during training.
    tensorboard : bool
        Whether to use tensorboard.
    writer : 
        Write tensorboard arguments to the summary.

    """
    def __init__(self, model, criterion, optimizer, device, tasktype="regression",
                    clip_grad_norm=True, max_norm=0.1, lr_scheduler=None,
                    evaluate_training=True, metric_sample=5000, item_sample=100000, 
                    patience=10,
                    resume=False, resume_dirs=None,
                    use_tensorbord=False, tensorbord_args={}):
        super().__init__()
        self.model = model.to(device)
        self.device = device
        self.criterion = criterion
        self.optimizer = optimizer
        self.patience = patience
        self.clip_grad_norm = clip_grad_norm
        self.max_norm = max_norm
        self.lr_scheduler = lr_scheduler
        self.selent = True
        # task related
        self.pred_prob = tasktype != "regression"
        self.evaluate_training = evaluate_training
        if not evaluate_training:
            logging.warning("setting trainer not to evaluate during training phase. \n \
                Note, the train/test/val accuracy will be set to 0.5, instead of real accuracy.")
        # sample tasks for calculate metric
        self.metric_sample = metric_sample
        self.item_sample = item_sample
        # Trainer will maintain logs information
        self.logs = {
            "train_batch_loss_list":[], 
            "val_batch_loss_list":[], 
            "test_batch_loss_list":[],

            "train_loss_list":[], 
            "val_loss_list":[], 
            "test_loss_list":[],

            "val_metric_list":[], 
            "test_metric_list":[], 
            "train_metric_list":[],
            "lrs":[],

            "best_val_loss":np.inf,
            "best_val_r":-np.inf,
            "best_val_epoch":0,
            }

        # TODO resume model training
        if resume and resume_dirs:
            pass
        
        # TODO capture error when loading tensorbord in lower version
        if use_tensorbord:
            from torch.utils.tensorboard import SummaryWriter
            self.tensorbord = True
            self.writer = SummaryWriter(**tensorbord_args)
        else:
            self.tensorbord = False

    def train_until_converge(self, train_loader, test_loader, validate_loader, EPOCH=100, resume=False, verbose_step=5):
        """
        Train until converge.
        
        Parameters
        ----------
        train_loader : 
            The data loader defined by using training dataset.
        test_loader : 
            The data loader defined by using testing dataset.
        validate_loader : 
            The data loader defined by using validation dataset. 
            validate_loader could be `None`, and will skip evaluate on val-set and early-stop.
        EPOCH : int
            An adjustable hyperparameter, The number of times to train the model until converge.
        resume : bool
            Whether to resume the model training. Default is False.
        verbose_step : int
            The number of steps to print the loss value. Default is 5.

        Attributes
        ----------
        tensorboard : bool
            Whether to use tensorboard.
        add_graph : 
            To iterate and visualize some samples in the training dataset.
        train_per_epoch : 
            Train each epoch.
        predict : 
            To predict based on the input data
        evaluate : 
            To evaluate the performance of the model.
        save_checkpoint : 
            save checkpoint.
        writer : 
            To add scalar or graph to the summary.
        load_best_model : 
            load the best model.
        
        """
        # graph
        if self.tensorbord:
            try:
                self.add_graph(next(iter(train_loader))[0][:2].shape)
            except BaseException as e:
                logging.warning("tensorbord cannot added graph")
                logging.warning(e)

        for epoch in range(EPOCH):
            # train
            train_batch_loss, train_loss = self.train_per_epoch(train_loader, epoch, verbose_step=20)

            # train metrics
            if self.evaluate_training:
                train_batch_loss, train_loss, train_pred_prob, train_target_prob = self.predict(train_loader)
                train_metric = self.evaluate(train_pred_prob, train_target_prob)
            else:
                train_metric = 0.5

            # test metrics
            if test_loader:
                test_batch_loss, test_loss, test_pred_prob, test_target_prob = self.predict(test_loader)
                test_metric = self.evaluate(test_pred_prob, test_target_prob) if self.evaluate_training else 0.5
            else:
                logging.info("The test_loader is None, \
                    and the test metrics will be set as test metrics.")
                test_batch_loss, test_loss, test_pred_prob, test_target_prob = val_batch_loss, val_loss, val_pred_prob, val_target_prob
                test_metric = val_metric

            # validation metrics
            if validate_loader:
                val_batch_loss, val_loss, val_pred_prob, val_target_prob = self.predict(validate_loader)
                val_metric = self.evaluate(val_pred_prob, val_target_prob) if self.evaluate_training else 0.5
            else:
                logging.info("The validate_loader is None, \
                    and the validation metrics will be set as test metrics.")
                val_batch_loss, val_loss, val_pred_prob, val_target_prob = test_batch_loss, test_loss, test_pred_prob, test_target_prob
                val_metric = test_metric

            # lr_scheduler
            _lr = self.optimizer.param_groups[0]['lr']

            # logs
            logging.info("Train\t Accuracy: %.4f\t Loss: %.4f\t\n" % (train_metric, train_loss))
            logging.info("Eval\t Accuracy: %.4f\t Loss: %.4f\t\n" % (val_metric, val_loss))
            
            self.logs['train_batch_loss_list'].extend(train_batch_loss)
            self.logs['val_batch_loss_list'].extend(val_batch_loss)
            self.logs['test_batch_loss_list'].extend(test_batch_loss)
            self.logs['train_loss_list'].append(train_loss)
            self.logs['val_loss_list'].append(val_loss)
            self.logs['test_loss_list'].append(test_loss)
            self.logs['train_metric_list'].append(train_metric)
            self.logs['val_metric_list'].append(val_metric)
            self.logs['test_metric_list'].append(test_metric)
            self.logs['lrs'].append(_lr)
            self.show_trainer_log()

            if self.tensorbord:
                self.writer.add_scalar('Loss/train', train_loss, epoch)
                self.writer.add_scalar('Loss/test', test_loss, epoch)
                self.writer.add_scalar('Accuracy/train', train_metric, epoch)
                self.writer.add_scalar('Accuracy/test', test_metric, epoch)

                try:
                    self.writer.add_histogram('Embedding.conv.bias', self.model.Embedding.conv.bias, epoch)
                    self.writer.add_histogram('Embedding.conv.weight', self.model.Embedding.conv.weight, epoch)
                except:
                    pass

            # update best
            if val_loss < self.logs["best_val_loss"]:
                self.logs["best_val_loss"] = val_loss
                self.logs["best_val_r"] = val_metric
                self.logs["best_val_epoch"] = epoch

                logging.info("Eval\t Best Eval Accuracy: %.4f\t Loss: %.4f\t at Epoch: %d\t lr: %.8f\n" % (val_metric, val_loss, epoch, _lr))
                logging.info("Eval\t Test Accuracy: %.4f\t Loss: %.4f\n" % (test_metric, test_loss))
                
                self.save_checkpoint(best=True)

            # or checkpoint
            elif epoch % 20 == 1:
                self.save_checkpoint()

            # early stop
            if epoch >= self.logs["best_val_epoch"] + self.patience: #patience_epoch:
                break
        
        # final
        fname = time.strftime("./Log/best_model@" + '%m%d_%H:%M:%S')
        shutil.copyfile("./Log/best_model.pth", fname+".params.pth")
        pickle.dump(self.logs, open(fname+".chekc_train_log.p", 'wb'))

        # final logs
        self.show_trainer_log(final=True)
        logging.info("\n"+self.model.__str__()+'\n')
        logging.info("Best Val Loss\t%.8f\t@epoch%d" % (self.logs['best_val_loss'], self.logs["best_val_epoch"]))
        logging.info("Best Val Metric\t%.8f\t@epoch%d" % (self.logs['best_val_r'], self.logs["best_val_epoch"]))
        logging.info("Best Test Metric\t%.8f\t@epoch%d" % (np.max(self.logs['test_metric_list']), np.argmax(self.logs['test_metric_list'])))
        logging.info("Best Test Loss\t%.8f\t@epoch%d" % (np.min(self.logs['test_loss_list']), np.argmin(self.logs['test_loss_list'])))

        self.load_best_model()

    def train_per_epoch(self, train_loader, epoch, verbose_step=5):
        """
        Train each epoch.

        Parameters
        ----------
        train_loader : 
            The data loader defined by using training dataset.
        epoch : int
            An adjustable hyperparameter, The number of times to train the model until converge.
        verbose_step : int
            The number of steps to print the loss value. Default is 5.
        
        Attributes
        ----------
        train_batch : 
            To train the batch fetched from the train_loader.
        
        Returns
        -------
        batch_losses : list
            The total loss values of all batches.
        average_loss : list
            The average of the losses of batches.

        """
        # train one epoch
        batch_losses = []
        self.model.train()
        for batch_idx, (data, target) in enumerate(train_loader):
            loss = self.train_batch(data, target)

            if batch_idx % verbose_step == 0:
                logging.info('Train Epoch: {} [{}/{} ({:.0f}%)]\tLoss: {:.6f}'.format(
                    epoch, batch_idx * len(data), len(train_loader.dataset),
                    100. * batch_idx / len(train_loader), loss))
            
            batch_losses.append(loss)
            average_loss = np.average(batch_losses)
        return batch_losses, average_loss
    
    def train_batch(self, data, target):
        """
        To meassure the loss of the batch.
        
        Parameters
        ----------
        data : numpy.ndarray
            The input data.
        target : numpy.ndarray
            True value that the user model was trying to predict.


        Attributes
        ----------
        model : 
            The container of parameters and methods.
        device : str
            An object representing the device on which the model is or will be allocated.
            ('cuda' or 'cpu').
        criterion : 
            The criterion the model aims to minimize.
        optimizer : 
            The algorithm that updates the weights of Parameter during the backward step.
        
        Returns
        -------
        loss : float
            The loss value.
        """
        self.model.train()
        if isinstance(data, list):
            data = [d.to(self.device) if isinstance(d, torch.Tensor) else d for d in data]
        else:
            data = data.to(self.device)
        if isinstance(target, list):
            target = [t.to(self.device) if isinstance(t, torch.Tensor) else t for t in target]
        else:
            target = target.to(self.device)
        output = self.model(data)
        loss = self.criterion(output, target)
        self.optimizer.zero_grad()
        loss.backward()
        if self.clip_grad_norm:
            torch.nn.utils.clip_grad_norm_(parameters=self.model.parameters(), 
                max_norm=self.max_norm, norm_type=2)
        self.optimizer.step()
        if self.lr_scheduler is not None:
            self.lr_scheduler.step()
        return loss.cpu().item()

    def predict(self, data_loader):
        """
        To predict based on the input data.

        parameters
        ----------
        data_loader : 
            The object to prepare the data for training, which retrieves batches of features and labels from the dataset iteratively.
        
        Attributes
        ----------
        model : 
            The container of parameters and methods.
        device : str
            An object representing the device on which the model is or will be allocated.
            ('cuda' or 'cpu').
        criterion : 
            The criterion the model aims to minimize.
        
        Returns
        -------
        batch_losses : list
            The total loss values of all batches.
        average_loss : list
            The average of the losses of batches.
        all_predictions : array
            All the predictions of the input data organized in arrays.
        all_targets : array
            All the targets organized in arrays.

        """
        batch_losses, all_predictions, all_targets = [], [], []
        self.model.eval()
        with torch.no_grad():
            for inputs, targets in data_loader:
                if isinstance(inputs, list):
                    inputs = [d.to(self.device) if isinstance(d, torch.Tensor) else d for d in inputs]
                else:
                    inputs = inputs.to(self.device)
                if isinstance(targets, list):
                    targets = [t.to(self.device) if isinstance(t, torch.Tensor) else t for t in targets]
                else:
                    targets = targets.to(self.device)
                # inputs, targets = inputs.to(self.device), targets.to(self.device)
                output = self.model(inputs)
                test_loss = self.criterion(output, targets)
                batch_losses.append(test_loss.cpu().item())
                
                if isinstance(targets, list):
                    all_targets.append([t.cpu().data.numpy() if isinstance(t, torch.Tensor) else t for t in targets])
                else:
                    all_targets.append(targets.cpu().data.numpy())
                if isinstance(output, (list, tuple)):
                    all_predictions.append([o.cpu().data.numpy() if isinstance(o, torch.Tensor) else o for o in output])
                else:
                    all_predictions.append(output.cpu().data.numpy())

        average_loss = np.average(batch_losses)
        if isinstance(targets, list):
            all_predictions = [np.vstack([batch[i] for batch in all_predictions]) 
                                for i in range(len(targets))]
            all_targets = [np.vstack([batch[i] for batch in all_targets]) 
                                for i in range(len(targets))]
        else:
            all_predictions = np.vstack(all_predictions)
            all_targets = np.vstack(all_targets)

        return batch_losses, average_loss, all_predictions, all_targets

    def evaluate(self, predict_prob, target_prob, sample_tasks=True, sample_items=True):
        """
        To evaluate the performance of the model.

        Parameters
        ----------
        predict_prob : list
            prediction probability
        target_prob : list
            target probability
        sample_tasks : bool
            Whether to sample tasks. Default is TRUE.
        sample_items : bool
            Whether to sample items. Default is TRUE.
        
        Attributes
        ----------
        metric_sample : int
            The number of metric sample.
        item_sample : int
            The number of item sample.
        pred_prob : bool
            Whether to calculate AUC-ROC curve.
        
        Returns
        -------
        metric : float
            The arithmetic mean of the auc/correlation, which is calculated by target probability and prediction probability.
        """
        item_size, output_size = target_prob.shape

        if self.metric_sample < output_size and sample_tasks and self.selent:
            metric_sample_idx = random.sample(range(output_size), self.metric_sample)
            target_prob = target_prob[:, metric_sample_idx]
            predict_prob = predict_prob[:, metric_sample_idx]
        logging.info("sampled %d tasks for metric" % self.metric_sample)

        if self.item_sample < item_size and sample_items and self.selent:
            item_sample_idx = random.sample(range(item_size), self.item_sample)
            target_prob = target_prob[item_sample_idx, :]
            predict_prob = predict_prob[item_sample_idx, :]
        logging.info("sampled %d items for metric" % (self.item_sample//10))

        if self.pred_prob:
            fpr, tpr, roc_auc = calculate_roc(target_prob, predict_prob)
            roc_l = [roc_auc[k] for k in roc_auc.keys() if roc_auc[k] >=0 and k not in ["macro", "micro"]]
            metric = np.mean(roc_l)
        
        else:
            correlation, pvalue = calculate_correlation(target_prob, predict_prob)
            correlation_l = [correlation[k] for k in correlation.keys() if k not in ["macro", "micro"]] # dict keys ordered by default in py3.7+
            metric = np.mean(correlation_l)

        return metric

    def add_graph(self, input_shape=(2, 4, 1000)):
        """
        add graph.

        Parameters
        ----------
        input_shape : tuple

        """
        self.writer.add_graph(self.model, torch.rand(input_shape).to(self.device))

    def load_best_model(self, fpath="./Log/best_model.pth"):
        """
        load the best model.

        Parameters
        ----------
        fpath : str
            The file path of the best model.

        """
        self.model.load(fpath)
        return self.model

    def get_current_model(self):
        """
        get the current model.
        """
        return self.model

    def save_checkpoint(self, best=False):
        """
        save checkpoint.

        Parameters
        ----------
        best : bool
            Whether the current model is the best model. Default is False.

        """
        os.makedirs("./Log", exist_ok=True)
        self.model.save("./Log/chekc_model.pth")
        torch.save(self.model, "./Log/chekc_model.p")
        pickle.dump(self.logs, open("./Log/chekc_train_log.p", 'wb'))

        if best:
            self.model.save("./Log/best_model.pth")
            torch.save(self.model, "./Log/best_model.p")
    
    def show_trainer_log(self, final=False):
        show_train_log(self.logs['train_loss_list'], 
                    self.logs['train_metric_list'], 
                    self.logs['val_loss_list'], 
                    self.logs['val_metric_list'],
                    self.logs['test_loss_list'], 
                    self.logs['test_metric_list'], 
                    self.logs['lrs'], output_fname='current_logs.pdf')
        show_train_log(loss_train=self.logs['train_loss_list'], 
                    loss_val=self.logs['val_loss_list'], 
                    loss_test=self.logs['test_loss_list'], 
                    output_fname='current_logs_loss.pdf')
        show_train_log(loss_val=self.logs['val_batch_loss_list'], output_fname='current_batch_loss_val.pdf')
        show_train_log(loss_test=self.logs['test_batch_loss_list'], output_fname='current_batch_loss_test.pdf')
        show_train_log(loss_train=self.logs['train_batch_loss_list'], output_fname='current_batch_loss_train.pdf')
        show_train_log(loss_val=self.logs['val_loss_list'], output_fname='current_loss_val.pdf')
        show_train_log(loss_test=self.logs['test_loss_list'], output_fname='current_loss_test.pdf')
        show_train_log(loss_train=self.logs['train_loss_list'], output_fname='current_loss_train.pdf')
        show_train_log(acc_val=self.logs['train_metric_list'], output_fname='current_acc_train.pdf')
        show_train_log(acc_val=self.logs['val_metric_list'], output_fname='current_acc_val.pdf')
        show_train_log(acc_test=self.logs['test_metric_list'], output_fname='current_acc_test.pdf')
        show_train_log(lrs=self.logs['lrs'], output_fname='current_lrs.pdf')

        if final:
            show_train_log(loss_val=self.logs['val_loss_list'], output_fname='final_loss_val.pdf')
            show_train_log(loss_train=self.logs['train_loss_list'], output_fname='final_loss_train.pdf')
            show_train_log(acc_val=self.logs['train_metric_list'], output_fname='final_acc_train.pdf')
            show_train_log(acc_val=self.logs['val_metric_list'], output_fname='final_acc_val.pdf')
            show_train_log(acc_test=self.logs['test_metric_list'], output_fname='final_acc_test.pdf')
            show_train_log(lrs=self.logs['lrs'], output_fname='final_lrs.pdf')


def show_train_log(loss_train=None, loss_val=None, acc_val=None, 
                   loss_test=None, acc_test=None,
                   lrs=None,
                    fig_size=(12,8),
                    save=True,
                    output_dir='Figures',
                    output_fname="Training_loss_log.pdf",
                    style="seaborn-colorblind",
                    fig_title="Training Log",
                    dpi=500):
    """function show train log in NvTK.

    Parameters
    ----------
    loss_train : list
        traing loss
    loss_val : list
        validation loss
    kernel_size : int, optional
        Size of the convolving kernel

    """
    
    os.makedirs(output_dir, exist_ok=True)
    
    import matplotlib
    backend = matplotlib.get_backend()
    if "inline" not in backend:
        matplotlib.use("PDF")
    import matplotlib.pyplot as plt
    plt.style.use(style)
    plt.figure()

    if loss_train:
        plt.plot(range(1, len(loss_train)+1), loss_train, 'b', label='Training Loss')
    if loss_val:
        plt.plot(range(1, len(loss_val)+1), loss_val, 'r', label='Validation Loss')
    if loss_test:
        plt.plot(range(1, len(loss_test)+1), loss_test, 'black', label='Test Loss')
    rate = 1
    if acc_val:
        plt.plot(range(1, len(acc_val)+1), list(map(lambda x: x*rate, acc_val)), 'g', label=str(rate)+'X Validation Accuracy')
    if acc_test:
        plt.plot(range(1, len(acc_test)+1), list(map(lambda x: x*rate, acc_test)), 'purple', label=str(rate)+'X Test Accuracy')
    if lrs:
        rate = int(1/(np.median(lrs)+1e-6))
        plt.plot(range(1, len(lrs)+1), list(map(lambda x: x*rate, lrs)), 'y--', alpha=0.2, label=str(rate)+'X Learning Rates')
    plt.title(fig_title)
    plt.legend()
    if save:
        plt.savefig(os.path.join(output_dir, output_fname),
                    format="pdf", dpi=dpi)
    else:
        plt.show()
    plt.close()


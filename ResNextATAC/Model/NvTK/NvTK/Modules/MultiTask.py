import logging
from copy import copy

import torch
from torch import nn

import networkx as nx

LAYER_KEY = 'layers'
NAME_KEY = 'name'
ANCHOR_KEY = 'anchor_layer'
LOSS_KEY = 'loss'
LOSS_REG_KEY = 'loss_weight'
AUTO_WEIGHT_KEY = 'auto'
WEIGHT_INIT_KEY = 'loss_init_val'
MISSING_WEIGHT_MSG = "Expect {0} for task {1} but none provided."

class MTLModel(nn.Module):
    """
    A torch.nn.Module built from a set of shared and task specific layers
    Attributes
    ----------
    g : networkx.Graph
        The meta-computation graph
    task_layers : list
        A list which holds the layers for which to build the computation graph
    output_tasks : list
        A list which holds the tasks for which the output should be returned
    layer_names : list
        A list of the names of each layer
    losses : dict
        A dictionary which maps the name of a layer to its loss function
    loss_weights : dict
        A dictionary which maps the name of a layer to the weight of its loss
        function
    """

    def __init__(self, task_layers, output_tasks):
        super(MTLModel, self).__init__()
        self.task_layers = task_layers
        self.output_tasks = output_tasks
        self.layer_names = [t[NAME_KEY] for t in task_layers]

        self._initialize_graph()

        self._initialize_losses()
        self._initialize_loss_weights()

    def _initialize_losses(self):
        self.losses = {task[NAME_KEY]: task[LOSS_KEY] for task in self.task_layers if LOSS_KEY in task.keys()}

    def _initialize_loss_weights(self):
        self.loss_weights = {}
        for task in self.task_layers:
            self._set_loss_weight(task)

    def _set_loss_weight(self, task):
        task_name = task[NAME_KEY]
        if LOSS_REG_KEY in task.keys():
            if task[LOSS_REG_KEY] == AUTO_WEIGHT_KEY:
                assert WEIGHT_INIT_KEY in task.keys(),\
                        MISSING_WEIGHT_MSG.format(WEIGHT_INIT_KEY, task_name)
                loss_weight = task[WEIGHT_INIT_KEY]
                loss_name = f'{task_name}_loss'
                loss_weight = torch.nn.Parameter(torch.full((1,), loss_weight))
                setattr(self, loss_name, loss_weight)
                self.loss_weights[task_name] = getattr(self, loss_name)
            else:
                self.loss_weights[task_name] = task[LOSS_REG_KEY]
    
    def _initialize_graph(self):
        self.g = nx.DiGraph()
        self.g.add_node('root')
        self._build_graph()

    def _bfs_forward(self, start_node):
        ''' Here we iteratore through the graph in a BFS-fashion starting from
        `start_node`, typically this is the `root` node. This node is skipped
        and we pass the input data and resulting outputs from all layers foward.
        '''
        visited = {node: False for node in self.layer_names}

        # First node is visited
        queue = [start_node]
        visited[start_node] = True

        while queue:
            node = queue.pop(0)
            if node != start_node:
                input_nodes = self.g.predecessors(node)
                if logging.getLogger().level == logging.DEBUG:
                    l = copy(input_nodes)
                    print(f"Feeding output from {list(l)} into {node}")
                cur_layer = getattr(self, node)

                # Get the output from the layers that serve as input
                output_pre_layers = []
                output_complete = True
                for n in input_nodes:
                    # If an output is not ready yet, because that node has not
                    # been computed, we put the current node back into the queue
                    if n not in self.outputs.keys():
                        if logging.getLogger().level == logging.DEBUG:
                            print(f"No output for layer {n} yet")
                        output_complete = False
                        break
                    else:
                        output_pre_layers.append(self.outputs[n])

                if not output_complete:
                    if logging.getLogger().level == logging.DEBUG:
                        print(f"Putting {node} back into the queue.")
                    queue.append(node)
                else:
                    cur_output = cur_layer(*output_pre_layers)
                    self.outputs[node] = cur_output

            for i in self.g.successors(node):
                if visited[i] == False:
                    queue.append(i)
                    visited[i] = True

        losses, loss_weights = self._get_losses()
        return [self.outputs[t] for t in self.output_tasks], losses, loss_weights

    def _get_losses(self):
        losses = []
        loss_weights = []
        for t in self.output_tasks:
            losses.append(self.losses.get(t))
            loss_weights.append(self.loss_weights.get(t))
        return losses, loss_weights

    def _build_graph(self):
        for layer in self.task_layers:
            self._add_layer(layer)
            self._add_to_graph(layer)

    def _add_to_graph(self, layer):
        layer_name = layer[NAME_KEY]
        self._add_node(layer_name)

        if 'anchor_layer' not in layer.keys():
            # If there is no anchor layer, we expect it to be a layer which
            # receives data inputs and is hence connected to the root node
            self.g.add_edge('root', layer_name)
        else:
            anchor_layer = layer[ANCHOR_KEY]
            if isinstance(anchor_layer, list):
                for a_l_name in anchor_layer:
                    self._add_node(a_l_name)
                    self.g.add_edge(a_l_name, layer_name)
            else:
                self._add_node(anchor_layer)
                self.g.add_edge(anchor_layer, layer_name)

    def _add_node(self, layer):
        if isinstance(layer, str):
            layer_name = layer
            self.g.add_node(layer_name)
        else:
            layer_name = layer[NAME_KEY]
            self.g.add_node(layer_name)
            if 'anchor_layer' not in layer.keys():
                self.g.add_edge('root', layer_name)
    
    def _add_layer(self, layer):
        layer_modules = layer[LAYER_KEY]
        layer_name_main = layer[NAME_KEY]
        setattr(self, layer_name_main, layer_modules)

    def forward(self, input):
        self.outputs = {'root': input}
        return self._bfs_forward('root')


class ConcatenateTask(nn.Module):
    def forward(self, x):
        return torch.cat(x, dim=1)


class MultiTaskWrapper(nn.Module):
    def __init__(self, sub_models, concat_dim=1):
        """
        The multi-model wrapper class can be used to concatenate the
        outputs of multiple models along a pre-specified axis. The wrapper
        can be used to load and run multiple trained models during prediction
        functions. This class should not be used for training. 

        Parameters
        ----------
        sub_models : list(torch.nn.Module)
            The 'sub-models' that are used in this multi-model wrapper class.
        concat_dim : int, optional
            Default is 1. The dimension along which to concatenate the models'
            predictions.
        """
        super().__init__()
        self.sub_models = sub_models
        self._concat_dim = concat_dim

    def cuda(self):
        for sm in self.sub_models:
            sm.cuda()

    def eval(self):
        for sm in self.sub_models:
            sm.eval()

    def forward(self, x):
        return torch.cat([sm(x) for sm in self.sub_models], self._concat_dim)


class MTLTrainer(object):
    def __init__(self, model, criterion, optimizer, device):
        super().__init__()
        # train one epoch
        self.model = model.to(device)
        self.device = device
        self.criterion = criterion
        self.optimizer = optimizer

    def train_per_epoch(self, train_loader, epoch, verbose_step=5):
        batch_losses = []
        self.model.train()
        for batch_idx, (data, target) in enumerate(train_loader):
            data, target = data.to(self.device), target.to(self.device)
            self.optimizer.zero_grad()

            # Our model will return a list of predictions (from the layers specified in `output_tasks`),
            # loss functions, and regularization parameters (as defined in the tasks variable)
            y_hat, l_funcs, l_weights = self.model(data)

            loss = 0
            idx = 0
            # We can now iterate over the tasks and accumulate the losses
            for i, cnt in enumerate(anno_cnt.values):
                loss += l_weights[i] * l_funcs[i](y_hat[i], target[:,idx:idx+cnt])
                idx += cnt

            # self.criterion(output, target)
            loss.backward()
            self.optimizer.step()
            
            if verbose_step:
                if batch_idx % verbose_step == 0:
                    logging.info('Train Epoch: {} [{}/{} ({:.0f}%)]\tLoss: {:.6f}'.format(
                        epoch, batch_idx * len(data), len(train_loader.dataset),
                        100. * batch_idx / len(train_loader), loss.data))
            batch_losses.append(loss.cpu().item())
            average_loss = np.average(batch_losses)

        return average_loss

    def evaluate(self, data):
        eval_loader = data
        batch_losses, all_predictions, all_targets = [], [], []
        self.model.eval()
        with torch.no_grad():
            for inputs, targets in eval_loader:
                inputs, targets = inputs.to(self.device), targets.to(self.device)

                y_hat, l_funcs, l_weights = self.model(inputs)

                loss = 0
                idx = 0
                # We can now iterate over the tasks and accumulate the losses
                for i, cnt in enumerate(anno_cnt.values):
                    loss += l_weights[i] * l_funcs[i](y_hat[i], targets[:,idx:idx+cnt])
                    idx += cnt

                test_loss = loss
                batch_losses.append(test_loss.cpu().item())

                output = torch.cat(y_hat, dim=1)
                all_predictions.append(output.cpu().data.numpy())
                all_targets.append(targets.cpu().data.numpy())
                
        average_loss = np.average(batch_losses)
        all_predictions = np.vstack(all_predictions)
        all_targets = np.vstack(all_targets)
        return average_loss, all_predictions, all_targets

    def get_current_model(self):
        return self.model


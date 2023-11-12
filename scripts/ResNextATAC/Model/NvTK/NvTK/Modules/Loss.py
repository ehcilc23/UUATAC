import torch
from torch import nn

class BCEFocalLoss(nn.Module):
    def __init__(self, gamma=2, weight=None, reduction='mean'):
        super().__init__()
        self.gamma = gamma
        self.weight = weight
        self.reduction = reduction
 
    def forward(self, inputs, targets):
        BCE_loss = F.binary_cross_entropy(inputs, targets, weight=self.weight, reduction='none')
        probs = torch.exp(-BCE_loss)
        focal_weight = torch.pow(1-probs, self.gamma)
        
        loss = focal_weight * BCE_loss
            
        if self.reduction == 'mean':
            loss = torch.mean(loss)
        elif self.reduction == 'sum':
            loss = torch.sum(loss)
            
        return loss


class MTLoss(nn.Module):
    def __init__(self, lamda=1e-8):
        super().__init__()
        self.lamda = lamda
        if pred_prob:
            self.loss_fn = nn.BCELoss()
        else:
            self.loss_fn = nn.MSELoss()

        if use_focal_loss:
            self.loss_fn = BCEFocalLoss()

    def forward(self, pred, target):
        L1_loss = 0
        for param in model.parameters():
            L1_loss += torch.sum(torch.abs(param))
        logging.debug(L1_loss)

        loss = self.loss_fn(pred, target)
        logging.debug(loss)

        MTloss = loss + self.lamda * L1_loss
        return MTloss


class scMTLoss(nn.Module):
    def __init__(self, annotation):
        super().__init__()
        if pred_prob:
            self.loss_fn = nn.BCELoss(reduction="none")
        else:
            self.loss_fn = nn.MSELoss(reduction="none")
        self.masks = [torch.ByteTensor((annotation == factor).astype(int)) for factor in np.unique(annotation)]

    def forward(self, pred, target):
        loss = self.loss_fn(pred, target)
        logging.debug(loss)
        
        factor_loss = []
        for mask in self.masks:
            factor_loss.append(loss[:, mask].mean())

        MTloss = torch.stack(factor_loss).mean()
        return MTloss


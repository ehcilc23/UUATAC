import numpy as np
from NvTK.Evaluator import *

prediction = np.random.rand((10, 4))

target = np.random.choice(range(4), 10)

target = onehot_encode(target)
assert target.shape == (10, 4)

fpr, tpr, roc_auc = calculate_roc(target, prediction)
visualize_roc_curves(target, prediction)

p, r, ap = calculate_pr(target, prediction)
visualize_pr_curves(target, prediction)

target = np.random.rand((10, 4))
c, p = calculate_correlation(target, prediction, method="pearson")
c, p = calculate_correlation(target, prediction, method="spearman")
c, p = calculate_correlation(target, prediction, method="kendall")

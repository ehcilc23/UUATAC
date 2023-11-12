"""Evaluate Metrics provided in NvTK."""

import os, logging, itertools

import numpy as np
import pandas as pd

from sklearn.metrics import auc, roc_curve, precision_recall_curve, average_precision_score
from scipy.stats import pearsonr, spearmanr, kendalltau


# class Metric(object):
#     def __init__(self, task_name, tasktype=None, metrics={"auroc":calculate_roc}):
#         super().__init__()
#         self.task_name = task_name
#         self.results = {}
#         self.metrics = {}
#         if tasktype in ["regression"]:
#             self.metrics["pcc"] = calculate_correlation
#             self.metrics["scc"] = [calculate_correlation, {"method":"spearman"}]
#         elif tasktype in ["classification", "binary_classification"]:
#             self.metrics["auroc"] = calculate_roc
#             self.metrics["scc"] = calculate_pr
#         for k, v in metrics.items():
#             self.metrics[k] = v
            
#     def calculate(self, target, prediction):
#         results = {}
#         for k, v in self.metrics.items():
#             if isinstance(v, list):
#                 metric, args = v
#                 results[k] = metric(target, prediction, **args)
#             else:
#                 metric = v
#                 results[k] = metric(target, prediction)
        
#         self.results = results
#         return results

#     def write(self, fname):
#         pd.DataFrame(self.results, index=self.task_name).T.to_csv(fname)

#     def draw(self):
#         pass
        

def onehot_encode(label):
    from sklearn.preprocessing import label_binarize
    return label_binarize(label, classes=range(np.max(label)+1))

# TODO 
def map_prob2label(y_pred_prob, map_fn=np.argmax):
    assert isinstance(y_pred_prob, np.ndarray)
    return np.array(list(map(map_fn, y_pred_prob)))


def calculate_roc(target, prediction):
    # assert len(np.shape(prediction))>1, "Input should be y_prediction_Probability"
    if len(np.shape(target)) == 1:
        target = onehot_encode(target)
    fpr, tpr, roc_auc = {}, {}, {} # orderedDict after python3.8
    n_classes = target.shape[-1]
    for index in range(n_classes):
        feature_targets = target[:, index]
        feature_preds = prediction[:, index]
        if len(np.unique(feature_targets)) > 1:
            fpr[index], tpr[index], _ = roc_curve(feature_targets, feature_preds)
            roc_auc[index] = auc(fpr[index], tpr[index])
        else:
            logging.warning("roc value was underestimated!")
            roc_auc[index] = 0
    fpr['micro'], tpr['micro'], _ = roc_curve(target.ravel(), prediction.ravel())
    roc_auc['micro'] = auc(fpr['micro'], tpr['micro'])
    # First aggregate all false positive rates
    all_fpr = np.unique(np.concatenate([fpr.get(i, [0]) for i in range(n_classes)]))
    # Then interpolate all ROC curves at this points
    mean_tpr = np.zeros_like(all_fpr)
    for i in range(n_classes):
        mean_tpr += np.interp(all_fpr, fpr.get(i, [0]), tpr.get(i, [0]))
    # Finally average it and compute AUC
    mean_tpr /= n_classes
    fpr["macro"] = all_fpr
    tpr["macro"] = mean_tpr
    roc_auc["macro"] = auc(fpr["macro"], tpr["macro"])
    return fpr, tpr, roc_auc


def show_auc_curve(fpr, tpr, roc_auc,
                    fig_size=(10,8),
                    save=True,
                    output_dir='Figures',
                    output_fname='roc_curves.pdf',
                    style="seaborn-colorblind",
                    fig_title="Feature ROC curves",
                    dpi=500):
    import matplotlib
    backend = matplotlib.get_backend()
    if "inline" not in backend:
        matplotlib.use("PDF")
    import matplotlib.pyplot as plt
    plt.style.use(style)
    plt.figure()
    # n_classes
    n_classes = len(roc_auc) - 2
    # Plot all ROC curves
    plt.figure(figsize=fig_size)
    lw = 1
    plt.plot(fpr["micro"], tpr["micro"],
             label='micro-average ROC curve (area = {0:0.2f})'
                   ''.format(roc_auc["micro"]),
             color='deeppink', linestyle=':', linewidth=4)
    plt.plot(fpr["macro"], tpr["macro"],
             label='macro-average ROC curve (area = {0:0.2f})'
                   ''.format(roc_auc["macro"]),
             color='navy', linestyle=':', linewidth=4)
    # colors = ['aqua', 'darkorange', 'cornflowerblue', 'red', 'blue', 'green', 'grey', 'black', 'yellow', 'purple', 'brown', 'darkblue', 'darkred', 'gold', 'orange', 'pink', 'violet', 'turquoise', 'tomato']
    colors = ["grey"]
    for i, color in zip(range(n_classes), itertools.cycle(colors)):
        plt.plot(fpr.get(i, [0]), tpr.get(i, [0]), color=color, lw=lw,
                #  label='ROC curve of class {0} (area = {1:0.2f})'
                #  ''.format(i, roc_auc.get(i, 0))
                 )
    plt.plot([0, 1], [0, 1], 'k--', lw=lw)
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title(fig_title)
    plt.legend(loc="lower right")
    if save:
        plt.savefig(os.path.join(output_dir, output_fname),
                    format="pdf", dpi=dpi)
    else:
        plt.show()
    plt.close()


def visualize_roc_curves(target, prediction,
                         fig_size=(10,8),
                         save=True,
                         output_dir='Figures',
                         output_fname='roc_curves.pdf',
                         style="seaborn-colorblind",
                         fig_title="Feature ROC curves",
                         dpi=500):
    """
    Output the ROC curves for each feature predicted by a model
    as an SVG.

    Parameters
    ----------
    prediction : numpy.ndarray
        Value predicted by user model.
    target : numpy.ndarray
        True value that the user model was trying to predict.
    output_dir : str
        The path to the directory to output the figures. Directories that
        do not currently exist will be automatically created.
    style : str, optional
        Default is "seaborn-colorblind". Specify a style available in
        `matplotlib.pyplot.style.available` to use.
    fig_title : str, optional
        Default is "Feature ROC curves". Set the figure title.
    dpi : int, optional
        Default is 500. Specify dots per inch (resolution) of the figure.

    Returns
    -------
    None
        Outputs the figure in `output_dir`.

    """
#     assert len(np.shape(prediction))>1, "Input should be y_prediction_Probability"
    if len(np.shape(target)) == 1:
        target = onehot_encode(target)
    os.makedirs(output_dir, exist_ok=True)
    # calculate_roc(target, prediction)
    fpr, tpr, roc_auc = calculate_roc(target, prediction)
    show_auc_curve(fpr, tpr, roc_auc,
                    fig_size=fig_size,
                    save=save,
                    output_dir=output_dir,
                    output_fname=output_fname,
                    style=style,
                    fig_title=fig_title,
                    dpi=dpi)


def calculate_pr(target, prediction):
    assert len(np.shape(prediction))>1, "Input should be y_prediction_Probability"
    if len(np.shape(target)) == 1:
        target = onehot_encode(target)
    # For each class
    precision, recall, average_precision = {}, {}, {}
    n_classes = target.shape[-1]
    for index in range(n_classes):
        feature_targets = target[:, index]
        feature_preds = prediction[:, index]
        if len(np.unique(feature_targets)) > 1:
            precision[index], recall[index], _ = precision_recall_curve(feature_targets, feature_preds)
            average_precision[index] = average_precision_score(feature_targets, feature_preds)
        else:
            logging.warning("pr value was underestimated!")
            average_precision[index] = 0
    # A "micro-average": quantifying score on all classes jointly
    precision["micro"], recall["micro"], _ = precision_recall_curve(target.ravel(), prediction.ravel())
    average_precision["micro"] = average_precision_score(target.ravel(), prediction.ravel(), average="micro")
    return precision, recall, average_precision


def show_pr_curve(precision, recall, average_precision,
                    fig_size=(10,8),
                    save=True,
                    output_dir='Figures',
                    output_fname='pr_curves.pdf',
                    style="seaborn-colorblind",
                    fig_title="Feature PR curves",
                    dpi=500):
    import matplotlib
    backend = matplotlib.get_backend()
    if "inline" not in backend:
        matplotlib.use("PDF")
    import matplotlib.pyplot as plt
    plt.style.use(style)
    plt.figure()
    # n_classes
    n_classes = len(precision) - 1
    # Plot all ROC curves
    plt.figure(figsize=fig_size)
    lw = 2 

    # setup plot details
    colors = ['aqua', 'darkorange', 'cornflowerblue', 'red', 'blue', 'green', 'grey', 'black', 'yellow', 'purple']

    plt.figure(figsize=fig_size)
    f_scores = np.linspace(0.2, 0.8, num=4)
    lines = []
    labels = []
    for f_score in f_scores:
        x = np.linspace(0.01, 1)
        y = f_score * x / (2 * x - f_score)
        l, = plt.plot(x[y >= 0], y[y >= 0], color='gray', alpha=0.2)
        plt.annotate('f1={0:0.1f}'.format(f_score), xy=(0.9, y[45] + 0.02))

    lines.append(l)
    labels.append('iso-f1 curves')
    l, = plt.plot(recall["micro"], precision["micro"], color='gold', lw=lw)
    lines.append(l)
    labels.append('micro-average Precision-recall (area = {0:0.2f})'
                ''.format(average_precision["micro"]))

    for i, color in zip(range(n_classes), colors):
        l, = plt.plot(recall[i], precision[i], color=color, lw=lw)
        lines.append(l)
        labels.append('Precision-recall for class {0} (area = {1:0.2f})'
                    ''.format(i, average_precision[i]))

    fig = plt.gcf()
    fig.subplots_adjust(bottom=0.25)
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('Recall')
    plt.ylabel('Precision')
    plt.title(fig_title)
    plt.legend(lines, labels, loc=(0, -.38), prop=dict(size=14))
    if save:
        plt.savefig(os.path.join(output_dir, output_fname),
                    format="pdf", dpi=dpi)
    else:
        plt.show()
    plt.close()


def visualize_pr_curves(target, prediction,
                        fig_size=(10,8),
                        save=True,
                        output_dir='Figures',
                        output_fname='pr_curves.pdf',
                        style="seaborn-colorblind",
                        fig_title="Feature precision-recall curves",
                        dpi=500):
    """
    Output the precision-recall (PR) curves for each feature predicted by
    a model as an SVG.

    Parameters
    ----------
    prediction : numpy.ndarray
        Value predicted by user model.
    target : numpy.ndarray
        True value that the user model was trying to predict.
    output_dir : str
        The path to the directory to output the figures. Directories that
        do not currently exist will be automatically created.
    report_gt_feature_n_positives : int, optional
        Default is 50. Do not visualize an PR curve for a feature with
        less than 50 positive examples in `target`.
    style : str, optional
        Default is "seaborn-colorblind". Specify a style available in
        `matplotlib.pyplot.style.available` to use.
    fig_title : str, optional
        Default is "Feature precision-recall curves". Set the figure title.
    dpi : int, optional
        Default is 500. Specify dots per inch (resolution) of the figure.

    Returns
    -------
    None
        Outputs the figure in `output_dir`.

    """
    os.makedirs(output_dir, exist_ok=True)
    assert len(np.shape(prediction))>1, "Input should be y_prediction_Probability"
    if len(np.shape(target)) == 1:
        target = onehot_encode(target)
    os.makedirs(output_dir, exist_ok=True)
    # calculate_roc(target, prediction)
    precision, recall, average_precision = calculate_pr(target, prediction)
    show_pr_curve(precision, recall, average_precision,
                    fig_size=fig_size,
                    save=save,
                    output_dir=output_dir,
                    output_fname=output_fname,
                    style=style,
                    fig_title=fig_title,
                    dpi=dpi)


def calculate_correlation(target, prediction, method="pearson"):
    if method == "pearson":
        correlation_fn = pearsonr
    elif method == "spearman":
        correlation_fn = spearmanr
    elif method == "kendall":
        correlation_fn = kendalltau
    # assert len(np.shape(prediction))>1, "Input should be y_prediction_Probability"
    if len(np.shape(target)) == 1:
        target = onehot_encode(target)
    correlation, pvalue = {}, {} # orderedDict after python3.8
    n_classes = target.shape[-1]
    for index in range(n_classes):
        feature_targets = target[:, index]
        feature_preds = prediction[:, index]
        if len(np.unique(feature_targets)) > 1:
            correlation[index], pvalue[index] = correlation_fn(feature_targets, feature_preds)
    return correlation, pvalue


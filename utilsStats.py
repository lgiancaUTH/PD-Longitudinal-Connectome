from scipy.stats import mannwhitneyu
from sklearn.metrics import roc_curve, auc, accuracy_score, confusion_matrix, roc_auc_score, cohen_kappa_score

import numpy as np


def sigTestAUC(data1, data2, disp='long'):
    '''
    return a string with AUC and significance based on the Mann Whitney test
    disp= short|long|auc
    '''
    u, p_value = mannwhitneyu(data1, data2, alternative='two-sided')
    # p_value *= 2 # no longer required

    p_val_str = ''
    pValStars = ''
    if (p_value <= 0.001):
        p_val_str = '***p<0.001'
        pValStars = '***'
    elif (p_value <= 0.01):
        p_val_str = '**p<0.01'
        pValStars = '**'
    elif (p_value <= 0.05):
        p_val_str = '*p<0.05'
        pValStars = '*'
    else:
        p_val_str = 'not sig. p={:0.2f}'.format(p_value)
        pValStars = ''

    aucVal = 1 - u / (len(data1) * len(data2))

    if disp == 'short':
        strOut = '{:0.2f}{:}'.format(aucVal, pValStars)
    elif disp == 'long':
        strOut = '{:0.2f} ({:})'.format(aucVal, p_val_str)
    else:
        strOut = '{:0.2f}'.format(aucVal)

    return strOut


def findCutoffPnt3(dataPos, dataNeg):
    """
    Find cutoff point minimizing the distance to Sens 1, spec 1 and calculate statistics (with kappa, code based on findCutoffPnt2).
    format confMat:
     array([[TN, FP],
            [ FN, TP]]))
    :param dataPos:
    :param dataNeg:
    :param dataPosGtNeg:
    :return: acc,sens,spec,roc_auc, cutoffTh, confusionMat, kappa
    """

    dataAll = np.concatenate((dataPos, dataNeg))
    lblArr = np.zeros(len(dataAll), dtype=bool)
    lblArr[0:len(dataPos)] = True

    fpr, tpr, thresholds = roc_curve(lblArr, dataAll, pos_label=True)
    roc_auc = auc(fpr, tpr)

    # invert comparison if (ROC<0.5) required
    if roc_auc < 0.5:
        lblArr = ~lblArr
        fpr, tpr, thresholds = roc_curve(lblArr, dataAll, pos_label=True)
        roc_auc = auc(fpr, tpr)
        print 'inverting labels'

    # calculate best cut-off based on distance to top corner of ROC curve
    distArr = np.sqrt(np.power(fpr, 2) + np.power((1 - tpr), 2))
    cutoffIdx = np.argsort(distArr)[0]
    cutoffTh = thresholds[cutoffIdx]

    lblOut = dataAll >= cutoffTh

    acc = accuracy_score(lblArr, lblOut)
    sens = tpr[cutoffIdx]
    spec = 1 - fpr[cutoffIdx]
    cfMat = confusion_matrix(lblArr, lblOut)

    kappa = cohen_kappa_score(lblOut, lblArr)

    return (acc, sens, spec, roc_auc, cutoffTh, cfMat, kappa)
# base directory
import sys
sys.path.append("/data/gdrive/datasets/ppmi/src")

#imports
import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt
import seaborn

import utils
import utilsStats

from sklearn.model_selection import KFold, StratifiedKFold
from sklearn.covariance import GraphLassoCV
from scipy.interpolate import interp1d
import sklearn.linear_model as lm
import sklearn.metrics as met
import sklearn.svm as svm
import sklearn.ensemble as ens
import sklearn.preprocessing as pre
import copy

from sklearn.manifold import TSNE
from sklearn.decomposition import PCA


# Load configuration
CONF_FILE = 'data/configuration.json'
CONF = utils.readConf(CONF_FILE)


# %qtconsole
#Freesurfer base dir
#FS_BASE_DIR = '/home/lgianca/freesurfer/subjects/'
FS_BASE_DIR = '/data/freesurfer_subj/'


def load_stats(statsfile):
    """
    load a freesurfer stats file
    """

    try:
        f = open(statsfile)
    except:
        print 'problem opening', statsfile
        return False

    lines = [i.strip() for i in f.readlines()]
    f.close()

    stats = {}
    for l in lines:
        if l.find('# ColHeaders') > -1:
            l_s = l.split()
            # print l_s
            colhdrs = l_s[2:]
            for colhdr in colhdrs:
                stats[colhdr] = []
        elif l.find('#') == 0:
            continue
        else:
            l_s = l.split()
            for c in range(len(colhdrs)):
                stats[colhdrs[c]].append(l_s[c])

    statsFr = pd.DataFrame(stats)
    statsFr = statsFr.set_index('SegId')
    return statsFr


def loadSubj(imgInfoFr, numId, typeIn=0, colName='fsName'):
    """
    load subject statistics
    numId: 0-based id
    typeIn = 0: subcortical, 1 including WM, 2 only WM
    """
    fileInAseg = FS_BASE_DIR + imgInfoFr[colName].iloc[numId] + '/stats/aseg.stats'
    asegStats = load_stats(fileInAseg)
    fileInWmparc = FS_BASE_DIR + imgInfoFr[colName].iloc[numId] + '/stats/wmparc.stats'
    wmparcStats = load_stats(fileInWmparc)

    if typeIn == 0:
        return asegStats
    elif typeIn == 1:
        return asegStats.append(wmparcStats)
    elif typeIn == 2:
        return wmparcStats


def loadConnectome(imgInfoFr, numId,
                   connName='0_connectomeLocal.csv',
                   baseDir='/data/bigdata/PPMIcomp/', colName='fsName'):
    """
    load connectome as a np matrix
    numId: 0-based id
    """

    connArr = np.loadtxt(baseDir + imgInfoFr[colName].iloc[numId] + '/' + connName)
    if len(connArr[np.isnan(connArr)]) > 0:
        print 'Nan values in: ', imgInfoFr[colName].iloc[numId]

    if connArr.sum()==0:
        print 'All zeros in ', imgInfoFr[colName].iloc[numId]

    return connArr


def compConnVec(connArr):
    """
    Convert connectome Mat to feat vectorX0
    """
    lblNum = connArr.shape[0]
    # find upper indices
    indConn = np.triu_indices(lblNum, 0)
    indConnNum = len(indConn[0])
    # linearize
    X = connArr[indConn].tolist()

    return X


def compVec2Conn(connVec):
    """
    Convert  feat vector to connectome Mat
    """
    connNum = len(connVec)
    # use inverse of gauss method to compute the original size
    l = (np.sqrt(1. + 8 * connNum) - 1) / 2.
    l = int(l)

    # create connectome mat
    connArr = np.zeros((l, l))
    # find upper indices
    indConn = np.triu_indices(l, 0)

    connArr[indConn[0], indConn[1]] = connVec

    return connArr


def rocBootstrap(cntArr, pdArr, bootstrapsNumIn=1000):
    '''
    ROC bootstrap
    '''
    n_bootstraps = bootstrapsNumIn
    rng_seed = 42  # control reproducibility
    bootstrapped_scores = []
    bootRocLst = []
    y_true = np.append([0] * len(cntArr), [1] * len(pdArr))
    y_pred = np.append(cntArr, pdArr)
    rng = np.random.RandomState(rng_seed)
    # create fpr grid for interpolation
    fprGridVec = np.linspace(0, 1, 100, endpoint=True)
    # matrix containing all tpr corresponding to fprGridVec
    tprGridMat = np.zeros((len(fprGridVec), n_bootstraps))
    for i in range(n_bootstraps):
        # bootstrap by sampling with replacement on the prediction indices
        indices = rng.randint(0, len(y_pred) - 1, len(y_pred))
        if len(np.unique(y_true[indices])) < 2:
            # We need at least one positive and one negative sample for ROC AUC
            # to be defined: reject the sample
            continue

        score = met.roc_auc_score(y_true[indices], y_pred[indices])
        tmpFpr, tmpTpr, _ = met.roc_curve(y_true[indices], y_pred[indices])
        tmpFpr = np.concatenate(([0], tmpFpr, [1]))
        tmpTpr = np.concatenate(([0], tmpTpr, [1]))

        # interpolate for comparable ROCs
        fInter = interp1d(tmpFpr, tmpTpr, kind='nearest')
        tprGridMat[:, i] = fInter(fprGridVec)

        bootstrapped_scores.append(score)
        bootRocLst.append([tmpFpr, tmpTpr])
        # print("Bootstrap #{} ROC area: {:0.3f}".format(i + 1, score))

    # confidence interval for AUCs
    sorted_scores = np.array(bootstrapped_scores)
    sorted_scores.sort()
    confidence_lower = sorted_scores[int(0.05 * len(sorted_scores))]
    confidence_upper = sorted_scores[int(0.95 * len(sorted_scores))]

    return (confidence_lower, confidence_upper, fprGridVec, tprGridMat)


def plotRocAndConf(fprGridVec, tprGridMat, labelIn=''):
    n_bootstraps = tprGridMat.shape[1]

    # confidence interval for ROC
    tprGridMatS = np.sort(tprGridMat, axis=1)
    tprLow05 = tprGridMatS[:, int(0.05 * n_bootstraps)]
    tprTop95 = tprGridMatS[:, int(0.95 * n_bootstraps)]
    tprMean = np.mean(tprGridMat, axis=1)

    plt.hold(True)
    ax = plt.gca()  # kwargs.pop('ax', plt.gca())
    base_line, = ax.plot(fprGridVec, tprMean, '-', linewidth=4, label=labelIn)
    ax.fill_between(fprGridVec, tprLow05, tprTop95, facecolor=base_line.get_color(), alpha=0.2)
    
    
    
def getFSnodes(connFlIn=None):
    """
    Return nodes indexes and names according to connFlIn
    default file: /data/repo/mrtrix3/src/connectome/tables/fs_default.txt
    """
    
    if connFlIn is None:
        #connFlIn = '/data/repo/mrtrix3/src/connectome/tables/fs_default.txt'
        connFlIn = './fs_default.txt'

    if connFlIn is 'PD25':
        connFlIn = './PD25.txt'
    connFr = pd.read_csv( connFlIn, sep="\s+", skiprows=3, header=None, names=['areaShort','area','x','y','z','c']  )
    
    return connFr

def plotTSNE(Xin, y):
    tnseMod = TSNE(n_components=2, random_state=540, perplexity=30)
    np.set_printoptions(suppress=True)

    #- scale data
    scaler = pre.RobustScaler( with_centering=True, with_scaling=False, quantile_range=(25.0, 75.0) )
    Xscl = scaler.fit_transform( Xin )
    #-

    # PCA
    Xscl2 = PCA( ).fit_transform( Xscl )

    resTsne = tnseMod.fit_transform(Xscl2)
    q = resTsne
    color = y
    plt.figure( figsize=(12,12) )
    plt.scatter( q[:, 0], q[:, 1], c=color, s=200,cmap=plt.cm.Spectral)
    plt.colorbar()
    
def makeSym(matIn):
    """
    Make symmetric matrix
    """
    # make sure the diagonal is computed only once
    matOut = matIn.T + np.triu( matIn, 1 )
    
    return matOut
    
    
    

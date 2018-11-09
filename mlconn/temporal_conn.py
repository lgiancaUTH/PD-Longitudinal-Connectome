
#imports
import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt
import seaborn

import utils
import utilsStats

from sklearn.model_selection import KFold, StratifiedKFold, GridSearchCV, RandomizedSearchCV
from scipy.stats import randint as sp_randint
from sklearn.covariance import GraphLassoCV
from scipy.interpolate import interp1d
import sklearn.linear_model as lm
# from sklearn.metrics import roc_auc
import sklearn.svm as svm
import sklearn.ensemble as ens
import sklearn.preprocessing as pre
import copy



# Load configuration
CONF_FILE = 'data/configuration.json'
CONF = utils.readConf(CONF_FILE)

def genImgInfo( fileFullInfo, getLastBool=True ):
    firstLast = 'first'
    if getLastBool == True:
        firstLast = 'last'
    elif getLastBool is 'Second':
        firstLast = 'second'

    # load/generate image info for prodromal
    imgInfoFr = utils.getImgsInfoFr(fileFullInfo, 'Prodromal', 'MRI', 'GRAPPA', getLast=getLastBool)

    #=== Experiment Last MRI Prodromal matching PD and CNT
    # load/generate image info

    prodrFr = utils.getImgsInfoFr(fileFullInfo, 'Prodromal', 'MRI', 'GRAPPA', getLast=getLastBool)
    # match PD
    pdFr = utils.getImgsInfoFr(fileFullInfo, 'PD', 'MRI', 'GRAPPA', getLast=getLastBool)
    newPdFr = utils.matchSubj(prodrFr, pdFr, 2)
    # match Controls
    cntFr = utils.getImgsInfoFr(fileFullInfo, 'Control', 'MRI', 'GRAPPA', getLast=getLastBool)
    newCntFr = utils.matchSubj(prodrFr, cntFr, 2)
    # join files
    maCntPdFr =  newPdFr.append(newCntFr)
    imgInfoFr = imgInfoFr.append( maCntPdFr )

    # add freesurfer name ppmi_[id]_last
    imgInfoFr['fsName'] = imgInfoFr['Subject'].apply(lambda x:  'ppmi_' + str(int(x)) + '_' + firstLast )
    # format acquisition date
    imgInfoFr['adate'] = pd.to_datetime( imgInfoFr['Acq Date'], format="%m/%d/%Y" )
    
    return imgInfoFr


def genImgInfoPD( fileFullInfo, getLastBool=True ):
    firstLast = 'first'
    if getLastBool:
        firstLast = 'last'
    # load/generate image info for prodromal
    # imgInfoFr = utils.getImgsInfoFr(fileFullInfo, 'Prodromal', 'MRI', 'GRAPPA', getLast=getLastBool)

    # === Experiment Last MRI Prodromal matching PD and CNT
    # load/generate image info

    prodrFr = utils.getImgsInfoFr(fileFullInfo, 'Prodromal', 'MRI', 'GRAPPA', getLast=getLastBool)
    # match PD
    pdFr = utils.getImgsInfoFr(fileFullInfo, 'PD', 'MRI', 'GRAPPA', getLast=getLastBool)
    newPdFr = utils.matchSubj(prodrFr, pdFr, 4)
    # match Controls
    cntFr = utils.getImgsInfoFr(fileFullInfo, 'Control', 'MRI', 'GRAPPA', getLast=getLastBool)
    newCntFr = utils.matchSubj(prodrFr, cntFr, 3)
    # join files
    imgInfoFr = prodrFr.append(newPdFr)
    imgInfoFr = imgInfoFr.append(newCntFr)

    # ===
    # get all DTI
    dtiInfoFr = utils.getImgsInfoFr(fileFullInfo, None, 'DTI', 'DTI', getLast=getLastBool)

    # match file names based on imgInfoFr. Set main fName for DTI
    imgInfoFr2 = imgInfoFr.set_index(['Subject', 'Visit'], verify_integrity=True) \
        .join(dtiInfoFr.set_index(['Subject', 'Visit'], verify_integrity=True)['fName'], lsuffix='_T1')
    imgInfoFr2 = imgInfoFr2.reset_index()

    # check if there are unmatched label between first DTI and first MRI
    # (probably because the first visit had only MRI)
    unmatchedLbl = imgInfoFr2.fName.isnull()
    if np.any(unmatchedLbl):
        # remove volumes if unmatched
        imgInfoFr2 = imgInfoFr2[~imgInfoFr2.fName.isnull()]
        print 'removing ', np.sum(unmatchedLbl), ' volumes because DTI unmatched!'

    # add freesurfer name ppmi_[id]_last
    imgInfoFr2['fsName'] = imgInfoFr2['Subject'].apply(lambda x:  'ppmi_' + str(int(x)) + '_' + firstLast )
    # format acquisition date
    imgInfoFr2['adate'] = pd.to_datetime( imgInfoFr2['Acq Date'], format="%m/%d/%Y" )

    return imgInfoFr2


def getImgInfoFL(fileFullInfoIn=None):
    """
    Get Dataframe with first and last visit
    """
    
    # select first or last temporal info
    fileFullInfo = fileFullInfoIn
    if fileFullInfo is None:
        fileFullInfo = 'full_DTI_only_3_17_2017.csv'

    imgInfoFirstFr = genImgInfo( fileFullInfo, getLastBool=False ) #
    imgInfoLastFr = genImgInfo( fileFullInfo, getLastBool='Second') # or 'Second' or True

    imgInfoFr = imgInfoFirstFr.set_index('Subject')[ ['Visit', 'fsName', 'Age', 'Group', 'Sex', 'adate'] ].join(imgInfoLastFr.set_index('Subject')[ ['Visit', 'fsName', 'Age', 'adate'] ], rsuffix='_lst' )
    # calculate difference in days
    imgInfoFr['daysDiff'] = (imgInfoFr.adate_lst - imgInfoFr.adate).dt.days
    
    return imgInfoFr



    
    
def runMLord( X, y, Xout, paramIn=None ):
    """
    Run machine learning models.
    Ensures that the output will follow the initial sample ordering
    :param X:
    :param y:
    :param Xout:
    :param **paramIn: if set, paramIn['N_SPLITS'] -> number of splits
                              paramIn['NUM_RAND_WEIGHTS'] -> number of random iteration for permutation test
    
    :return: (model dictionary, gsDic)
    """

    N_SPLITS = 21 
    RND_SEED= 6543215468

    #estimate random weights if > 0. Number represent the iterations per fold
    NUM_RAND_WEIGHTS = 65

    SDA_COMP = None # None equals to n_classes-1
    SDA_FEAT = 0

    # set parameters from paramIn
    if paramIn is not None:
        if 'N_SPLITS' in paramIn.keys():
            N_SPLITS = paramIn['N_SPLITS']
            print 'N_SPLITS', N_SPLITS
        if 'NUM_RAND_WEIGHTS' in paramIn.keys():
            NUM_RAND_WEIGHTS = paramIn['NUM_RAND_WEIGHTS']
            print 'NUM_RAND_WEIGHTS', NUM_RAND_WEIGHTS

    #--- Validation
    skf = StratifiedKFold(n_splits=N_SPLITS, random_state=RND_SEED)

    #MODEL DICTIONARY
    mdlDic = {}
    # liblinear coordinate descent
    mdlDic[0] = {'name': 'Logistic Regression (L1 reg.)', 'auc': [], 'scores': [], 'featWeight': [], 'rndFeatWeights': [], 'scoresOut': [], \
                 'model': lm.LogisticRegression(penalty='l1'), 'y': []}
    # liblinear coordinate descent
    mdlDic[1] = {'name': 'Logistic Regression (L2 reg.)', 'auc': [], 'scores': [], 'featWeight': [], 'rndFeatWeights': [], 'scoresOut': [],\
                 'model': lm.LogisticRegression(penalty='l2'), 'y': []}
    # liblinear coordinate descent
    mdlDic[2] = {'name': 'Elastic Net', 'auc': [], 'scores': [], 'featWeight': [], 'rndFeatWeights': [], 'scoresOut': [],\
                 'model': lm.ElasticNet(alpha=1, l1_ratio=0.005, fit_intercept=True, normalize=False), 'y': []}

    mdlDic[3] = {'name': 'Linear SVM', 'auc': [], 'scores': [], 'featWeight': [], 'rndFeatWeights': [], 'scoresOut': [],\
                 'model': svm.SVC(kernel="linear", C=0.0001, probability=True), 'y': []} 

    mdlDic[4] = {'name': 'Random Forest Classifier', 'auc': [], 'scores': [], 'featWeight': [], 'rndFeatWeights': [], 'scoresOut': [], \
                 'model': ens.RandomForestClassifier(max_depth=100, n_estimators=100), 'y': []}


    for train_index, test_index in skf.split(range(len(y)), y):
        # split
        trainX = X[train_index, :]
        trainY = y[train_index]
        testX = X[test_index, :]
        testY = y[test_index]

        # standardize #just for opt purposes
    #     scaler = pre.StandardScaler( with_mean=True, with_std=True )
        scaler = pre.RobustScaler( with_centering=True, with_scaling=True, quantile_range=(25.0, 75.0) )
        scaler.fit( trainX )
        trainX = scaler.transform( trainX )
        testX = scaler.transform( testX )
        testXout = scaler.transform( Xout ) # out of cross validation

        for mId in mdlDic:

            mdlDic[mId]['model'].fit(trainX, trainY) #we get the parameters for the model

            # detect probabilities (using predict_proba when available)
            p = None
            pOut = None
            if 'predict_proba' in dir(mdlDic[mId]['model']):
                p = mdlDic[mId]['model'].predict_proba(testX)[:,1]
                pOut = mdlDic[mId]['model'].predict_proba(testXout)[:,1]
            else:
                p = mdlDic[mId]['model'].predict(testX)
                pOut =  mdlDic[mId]['model'].predict(testXout)

            # store coefficients (if available)
            if 'coef_' in dir(mdlDic[mId]['model']):
                if not mId==6:
                    mdlDic[mId]['featWeight'].append( mdlDic[mId]['model'].coef_.flatten() )
                #print np.array(mdlDic[mId]['featWeight'])

                #estimate
                if NUM_RAND_WEIGHTS: 
                    # generate a copy of the model (not reference!)
                    tmpModel = copy.deepcopy(mdlDic[mId]['model'])
                    for idxEst in range(NUM_RAND_WEIGHTS): #this is used later on for the permutation test
                        # random labels
                        yRand = np.random.randint( 0, 2, trainY.shape[0] )
                        # fit
                        tmpModel.fit( trainX, yRand )
                        # store weights
                        if not mId == 6:
                            mdlDic[mId]['rndFeatWeights'].append( tmpModel.coef_.flatten() ) #flatten-> collapse into 1D array


            # init arrays
            if mdlDic[mId]['scores'] == []:
                mdlDic[mId]['scores'] = np.ones(len(y))*-1
            if mdlDic[mId]['y'] == []:
                mdlDic[mId]['y'] = np.ones(len(y))*-1
                
            # store prediction and Ys to the right positions (if a smaple is tested twice it will be replaced)
            mdlDic[mId]['scores'][test_index] = p
            mdlDic[mId]['y'][test_index] = y[test_index]
            
            # stack on 0-axis since they are the predictions on the same samples
            if mdlDic[mId]['scoresOut'] == []:
                mdlDic[mId]['scoresOut'] = pOut
            else:
                mdlDic[mId]['scoresOut'] = np.vstack( (mdlDic[mId]['scoresOut'], pOut) )

        


    return mdlDic
    #---


def printMLperf( mdlDic, yCvArr, stdInfo=True, latexInfo=False ):
    """
    Print performance metric 
    :param mdlDic:
    :param yCvArr:
    :param stdInfo: set true to plot standard info
    :param latexInfo: set true to format it as a table
    :return:
    """
    for mId in mdlDic:
        print '-' * 20, mdlDic[mId]['name']

        # global AUC/significance
        scoresArr = mdlDic[mId]['scores']
        # stats
        aucStr= utilsStats.sigTestAUC( scoresArr[yCvArr==0], scoresArr[yCvArr==1], disp='long' )
        (acc, sens, spec, roc_auc, cutoffTh, cfMat, kappa ) = utilsStats.findCutoffPnt3(scoresArr[yCvArr==1], scoresArr[yCvArr==0])
        # balanced accuracy
        balAcc = (sens+spec)/2
        # store aucStr
        mdlDic[mId]['aucStr'] = aucStr

        print aucStr
        if stdInfo:
            print 'roc_auc_t: {:0.3f}, sens: {:0.3f}, spec: {:0.3f}, acc: {:0.3f}, kappa: {:0.3f}, bal. acc: {:0.3f}, cutoffTh: {:0.3f}'.format(roc_auc, sens, spec, acc, kappa, balAcc, cutoffTh )
        if latexInfo:
            print mdlDic[mId]['name']+' & {:0.3f} & {:0.3f} & {:0.3f} & {:0.3f} & {:0.3f} \\ \hline'.format(roc_auc, sens, spec, acc, kappa )

        
def printMLperfOutScores( mdlDic, yCvArr, stdInfo=True, latexInfo=False ):
    """
    Print performance metric comparing class 0 with scores left out  mdlDic[x]['scoresOut']
    :param mdlDic:
    :param yCvArr:
    :param stdInfo: set true to plot standard info
    :param latexInfo: set true to format it as a table
    :return:
    """
    for mId in mdlDic:
        print '-' * 20, mdlDic[mId]['name']


	# add out of cross validation predictions
	scoresOutArr = np.median( mdlDic[mId]['scoresOut'] , 0)

        # global AUC/significance
        scoresArr = mdlDic[mId]['scores']
        # stats
        aucStr= utilsStats.sigTestAUC( scoresArr[yCvArr==0], scoresOutArr, disp='long' )
        (acc, sens, spec, roc_auc, cutoffTh, cfMat, kappa ) = utilsStats.findCutoffPnt3(scoresOutArr, scoresArr[yCvArr==0])
        # balanced accuracy
        balAcc = (sens+spec)/2
        
        # store aucStr
        mdlDic[mId]['aucStr'] = aucStr

        print aucStr
        if stdInfo:
            print 'roc_auc_t: {:0.3f}, sens: {:0.3f}, spec: {:0.3f}, acc: {:0.3f}, kappa: {:0.3f}, bal. acc: {:0.3f}, cutoffTh: {:0.3f}'.format(roc_auc, sens, spec, acc, kappa, balAcc, cutoffTh )
        if latexInfo:
            print mdlDic[mId]['name']+' & {:0.3f} & {:0.3f} & {:0.3f} & {:0.3f} & {:0.3f} \\ \hline'.format(roc_auc, sens, spec, acc, kappa )

        
        
        
        
# if __name__ == "__main__":
#    #GT
#    imgInfoFr = getImgInfoFL('PD25')
#
#    print( len(imgInfoFr) )
#    pass

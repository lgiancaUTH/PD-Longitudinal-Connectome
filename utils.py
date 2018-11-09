"""
Common utilities to analyse PPMI data

"""

import numpy as np
import os
import re
import pandas as pd
import time
import json

# DICOM extension use lowercase
DICOM_EXT = 'dcm'


def readConf(confFile):
    """
    Read configuration
    :param confFile: file name
    :return: configuration dictionary
    """
    config = None
    with open(confFile, 'r') as f:
        config = json.load(f)

    return config


# Load configuration
CONF_FILE = 'data/configuration.json'
CONF = readConf(CONF_FILE)



def renameDICOMFiles(dirIn, tmpDir='/tmp'):
    """
    Rename DICOM files to allow for alphabetical ordering (using a tmp directory). It assumes ".dcm" extension set by (DICOM_EXT)
    :param dirIn:directory containing the files
    :param tmpDir: tmp directory for processing 
    
    :return: none
    """
    filesLst = os.listdir(dirIn)
    print 'Renaming DICOM files...'
    for fName in filesLst:
        if fName[-3:].lower() != DICOM_EXT:
            continue

        p = re.compile('^\S+_raw_\d+_(\d+).+$')
        m = p.match(fName)
        # make sure the pattern matches
        assert( m is not None )
        # new file starting with number
        fNameNew = m.group(1) + '.' + DICOM_EXT
        # rename file to new directory
        os.system('cp ' + dirIn + '/' + fName + ' ' + tmpDir + '/' + fNameNew)


    print '     Renamed'

    pass


def convertToMif(dirIn, outFile='outDti.mif'):
    """
    Convert to Mif (nifty format) 
    :rtype: object
    :param dirIn: directory with dcm files
    :return: 
    """

    outFile[len(outFile)-6:len(outFile)]
    if not outFile[len(outFile)-6:len(outFile)] in ('T1.mif','T2.mif'): #to convert to DTI.mif
        print 'ConvertToMif start... DTI'
        # firstFile = os.listdir(dirIn)[0]
        # get gradients
        os.system( 'mrinfo ' + dirIn  + ' -export_grad_mrtrix '+dirIn+'grad.b' + ' -quiet')
        # mrinfo S179002 -export_grad_mrtrix grad.b
        # convert to Nifty
        os.system('mrconvert ' + dirIn + '[].'+DICOM_EXT + ' -grad '+dirIn+'grad.b ' + outFile + ' -force -quiet')
        # mrconvert S179002_2/PPMI_3102_[]_S179002_I353459.dcm -grad grad.b comb3.mif
    else: #to convert to T1.mif
        print 'ConvertToMif start... ' + outFile[len(outFile)-6:len(outFile)-4]
        os.system('mrinfo ' + dirIn  + ' -quiet')
        # convert to Nifty
        os.system('mrconvert ' + dirIn + ' ' + outFile + ' -quiet')

    print '     Converted'




class Param:
    """
    Class describing the parameters for tractography, registration connectome
    """

    def __init__(self):
        self.pName = '0'
        # Tractography parameters
        self.trkgenCmd = '-maxlength 250 -number 10M -cutoff 0.06'
        self.tcksiftCmd = '-term_number 1M'
        self.fslSkullRm = '-f 0.2'

def genTractography( dirIn, dirOut, dirFS, pIn, check=False ):
    """
    Generate tractography on dirOut

    :param dirIn: PPMI directory with DTI info 
    :param dirOut: output directory
    :param dirFS: Freesufer directory containing the computed information for the subject
    :param pIn: parameters as defined by Param class
    :param check: set to True to check if connectome exist
    :return: None (if check is true returns True, False depending if the connectome has been computed)
    """


    #==== Set file names
    dtiVolFname = dirOut + pIn.pName + '_outDti.mif' #input DWI
    dtiVolFnameNii = dirOut + pIn.pName + '_outDti.nii.gz'
    maskFname = dirOut + pIn.pName + '_outDti_mask.mif' #brain mask
    dtiRespFname = dirOut + pIn.pName + '_outDti_resp.txt' #response function
    dtiFodFname = dirOut + pIn.pName + '_outDti_fod.mif' #fibre orientation distribution estimation
    dtiTrkFname = dirOut + pIn.pName + '_outDti_trk.tck' #tractography file
    dtiSiftTrkFname = dirOut + pIn.pName + '_outDti_sift_trk.tck' #tractography file with reduced streamline count-> more biologically meaningful estimates
    dtiSiftTrkvisFname = dirOut + pIn.pName + '_outDti_sift_trk.trk'
    dirTmp = dirOut + 'tmp/'
    mriFnameNii = dirOut + 'brain.nii.gz'
    dtiVolB0FnameNii = dirOut + pIn.pName + '_outDti-b0.nii.gz' # DTI b0
    dtiVolB0nskFnameNii = dirOut + pIn.pName + '_outDti-b0_sk.nii.gz' # no skull
    dtiVolB0regFnameNii = dirOut + pIn.pName + '_outDti-b0-reg.nii.gz'  # name WRONG! it is the MRI registered to DTI b0 space
    regMatMri2b0 = dirOut + pIn.pName + '_brain-to-dti_dec.mat'  # registration matrix
    regMatMri2b0f = dirOut + pIn.pName + '_brain-to-dti_dec2.mat' # registration matrix (fixed)
    outNodesVol = dirOut + 'nodes-aparc+aseg.nii.gz'
    outNodesVol2Dti = dirOut +  pIn.pName + '_nodes-aparc+aseg-regDti.nii.gz'
    outConnectome = dirOut +  pIn.pName + '_connectome.csv'
    # ====

    # check for connectome existance
    if check:
        if os.path.exists(outConnectome):
            return True
        else:
            return False
    # stop execution if check parameter is True


    # create tmp dir (and parent directory) if needed
    if not os.path.exists(dirTmp):
        os.makedirs(dirTmp)

    renameDICOMFiles(dirIn, dirTmp)
    if not os.path.exists(dtiVolFname):
        convertToMif(dirTmp, dtiVolFname)
    # clean up tmp dir
    os.system('rm ' + dirTmp + '/*.' + DICOM_EXT)
    os.system('rm ' + dirTmp + '/grad.b')

    #============================ Tractography
    os.system('dwi2mask ' + dtiVolFname + ' ' + maskFname + ' -force')
    os.system('/data/repo/mrtrix3/scripts/dwi2response tournier ' + dtiVolFname + ' ' + dtiRespFname + ' -force')
    os.system('dwi2fod csd ' + dtiVolFname + ' ' + dtiRespFname + ' ' + dtiFodFname + ' -mask ' + maskFname + ' -force')
    os.system(
        'tckgen ' + dtiFodFname + ' ' + dtiTrkFname + ' -seed_image ' + maskFname + ' -mask ' + maskFname + ' '
        + pIn.trkgenCmd + ' -force')
    os.system('tcksift ' + dtiTrkFname + ' ' + dtiFodFname + ' ' + dtiSiftTrkFname + ' ' + pIn.tcksiftCmd + ' -force')
    #============================

    #============================ Registration
    # convert MRI file from Freesurfer dir
    os.system('mrconvert ' + dirFS + '/brain.mgz '+ mriFnameNii + ' -force')
    # convert dti
    os.system('mrconvert ' + dtiVolFname + ' ' + dtiVolFnameNii + ' -force')
    # get b0 only
    os.system('fsl5.0-fslroi '+dtiVolFnameNii+' ' + dtiVolB0FnameNii + ' 0 1')
    # remove skull
    os.system('fsl5.0-bet ' + dtiVolB0FnameNii + ' ' + dtiVolB0nskFnameNii + ' ' + pIn.fslSkullRm)
    # register t1 to DTI
    os.system( 'fsl5.0-flirt -in '+mriFnameNii+' -ref '+dtiVolB0nskFnameNii+' -o '+dtiVolB0regFnameNii+
               ' -omat '+regMatMri2b0+' -usesqform' )
    # fix fsl registration matrix
    #os.system('/home/lgianca/bin/cnv2dec.sh '+regMatMri2b0+' > ' + regMatMri2b0f)
    os.system('/home/openanogales/bin/cnv2dec.sh '+regMatMri2b0+' > ' + regMatMri2b0f)
    # convert labelled data (from aparc+aseg)
    os.system(  'labelconvert '+dirFS+'/aparc+aseg.mgz /usr/local/freesurfer/FreeSurferColorLUT.txt ' +
                '/home/openanogales/ppmi/src/fs_default.txt ' + outNodesVol )
    # apply registration to labels
    os.system( 'fsl5.0-flirt -in '+outNodesVol+' -applyxfm -init '+regMatMri2b0f+ ' -out '+outNodesVol2Dti+
               ' -paddingsize 0.0 -interp nearestneighbour -ref ' + dtiVolB0nskFnameNii )
    # map streamlines to the parcellated image to produce a connectome
    os.system( 'tck2connectome '+dtiSiftTrkFname+' '+outNodesVol2Dti+' ' + outConnectome )
    #============================

    # convert tracts to trackvis format
    cnvMR2trackvis(dtiSiftTrkFname, dtiVolB0FnameNii, dtiSiftTrkvisFname)

    pass


def regATLAS( dirInT1, dirInT2, dirInDTI, dirOut, dirInAux, pIn, check=False, connectome='global' ):
    """
    Registrate ATLAS to image on dirOut

    :param dirInT1: PPMI directory with image T1 info
    :param dirInT2: PPMI directory with image T2 info
    :param dirInDTI: PPMI directory with image DTI info
    :param dirOut: output directory
    :param check: set to True to check if connectome exist
    :param connectome: decide which connectome to use, global or local - only of the PD area.
    :return: None (if check is true returns True, False depending if the connectome has been computed)
    """

    # ==== Set file names
    dirTmp = dirOut + 'tmp/'
    dtiSiftTrkFnameAux = dirInAux + pIn.pName + '_outDti_sift_trk.tck'
    dtiVolFnameNiiAux = dirInAux +  pIn.pName + '_outDti.nii.gz'

    dtiVolFname = dirOut + pIn.pName + '_outDti.mif'  # input DWI
    dtiVolFnameNii = dirOut + pIn.pName + '_outDti.nii.gz'
    maskFname = dirOut + pIn.pName + '_outDti_mask.mif'  # brain mask
    maskFnameNii = dirOut + pIn.pName + '_outDti_mask.nii.gz'
    dtiRespFname = dirOut + pIn.pName + '_outDti_resp.txt'  # response function
    dtiFodFname = dirOut + pIn.pName + '_outDti_fod.mif'  # fibre orientation distribution estimation
    dtiTrkFname = dirOut + pIn.pName + '_outDti_trk.tck'  # tractography file
    dtiTrkLocalFname = dirOut + pIn.pName + '_outDtiLocal_trk.tck'
    dtiSiftTrkFname = dirOut + pIn.pName + '_outDti_sift_trk.tck'  # tractography file with reduced streamline count-> more biologically meaningful estimates
    dtiSiftTrkLocalFname = dirOut + '2_outDti_sift_trk_local.tck' #local tractography

    dtiSiftTrkLocalFnameTrim = dirOut + '2_outDti_sift_trk_local_trim70.tck'  # local tractography trim
    dtiSiftTrkLocalFnameTrim_aux = dirOut + '2_outDti_sift_trk_local_trim70_aux.tck'
    dtiSiftTrkLocalFnameShiftTrim = dirOut + '2_outDti_sift_trk_local_ShiftTrim.tck'

    dtiSiftTrkLocalFnameShiftTrim2 = dirOut + '2_outDti_sift_trk_local_ShiftTrim3.tck'

    TrkLocal_samples = dirOut + pIn.pName + '_trk_local_samples.tck'  # local tractography
    dtiSiftTrkvisFname = dirOut + pIn.pName + '_outDti_sift_trk.trk'
    dtiSiftTrkvisLocalFname = dirOut + pIn.pName + '_outDti_sift_trk_local.trk'
    dtiVolB0FnameNii = dirOut + pIn.pName + '_outDti-b0.nii.gz'  # DTI b0
    dtiVolB0nskFnameNii = dirOut + pIn.pName + '_outDti-b0_sk.nii.gz'  # no skull
    dtiVolB0MasknskFnameNii = dirOut + pIn.pName + '_outDti-b0_sk_mask.nii.gz'  # no skull mask
    regMatMri2b0 = dirOut + pIn.pName + '_brain-to-dti_dec.mat'  # registration matrix
    regMatMri2b0f = dirOut + pIn.pName + '_brain-to-dti_dec2.mat'  # registration matrix (fixed)

    dtiTensorVolFname = dirOut + pIn.pName + '_outDtensor.mif' #tensor map
    FaMapFname = dirOut + pIn.pName + '_outFA.mif' #FA map
    RdMapFname = dirOut + pIn.pName + '_outRD.mif' #RD map
    AdcMapFname = dirOut + pIn.pName + '_outADC.mif' #MD map
    AdMapFname = dirOut + pIn.pName + '_outAD.mif' #AD map
    basenameFSL = dirOut + pIn.pName + '_outFSL' #out put basename for all fsl files creates when computing FA and MD maps
    FaMapFnameFSLNii = dirOut + pIn.pName + '_outFSL_FA.nii.gz' #FA map computed with FSL
    FaMapFnameFSL = dirOut + pIn.pName + '_outFSL_FA.mif'
    AdcMapFnameFSLNii = dirOut + pIn.pName + '_outFSL_MD.nii.gz'  #MD map computed with FSL
    AdcMapFnameFSL = dirOut + pIn.pName + '_outFSL_MD.mif'
    MeanFaMapFname = dirOut + '2_outmeanFA.csv'  # mean track FA map
    MeanRdMapFname = dirOut + '2_outmeanRD.csv'  # mean track RD map
    MeanAdcMapFname = dirOut + '2_outmeanADC.csv'  # mean track MD map
    MeanAdMapFname = dirOut + '2_outmeanAD.csv'  # mean track AD map
    MeanFaMapFnameFSL = dirOut + '2_outmeanFSL_FA.csv'  # mean track FA map FSL
    MeanAdcMapFnameFSL = dirOut + '2_outmeanFSL_MD.csv'  # mean track MD map FSL

    MeanFaMapFnameTrim = dirOut + '2_outmeanFA_trim70.csv'  # mean track FA map trim
    MeanRdMapFnameTrim = dirOut + '2_outmeanRD_trim70.csv'  # mean track RD map trim
    MeanAdcMapFnameTrim = dirOut + '2_outmeanADC_trim70.csv'  # mean track MD map trim
    MeanAdMapFnameTrim = dirOut + '2_outmeanAD_trim70.csv'  # mean track AD map trim
    MeanFaMapFnameFSLTrim = dirOut + '2_outmeanFSL_FA_trim70.csv'  # mean track FA map trim FSL
    MeanAdcMapFnameFSLTrim = dirOut + '2_outmeanFSL_MD_trim70.csv'  # mean track FA map trim FSL
    MeanFaMapFnameFSLShiftTrim = dirOut + '2_outmeanFSL_FA_ShiftTrim.csv'  # mean track FA map trim FSL
    MeanAdcMapFnameFSLShiftTrim = dirOut + '2_outmeanFSL_MD_ShiftTrim.csv'  # mean track FA map trim FSL

    MeanFaMapFnameFSLShiftTrim2 = dirOut + '2_outmeanFSL_FA_ShiftTrim3.csv'  # mean track FA map trim FSL
    MeanAdcMapFnameFSLShiftTrim2 = dirOut + '2_outmeanFSL_MD_ShiftTrim3.csv'  # mean track FA map trim FS

    vecMapFname = dirOut + pIn.pName + '_outVec.mif' #vector
    vecMapFnameNii = dirOut + pIn.pName + '_outVec.nii.gz' #for visualization with RGB.


    outConnectome = dirOut + pIn.pName + '_connectome.csv' # connectome
    outConnectomeLocal = dirOut + '2_connectomeLocal.csv'  # Local connectome
    outNodes2fnirtVolConnectome = dirOut + pIn.pName + '_outNodes-2fnirt-Volconnectome.csv'  # Volume Local connectome
    outNormVolConnectomeLocal = dirOut + '2_NormVolconnectomeLocal.csv'  # Local connectome normalize with the inverse volume of the two nodes
    outFaConnectomeLocal = dirOut + '1_FaconnectomeLocal.csv'  # Local FA connectome
    outRdConnectomeLocal = dirOut + '1_RdconnectomeLocal.csv'  # Local RD connectome
    outAdcConnectomeLocal = dirOut + '1_AdcconnectomeLocal.csv'  # Local ADC connectome
    outAdConnectomeLocal = dirOut + '1_AdconnectomeLocal.csv'  # Local AD connectome
    outFaConnectomeLocalFSL = dirOut + '2_FaconnectomeLocalFSL.csv'  # Local FA connectome FSL
    outAdcConnectomeLocalFSL = dirOut + '2_AdcconnectomeLocalFSL.csv'  # Local ADC connectome FSL

    outConnectomeLocalTrim = dirOut + '2_connectomeLocal_trim70.csv'  # Local connectome trim
    outConnectomeLocalShiftTrim = dirOut + '2_connectomeLocal_ShiftTrim.csv'
    outConnectomeLocalShiftTrim2 = dirOut + '2_connectomeLocal_ShiftTrim3.csv'
    outNormVolConnectomeLocalTrim = dirOut + '2_NormVolconnectomeLocal_trim70.csv'  # Local connectome trim normalize with the inverse volume of the two nodes
    outNormVolConnectomeLocalShiftTrim = dirOut + '2_NormVolconnectomeLocal_ShiftTrim.csv'
    outNormVolConnectomeLocalShiftTrim2 = dirOut + '2_NormVolconnectomeLocal_ShiftTrim3.csv'

    outFaConnectomeLocalTrim = dirOut + '2_FaconnectomeLocal_trim70.csv'  # Local FA connectome trim
    outRdConnectomeLocalTrim = dirOut + '2_RdconnectomeLocal_trim70.csv'  # Local RD connectome trim
    outAdcConnectomeLocalTrim = dirOut + '2_AdcconnectomeLocal_trim70.csv'  # Local ADC connectome trim
    outAdConnectomeLocalTrim = dirOut + '2_AdconnectomeLocal_trim70.csv'  # Local AD connectome trim
    outFaConnectomeLocalFSLTrim = dirOut + '2_FaconnectomeLocalFSL_trim70.csv'  # Local FA connectome trim FSL
    outAdcConnectomeLocalFSLTrim = dirOut + '2_AdcconnectomeLocalFSL_trim70.csv'  # Local MD connectome trim FSL
    outFaConnectomeLocalFSLShiftTrim = dirOut + '2_FaconnectomeLocalFSL_ShiftTrim.csv'
    outAdcConnectomeLocalFSLShiftTrim = dirOut + '2_AdcconnectomeLocalFSL_ShiftTrim.csv'

    outFaConnectomeLocalFSLShiftTrim2 = dirOut + '2_FaconnectomeLocalFSL_ShiftTrim3.csv'
    outAdcConnectomeLocalFSLShiftTrim2 = dirOut + '2_AdcconnectomeLocalFSL_ShiftTrim3.csv'

    dirATLAS = '/data/PD25/mni_PD25_20170213_nifti/' #atlas directory
    T1VolFname = dirOut + pIn.pName + '_outT1.mif' #input T1 GRAPPA volum [MIF]
    T1VolFnameNii = dirOut + pIn.pName + '_outT1.nii.gz' #input T1 GRAPPA volum [NIFTI]
    T1VolnskFnameNii = dirOut + pIn.pName + '_outT1_sk.nii.gz' #input T1 GRAPPA volum without skull [NIFTI] - Final one.
    T1VolnskFnameNii2 = dirOut + pIn.pName + '_outT1_sk_2.nii.gz'  # input T1 GRAPPA volum without skull [NIFTI] - Intermidiate. Appling T2 brain mask.
    T1VolnskFnameDilNii = dirOut + pIn.pName + '_outT1_sk_dil.nii.gz'  # input T1 GRAPPA volum with some skull because it the mask was dilated [NIFTI]
    T1VolnskFnameMaskNii = dirOut + pIn.pName + '_outT1_sk_mask.nii.gz'  # input T1 GRAPPA mask without skull [NIFTI]
    T1VolregDTIb0FnameNii = dirOut + pIn.pName + '_outT1-reg-DTI-b0.nii.gz'  # it is the MRI (T1) flirt registered to DTI b0 space
    T1VolregFnirtWpDTIb0FnameNii = dirOut + pIn.pName + '_outT1-reg-fnirt-wp-T1.nii.gz' # it is the wp T1 fnirt registered to DTI b0 space
    T1VolregFnirtDTIb0FnameNii = dirOut + pIn.pName + '_outT1-reg-fnirt-DTI.nii.gz'  # it is the T1 fnirt registered to DTI b0 space.
    T1VolATLAS = dirATLAS + 'PD25-T1MPRAGE-template-1mm.nii.gz' #ATLAS volume with 1x1mm resolution
    maskVolATLAS = dirATLAS + 'PD25-atlas-mask-1mm.nii.gz' #ATLAS mask with 1x1mm resolution
    T1VolATLASnsk = dirATLAS + 'PD25-T1MPRAGE-template-1mm_sk.nii.gz' #ATLAS T1 volume without skull
    segVolATLAS = dirATLAS + 'PD25-subcortical-1mm.nii.gz' #ATLAS 8 subcortical structure segmentation with 1x1mm resolution
    outNodesVolregDTI = dirOut + pIn.pName + '_outNodes-reg-DTI.nii.gz'  # it is the segmentation registered to the DTI data.
    outNodes2fnirtVolregDTI = dirOut + pIn.pName + '_outNodes-2fnirt-reg-DTI.nii.gz'
    outNodes2fnirtHistregDTI = dirOut + pIn.pName + '_outNodes-hist-2fnirt-reg-DTI.txt'
    outNodesMaskVolregDTI = dirOut + pIn.pName + '_outNodesMask-reg-DTI.nii.gz'  # it is the segmentation Mask registered to the DTI data.
    outNodesMask2fnirtVolregDTI = dirOut + pIn.pName + '_outNodesMask-2fnirt-reg-DTI.nii.gz'
    T1ATLASVolregFlirtT1VolnskFnameMat = dirOut + pIn.pName + '_outATLAS-reg-flirt-T1.mat'  # it is the ATLAS (T1) registered to T1 with flirt. Registration matrix (fixed)
    T1ATLASVolregFlirtT1VolnskFnameMatf = dirOut + pIn.pName + '_outATLAS-reg-flirt-T1_2.mat'  # registration matrix (fixed)
    T1ATLASVolregFlirtT1VolnskFnameNii = dirOut + pIn.pName + '_outATLAS-reg-flirt-T1.nii.gz'  # it is the ATLAS (T1) registered to T1 with flirt.
    outNodesVolregFnirtT1VolnskFnameNii = dirOut + pIn.pName + '_outsegNodes-reg-fnirt-T1.nii.gz' # it is the segmented ATLAS registered to T1 with fnirt.
    T1ATLASVolregFnirtWpT1VolnskFnameNii = dirOut + pIn.pName + '_outATLAS-reg-fnirt-wp-T1.nii.gz'  # it is the warp coefficients of the ATLAS to T1 transformation
    T1ATLASVolregFnirtT1VolnskFnameNii = dirOut + pIn.pName + '_outATLAS-reg-fnirt-T1.nii.gz' # it is the T1 ATLAS fnirt registered to T1 space.
    T1ATLASVolnskregFnirtT1VolnskFnameNii = dirOut + pIn.pName + '_outATLASnsk-reg-fnirt-T1.nii.gz'  # it is the T1 ATLAS nsk fnirt registered to T1 space.
    T1ATLASVolregDTI = dirOut + pIn.pName + '_outATLAS-reg-DTI.nii.gz' #it is the ATLAS registered to the DTI data.
    T1ATLASVol2fnirtregDTI = dirOut + pIn.pName + '_outATLAS-2fnirt-reg-DTI.nii.gz'  # it is the ATLAS registered to the DTI data two fnirt steps.

    T2VolFname = dirOut + pIn.pName + '_outT2.mif' #input T2 volum [MIF]
    T2VolFnameNii = dirOut + pIn.pName + '_outT2.nii.gz'  # input T2 volum [NIFTI]
    T2VolnskFnameNii = dirOut + pIn.pName + '_outT2_sk.nii.gz'  # input T2 volum without skull [NIFTI]
    T2VolnskFnameMaskNii = dirOut + pIn.pName + '_outT2_sk_mask.nii.gz'  # input T2 mask without skull [NIFTI]
    outT2VolnskFnameMask2T1Nii = dirOut + pIn.pName + '_outT2_sk_mask_T1.nii.gz' #T2 brain mask applied ot the T1 space [NIFT]
    outT2VolnskFnameMask2T1DilNii = dirOut + pIn.pName + '_outT2_sk_mask_T1_dil.nii.gz'  # T2 brain mask applied to the T1 space - Dilated. [NIFT]
    outT2VolnskFnameMask2T1Mat = dirOut + pIn.pName + '_outT2_sk_mask_T1.mat' #T2 brain mask registration matrix applied ot the T1 space [NIFT]


    # check for registration existance
    if check:
        if connectome is 'global':
            if os.path.exists(outConnectome):
                return True
            else:
                return False
        else:
            if os.path.exists(outConnectomeLocal):
                return True
            else:
                return False



    # create tmp dir (and parent directory) if needed
    if not os.path.exists(dirTmp):
        os.makedirs(dirTmp)

    # ============================ Tractography # Only if it's the global one, if it's the local one we have to do it
    #after registering the ATLAS to the DTI.
    if connectome is 'global':
        print 'Start Global tractography...'
        if not os.path.exists(dtiSiftTrkFname): #First: we check if it was computed before.
            if not os.path.exists(dtiSiftTrkFnameAux): #same code than in genTractography()
                renameDICOMFiles(dirInDTI, dirTmp)
                if not os.path.exists(dtiVolFname):
                    convertToMif(dirTmp, dtiVolFname) #first DTI to MIF
                # clean up tmp dir
                os.system('rm ' + dirTmp + '/*.' + DICOM_EXT)
                os.system('rm ' + dirTmp + '/grad.b')

                os.system('dwi2mask ' + dtiVolFname + ' ' + maskFname + ' -force') #get brain mask DWI
                if not os.path.exists(dtiRespFname): #-force would give me problems....
                    os.system('/home/openanogales/mrtrix3/scripts/dwi2response tournier ' + dtiVolFname + ' ' + dtiRespFname ) #set tournier as the alg to compute the response function
                os.system('dwi2fod csd ' + dtiVolFname + ' ' + dtiRespFname + ' ' + dtiFodFname + ' -mask ' + maskFname + ' -force') #estimate fiber orientation distribution (FOD)
                os.system('tckgen ' + dtiFodFname + ' ' + dtiTrkFname + ' -seed_image ' + maskFname + ' -mask ' + maskFname + ' ' + pIn.trkgenCmd + ' -force') #get the streamlines tractography with the default algorithm (iFOD which is probabilistic).
                os.system('tcksift ' + dtiTrkFname + ' ' + dtiFodFname + ' ' + dtiSiftTrkFname + ' ' + pIn.tcksiftCmd + ' -force') #filter the fiber trackig dataset. Step suggested by MRtrix3
            else:
                os.system('cp ' + dtiSiftTrkFnameAux + ' ' + dtiSiftTrkFname) # bring the filtered tracktography to the current working folder.
        print '     Completed.'


    # ============================ Registrations
    # convert MRI files (T1 and T2) from PPMI DICOM dataset to NIFTI
    # T1
    if not os.path.exists(T1VolFnameNii):
        renameDICOMFiles(dirInT1, dirTmp)
        if not os.path.exists(T1VolFname):
            convertToMif(dirTmp, T1VolFname) #first T1 to MIF
        # clean up tmp dir
        os.system('rm ' + dirTmp + '/*.' + DICOM_EXT)
        os.system('mrconvert ' + T1VolFname + ' ' + T1VolFnameNii + ' -force -quiet') #then to NIFTI

    # T2
    if not os.path.exists(T2VolFnameNii):
        renameDICOMFiles(dirInT2, dirTmp)
        if not os.path.exists(T2VolFname):
            convertToMif(dirTmp, T2VolFname)  # first T2 to MIF
        # clean up tmp dir
        os.system('rm ' + dirTmp + '/*.' + DICOM_EXT)
        os.system('mrconvert ' + T2VolFname + ' ' + T2VolFnameNii + ' -force -quiet')  # then to NIFTI


    # extract brain from T1 data. Experimentally, I found that extracting the brain from T1 data somestimes does not work.
    # that's why 1) we check if the brain is extracted, 2) manually decide whether we extract it from the T1 or T2,
    # 3) if from T2, we extract the brain, register the mask to the T1, apply registration, dilate the mask,
    # apply registered mask to the T1 brain, extract the brain from this masked T1 volume.
    if not os.path.exists(T1VolnskFnameNii):
        T1Bool = False #Decide whether we use the T2 to extract the brain from T1 or not.
        if T1Bool:
            print 'Extract brain from T1...'
            os.system('fsl5.0-bet ' + T1VolFnameNii + ' ' + T1VolnskFnameNii + ' -m -R ') # remove skull of T1. The mask is given an automatic name. -R to compute the center of grabitiy several times.
            print '     Completed!.'
        else:
            print 'Extract brain from T2...'
            os.system('fsl5.0-bet ' + T2VolFnameNii + ' ' + T2VolnskFnameNii + ' -m -R ')  # remove skull of T2.
            os.system('fsl5.0-flirt -in ' + T2VolnskFnameMaskNii + ' -ref ' + T1VolFnameNii + ' -out ' + outT2VolnskFnameMask2T1Nii
                      + ' -omat ' + outT2VolnskFnameMask2T1Mat + ' -searchrx 0 0 -searchry 0 0 -searchrz 0 0 '
                      + '-2D -dof 12 -interp nearestneighbour') #Register mask to T1. We need to specify that images are alligned and that the interpolation should be done with nearestneighbour.
            os.system('fsl5.0-fslmaths ' + T1VolFnameNii + ' -mas ' + outT2VolnskFnameMask2T1Nii + ' ' + T1VolnskFnameNii2) #apply the mask to the T1 volume
            os.system('fsl5.0-fslmaths -dt char ' + outT2VolnskFnameMask2T1Nii + ' -kernel boxv 5 -dilD ' + outT2VolnskFnameMask2T1DilNii)  # dilate the mask. It is equivalent to make it bigger.
            os.system('fsl5.0-fslmaths ' + T1VolFnameNii + ' -mas ' + outT2VolnskFnameMask2T1DilNii + ' ' + T1VolnskFnameDilNii) #apply the dilated mask to the T1 volume
            os.system('fsl5.0-bet ' + T1VolnskFnameDilNii + ' ' + T1VolnskFnameNii + ' -m -R -f 0.3') #extract the brain from the T1.
            print '     Completed!.'


    # ========Register ATLAS to T1 image, no-linear registration
    if not os.path.exists(T1ATLASVolregFnirtWpT1VolnskFnameNii):
        print 'Start non-linear registration ATLAS to T1... '
        #1) We need to have the mask of the brain for the ATLAS and the T1 images.
        if not os.path.exists(T1VolATLASnsk):
            os.system('fsl5.0-fslmaths ' + T1VolATLAS + ' -mas ' + maskVolATLAS + ' ' + T1VolATLASnsk)
        #2) We perform a linear affine registration with flirt to have a good enough starting point in the non-linear registration.
        #       This step is really important to avoid local minima, if this registration fails, the non-linear registration is likely to fail as well.
        print '     Computing seed...'
        os.system('fsl5.0-flirt -in ' + T1VolATLASnsk + ' -ref ' + T1VolnskFnameNii + ' -omat ' + T1ATLASVolregFlirtT1VolnskFnameMat +
                   ' -o ' + T1ATLASVolregFlirtT1VolnskFnameNii + ' -usesqform')
        os.system('/home/openanogales/bin/cnv2dec.sh ' + T1ATLASVolregFlirtT1VolnskFnameMat + ' > ' + T1ATLASVolregFlirtT1VolnskFnameMatf)  # fix fsl registration matrix
        #3) We perform the non-linear registration as specified in the config file.
        print '     Starting fnirt...'
        t = time.time()
        os.system('fsl5.0-fnirt --ref=' + T1VolFnameNii + ' --refmask=' + T1VolnskFnameMaskNii + ' --aff=' + T1ATLASVolregFlirtT1VolnskFnameMatf +
                  ' --cout=' + T1ATLASVolregFnirtWpT1VolnskFnameNii + ' --config=/home/openanogales/ppmi/src/PD25_2_T1_1mm.cnf')
        elapsed = time.time() - t
        print '     Time: ' + str(elapsed/60) + ' min.'
        os.system('fsl5.0-applywarp --ref=' + T1VolFnameNii + ' --in=' + T1VolATLAS + ' --warp=' + T1ATLASVolregFnirtWpT1VolnskFnameNii +
                  ' --out=' + T1ATLASVolregFnirtT1VolnskFnameNii)
        os.system('fsl5.0-applywarp --ref=' + T1VolFnameNii + ' --in=' + segVolATLAS + ' --warp=' + T1ATLASVolregFnirtWpT1VolnskFnameNii +
                  ' --out=' + outNodesVolregFnirtT1VolnskFnameNii)

    print '     Completed.'


    # ========Register T1 to DTI image. Linear - rigid registration & non linear.
    print 'Start T1 to DTI registration....'
    if not os.path.exists(dtiVolFnameNii):
        if not os.path.exists(dtiVolFnameNiiAux):
            renameDICOMFiles(dirInDTI, dirTmp)
            convertToMif(dirTmp, dtiVolFname)  # first DTI to MIF

            # clean up tmp dir
            os.system('rm ' + dirTmp + '/*.' + DICOM_EXT)
            os.system('rm ' + dirTmp + '/grad.b')

            os.system('mrconvert ' + dtiVolFname + ' ' + dtiVolFnameNii + ' -force -quiet') #second MIF to NIFTI
        else:
            os.system('cp ' + dtiVolFnameNiiAux + ' ' + dtiVolFnameNii)

    # compute FA map, and vec for visualization with MRTRIX
    if not os.path.exists(AdMapFname):
        os.system('dwi2tensor ' + dtiVolFname + ' ' + dtiTensorVolFname + ' -quiet -force')
        os.system('tensor2metric ' + dtiTensorVolFname + ' -fa ' + FaMapFname + ' -vec ' + vecMapFname + ' -quiet  -force')
        os.system('tensor2metric ' + dtiTensorVolFname + ' -rd ' + RdMapFname + ' -quiet  -force')
        os.system('tensor2metric ' + dtiTensorVolFname + ' -adc ' + AdcMapFname + ' -quiet  -force')
        os.system('tensor2metric ' + dtiTensorVolFname + ' -ad ' + AdMapFname  + ' -quiet  -force')
        # os.system('mrconvert ' + FaMapFname + ' ' + FaMapFnameNii + ' -force -quiet  -force')  # second MIF to NIFTI
        # os.system('mrconvert ' + RdMapFname + ' ' + RdMapFnameNii + ' -force -quiet  -force')  # second MIF to NIFTI
        # os.system('mrconvert ' + AdcMapFname + ' ' + AdcMapFnameNii + ' -force -quiet  -force')  # second MIF to NIFTI
        # os.system('mrconvert ' + AdMapFname + ' ' + AdMapFnameNii + ' -force -quiet  -force')  # second MIF to NIFTI
        os.system('mrconvert ' + vecMapFname + ' ' + vecMapFnameNii + ' -force -quiet  -force')  # second MIF to NIFTI

    # compute FA and MD map, and vec for visualization with FSL, RD could also be computed but we need to use the eigenvalues
    if not os.path.exists(FaMapFnameFSL):
        print 'Computing FA and MD with fsl...'
        renameDICOMFiles(dirInDTI, dirTmp)
        os.system('mrinfo ' + dirTmp + ' -export_grad_fsl ' + dirTmp + 'bvec.bvec' + ' ' + dirTmp + 'bval.bval' + ' -quiet')
        if not os.path.exists(maskFname):
            renameDICOMFiles(dirInDTI, dirTmp)
            if not os.path.exists(dtiVolFname):
                convertToMif(dirTmp, dtiVolFname)  # first DTI to MIF
            os.system('dwi2mask ' + dtiVolFname + ' ' + maskFname + ' -force')  # get brain mask DWI

        os.system('mrconvert ' + maskFname + ' ' + maskFnameNii + ' -force -quiet')  # then to NIFTI
        os.system('fsl5.0-dtifit -k ' + dtiVolFnameNii + ' -o ' + basenameFSL + ' -m ' + maskFnameNii + ' -r ' + dirTmp + 'bvec.bvec -b ' +  dirTmp + 'bval.bval')
        # clean up tmp dir
        os.system('rm ' + dirTmp + '/*.' + DICOM_EXT)
        os.system('rm ' + dirTmp + '/grad.b')
        os.system('mrconvert ' + FaMapFnameFSLNii + ' ' + FaMapFnameFSL + ' -quiet')
        os.system('mrconvert ' + AdcMapFnameFSLNii + ' ' + AdcMapFnameFSL + ' -quiet')
        print '     Completed!.'


    if not os.path.exists(regMatMri2b0f):
        os.system('fsl5.0-fslroi ' + dtiVolFnameNii + ' ' + dtiVolB0FnameNii + ' 0 1') # get b0 only
        os.system('fsl5.0-bet ' + dtiVolB0FnameNii + ' ' + dtiVolB0nskFnameNii + ' -R -m ' + pIn.fslSkullRm) # remove skull of DTI
        os.system('fsl5.0-flirt -in ' + T1VolnskFnameNii + ' -ref ' + dtiVolB0nskFnameNii + ' -o ' + T1VolregDTIb0FnameNii + ' -omat ' + regMatMri2b0 + ' -usesqform -cost normmi -searchcost normmi ') # register T1 to DTI. Current one seems best one.
                  #' -omat ' + regMatMri2b0 + ' -usesqform -cost mutualinfo -searchcost mutualinfo -schedule /usr/share/fsl/5.0/etc/flirtsch/sch2D_6dof') #register T1 to DTI. Current one seems best one.
        # os.system('/home/lgianca/bin/cnv2dec.sh '+regMatMri2b0+' > ' + regMatMri2b0f)
        os.system('/home/openanogales/bin/cnv2dec.sh ' + regMatMri2b0 + ' > ' + regMatMri2b0f) # fix fsl registration matrix

    if not os.path.exists(T1VolregFnirtWpDTIb0FnameNii):
        print '     Starting fnirt...'
        t = time.time()
        os.system('fsl5.0-fnirt --in=' + T1VolFnameNii + ' --inmask=' + T1VolnskFnameMaskNii + ' --ref=' + dtiVolB0FnameNii + ' --refmask=' + dtiVolB0MasknskFnameNii + ' --aff=' + regMatMri2b0f +
            ' --cout=' +  T1VolregFnirtWpDTIb0FnameNii)
        # os.system('fsl5.0-fnirt --in=' + T1VolFnameNii + ' --inmask=' + T1VolnskFnameMaskNii + ' --ref=' + dtiVolB0FnameNii + ' --refmask=' + dtiVolB0MasknskFnameNii + ' --aff=' + regMatMri2b0f +
        #     ' --cout=' +  T1VolregFnirtWpDTIb0FnameNii + ' --intmod=global_non_linear_with_bias --regmod=membrane_energy') #sometimes is better
        os.system('fsl5.0-applywarp --ref=' + dtiVolB0FnameNii + ' --in=' + T1VolFnameNii + ' --warp=' + T1VolregFnirtWpDTIb0FnameNii +
            ' --out=' + T1VolregFnirtDTIb0FnameNii)
        elapsed = time.time() - t
        print '     Time: ' + str(elapsed / 60) + ' min.'
    print '     Completed.'


    # ========Apply transformation from the PD25 ATLAS to the DTI vol.
    print 'Start ATLAS to DTI transformation...'
    if not os.path.exists(T1ATLASVol2fnirtregDTI):
        os.system('fsl5.0-applywarp --ref=' + dtiVolB0nskFnameNii + ' --in=' + T1VolATLASnsk + ' --warp=' + T1ATLASVolregFnirtWpT1VolnskFnameNii
                + ' --postmat=' + regMatMri2b0f + ' --out=' + T1ATLASVolregDTI) #apply transformation to the T1 ATLAS.
        os.system('fsl5.0-applywarp --ref=' + T1VolFnameNii + ' --in=' + T1VolATLASnsk + ' --warp=' + T1ATLASVolregFnirtWpT1VolnskFnameNii +
            ' --out=' + T1ATLASVolnskregFnirtT1VolnskFnameNii)
        os.system('fsl5.0-applywarp --ref=' + dtiVolB0nskFnameNii + ' --in=' + T1ATLASVolnskregFnirtT1VolnskFnameNii + ' --warp=' + T1VolregFnirtWpDTIb0FnameNii
                + ' --out=' + T1ATLASVol2fnirtregDTI) #apply second fnirt transformation to the T1 ATLAS.

    if not os.path.exists(outNodesMask2fnirtVolregDTI):
        os.system('fsl5.0-applywarp --ref=' + dtiVolB0nskFnameNii + ' --in=' + segVolATLAS + ' --warp=' + T1ATLASVolregFnirtWpT1VolnskFnameNii
            + ' --postmat=' + regMatMri2b0f + ' --interp=nn --out=' + outNodesVolregDTI) #apply transformation to the segmentation.
        os.system('fsl5.0-fslmaths ' + outNodesVolregDTI + ' -div ' + outNodesVolregDTI + ' ' + outNodesMaskVolregDTI)  # get nodes mask

        os.system('fsl5.0-applywarp --ref=' + T1VolFnameNii + ' --in=' + segVolATLAS + ' --warp=' + T1ATLASVolregFnirtWpT1VolnskFnameNii +
            ' --interp=nn --out=' + outNodesVolregFnirtT1VolnskFnameNii)
        os.system('fsl5.0-applywarp --ref=' + dtiVolB0nskFnameNii + ' --in=' + outNodesVolregFnirtT1VolnskFnameNii + ' --warp=' + T1VolregFnirtWpDTIb0FnameNii
            + ' --interp=nn --out=' + outNodes2fnirtVolregDTI) #apply second fnirt transformation to the T1 segmentation
        os.system('fsl5.0-fslmaths ' + outNodes2fnirtVolregDTI + ' -div ' + outNodes2fnirtVolregDTI + ' ' + outNodesMask2fnirtVolregDTI)  # get nodes mask

        print '     Completed.'

    if connectome is 'global':
        # map streamlines to the parcellated image to produce a connectome
        os.system('tck2connectome ' + dtiSiftTrkFname + ' ' + outNodesVolregDTI + ' ' + outConnectome + ' -quiet')
        # ============================

        # convert tracts to trackvis format
        cnvMR2trackvis(dtiSiftTrkFname, dtiVolB0FnameNii, dtiSiftTrkvisFname)
        pass



    # ============================ Local Tractography
    # if False:
    if not os.path.exists(MeanAdMapFname):
        print 'Start Local tractography...'
        if not os.path.exists(dtiFodFname): #check if we computed before the needed files.
            renameDICOMFiles(dirInDTI, dirTmp)
            if not os.path.exists(dtiVolFname):
                convertToMif(dirTmp, dtiVolFname)  # first DTI to MIF
            # clean up tmp dir
            os.system('rm ' + dirTmp + '/*.' + DICOM_EXT)
            os.system('rm ' + dirTmp + '/grad.b')

            os.system('dwi2mask ' + dtiVolFname + ' ' + maskFname + ' -force')  # get brain mask DWI
            if not os.path.exists(dtiRespFname):  # -force would give me problems....
                os.system('/home/openanogales/mrtrix3/scripts/dwi2response tournier ' + dtiVolFname + ' ' + dtiRespFname)  # set tournier as the alg to compute the response function
            os.system('dwi2fod csd ' + dtiVolFname + ' ' + dtiRespFname + ' ' + dtiFodFname + ' -mask ' + maskFname + ' -force')

        os.system('tckgen ' + dtiFodFname + ' ' + dtiTrkLocalFname + ' -seed_image ' + outNodesMask2fnirtVolregDTI + ' -mask ' + maskFname + ' ' + pIn.trkgenCmd + ' -force')  # get the streamlines tractography with the default algorithm (iFOD which is probabilistic).
        os.system('tcksift ' + dtiTrkLocalFname + ' ' + dtiFodFname + ' ' + dtiSiftTrkLocalFname + ' ' + pIn.tcksiftCmd + ' -force')

        # map streamlines to the parcellated image to produce a connectome
        os.system('tck2connectome ' + dtiSiftTrkLocalFname + ' ' + outNodes2fnirtVolregDTI + ' ' + outConnectomeLocal + ' -quiet -force')

        # map parametric maps to streamlines
        # os.system('tckresample ' + dtiSiftTrkLocalFname + ' -line 20 1,2,3 4,5,10 ' + TrkLocal_samples)
        # Example
        # http://mrtrix.readthedocs.io/en/latest/troubleshooting/FAQ.html?highlight=tcksample
        os.system('tcksample ' + dtiSiftTrkLocalFname + ' ' + FaMapFname + ' ' + MeanFaMapFname + ' -stat_tck mean -quiet')
        os.system('tcksample ' + dtiSiftTrkLocalFname + ' ' + RdMapFname + ' ' + MeanRdMapFname + ' -stat_tck mean -quiet')
        os.system('tcksample ' + dtiSiftTrkLocalFname + ' ' + AdcMapFname + ' ' + MeanAdcMapFname + ' -stat_tck mean -quiet')
        os.system('tcksample ' + dtiSiftTrkLocalFname + ' ' + AdMapFname + ' ' + MeanAdMapFname + ' -stat_tck mean -quiet')

        os.system('tck2connectome ' + dtiSiftTrkLocalFname + ' ' + outNodes2fnirtVolregDTI + ' ' + outFaConnectomeLocal + ' -scale_file ' + MeanFaMapFname + ' -stat_edge mean -quiet')
        os.system('tck2connectome ' + dtiSiftTrkLocalFname + ' ' + outNodes2fnirtVolregDTI + ' ' + outRdConnectomeLocal + ' -scale_file ' + MeanRdMapFname + ' -stat_edge mean -quiet')
        os.system('tck2connectome ' + dtiSiftTrkLocalFname + ' ' + outNodes2fnirtVolregDTI + ' ' + outAdcConnectomeLocal + ' -scale_file ' + MeanAdcMapFname + ' -stat_edge mean -quiet')
        os.system('tck2connectome ' + dtiSiftTrkLocalFname + ' ' + outNodes2fnirtVolregDTI + ' ' + outAdConnectomeLocal + ' -scale_file ' + MeanAdMapFname + ' -stat_edge mean -quiet')
        # ============================

        # convert tracts to trackvis format
        # cnvMR2trackvis(dtiSiftTrkLocalFname, dtiVolB0FnameNii, dtiSiftTrkvisLocalFname)

    # compute the connectome for the quantitative maps computed with MRTRIX.
    # if not False:
    if not os.path.exists(dtiSiftTrkLocalFnameTrim):
        print 'triming tracktography...'
        os.system('tckedit -maxlength 70 ' + dtiSiftTrkLocalFname + ' ' + dtiSiftTrkLocalFnameTrim_aux + ' -force')
        os.system('tckedit -number 100k ' + dtiSiftTrkLocalFnameTrim_aux + ' ' + dtiSiftTrkLocalFnameTrim + ' -force')
        os.system('rm ' + dtiSiftTrkLocalFnameTrim_aux)
        os.system('tck2connectome ' + dtiSiftTrkLocalFnameTrim + ' ' + outNodes2fnirtVolregDTI + ' ' + outConnectomeLocalTrim + ' -quiet -force')
        os.system('tcksample ' + dtiSiftTrkLocalFnameTrim + ' ' + FaMapFname + ' ' + MeanFaMapFnameTrim + ' -stat_tck mean -precise -use_tdi_fraction -quiet -force')
        os.system('tcksample ' + dtiSiftTrkLocalFnameTrim + ' ' + RdMapFname + ' ' + MeanRdMapFnameTrim + ' -stat_tck mean -precise -use_tdi_fraction -quiet -force')
        os.system('tcksample ' + dtiSiftTrkLocalFnameTrim + ' ' + AdcMapFname + ' ' + MeanAdcMapFnameTrim + ' -stat_tck mean -precise -use_tdi_fraction -quiet -force')
        os.system('tcksample ' + dtiSiftTrkLocalFnameTrim + ' ' + AdMapFname + ' ' + MeanAdMapFnameTrim + ' -stat_tck mean -precise -use_tdi_fraction -quiet -force')

        os.system('tck2connectome ' + dtiSiftTrkLocalFnameTrim + ' ' + outNodes2fnirtVolregDTI + ' ' + outFaConnectomeLocalTrim + ' -scale_file ' + MeanFaMapFnameTrim + ' -stat_edge mean -quiet -force')
        os.system('tck2connectome ' + dtiSiftTrkLocalFnameTrim + ' ' + outNodes2fnirtVolregDTI + ' ' + outRdConnectomeLocalTrim + ' -scale_file ' + MeanRdMapFnameTrim + ' -stat_edge mean -quiet -force')
        os.system('tck2connectome ' + dtiSiftTrkLocalFnameTrim + ' ' + outNodes2fnirtVolregDTI + ' ' + outAdcConnectomeLocalTrim + ' -scale_file ' + MeanAdcMapFnameTrim + ' -stat_edge mean -quiet -force')
        os.system('tck2connectome ' + dtiSiftTrkLocalFnameTrim + ' ' + outNodes2fnirtVolregDTI + ' ' + outAdConnectomeLocalTrim + ' -scale_file ' + MeanAdMapFnameTrim + ' -stat_edge mean -quiet -force')

    # compute the connectome for the FA and MD from the FSL - FA and MD maps
    if not os.path.exists(outFaConnectomeLocalFSLTrim):
        print 'Connectomes for FSL...'
        # os.system('tcksample ' + dtiSiftTrkLocalFname + ' ' + FaMapFnameFSL + ' ' + MeanFaMapFnameFSL + ' -stat_tck mean -precise -use_tdi_fraction -quiet -force')
        # os.system('tcksample ' + dtiSiftTrkLocalFname + ' ' + AdcMapFnameFSL + ' ' + MeanAdcMapFnameFSL + ' -stat_tck mean -precise -use_tdi_fraction -quiet -force')
        #
        # os.system('tck2connectome ' + dtiSiftTrkLocalFname + ' ' + outNodes2fnirtVolregDTI + ' ' + outFaConnectomeLocalFSL + ' -scale_file ' + MeanFaMapFnameFSL + ' -stat_edge mean -quiet')
        # os.system('tck2connectome ' + dtiSiftTrkLocalFname + ' ' + outNodes2fnirtVolregDTI + ' ' + outAdcConnectomeLocalFSL + ' -scale_file ' + MeanAdcMapFnameFSL + ' -stat_edge mean -quiet')

        os.system('tcksample ' + dtiSiftTrkLocalFnameTrim + ' ' + FaMapFnameFSL + ' ' + MeanFaMapFnameFSLTrim + ' -stat_tck mean -precise -use_tdi_fraction -quiet -force')
        os.system('tcksample ' + dtiSiftTrkLocalFnameTrim + ' ' + AdcMapFnameFSL + ' ' + MeanAdcMapFnameFSLTrim + ' -stat_tck mean -precise -use_tdi_fraction -quiet -force')

        os.system('tck2connectome ' + dtiSiftTrkLocalFnameTrim + ' ' + outNodes2fnirtVolregDTI + ' ' + outFaConnectomeLocalFSLTrim + ' -scale_file ' + MeanFaMapFnameFSLTrim + ' -stat_edge mean -quiet')
        os.system('tck2connectome ' + dtiSiftTrkLocalFnameTrim + ' ' + outNodes2fnirtVolregDTI + ' ' + outAdcConnectomeLocalFSLTrim + ' -scale_file ' + MeanAdcMapFnameFSLTrim + ' -stat_edge mean -quiet')


    if not os.path.exists(outAdcConnectomeLocalFSLShiftTrim):
        os.system('tcksift ' + dtiSiftTrkLocalFname + ' ' + dtiFodFname + ' ' + dtiSiftTrkLocalFnameShiftTrim + ' -term_number 200K -force')
        os.system('tck2connectome ' + dtiSiftTrkLocalFnameShiftTrim + ' ' + outNodes2fnirtVolregDTI + ' ' + outConnectomeLocalShiftTrim + ' -quiet -force')
        os.system('tcksample ' + dtiSiftTrkLocalFnameShiftTrim + ' ' + FaMapFnameFSL + ' ' + MeanFaMapFnameFSLShiftTrim + ' -stat_tck mean -precise -use_tdi_fraction -quiet -force')
        os.system('tcksample ' + dtiSiftTrkLocalFnameShiftTrim + ' ' + AdcMapFnameFSL + ' ' + MeanAdcMapFnameFSLShiftTrim + ' -stat_tck mean -precise -use_tdi_fraction -quiet -force')

        os.system('tck2connectome ' + dtiSiftTrkLocalFnameShiftTrim + ' ' + outNodes2fnirtVolregDTI + ' ' + outFaConnectomeLocalFSLShiftTrim + ' -scale_file ' + MeanFaMapFnameFSLShiftTrim + ' -stat_edge mean -quiet')
        os.system('tck2connectome ' + dtiSiftTrkLocalFnameShiftTrim + ' ' + outNodes2fnirtVolregDTI + ' ' + outAdcConnectomeLocalFSLShiftTrim + ' -scale_file ' + MeanAdcMapFnameFSLShiftTrim + ' -stat_edge mean -quiet')

    if not os.path.exists(outAdcConnectomeLocalFSLShiftTrim2):
        os.system('tckgen ' + dtiFodFname + ' ' + dtiTrkLocalFname + ' -seed_image ' + outNodesMask2fnirtVolregDTI + ' -mask ' + maskFname + ' ' + pIn.trkgenCmd + ' -force')

        os.system('tcksift ' + dtiTrkLocalFname + ' ' + dtiFodFname + ' ' + dtiSiftTrkLocalFnameShiftTrim2 + ' -term_number 200K -force')
        os.system('tck2connectome ' + dtiSiftTrkLocalFnameShiftTrim2 + ' ' + outNodes2fnirtVolregDTI + ' ' + outConnectomeLocalShiftTrim2 + ' -quiet -force')
        os.system('tcksample ' + dtiSiftTrkLocalFnameShiftTrim2 + ' ' + FaMapFnameFSL + ' ' + MeanFaMapFnameFSLShiftTrim2 + ' -stat_tck mean -precise -use_tdi_fraction -quiet -force')
        os.system('tcksample ' + dtiSiftTrkLocalFnameShiftTrim2 + ' ' + AdcMapFnameFSL + ' ' + MeanAdcMapFnameFSLShiftTrim2 + ' -stat_tck mean -precise -use_tdi_fraction -quiet -force')

        os.system('tck2connectome ' + dtiSiftTrkLocalFnameShiftTrim2 + ' ' + outNodes2fnirtVolregDTI + ' ' + outFaConnectomeLocalFSLShiftTrim2 + ' -scale_file ' + MeanFaMapFnameFSLShiftTrim2 + ' -stat_edge mean -quiet')
        os.system('tck2connectome ' + dtiSiftTrkLocalFnameShiftTrim2 + ' ' + outNodes2fnirtVolregDTI + ' ' + outAdcConnectomeLocalFSLShiftTrim2 + ' -scale_file ' + MeanAdcMapFnameFSLShiftTrim2 + ' -stat_edge mean -quiet')


    # compute connectome normalized by the volume of both nodes. Modify outConnectomeLocalInv and dtiSiftTrkLocalFname according to the connectome we want
    if not os.path.exists(outNormVolConnectomeLocalShiftTrim2):
        # get vol connectome
        if not os.path.exists(outNodes2fnirtVolConnectome):
            os.system('fsl5.0-fslstats ' + outNodes2fnirtVolregDTI + ' -H 17 0 16 >> ' + outNodes2fnirtHistregDTI) #save histogram of nodes - from nodes 0 to 16
            histArr = np.loadtxt(outNodes2fnirtHistregDTI)
            histArr = histArr[1:17] #remove background (bin 0)
            connVolArr=np.diag(histArr)
            for x in range(0, 16): #make connectome
                for y in range(x, 16):
                    if x is not y:
                        connVolArr[x, y] = histArr[x] + histArr[y]

            np.savetxt(outNodes2fnirtVolConnectome, connVolArr, fmt='%.0f')
        else:
            connVolArr = np.loadtxt(outNodes2fnirtVolConnectome)

        #normalize streamlines connectome
        connVolArr[connVolArr == 0] = 1 #to avoid nans
        connArr = np.loadtxt(outConnectomeLocalShiftTrim2)
        connNormArr = np.divide(connArr,connVolArr)
        np.savetxt(outNormVolConnectomeLocalShiftTrim2, connNormArr)

        # the following command line should've done all the above... however it does not. That's why I implement the rest of the code
        # os.system('tck2connectome ' + dtiSiftTrkLocalFname + ' ' + outNodes2fnirtVolregDTI + ' ' + outConnectomeLocalInv + ' -scale_invnodevol -quiet -force')





    pass



def multipleConnectomes( dirOut, pIn, nCon = 5):
    """
    Compute n number of tracktographies and connectomes.

    :param dirInT1: PPMI directory with image T1 info
    :param dirInT2: PPMI directory with image T2 info
    :param dirInDTI: PPMI directory with image DTI info
    :param dirOut: output directory
    :param check: set to True to check if connectome exist
    :param connectome: decide which connectome to use, global or local - only of the PD area.
    :return: None (if check is true returns True, False depending if the connectome has been computed)
    """

    # ==== Set file names
    maskFname = dirOut + pIn.pName + '_outDti_mask.mif'  # brain mask
    dtiFodFname = dirOut + pIn.pName + '_outDti_fod.mif'  # fibre orientation distribution estimation
    dtiTrkLocalFname = dirOut + pIn.pName + '_outDtiLocal_trk.tck'

    FaMapFnameFSL = dirOut + pIn.pName + '_outFSL_FA.mif'
    AdcMapFnameFSL = dirOut + pIn.pName + '_outFSL_MD.mif'

    outNodes2fnirtVolregDTI = dirOut + pIn.pName + '_outNodes-2fnirt-reg-DTI.nii.gz'
    outNodesMask2fnirtVolregDTI = dirOut + pIn.pName + '_outNodesMask-2fnirt-reg-DTI.nii.gz'

    MeanFaMapFnameFSLShiftTrim = dirOut + '2_outmeanFSL_FA_ShiftTrim_aux.csv'  # mean track FA map trim FSL
    MeanAdcMapFnameFSLShiftTrim = dirOut + '2_outmeanFSL_MD_ShiftTrim_aux.csv'  # mean track FA map trim FSL


    for x in range(4,nCon+4):
        print 'Starting connectome ', x-4
        dtiSiftTrkLocalFnameShiftTrim = dirOut + '2_outDti_sift_trk_local_ShiftTrim_' + str(x) + '.tck'
        dtiSiftTrkLocalFnameShiftTrimAux = dirOut + '2_outDti_sift_trk_local_ShiftTrimAux_' + str(x) + '.tck'
        outConnectomeLocalShiftTrim = dirOut + '2_connectomeLocal_ShiftTrim_' + str(x) + '.csv'
        outFaConnectomeLocalFSLShiftTrim = dirOut + '2_FaconnectomeLocalFSL_ShiftTrim_' + str(x) + '.csv'
        outAdcConnectomeLocalFSLShiftTrim = dirOut + '2_AdcconnectomeLocalFSL_ShiftTrim_' + str(x) + '.csv'

        os.system('tckgen ' + dtiFodFname + ' ' + dtiTrkLocalFname + ' -seed_image ' + outNodesMask2fnirtVolregDTI + ' -mask ' + maskFname + ' ' + pIn.trkgenCmd + ' -force')  # get the streamlines tractography with the default algorithm (iFOD which is probabilistic).
        # 1) or 2)
        # 1 -
        # os.system('tcksift ' + dtiTrkLocalFname + ' ' + dtiFodFname + ' ' + dtiSiftTrkLocalFnameShiftTrim + ' -term_number 200K -force')
        # 2 -
        os.system('tcksift ' + dtiTrkLocalFname + ' ' + dtiFodFname + ' ' + dtiSiftTrkLocalFnameShiftTrimAux + ' -term_number 350K -force')
        os.system('tcksift ' + dtiSiftTrkLocalFnameShiftTrimAux + ' ' + dtiFodFname + ' ' + dtiSiftTrkLocalFnameShiftTrim + ' -term_number 200K -force')

        os.system('tck2connectome ' + dtiSiftTrkLocalFnameShiftTrim + ' ' + outNodes2fnirtVolregDTI + ' ' + outConnectomeLocalShiftTrim + ' -quiet -force')
        os.system('tcksample ' + dtiSiftTrkLocalFnameShiftTrim + ' ' + FaMapFnameFSL + ' ' + MeanFaMapFnameFSLShiftTrim + ' -stat_tck mean -precise -use_tdi_fraction -quiet -force')
        os.system('tcksample ' + dtiSiftTrkLocalFnameShiftTrim + ' ' + AdcMapFnameFSL + ' ' + MeanAdcMapFnameFSLShiftTrim + ' -stat_tck mean -precise -use_tdi_fraction -quiet -force')

        os.system('tck2connectome ' + dtiSiftTrkLocalFnameShiftTrim + ' ' + outNodes2fnirtVolregDTI + ' ' + outFaConnectomeLocalFSLShiftTrim + ' -scale_file ' + MeanFaMapFnameFSLShiftTrim + ' -stat_edge mean -quiet')
        os.system('tck2connectome ' + dtiSiftTrkLocalFnameShiftTrim + ' ' + outNodes2fnirtVolregDTI + ' ' + outAdcConnectomeLocalFSLShiftTrim + ' -scale_file ' + MeanAdcMapFnameFSLShiftTrim + ' -stat_edge mean -quiet')

        os.system('rm ' + MeanFaMapFnameFSLShiftTrim + ' ' + MeanAdcMapFnameFSLShiftTrim + ' ' )

    pass


def cnvMR2trackvis(fileIn, volIn, fileOut):
    """
    Convert Mrtrix file to trackvis
    :param fileIn: 
    :param volIn: 
    :param fileOut: 
    :return: 
    """
    import nipype.interfaces.mrtrix as mrt
    tck2trk = mrt.MRTrix2TrackVis()
    tck2trk.inputs.in_file = fileIn
    tck2trk.inputs.image_file = volIn
    tck2trk.inputs.out_filename = fileOut
    tck2trk.run()




def getImgsInfoFr(baseFile='full_DTI_only_3_17_2017.csv', groupIn='Prodromal', modalityIn='MRI', modDescrIn='GRAPPA', getLast=True ):
    """
    Get img location according to parameters
    :param baseFile: csv file from PPMI data 
    :param groupIn: subject group (i.e. Prodromal) (set None for all)
    :param modalityIn: general modality (i.e. MRI) 
    :param modDescrIn: string description of the modality (i.e. GRAPPA) check nb/Find Subjects for Pipeline for examples
    :param getLast: if True get the last volume if not the first, if Second get the first acq after 300 days.
    :return: DataFrame with subject info and img location list
    """
    # only DTI
    info2Fr = pd.read_csv( baseFile )

    # filter modality
    info2Fr = info2Fr[info2Fr.Modality == modalityIn]
    # filter cohort
    if groupIn is not None:
        info2Fr = info2Fr[info2Fr.Group == groupIn]

    info2Fr = info2Fr.sort_values(['Subject', 'Visit'])
    # keep grappa only
    info2Fr = info2Fr[info2Fr.Description.str.contains(modDescrIn)]
    # get last visit
    subjArr = info2Fr.groupby('Subject').first().index.values

    # create dataframe for MRI at first/last visit  with data to analyze
    mriLastVisitFr = pd.DataFrame()
    for subjId in subjArr:
        #- select first or last sample
        tmpSid = -1
        if getLast == False:
            tmpSid = 0
        elif getLast == 'Second':
            tmpSid = 2

        # if subjId in (3814,60043):
        #     print 'there'

        #OSCAR:
        if tmpSid is 2: #get second acquistion in time.
            # keysSbInt = info2Fr[info2Fr.Subject == subjId]['Acq Date'].keys()
            dateDt = pd.to_datetime(info2Fr[info2Fr.Subject == subjId]['Acq Date'])
            dateSortDt = dateDt.sort_values()
            keysSbInt = dateSortDt.keys()
            zeroTm = pd.to_datetime(info2Fr[info2Fr.Subject == subjId]['Acq Date'][keysSbInt[0]], format="%m/%d/%Y" )

            for idx, val in enumerate(keysSbInt):
                # compute time difference
                nthTm = pd.to_datetime(info2Fr[info2Fr.Subject == subjId]['Acq Date'][val], format="%m/%d/%Y")
                # print subjId, (nthTm - zeroTm).days
                if (nthTm - zeroTm).days > 270:  # we stop on the first acq that is more than 270 days apart from the first acq. (270 because that's the case of a PD)
                    tmpSr = info2Fr[info2Fr.Subject == subjId].loc[val]
                    break
                else: #if it does not find any acq. after 300 we just return the first one.
                    tmpSr = info2Fr[info2Fr.Subject == subjId].iloc[0]
        else:  #get last acquisition regarding of the order
            tmpSr = info2Fr[info2Fr.Subject == subjId].iloc[tmpSid]


        # Luca:
        # tmpSr = info2Fr[info2Fr.Subject == subjId].iloc[tmpSid]

        #-
        # convert description field to valid dir name
        descrStr = str(tmpSr['Description'])
        for c in ' *()/':
            descrStr = descrStr.replace(c, "_")
        # for c in '*()/': #apparently 3221 has a bug
        #     descrStr = descrStr.replace(c, "_")
        fName = str(tmpSr['Subject']) + '/' + descrStr + '/'
        # generate first dir
        fName = str(tmpSr['Subject']) + '/' + descrStr + '/'
        # find all dirs (each is a different acquisition)
        lstFiles = os.listdir(CONF['DIR_PPMI_IMGS']  + '/' + fName)
        # sort by reverse order, possible because the format is 2016-03-09_10_33_34.0
        if getLast in (True, 'Second'):
            lstFiles.sort(reverse=True)
        else:
            lstFiles.sort(reverse=False)
        # add date
        fName += lstFiles[0]
        # add single directory
        fName += '/' + os.listdir(CONF['DIR_PPMI_IMGS'] + '/' + fName)[0]

        tmpSr['fName'] = fName

        # add to dataframe
        mriLastVisitFr = mriLastVisitFr.append(tmpSr)

    # convert column types
    mriLastVisitFr['Subject'] = mriLastVisitFr['Subject'].astype(int)
    mriLastVisitFr['Visit'] = mriLastVisitFr['Visit'].astype(int)
    mriLastVisitFr['Image Data ID'] = mriLastVisitFr['Image Data ID'].astype(int)

    return mriLastVisitFr


def matchSubj(smallFr, largeFr, numToMatch=2):
    """
    Match subjects according to Age and Sex from the smallFr to largeFr
    :param smallFr: Dataframe generated by getImgsInfoFr
    :param largeFr: Dataframe generated by getImgsInfoFr
    :param numToMatch: num of subjects for each of the smallFr
    :return: pruned largeFr
    """
    # init
    largeFr2 = largeFr.copy()
    largeFr2['sel'] = 0 # selected flag
    outFr = pd.DataFrame()
    for id, row in smallFr.iterrows():
        # find closest ids first from sex and then age (that were not previously selected)
        closestIds = \
            np.abs(row['Age'] - largeFr2[(largeFr2['Sex'] == row['Sex']) & (largeFr2['sel'] == 0)]['Age'])\
                .sort_values().iloc[0:numToMatch].index.values
        assert(len(closestIds) == numToMatch)
        # set as selected
        largeFr2.loc[closestIds, 'sel'] = 1
        # append found (from original Frame, avoiding "sel" column)
        outFr = outFr.append( largeFr.loc[closestIds] )

    return outFr



if __name__ == "__main__":

    # set  parameters
    trkParam = Param()
    # directories
    dirIn = '/data/bigdata/PPMI/ppmi_imgs/full/PPMI/14426/DTI_gated/2016-03-17_13_56_22.0/S405055/'
    dirOut = '/data/bigdata/PPMIcomp/ppmi_14426_last/'
    dirFS = '/data/freesurfer_subj/ppmi_14426_last/mri/'
    genTractography(dirIn, dirOut, dirFS, trkParam)

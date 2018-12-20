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




if __name__ == "__main__":

    # set  parameters
    trkParam = Param()
    # directories
    dirIn = '/data/bigdata/PPMI/ppmi_imgs/full/PPMI/14426/DTI_gated/2016-03-17_13_56_22.0/S405055/'
    dirOut = '/data/bigdata/PPMIcomp/ppmi_14426_last/'
    dirFS = '/data/freesurfer_subj/ppmi_14426_last/mri/'
    genTractography(dirIn, dirOut, dirFS, trkParam)

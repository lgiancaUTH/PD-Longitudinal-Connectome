
# Longitudinal Connectomes Analysis

This repository contains supporting code for the following manuscript
> Óscar Peña-Nogales, Timothy M. Ellmore, Rodrigo de Luis-García, Jessika Suescun, Mya C. Schiess and Luca Giancardo. Longitudinal Connectomes as a Candidate Progression Marker for Prodromal Parkinson's Disease. (under review)


In this work, we aim to identify and quantify a longitudinal degenerative Parkinson's disease pattern from the diffusion magnetic resonance imaging connectivity information that can be found in a de novo early Parkinson's disease cohort (n=21) and in a cohort at high risk of being in the Parkinson's disease prodromal phase (n=16) that was not present in cohort matched Controls (n=30). This progression pattern was numerically quantified with a longitudinal brain connectome progression score. This score is generated by an interpretable machine learning algorithm trained, with cross-validation, on the longitudinal connectivity information of Parkinson's disease and Control groups computed on a nigrostriatal pathway-specific parcellation atlas. 


The second half computes statistical tests and comparison with the clinical metrics. I order to run this analysis you will have to request access to the PPMI dataset and change the following parameters in the configuration.json file

> "DIR_PPMI_IMGS": "raw images location", (optional)

> "DIR_PPMI_CLIN_INFO": "clinical PPMI data location"

## No installation
You can visualize the analysis directly by clicking on [longitudinal-connectomes-analysis.ipynb](https://github.com/lgiancaUTH/PD-Longitudinal-Connectome/blob/master/longitudinal-connectomes-analysis.ipynb)

## Installation
This code has been tested with Python 2.7 and conda/pip 

```
git clone https://github.com/lgiancaUTH/PD-Longitudinal-Connectome.git
cd PD-Longitudinal-Connectome
```
Optionally, create a (conda) virtual environment and activate it
```
conda create -n valve-env python=2.7 

conda activate pd-long-env
```

and finally, install the required dependencies
```
pip install -r requirements.txt
```

if rpy2 fails to compile, you can install the precompiled packages with conda
```
conda install rpy2==2.8.5
```


## Run
Run the jupyter notebook server  
```
jupyter notebook
```
and open [longitudinal-connectomes-analysis.ipynb](https://github.com/lgiancaUTH/PD-Longitudinal-Connectome/blob/master/longitudinal-connectomes-analysis.ipynb)

In order to run the full analysis, you will have to request the PPMI datasets at www.ppmi-info.org/data and change the [data/configuration.json] file accordingly.


# Parkinson’s Progression Markers Initiative (PPMI)

Data used in the preparation of this article were obtained from the Parkinson’s Progression Markers Initiative (PPMI) database (www.ppmi-info.org/data); for up to date information on the study, visit www.ppmi-info.org.

The Parkinson’s Progression Markers Initiative—a public-private partnership—is funded by the Michael J. Fox Foundation for Parkinson’s Research and funding partners including AbbVie, Avid, Biogen, BioLegend Bristol-Myers Squibb, GE Healthcare, Genentech, GlaxoSmithKline, Lilly, Lundbeck, Merck, Meso Scale Discovery, Pfizer, Piramal, Roche, Sanofi Genzyme, Servier, Takeda, Teva, UCB, and Golub Capital
~                                                                                                                                        



## Other information
This code is free to use for any non-commercial purposes, provided that the original publication is cited. 

We refer to the original publication for additional infomation and acknowledgements.


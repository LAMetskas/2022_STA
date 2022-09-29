# 2022_STA

This script repository is designed to assist tomography data processing. The intended audience are users who are comfortable running existing software packages but are looking to either automate the process or improve outputs. The authors accept no responsibility for incorrect results stemming from over-fitting or other processing errors potentially enabled by these methods, or for silent bugs and errors caused by use of versions other than those for which the wrappers and scripts were written. All code is supplied as-is under share-and-share-alike licenses as specified in the individual scripts.

These scripts are designed for specific software packages and versions as follows:

SerialEM_CollectionScript.txt is a script for data collection on a Titan Krios (K3 camera) running SerialEM version 3.8.
IMOD_FrameAlignment.sh is a wrapper for the alignframes function in IMOD version 4.9.8
ctffind4_CorrectOutput.m is a Matlab R2017b script to detect and correct failures in outputs from ctffind4 v. 4.1.13.
novaCTF_runWBPandSIRT.sh is a wrapper for novaCTF, October 2017 updated release.
Dynamo_BulkInputOfIMODModels.m is a Matlab R2017b script for Dynamo v 1.1.319, and inputs model2point outputs from IMOD 4.9.8.

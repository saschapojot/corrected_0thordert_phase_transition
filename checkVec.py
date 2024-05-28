import xml.etree.ElementTree as ET
import numpy as np
# from datetime import datetime
import statsmodels.api as sm
from scipy import stats
import glob
import sys
import re
from copy import deepcopy
import warnings
from scipy.stats import ks_2samp
import pickle


#This script checks if a vector reaches equilibrium
sigWrongArg="wrong number of arguments"
sigEq="equilibrium"
sigContinue="continue"
sigStop="stop"

if (len(sys.argv)!=2):
    print(sigWrongArg)
    exit()


xmlFilesPath=str(sys.argv[1])
#fetch files in the directory
inXMLFileNames=[]
startVals=[]

for file in glob.glob(xmlFilesPath+"/*.pkl"):
    inXMLFileNames.append(file)
    matchStart=re.search(r"loopStart(-?\d+(\.\d+)?)loopEnd",file)
    if matchStart:
        startVals.append(matchStart.group(1))


def str2int(valList):
    ret = [int(strTmp) for strTmp in valList]
    return ret


startVals = str2int(startVals)

start_inds = np.argsort(startVals)

#sort files by starting value
inXMLFileNames=[inXMLFileNames[ind] for ind in start_inds]

#ensure the file number is a multiple of 3
if len(inXMLFileNames)%3==0:
    xmlFileToBeParsed=deepcopy(inXMLFileNames)
elif len(inXMLFileNames)%3==1:
    xmlFileToBeParsed=deepcopy(inXMLFileNames[1:])
else:
    xmlFileToBeParsed=deepcopy(inXMLFileNames[2:])

lastFileNum=1500 if len(xmlFileToBeParsed)>2250 else int(len(xmlFileToBeParsed)/3*2)

xmlFileToBeParsed=xmlFileToBeParsed[-lastFileNum:]
startingFileInd=len(xmlFileToBeParsed)-lastFileNum
def parse1File(fileName):
    """

    :param fileName: xml file name
    :return: the values in the vector
    """

    with open(fileName,"rb") as fptr:
        vec=list((pickle.load(fptr)))
        return vec


#combine all vectors
vecValsCombined=parse1File(xmlFileToBeParsed[0])

for file in xmlFileToBeParsed[1:]:
    vecValsCombined+=parse1File(file)


vecValsCombined=np.array(vecValsCombined)

#all values equal
meanU=np.mean(vecValsCombined)

diff=np.linalg.norm(vecValsCombined-meanU,ord=1)/len(vecValsCombined)


if diff<1e-7:
    print(sigStop+" same"+", fileNum="+str(lastFileNum))
    exit()


#check if the whole vector has the same value
with warnings.catch_warnings():
    warnings.filterwarnings("error")
    try:
        vecAutc=sm.tsa.acf(vecValsCombined)
    except Warning as w:
        print(sigStop+" same"+", fileNum="+str(lastFileNum))
        exit()


halfLength=int(len(vecValsCombined)/2)
part0=vecValsCombined[:halfLength]
part1=vecValsCombined[halfLength:]

same0=False
same1=False

#check if the part0 has the same value
with warnings.catch_warnings():
    warnings.filterwarnings("error")
    try:
        autc0=sm.tsa.acf(part0)
    except Warning as w:
        same0=True


#check if the part1 has the same value
with warnings.catch_warnings():
    warnings.filterwarnings("error")
    try:
        autc1=sm.tsa.acf(part1)
    except Warning as w:
        same1=True


if same0 and same1:
    print(sigStop+" same"+", fileNum="+str(lastFileNum))
    exit()


elif same0==True and same1==False:
    # print("f0 True f1 False")
    print(sigContinue)
    exit()

elif same0==False and same1==True:
    # print("f0 False f1 True")
    print(sigContinue)
    exit()

def Jackknife(vec):
    """

    :param vec:
    :return: the mean and half length  of 0.95 confidence interval computed by Jackkknife
    """

    psMean=np.mean(vec)

    psVar=np.var(vec,ddof=1)

    n=len(vec)

    hfLen=1.96*np.sqrt(psVar/n)
    return psMean,hfLen


#computation of auto-correlation
NLags=int(np.ceil(len(vecValsCombined)*5/6))
acfOfVec=sm.tsa.acf(vecValsCombined,nlags=NLags)
eps=1e-3

lagVal=0
if np.min(np.abs(acfOfVec))>eps:
    print("high correlation")
    # print(np.min(np.abs(acfOfVec)))

    print(sigContinue)
    exit()

else:
    lagVal=np.where(np.abs(acfOfVec)<=eps)[0][0]
    vecValsSelected=vecValsCombined[::lagVal]
    lengthTmp=len(vecValsSelected)
    if lengthTmp%2==1:
        lengthTmp-=1
    vecValsToCompute=vecValsSelected[:lengthTmp]
    lenPart=int(len(vecValsToCompute)/2)
    selectedFromPart0=vecValsToCompute[:lenPart]
    selectedFromPart1=vecValsToCompute[lenPart:]

    result = ks_2samp(selectedFromPart0, selectedFromPart1)
    if result.pvalue>0.1:
        print(sigEq+" ,lag="+str(lagVal)+", fileNum="+str(lastFileNum)+", startingFileInd="+str(startingFileInd))
        exit()
    else:
        print(sigContinue)
        exit()
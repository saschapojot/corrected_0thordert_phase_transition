import pandas as pd

import numpy as np
import glob
import sys
import re
import matplotlib.pyplot as plt
from datetime import datetime

import pandas as pd

pathData="../../version1Data/1d/funcLJPotPBC/row0"
TVals=[]
TFileNames=[]

for TFile in glob.glob(pathData+"/T*"):
    TFileNames.append(TFile)

    matchT=re.search(r"T(\d+(\.\d+)?)",TFile)

    TVals.append(float(matchT.group(1)))

#sort T files

sortedInds=np.argsort(TVals)
sortedTVals=[TVals[ind] for ind in sortedInds]
sortedTFiles=[TFileNames[ind] for ind in sortedInds]


def readGAA(oneTFile):
    """

    :param oneTFile:
    :return: fA, plot of fA
    """
    matchT=re.search(r"T(\d+(\.\d+)?)",oneTFile)
    TVal=float(matchT.group(1))
    GAAFile=oneTFile+"/GAA.csv"
    GAAData=pd.read_csv(GAAFile,header=None)
    GAAMat=GAAData.to_numpy()
    N,_=GAAMat.shape
    fA=[]
    for i in range(0,N):
        supDiag=np.diagonal(GAAMat,offset=i)
        fA.append(np.sum(supDiag)/supDiag.size)

    outfAFile=oneTFile+"/fA.csv"
    fA=np.array(fA)

    np.savetxt(outfAFile, fA, delimiter=',')

def readGBB(oneTFile):
    """

    :param oneTFile:
    :return: fB, plot of fB
    """
    matchT=re.search(r"T(\d+(\.\d+)?)",oneTFile)
    TVal=float(matchT.group(1))
    GBBFile=oneTFile+"/GBB.csv"
    GBBData=pd.read_csv(GBBFile,header=None)
    GBBMat=GBBData.to_numpy()
    N,_=GBBMat.shape
    fB=[]
    for i in range(0,N):
        supDiag=np.diagonal(GBBMat,offset=i)
        fB.append(np.sum(supDiag)/supDiag.size)

    outfBFile=oneTFile+"/fB.csv"
    fB=np.array(fB)

    np.savetxt(outfBFile, fB, delimiter=',')

def readGAB(oneTFile):
    """
    :param oneTFile:
    :return: fAB,
    """
    matchT=re.search(r"T(\d+(\.\d+)?)",oneTFile)
    TVal=float(matchT.group(1))
    GABFile=oneTFile+"/GAB.csv"
    GABData=pd.read_csv(GABFile,header=None)
    GABMat=GABData.to_numpy()
    N,_=GABMat.shape
    fAB=[]
    for i in range(-N+1,N):
        offDiag=np.diagonal(GABMat,offset=i)
        fAB.append(np.sum(offDiag)/len(offDiag))
    outfABFile=oneTFile+"/fAB.csv"
    fAB=np.array(fAB)

    np.savetxt(outfABFile, fAB, delimiter=',')


for oneTFile in sortedTFiles:
    readGAA(oneTFile)
    readGBB(oneTFile)
    readGAB(oneTFile)

#compare GAA
fAArray=[]
for oneTFile in sortedTFiles:
    pathfATmp=oneTFile+"/fA.csv"
    dfTmp=pd.read_csv(pathfATmp,header=None)
    arrTmp=dfTmp.to_numpy()

    fAArray.append(arrTmp.T[0,:])
#compare GBB
fBArray=[]
for oneTFile in sortedTFiles:
    pathfBTmp=oneTFile+"/fB.csv"
    dfTmp=pd.read_csv(pathfBTmp,header=None)
    arrTmp=dfTmp.to_numpy()

    fBArray.append(arrTmp.T[0,:])
#compare GAB
fABArray=[]
for oneTFile in sortedTFiles:
    pathfABTmp=oneTFile+"/fAB.csv"
    dfTmp=pd.read_csv(pathfABTmp,header=None)
    arrTmp=dfTmp.to_numpy()
    fABArray.append(arrTmp.T[0,:])

#each row is GAA for one temperature, each column is the GAA with the same distance under different tempeatures
fAArray=np.array(fAArray)
#each row is GAB for one temperature, each column is the GAA with the same distance under different tempeatures

fBArray=np.array(fBArray)
#each row is GAB for one temperature, each column is the GAA with the same distance under different tempeatures
#for fABArray, indices 0 to N-2 are negative distances,
#index N-1 is 0-distance
#indices N to 2N-2 are positive distances
fABArray=np.array(fABArray)

N=fAArray[0,:].size



#variance of A
# j=0
plt.figure()
varATmp=fAArray[:,0]
plt.plot(sortedTVals,varATmp,color="blue")
plt.scatter(sortedTVals,varATmp,color="red")
plt.xlabel("T")
plt.ylabel("variance of $x_{A}$")


plt.savefig(pathData+"/varA.pdf")
plt.close()

#variance of B
plt.figure()

varBTmp=fBArray[:,0]
plt.plot(sortedTVals,varBTmp,color="blue")
plt.scatter(sortedTVals,varBTmp,color="red")
plt.xlabel("T")
plt.ylabel("variance of $x_{B}$")


plt.savefig(pathData+"/varB.pdf")
plt.close()

def plotGAA(j):
    """

    :param j: distance
    :return: correlation function of xA
    """
    colATmp=fAArray[:,j]
    plt.plot(sortedTVals,colATmp,color="blue")
    plt.scatter(sortedTVals,colATmp,color="red")
    plt.xlabel("T")
    plt.ylabel("$G^{AA}($"+str(j)+"$)$")


    plt.savefig(pathData+"/GAA"+str(j)+".pdf")

    plt.close()

def plotrAA(j):
    """

    :param j: distance
    :return: correlation coefficient of xA
    """
    colATmp=fAArray[:,j]
    colA0=fAArray[:,0]
    rhoTmp=colATmp/colA0
    plt.plot(sortedTVals,rhoTmp,color="green")
    plt.scatter(sortedTVals,rhoTmp,color="magenta")
    plt.xlabel("T")
    plt.ylabel("$\\rho^{AA}($"+str(j)+"$)$")


    plt.savefig(pathData+"/rhoAA"+str(j)+".pdf")

    plt.close()


def plotGBB(j):
    """

    :param j: distance
    :return: correlation function of xB
    """
    colBTmp=fBArray[:,j]
    plt.plot(sortedTVals,colBTmp,color="blue")
    plt.scatter(sortedTVals,colBTmp,color="red")
    plt.xlabel("T")
    plt.ylabel("$G^{BB}($"+str(j)+"$)$")


    plt.savefig(pathData+"/GBB"+str(j)+".pdf")

    plt.close()

def plotrBB(j):
    """

    :param j:
    :return: correlation coefficient of xB
    """
    colBTmp=fBArray[:,j]
    colB0=fBArray[:,0]
    rhoTmp=colBTmp/colB0
    plt.plot(sortedTVals,rhoTmp,color="green")
    plt.scatter(sortedTVals,rhoTmp,color="magenta")
    plt.xlabel("T")
    plt.ylabel("$\\rho^{BB}($"+str(j)+"$)$")


    plt.savefig(pathData+"/rhoBB"+str(j)+".pdf")

    plt.close()


def plotGAB(j):
    """

    :param j: distance, -(N-1), ..., -1, 0, 1,...,N-1
    :return: correlation function of xA, xB
    """
    ind=j+N-1
    colABTmp=fABArray[:,ind]
    plt.plot(sortedTVals,colABTmp,color="blue")
    plt.scatter(sortedTVals,colABTmp,color="red")
    plt.xlabel("T")
    plt.ylabel("$G^{AB}($"+str(j)+"$)$")
    plt.savefig(pathData+"/GAB"+str(j)+".pdf")

    plt.close()


def plotrAB(j):
    """

    :param j: distance, -(N-1), ..., -1, 0, 1,...,N-1
    :return: correlation coefficient of xA, xB
    """
    ind=j+N-1
    colABTmp=fABArray[:,ind]
    colA0=fAArray[:,0]
    colB0=fBArray[:,0]
    sdAB=np.sqrt(colA0*colB0)
    rhoTmp=colABTmp/sdAB
    plt.plot(sortedTVals,rhoTmp,color="green")
    plt.scatter(sortedTVals,rhoTmp,color="magenta")
    plt.xlabel("T")
    plt.ylabel("$\\rho^{AB}($"+str(j)+"$)$")


    plt.savefig(pathData+"/rhoAB"+str(j)+".pdf")

    plt.close()



plotGAA(int(N/2))
plotrAA(int(N/2))

plotGBB(int(N/2))
plotrBB(int(N/2))

plotGAB(int(N/2))
plotrAB(int(N/2))

plotGAB(-int(N/2))
plotrAB(-int(N/2))
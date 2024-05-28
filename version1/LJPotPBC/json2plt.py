import numpy as np
import glob
import sys
import re
import matplotlib.pyplot as plt
from datetime import datetime
import json

#This script loads json data and plot


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


def pltU(oneTFile):
    """

    :param oneTFile: corresponds to one temperature
    :return: U plots
    """
    matchT=re.search(r"T(\d+(\.\d+)?)",oneTFile)
    TVal=float(matchT.group(1))
    UFilePath=oneTFile+"/jsonData/jsonU/UData.json"
    with open(UFilePath, 'r') as fptr:
        data = json.load(fptr)
    UVec=np.array(data["U"])
    meanU=np.mean(UVec)

    varU=np.var(UVec,ddof=1)
    sigmaU=np.sqrt(varU)

    nbins=500
    fig=plt.figure()
    axU=fig.add_subplot()
    (n0,_,_)=axU.hist(UVec,bins=nbins)
    meanU=np.round(meanU,4)
    sigmaU=np.round(sigmaU,4)

    axU.set_title("T="+str(np.round(TVal,3)))
    axU.set_xlabel("$U$")
    axU.set_ylabel("#")
    xPosUText=(np.max(UVec)-np.min(UVec))*1/2+np.min(UVec)
    yPosUText=np.max(n0)*2/3
    axU.text(xPosUText,yPosUText,"mean="+str(meanU)+"\nsd="+str(sigmaU))
    plt.axvline(x=meanU,color="red",label="mean")
    axU.text(meanU*1.1,0.5*np.max(n0),str(meanU)+"$\pm$"+str(sigmaU),color="red")
    axU.hlines(y=0,xmin=meanU-sigmaU,xmax=meanU+sigmaU,color="green",linewidth=15)

    plt.legend(loc="best")

    EHistOut="T"+str(TVal)+"UHist.png"
    plt.savefig(oneTFile+"/"+EHistOut)

    plt.close()

    ### test normal distribution for mean U
    #block mean
    USelectedAll=UVec
    def meanPerBlock(length):
        blockNum=int(np.floor(len(USelectedAll)/length))
        UMeanBlock=[]
        for blkNum in range(0,blockNum):
            blkU=USelectedAll[blkNum*length:(blkNum+1)*length]
            UMeanBlock.append(np.mean(blkU))
        return UMeanBlock
    fig=plt.figure(figsize=(20,20))
    fig.tight_layout(pad=5.0)
    lengthVals=[5,10,20,40]
    for i in range(0,len(lengthVals)):
        l=lengthVals[i]
        UMeanBlk=meanPerBlock(l)
        ax=fig.add_subplot(2,2,i+1)
        (n,_,_)=ax.hist(UMeanBlk,bins=100,color="aqua")
        xPosTextBlk=(np.max(UMeanBlk)-np.min(UMeanBlk))*1/7+np.min(UMeanBlk)
        yPosTextBlk=np.max(n)*3/4
        meanTmp=np.mean(UMeanBlk)
        meanTmp=np.round(meanTmp,3)
        sdTmp=np.sqrt(np.var(UMeanBlk))
        sdTmp=np.round(sdTmp,3)
        ax.set_title("L="+str(l))
        ax.text(xPosTextBlk,yPosTextBlk,"mean="+str(meanTmp)+", sd="+str(sdTmp))
    fig.suptitle("T="+str(TVal))
    plt.savefig(oneTFile+"/T"+str(TVal)+"UBlk.png")
    # plt.savefig(EBlkMeanDir+"/T"+str(TTmp)+"EBlk.png")
    plt.close()


def plt_xAxB(oneTFile):
    """

    :param oneTFile: corresponds to one temperature
    :return: plots of xA, xB
    """
    matchT=re.search(r"T(\d+(\.\d+)?)",oneTFile)
    TVal=float(matchT.group(1))
    pathAllForAUnitCell=[]
    jsonPath=oneTFile+"/jsonData/"
    for cellFile in glob.glob(jsonPath+"/jsonUnitCell*"):
        pathAllForAUnitCell.append(cellFile)
    #sort files

    unitCellIndex=[]
    for file in pathAllForAUnitCell:
        matchInd=re.search(r"jsonUnitCell(\d+)",file)
        unitCellIndex.append(int(matchInd.group(1)))

    sortedInd=np.argsort(unitCellIndex)
    sortedPathForAUnitCell=[pathAllForAUnitCell[j] for j in sortedInd]

    fig=plt.figure(figsize=(20,160))
    fig.tight_layout(pad=5.0)
    x_vertical_distance = 0.9
    xAMeanAll=[]
    xASdAll=[]
    xBMeanAll=[]
    xBSdAll=[]
    unitCellNum=len(sortedPathForAUnitCell)
    for j in range(0,unitCellNum):
        onePath=sortedPathForAUnitCell[j]
        jsonFile=onePath+"/xAxBData.json"
        with open(jsonFile, 'r') as fptr:
            dataAB = json.load(fptr)
        xAVec=np.array(dataAB["xA"])
        xBVec=np.array(dataAB["xB"])



        axx=fig.add_subplot(unitCellNum,1,j+1,sharex=axx if j != 0 else None)

        xAMean=np.mean(xAVec)
        xAMeanAll.append(xAMean)
        xBMean=np.mean(xBVec)
        xBMeanAll.append(xBMean)

        xAVar=np.var(xAVec,ddof=1)

        xBVar=np.var(xBVec,ddof=1)

        xASigma=np.sqrt(xAVar)
        xASdAll.append(xASigma)
        xBSigma=np.sqrt(xBVar)
        xBSdAll.append(xBSigma)

        nbins=500
        #plot A
        (nA,_,_)=axx.hist(xAVec,bins=nbins,color = "blue", ec="blue")
        xAMean=np.round(xAMean,4)
        xASigma=np.round(xASigma,4)

        axx.set_title("Unit cell "+str(j)+", T="+str(np.round(TVal,3)))

        plt.axvline(x=xAMean,color="red",label="mean A")
        axx.text(xAMean*1.1,0.5*np.max(nA),str(xAMean)+"$\pm$"+str(xASigma),color="red")

        #plot B
        nbins=500
        (nB,_,_)=axx.hist(xBVec,bins=nbins,color = "green", ec="green")
        plt.axvline(x=xBMean,color="magenta",label="mean B")
        axx.text(xBMean*1.1,0.5*np.max(nB),str(xBMean)+"$\pm$"+str(xBSigma),color="magenta")
        axx.set_ylabel("#")
        axx.set_xlabel("position")
    xHistOut="T"+str(TVal)+"xHist.pdf"
    plt.subplots_adjust(hspace=x_vertical_distance)
    plt.savefig(oneTFile+"/"+xHistOut)
    plt.close()

    #summary of distance between neighboring points

    smrDistFileName=oneTFile+"/distSummary.txt"
    #intracell
    vecDiffIntraCell=[xBMeanAll[i]-xAMeanAll[i] for i in range(0,len(xBMeanAll))]
    #intercell
    vecDiffInterCell=[xAMeanAll[i+1]-xBMeanAll[i] for i in range(0,len(xBMeanAll)-1)]

    fptrTxt=open(smrDistFileName,"w")
    fptrTxt.write("T="+str(TVal)+"\n")
    fptrTxt.write("Intracell: "+str(vecDiffIntraCell)+"\n")
    fptrTxt.write("Intercell: "+str(vecDiffInterCell)+"\n")
    fptrTxt.close()



    plt.figure(figsize=(12, 6))
    plt.ylim(-1, 1)


    for i in range(0,len(xAMeanAll)):
        plt.hlines(y=0,xmin=xAMeanAll[i]-xASdAll[i],xmax=xAMeanAll[i]+xASdAll[i],color="red",linewidth=2,alpha=0.2)
        plt.text(xAMeanAll[i],0.3,str(i)+"A",color="blue", ha='center')
        plt.text(xAMeanAll[i],-0.1,str(np.round(xAMeanAll[i],4)),color="blue", ha='center',fontsize=8)

    plt.scatter(xAMeanAll,[0]*len(xAMeanAll),color="blue",s=8,label="A")


    for i in range(0,len(xBMeanAll)):
        plt.hlines(y=0,xmin=xBMeanAll[i]-xBSdAll[i],xmax=xBMeanAll[i]+xBSdAll[i],color="magenta",linewidth=2,alpha=0.2)
        plt.text(xBMeanAll[i],0.1,str(i)+"B",color="green", ha='center')
        plt.text(xBMeanAll[i],-0.3,str(np.round(xBMeanAll[i],4)),color="green", ha='center',fontsize=8)

    plt.scatter(xBMeanAll,[0]*len(xBMeanAll),color="green",s=8,label="B")
    plt.legend(loc="best")
    plt.gca().get_yaxis().set_visible(False)
    plt.axhline(y=0, color='black', linewidth=0.5,alpha=0.3)
    gridOut="T"+str(TVal)+"grid.pdf"
    plt.title("T="+str(np.round(TVal,4)))
    plt.savefig(oneTFile+"/"+gridOut)

    plt.close()






tStatsStart=datetime.now()
for oneTFile in sortedTFiles:
    pltU(oneTFile)
    plt_xAxB(oneTFile)

tStatsEnd=datetime.now()
print("stats total time: ",tStatsEnd-tStatsStart)
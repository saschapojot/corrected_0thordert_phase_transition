import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pathlib import Path
import sys
from scipy.optimize import root

#python readCSV.py path rowNum
if len(sys.argv)!=3:
    print("wrong number of arguments")

fileName=sys.argv[1]
rowNum=int(sys.argv[2])

dfStr=pd.read_csv(fileName)

oneRow=dfStr.iloc[rowNum,:]

alpha1=float(oneRow.loc["alpha1"])
beta1=float(oneRow.loc["beta1"])
p1=float(oneRow.loc["p1"])
q1=float(oneRow.loc["q1"])


alpha2=float(oneRow.loc["alpha2"])
beta2=float(oneRow.loc["beta2"])
p2=float(oneRow.loc["p2"])
q2=float(oneRow.loc["q2"])

print("alpha1"+str(alpha1)+"beta1"+str(beta1)+"p1"+str(p1)+"q1"+str(q1)
      +"alpha2"+str(alpha2)+"beta2"+str(beta2)+"p2"+str(p2)+"q2"+str(q2))



def dV1(r):
    return -alpha1*p1*r**(-p1-1)+beta1*q1*r**(-q1-1)+4*r**3


def dV2(r):
    return -alpha2*p2*r**(-p2-1)+beta2*q2*r**(-q2-1)+4*r**3
def dV(r):
    return dV1(r)+dV2(r)

sol=root(dV,1,method="broyden2",tol=1e-9)

x=sol.x
print("x="+str(x))

# rValsAll=np.linspace(0.5,1.2,100)
# dV1ValsAll=dV1(rValsAll)
# dV2ValsAll=dV2(rValsAll)
# plt.figure()
# plt.plot(rValsAll,dV1ValsAll,color="black")
# plt.title("$dV_{1}$")
# plt.savefig("tmpdV1.png")
# plt.close()
#
# plt.figure()
# plt.plot(rValsAll,dV2ValsAll,color="blue")
# plt.title("$dV_{2}$")
# plt.savefig("tmpdV2.png")
# plt.close()

# sol1=root(dV1,1,method="broyden2",tol=1e-9)
#
# sol2=root(dV2,1,method="broyden2",tol=1e-9)
#
# print(sol1)
# print(sol2)
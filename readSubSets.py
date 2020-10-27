import glob
from netCDF4 import Dataset
import sys

argv=sys.argv
print(argv)
if len(argv)==2:
    iz=int(argv[1])
else:
    iz=10


def readZeroDeg(fname,nr1,nr2):
    f=Dataset(fname)
    zeroDeg=f['/NS/VER/binZeroDeg'][:,nr1:nr2]
    binClut=f['/NS/PRE/binClutterFreeBottom'][:,nr1:nr2]
    binSf=f['/NS/PRE/binRealSurface'][:,nr1:nr2]
    sfcRain=f['NS/SLV/precipRateNearSurface'][:,nr1:nr2]
    zFactKu=f['NS/PRE/zFactorMeasured'][:,nr1:nr2,:]
    zFactKa=f['MS/PRE/zFactorMeasured'][:,:,:]
    relFlag=f['NS/SRT/reliabFlag'][:,nr1:nr2]
    piaSRT=f['NS/SRT/pathAtten'][:,nr1:nr2]
    precipType=f['/NS/CSF/typePrecip'][:,nr1:nr2]
    precipType=(precipType/1e7).astype(int)
    stormT=f['NS/PRE/heightStormTop'][:,nr1:nr2]
    sfcType=f['/NS/PRE/landSurfaceType'][:,nr1:nr2]
    dm=f['NS/SLV/paramDSD'][:,nr1:nr2,:,1]
    return zeroDeg,sfcRain,binSf,zFactKu,binClut,piaSRT,relFlag,precipType,zFactKa,stormT,sfcType,dm
fs=glob.glob("SEAsiaCS/*HDF5")
zeroDegMinL=[]
from numpy import *
iprof=0
z1L=[]
zKuL=[]
zKaL=[]
addInfoL=[]
for f in fs[:]:
    nr1=12
    nr2=37
    zeroDeg,sfcRain,binSf,zFactKu,binClut,piaSRT,relFlag,pType,zFactKa,stormT,sfcType,dm=readZeroDeg(f,nr1,nr2)
    #stop
    a=nonzero(sfcRain>0)
    b=nonzero(zeroDeg[a]>0)
    for i1,j1 in zip(a[0][b],a[1][b]):
        if binSf[i1,j1]>157 and pType[i1,j1]==2:
            if binClut[i1,j1]-zeroDeg[i1,j1]<iz:
                continue
            else:
                iprof+=1
                #zprof1d=zeros(164,float)
                zKuL.append(zFactKu[i1,j1,:])
                zKaL.append(zFactKa[i1,j1,:])
                #print(zeroDeg[i1,j1]+8,binClut[i1,j1])
                zKu12=zFactKu[i1,j1,zeroDeg[i1,j1]:zeroDeg[i1,j1]+iz]
                zKu11=zFactKu[i1,j1,zeroDeg[i1,j1]-60:zeroDeg[i1,j1]]
                if zeroDeg[i1,j1]+iz>binClut[i1,j1]:
                    print("error")
                    stop
                z1prof1d=list(zKu11)
                z1prof1d.extend(zKu12)
                z1prof1d.append(zFactKu[i1,j1,binClut[i1,j1]])
                z1L.append(z1prof1d)
                if relFlag[i1,j1]==1 or relFlag[i1,j1]==2:
                    reliabF1=1
                else:
                    reliabF1=0
                addInfoL.append([(175-binSf[i1,j1])*0.125,(binSf[i1,j1]-binClut[i1,j1])*0.125,(175-zeroDeg[i1,j1])*0.125,\
                binSf[i1,j1],zeroDeg[i1,j1],binClut[i1,j1],piaSRT[i1,j1],dm[i1,j1,binClut[i1,j1]],\
                sfcRain[i1,j1],stormT[i1,j1],\
                sfcType[i1,j1],reliabF1])


import matplotlib.pyplot as plt
from numpy import *
z1L=array(z1L)
import pickle
pickle.dump([z1L,zKuL,zKaL,addInfoL],open('cvProfs%2.2i.pklz'%iz,'wb'))

z1L[z1L<0]=0
plt.plot(z1L.mean(axis=0))

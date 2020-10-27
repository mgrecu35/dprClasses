#addInfoL.append([(175-binSf[i1,j1])*0.125,(binSf[i1,j1]-binClut[i1,j1])*0.125,(175-zeroDeg[i1,j1])*0.125,\
#binSf[i1,j1],zeroDeg[i1,j1],binClut[i1,j1],piaSRT[i1,j1],sfcRain[i1,j1]])


import pickle
[z1L,zKuL,zKaL,addInfoL]=pickle.load(open('cvProfs.pklz','rb'))
n=len(zKuL)
zeta1L=[]
zeta2L=[]
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
addInfoL=np.array(addInfoL)
addInfoL[:,-5]/=(addInfoL[:,2]-addInfoL[:,0])
z1L[z1L<0]=0
from sklearn.cluster import KMeans
kmeans = KMeans(n_clusters=25, random_state=0).fit(z1L[:,:-1])
kmeansH = KMeans(n_clusters=3, random_state=0)
addInfoL=np.array(addInfoL)
#plt.suptitle("relative")
h1=60*0.125-np.arange(80)*0.125
plt.figure(figsize=(12,8))
piaClass=np.zeros((25),float)
piaClass2=np.zeros((25),float)
#ic=np.array([ 3, 15,  5,  1,  8, 12,  0, 13,  2,  7, 14, 11,  9,  6,  4, 10])

ic=range(25)
for i in range(25):
    a=np.nonzero(kmeans.labels_==ic[i])
    plt.subplot(5,5,i%25+1)
    for j in range(25):
        if i!=j:
            a1=np.nonzero(kmeans.labels_==j)
            plt.plot(z1L[a1[0],:-1].mean(axis=0),h1,color="gray",linewidth=0.5)
    plt.plot(z1L[a[0],:-1].mean(axis=0),h1)
    freq=len(a[0])/z1L.shape[0]
    a1=np.nonzero(addInfoL[a[0],-2]>0)
    b1=np.nonzero(addInfoL[a[0],-1][a1]==1)
    piaMean=addInfoL[a[0],-5][a1][b1].mean()
    a11=np.nonzero(addInfoL[a[0],-2]==0)
    b11=np.nonzero(addInfoL[a[0],-1][a11]==1)
    piaMean2=addInfoL[a[0],-5][a11][b11].mean()
    piaClass[i]=piaMean
    piaClass2[i]=piaMean2
    plt.title("Freq=%6.2f Land_fract=%6.2f\n piaM=(%6.2f,%6.2f)"%(freq,len(a1[0])/len(a[0]),piaMean2,piaMean))
    if i%5==0:
        plt.ylabel("Height (km)")
    if i>=20:
        plt.xlabel("dBZ")

plt.figure(figsize=(12,8))
ic=np.argsort(piaClass)

#ic=range(16)
kgains=[]
rEst=addInfoL[:,-5].copy()*0.0
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import DotProduct, WhiteKernel
kernel = DotProduct() + WhiteKernel()

for i in range(25):
    a=np.nonzero(kmeans.labels_==ic[i])
    plt.subplot(5,5,i%25+1)
    for j in range(25):
        if i!=j:
            a1=np.nonzero(kmeans.labels_==j)
            plt.plot(z1L[a1[0],:-1].mean(axis=0),h1,color="gray",linewidth=0.5)
    plt.plot(z1L[a[0],:-1].mean(axis=0),h1)
    covT=np.cov(z1L[a[0],:].T,(addInfoL[a[0],-5:-4]).T)
    nx=z1L.shape[1]
    covyy=covT[0:nx,0:nx]+np.eye(nx)*1.0
    covxy=covT[nx:,:nx]
    kgain=np.dot(covxy,np.linalg.pinv(covyy))
    zm=z1L[a[0],:].mean(axis=0)
    rm=(addInfoL[a[0],-5]).mean()
    rEstL=[]
    for ik in a[0]:
        r1=rm+np.dot(kgain,(z1L[ik,:]-zm))
        rEstL.append((r1[0]))
    rEst[a[0]]=rEstL
    X=z1L[a[0],:]
    y=addInfoL[a[0],-5]
    freq=len(a[0])/z1L.shape[0]
    sfcRainM=addInfoL[a[0],-4].mean()
    a1=np.nonzero(addInfoL[a[0],-2]>0)
    b1=np.nonzero(addInfoL[a[0],-1][a1]==1)
    piaMean=addInfoL[a[0],-5][a1][b1].mean()
    a11=np.nonzero(addInfoL[a[0],-2]==0)
    b11=np.nonzero(addInfoL[a[0],-1][a11]==1)
    piaMean2=addInfoL[a[0],-5][a11][b11].mean()
    piaClass[i]=piaMean
    piaClass2[i]=piaMean2
    plt.title("Freq=%6.2f Land_fract=%6.2f\n piaM=(%6.2f,%6.2f,%6.2f)"%(freq,len(a1[0])/len(a[0]),piaMean2,piaMean,sfcRainM))
    if i%5==0:
        plt.ylabel("Height (km)")
    if i>=20:
        plt.xlabel("dBZ")


plt.tight_layout()

using PyCall

include("readCMBTables.jl")

zKuJ,attKuJ,pwcJ,rrateJ,dmJ,zKuBBJ,attKuBBJ,pwcBBJ,rrateBBJ,dmBBJ,
zKuSJ,attKuSJ,pwcSJ,rrateSJ,dmSJ,nbinsj,
zKaJ,attKaJ,zKaBBJ,attKaBBJ,zKaSJ,attKaSJ,dNw,dNwBB,dNwSJ=readCMBTables()

pickle=pyimport("pickle")
fh=pybuiltin("open")("cvProfs20.pklz","rb")
np=pyimport("numpy")
z1L,zKuL,zKaL,addInfoL=pickle.load(fh)
addInfoL=np.array(addInfoL)
#addInfoL[:,-6]/=(addInfoL[:,2]-addInfoL[:,0])

h0=addInfoL[:,3]
hsfc=addInfoL[:,1]
fh.close()
fh=pybuiltin("open")("kFilterCvClasses20.pklz","rb")
zmL,kgainL,xL,kmeans,ic=pickle.load(fh)

using PyCall

include("readCMBTables.jl")
include("src_jl/psdInt_2.jl")
include("hbprof.jl")
newTables=readCMBTables()
KuJ,attKuJ,pwcJ,rrateJ,dmJ,zKuBBJ,attKuBBJ,pwcBBJ,rrateBBJ,dmBBJ,
zKuSJ,attKuSJ,pwcSJ,rrateSJ,dmSJ,nbinsj,
zKaJ,attKaJ,zKaBBJ,attKaBBJ,zKaSJ,attKaSJ,dNw,dNwBB,dNwSJ,
zKuHJ,attKuHJ,pwcHJ,rrateHJ,dmHJ,
zKaHJ,attKaHJ,pwcHJ2,rrateHJ2,dmHJ2=newTables
pickle=pyimport("pickle")
fh=pybuiltin("open")("cvProfs20.pklz","rb")
np=pyimport("numpy")
z1L,zKuL,zKaL,dmL,nwL,rateL,addInfoL=pickle.load(fh)
addInfoL=np.array(addInfoL)
#addInfoL[:,-6]/=(addInfoL[:,2]-addInfoL[:,0])

h0=addInfoL[:,3]
hsfc=addInfoL[:,1]
fh.close()
fh=pybuiltin("open")("kFilterCvClasses20.pklz","rb")
zmL,kgainL,xL,kmeans,ic=pickle.load(fh)
z1L[z1L.<0].=0.0;
labels=kmeans.predict(z1L[:,1:end-1]);

indx=findall(labels.==ic[end-3])
i=1

dr=0.25
dmRetL=[]
rmRetL=[]
piaL=[]
z1pL=[]
n1L=[]
dm1L=[]
piaBL=[]
rrate1DL=[]
piaKuL=[]
zKuL1=[]
zKuL2=[]
for i=1:1200
    z1p=zeros(176)
    bsfc,bzd,bcf=Int.(addInfoL[indx[i],4:6])
    piaKu=0.0
    btop=175-Int(trunc(addInfoL[indx[i],10]/125))
    dr=0.125
    zmax=maximum(zKuL[indx[i]][1:bzd+5])
    pia=max(addInfoL[indx[i],7],0.1)
    #piaMax=zmax-zKuL[indx[i]][bcf+1]+np.random.rand()*2
    piaMax=pia
    if addInfoL[indx[i],end]==1 || addInfoL[indx[i],end]==2
        pia=addInfoL[indx[i],7]
        ifull=0
        dm1d,dnw1d,piaB,rrate1D,zKuC,piaKu,zKuSim=iter_prof(btop,bzd,bcf,bsfc,indx,i,zKuL,zKaL,dr,piaMax,newTables,ifull)
        push!(dm1L,[dm1d[bcf+1],dm1d[bzd+1],addInfoL[indx[i],8]])
        push!(n1L,[dnw1d[bcf+1],dnw1d[bzd+1]])
        push!(piaBL,pia)
        push!(rmRetL,[addInfoL[indx[i],9],rrate1D[bcf+1]])
        push!(rrate1DL,rrate1D)
        push!(piaKuL,piaKu)
        push!(zKuL1,zKuL[indx[i]])
        push!(zKuL2,zKuSim)
    end
    #if addInfoL[indx[i],end]==1 || addInfoL[indx[i],end]==2
    #
        #dnv=zeros(88).+0.0
        #z13obs,z13c,eps,dm_hb,sfcRain_hb,piaHB=hb_cv(btop,bzd,bcf,bsfc,indx,i,zKuL,zKaL,dnv,dr,newTables)
        #push!(dmRetL,[addInfoL[indx[i],8],dm_hb])
        #push!(rmRetL,[addInfoL[indx[i],9],sfcRain_hb])
        #
    #end

end

rrate1DL=np.array(rrate1DL)
plt.pcolormesh(rrate1DL[1:200,169:-1:100]',cmap="jet")
zKuL1=np.array(zKuL1)
zKuL2=np.array(zKuL2)
iplot=1
if iplot==1
plt.figure()
plt.subplot(211)
plt.pcolormesh(zKuL1[1:200,165:-1:100]',cmap="jet",vmax=50,vmin=0)
plt.colorbar()
plt.subplot(212)
plt.pcolormesh(zKuL2[1:200,165:-1:100]',cmap="jet",vmax=50,vmin=0)
plt.colorbar()
end

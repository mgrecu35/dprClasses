using PyCall

include("readCMBTables.jl")
include("src_jl/psdInt_2.jl")
include("hbprof2columns.jl")
newTables=readCMBTables()
zKuJ,attKuJ,pwcJ,rrateJ,dmJ,zKuBBJ,attKuBBJ,pwcBBJ,rrateBBJ,dmBBJ,
zKuSJ,attKuSJ,pwcSJ,rrateSJ,dmSJ,nbinsj,
zKaJ,attKaJ,zKaBBJ,attKaBBJ,zKaSJ,attKaSJ,dNw,dNwBB,dNwSJ,
zKuHJ,attKuHJ,pwcHJ,rrateHJ,dmHJ,
zKaHJ,attKaHJ,pwcHJ2,rrateHJ2,dmHJ2=newTables
pickle=pyimport("pickle")
fh=pybuiltin("open")("cvProfs20.pklz","rb")
np=pyimport("numpy")
z1L,zKuLt,zKaLt,dmL,nwL,rateL,addInfoL=pickle.load(fh)
n=size(zKuLt)[1]
zKuL=[]
zKaL=[]
for i=1:n
    push!(zKuL,zKuLt[i][2,:])
    push!(zKaL,zKaLt[i][2,:])
end
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
nMemb=50
zKuSimEns=zeros(176,nMemb)
rrate1DEns=zeros(176,nMemb)
epsEns=zeros(1,nMemb)
piaEns=zeros(1,nMemb)
for i=1:300
    z1p=zeros(176)
    bsfc,bzd,bcf=Int.(addInfoL[indx[i],4:6])
    bzd=bzd+6
    piaKu=0.0
    btop=175-Int(trunc(addInfoL[indx[i],10]/125))
    dr=0.125
    zmax=maximum(zKuL[indx[i]][1:bzd+5])
    pia=max(addInfoL[indx[i],7],0.1)
    zKu1=copy(zKuL[indx[i]])
    zKu2=copy(zKuL[indx[i]])*0
    for k=bzd:bcf
        if(zKu1[k]>-10)
            zKu1[k]=zKuL[indx[i]][k]-25*(k-bzd)/(bcf-bzd)
            zKu2[k]=10*log10(10^(0.1*zKuL[indx[i]][k])-10^(0.1*zKu1[k])+1e-9)
        end
    end
    push!(zKuL1,zKu1)
    push!(zKuL2,zKu2)
end



#plt.plot(zKuSimF)
#plt.plot(zKuL[indx[i]])
#plt.ylim(0,45)
#rrate1DL=np.array(rrate1DL)
#plt.pcolormesh(rrate1DL[1:100,169:-1:100]',cmap="jet")
zKuL1=np.array(zKuL1)
zKuL2=np.array(zKuL2)
iplot=1
if iplot==1
plt.figure()
plt.subplot(211)
plt.pcolormesh(zKuL1[1:100,165:-1:100]',cmap="jet",vmax=50,vmin=0)
plt.colorbar()
plt.subplot(212)
plt.pcolormesh(zKuL2[1:100,165:-1:100]',cmap="jet",vmax=50,vmin=0)
plt.colorbar()
end
rrateSeto=zeros(100)
dmSeto=zeros(100)
zKuSeto=zeros(100)
for i=1:100
    rrateSeto[i]=i
    dmSeto[i]=(rrateSeto[i]/1.37)^(1/5.420)
    n11,n21=bisection(dmJ,dmSeto[i])
    dn=log10((rrateSeto[i]/rrateJ[n11]))
    zKuSeto[i]=zKuJ[n11]+10*dn
end
rrate1L=[]
rrate1dL=[]
sfcRainL=[]
piaL=[]
for i=1:1200
    dr=0.125
    bsfc,bzd,bcf=Int.(addInfoL[indx[i],4:6])
    piaKu=0.0
    btop=175-Int(trunc(addInfoL[indx[i],10]/125))
    dr=0.125
    zmax=maximum(zKuL[indx[i]][1:bzd+5])
    pia=max(addInfoL[indx[i],7],0.1)
    zKu1=copy(zKuL[indx[i]])
    zKu2=copy(zKuL[indx[i]])*0
    bzd+=0
    for k=bzd:bcf
        if(zKu1[k]>-10)
            zKu1[k]=zKuL[indx[i]][k]-0.1*(k-bzd)/(bcf-bzd)
            zKu2[k]=10*log10(10^(0.1*zKuL[indx[i]][k])-10^(0.1*zKu1[k])+1e-9)
        end
    end
    for k=1:bzd-1
        zKu1[k]+=0
    end

    dm1dd, dn1dd,piaKud, rrate1dd,rrate1ds,zKuCd,piaKu,piaKa=prof(btop,bzd,bcf,bsfc,zKu1,zKaL[indx[i]],
    dr,pia,newTables)
    zKu11=zKuLt[indx[i]][1,:]
    zKa11=zKaLt[indx[i]][1,:]
    dm1dd1, dn1dd1,piaKud1, rrate1dd1,rrate1ds1,zKuCd1,piaKu1,piaKa1=prof(btop,bzd,bcf,bsfc,zKu11,zKa11,
    dr,pia,newTables)
    zKu12=zKuLt[indx[i]][3,:]
    zKa12=zKaLt[indx[i]][3,:]
    dm1dd2, dn1dd2,piaKud2, rrate1dd2,rrate1ds2,zKuCd2,piaKu2,piaKa2=prof(btop,bzd,bcf,bsfc,zKu12,zKa12,
    dr,pia,newTables)
    appPIA_ku=-10*log10((10^(-0.1*piaKu1)+10^(-0.1*piaKu)+10^(-0.1*piaKu2))/3.0)
    appPIA_ka=-10*log10((10^(-0.1*piaKa1)+10^(-0.1*piaKa)+10^(-0.1*piaKa2))/3.0)
    nubf=-10*log10((10^(-0.1*piaKu1)+10^(-0.1*piaKu)+10^(-0.1*piaKu2))/3.0)/((piaKu+piaKu1+piaKu2)/3.0)
    nubf2=-10*log10((10^(-0.1*piaKa1)+10^(-0.1*piaKa)+10^(-0.1*piaKa2))/3.0)/((piaKa+piaKa1+piaKa2)/3.0)
    dm1d, dn1d,piaKu2, rrate1d,zKuC=profw(btop,bzd,bcf,bsfc,zKu2,dr,pia,newTables)
    if addInfoL[indx[i],end]==1
        piaU=pia/nubf*1.0
        dm1dd, dn1dd,piaKud, rrate1dd,rrate1ds,zKuCd,piaKu,piaKa=profPIA(btop,bzd,bcf,bsfc,zKu1,zKaL[indx[i]],
        dr,piaU,newTables)
        push!(rrate1L,rrate1ds[:])
        push!(rrate1dL,rrate1dd[:])
        push!(sfcRainL,[addInfoL[indx[i],9],rrate1dd[bcf+1]])
        push!(piaL,[piaKu,pia,nubf,nubf2,(piaKu+piaKu1+piaKu2)/3.0,(piaKa+piaKa1+piaKa2)/3.0,
        appPIA_ku,appPIA_ka])
    end
end
rrate1L=np.array(rrate1L)
rrate1dL=np.array(rrate1dL)
iplot=1
if iplot==1
plt.figure()
plt.subplot(111)
plt.pcolormesh(0.975*(rrate1dL)[1:500,165:-1:100]',cmap="jet",vmax=200,vmin=0)
plt.colorbar()
end
plt.figure();plt.plot(np.mean(rrate1dL,axis=0))
plt.plot(np.mean(rrate1L,axis=0))
plt.plot(np.mean(rrate1L,axis=0).+np.mean(rrate1dL,axis=0))

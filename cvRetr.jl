using PyCall

include("readCMBTables.jl")
include("src_jl/psdInt_2.jl")
include("hbprof.jl")
newTables=readCMBTables()
zKuJ,attKuJ,pwcJ,rrateJ,dmJ,zKuBBJ,attKuBBJ,pwcBBJ,rrateBBJ,dmBBJ,
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

indx=findall(labels.==ic[end])
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
    piaKu=0.0
    btop=175-Int(trunc(addInfoL[indx[i],10]/125))
    dr=0.125
    zmax=maximum(zKuL[indx[i]][1:bzd+5])
    pia=max(addInfoL[indx[i],7],0.1)
    zKu1=copy(zKuL[indx[i]])
    zKu2=copy(zKuL[indx[i]])*0
    for k=bzd:bcf
        if(zKu1[k]>-10)
            zKu1[k]=zKuL[indx[i]][k]-10*(k-bzd)/(bcf-bzd)
            zKu2[k]=10*log10(10^(0.1*zKuL[indx[i]][k])-10^(0.1*zKu1[k])+1e-9)
        end
    end
    push!(zKuL1,zKu1)
    push!(zKuL2,zKu2)
end

for i=1:-300
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
        for ien=1:nMemb
            dm1d,dnw1d,piaB,rrate1D,zKuC,piaKu,zKuSim,eps=iter_prof(btop,bzd,bcf,bsfc,indx,i,zKuL,
            zKaL,dr,piaMax,newTables,ifull)
            zKuSimEns[:,ien]=zKuSim
            zKuSimEns[bcf+1,ien]=10*piaKu
            rrate1DEns[:,ien]=rrate1D
            epsEns[1,ien]=eps
            piaEns[1,ien]=piaKu
        end
        zObs=zKuL[indx[i]][bzd+1:bcf+1]
        zObs[end]=10*pia
        println("$(pia) $(np.mean(piaEns))")
        nx=bcf-bzd+1
        covT=(np.cov(rrate1DEns[bzd+1:bcf+1,:],zKuSimEns[bzd+1:bcf+1,:]))
        covTpia=(np.cov(epsEns[1,:],zKuSimEns[bzd+1:bcf+1,:]))
        covXY=covT[1:nx,nx+1:end]
        #covXYpia=covTpia[1,2:end]
        covYY=covT[nx+1:end,nx+1:end].+2*np.eye(nx)
        covYY[end,end]=2.0
        kgain=np.dot(covXY,np.linalg.pinv(covYY))
        #kgainPia=np.dot(covXYpia,np.linalg.pinv(covYY))
        dZ=zObs.-np.mean(zKuSimEns[bzd+1:bcf+1,:],axis=1)
        dZ[zObs.<0].=0
        #dZ[end]=0
        dR=np.dot(kgain,dZ)
        #dEps=np.dot(kgainPia,dZ)
        rrate1D=np.mean(rrate1DEns,axis=1)
        rrate1D[bzd+1:bcf+1].=rrate1D[bzd+1:bcf+1].+0.85*dR
        rrate1D[rrate1D.<0.01].=0.01
        dr=0.125
        dm1d,dnw1d,piaB,rrate1D,zKuC,piaKu,zKuSimF=iter_profR(btop,bzd,bcf,bsfc,indx,i,zKuL,
        zKaL,dr,rrate1D,piaMax,exp(-0.5),newTables,ifull)

        push!(dm1L,[dm1d[bcf+1],dm1d[bzd+1],addInfoL[indx[i],8]])
        push!(n1L,[dnw1d[bcf+1],dnw1d[bzd+1]])
        push!(piaBL,pia)
        push!(rmRetL,[addInfoL[indx[i],9],rrate1D[bcf]])
        push!(rrate1DL,rrate1DEns[:,1])
        push!(piaKuL,piaKu)
        push!(zKuL1,zKuL[indx[i]])
        push!(zKuL2,zKuSimF)
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
i=1

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
    n11,n21=bisection(dmJ,dm)
    dn=log10((rrateSeto[i]/rrateJ[n11]))
    zKuSeto[i]=zKuJ[n11]*10^dn
end

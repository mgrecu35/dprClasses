
function back_prof(btop,bzd,bcf,bsfc,indx,i,zKuL,zKaL,dr,pia,newTables)
    KuJ,attKuJ,pwcJ,rrateJ,dmJ,zKuBBJ,attKuBBJ,pwcBBJ,rrateBBJ,dmBBJ,
    zKuSJ,attKuSJ,pwcSJ,rrateSJ,dmSJ,nbinsj,
    zKaJ,attKaJ,zKaBBJ,attKaBBJ,zKaSJ,attKaSJ,dNw,dNwBB,dNwSJ,
    zKuHJ,attKuHJ,pwcHJ,rrateHJ,dmHJ=newTables
    piaB=pia
    ztrue=zKuL[indx[i]][bcf+1]+piaB
    dm=0.9*(0.00113070*ztrue^2+0.0047*ztrue+0.4911)
    n1,n2=bisection(dmJ,dm)
    dn=(ztrue-zKuJ[n1])/10.0
    if dm>3.0
        dm=3
    end
    dm1d=zeros(176)
    dn1d=zeros(176)
    rrate1d=zeros(176)
    for k=bcf:-1:bzd
        ztrue=zKuL[indx[i]][k+1]+piaB
        dm=0.5*(0.00113070*ztrue^2+0.0047*ztrue+0.4911)
        if dm>3.0
            dm=3.0
        end
        n1,n2=bisection(dmJ,dm)
        dn=(ztrue-zKuJ[n1])/10.0
        attKu=attKuJ[n1]*10^dn

        piaB-=attKu*dr*2
        if piaB<0
            piaB=0
        end
        dm1d[k+1]=dm
        dn1d[k+1]=dn
        rrate1d[k+1]=rrateJ[n1]*10^dn
    end
    for k=bzd-1:-1:btop
        ztrue=(zKuL[indx[i]][k+1]-7)+piaB
        dm=0.5*(0.00113070*ztrue^2+0.0047*ztrue+0.4911)
        if dm>3.0
            dm=3.0
        end
        n1,n2=bisection(dmSJ,dm)
        dn=(ztrue-zKuSJ[n1])/10.0
        attKu=attKuSJ[n1]*10^dn

        piaB-=attKu*dr*2
        if piaB<0
            piaB=0
        end
        dm1d[k+1]=dm
        dn1d[k+1]=dn
        rrate1d[k+1]=rrateSJ[n1]*10^dn
    end
    #println(dm1d)
    return  dm1d, dn1d,piaB, rrate1d
end

function iter_prof(btop,bzd,bcf,bsfc,indx,i,zKuL,zKaL,dr,pia,newTables,ifull)
    zKuJ,attKuJ,pwcJ,rrateJ,dmJ,zKuBBJ,attKuBBJ,pwcBBJ,rrateBBJ,dmBBJ,
    zKuSJ,attKuSJ,pwcSJ,rrateSJ,dmSJ,nbinsj,
    zKaJ,attKaJ,zKaBBJ,attKaBBJ,zKaSJ,attKaSJ,dNw,dNwBB,dNwSJ,
    zKuHJ,attKuHJ,pwcHJ,rrateHJ,dmHJ,
    zKaHJ,attKaHJ,pwcHJ2,rrateHJ2,dmHJ2=newTables
    piaB=pia
    dm1d=zeros(176)
    dn1d=zeros(176)
    rrate1d=zeros(176)
    piaKu=0
    zKuC=copy(zKuL[indx[i]])
    zKuSim=zeros(176,4)
    piaKuV=zeros(4)
    vf=[2.,1,0.8,0.7]
    w=[1,1,2,2]
    dns=zeros(176)
    zKaSim=zeros(176)
    dZKa=zeros(176)
    for it=1:2
        piaKu=0.0
        piaKa=0.0
        for k=btop:bzd-1
            if zKuC[k+1]>10
                ztrue=zKuC[k+1]
                dm=(0.00113070*ztrue^2+0.0047*ztrue+0.4911)
                if dm>3.0
                    dm=3.0
                end
                if ztrue<40
                    ztrueS=ztrue
                    ztrueH=0.0
                else
                    ztrueS=40
                    ztrueH=10*log10(10^(0.1*ztrue)-10^(0.1*39.999))
                end
                n1,n2=bisection(zKuSJ,ztrueS-10*dns[k+1])
                attKuS=attKuSJ[n1]*10^dns[k+1]
                attKaS=attKuSJ[n1]*10^dns[k+1]
                n1H,n2H=bisection(zKuHJ,ztrueH)
                attKuH=attKuHJ[n1H]
                attKaH=attKaHJ[n1H]
                piaKu+=(attKuS+attKuH)*dr
                piaKa+=(attKaS+attKaH)*dr
                zKuC[k+1]=zKuL[indx[i]][k+1]+piaKu
                zKaSim[k+1]=10*log10(10^(0.1*zKaSJ[n1]+dns[k+1])+10^(0.1*zKaHJ[n1H]))-piaKa
                if(zKaL[indx[i]][k+1]>10)
                    dZKa[k+1]=zKaSim[k+1]-zKaL[indx[i]][k+1]
                end
                piaKu+=(attKuS+attKuH)*dr
                piaKa+=(attKaS+attKaH)*dr
                dm1d[k+1]=0.5*(dmSJ[n1]+dmHJ[n1H])
                dn1d[k+1]=0.0
                rrate1d[k+1]=rrateSJ[n1]*10^dns[k+1]+rrateHJ[n1H]
            end
            zKuSim[k+1,:].=zKuC[k+1]
        end
        println("$(it)=",np.mean(dZKa[btop+1:bzd]), " $(rrate1d[bzd]) $(dns[bzd]) $(dZKa[bzd])")
        dnbzd=dns[bzd]
        if it==1
            dns[btop+1:bzd].=dns[btop+1:bzd].-0.6*(dZKa[btop+1:bzd])
            dns[dns.>1].=1
            dns[dns.<-1].=-1
        end
        rate0=rrate1d[bzd]
        rv=0#np.random.randn()
        rv2=0
        piaKuV.=piaKu
        if ifull==1
            for k=bzd:bcf
                if zKuC[k+1]>0
                    ztrue=zKuC[k+1]
                    dm=0.9*(0.00113070*ztrue^2+0.0047*ztrue+0.4911)
                    if dm>3.0
                        dm=3.0
                    end
                    n1,n2=bisection(dmJ,dm)
                    dn=(ztrue-zKuJ[n1])/10.0
                    attKu=attKuJ[n1]*10^dn
                    piaKu+=attKu*dr
                    zKuC[k+1]=zKuL[indx[i]][k+1]+piaKu
                    piaKu+=attKu*dr
                    dm1d[k+1]=dm
                    dn1d[k+1]=dn
                    rrate1d[k+1]=rrateJ[n1]*10^dn
                end
            end
        else
            for k=bzd:bcf
                rv=0.75*rv+0.25*np.random.randn()
                rv2=0.75*rv2+0.25*np.random.randn()
                dn=dnbzd+rv2
                rate0=rate0#*exp(0.5*rv)
                rrate1d[k+1]=rate0*exp(0.5*rv)
                #R=1.370*eps^4.258*dm^5.420

                for isim=1:4
                    dm=(vf[isim]*rrate1d[k+1])^(1/5.420)
                    n1,n2=bisection(dmJ,dm)
                    dn=log10(rrateJ[n1]/(vf[isim]*rrate1d[k+1]))
                    n1,n2=bisection(rrateJ,vf[isim]*rrate1d[k+1]/10^dn)
                    attKu=attKuJ[n1]*10^dn
                    dm1d[k+1]=dmJ[n1]
                    piaKuV[isim]+=attKu*dr
                    zKuSim[k+1,isim]=zKuJ[n1]-piaKuV[isim]+10*dn
                    piaKuV[isim]+=attKu*dr
                end
            end
        end
    end
    w=w./sum(w)
    extMean=np.sum(w.*10 .^(-0.1*piaKuV))
    piaKu=-10*log10(extMean)
    nubf=piaKu/np.mean(piaKuV)
    zKuSim1=zeros(176)
    for isim=1:4
        zKuSim1=zKuSim1.+w[isim]*10 .^(0.1*zKuSim[:,isim])
    end
    zKuSim1=10*log10.(zKuSim1)
    #println(nubf,piaKuV,vf,w)
    return  dm1d, dn1d,piaKu, rrate1d,zKuC,piaKu,zKuSim1
end

function hb_cv(btop,bzd,bcf,bsfc,indx,i,zKuL,zKaL,dnv,dr,newTables)
    zKuJ,attKuJ,pwcJ,rrateJ,dmJ,
    zKuBBJ,attKuBBJ,pwcBBJ,rrateBBJ,dmBBJ,
    zKuSJ,attKuSJ,pwcSJ,rrateSJ,dmSJ,nbinsj,
    zKaJ,attKaJ,zKaBBJ,attKaBBJ,zKaSJ,attKaSJ,dNw,dNwBB,dNwSJ=newTables
    nodes=[0,0,0,0,0]
    nodes[1]=Int32(trunc(btop/2))
    nodes[3]=Int32(trunc(bzd/2))
    nodes[2]=nodes[3]-3
    nodes[4]=nodes[3]+2
    nodes[5]=Int32(trunc(bcf/2))
    if nodes[1]>nodes[2]
        nodes[1]=nodes[2]-3
    end
    z13obs=zeros(88)
    z35obs=zeros(88)
    dm=zeros(88)
    for k=1:88
        z13obs[k]=Cfloat(10.0*log10(1e-9+(10^(0.1*zKuL[indx[i]][2*k])+10^(0.1*zKuL[indx[i]][2*k-1]))/2))
        z35obs[k]=Cfloat(10.0*log10(1e-9+(10^(0.1*zKaL[indx[i]][2*k])+10^(0.1*zKaL[indx[i]][2*k-1]))/2))
    end
    z13c=copy(z13obs)
    iSurf=Int32(trunc(bsfc/2))
    zeta1d=zeros(88)
    beta=0.76
    if z13obs[nodes[5]+1]<0
        z13obs[nodes[5]+1]=0
    end
    piamax=maximum(z13obs[nodes[1]+1:nodes[5]+1])+3-z13obs[nodes[5]+1]
    println(piamax)
    zetaS=0.0
    q=0.2*log(10)
    eps=1
    piaKa=0.0
    zKaSim=zeros(88).-99
# log10(Nw)=0.279*Dm^2-2.1347*Dm+5.8102
# 0.2876 –2.1543 5.7352
#Dm=0.00093071*Zku^2+0.0027*Zku+0.5904
#Dm=0.00113070*Zku^2+0.0047*Zku+0.4911
    for it=1:2
        zetaS=0.0
        piaKa=0.0
        for i=nodes[1]:nodes[5]
            if i<=nodes[2]
                if z13obs[i+1]>0
                    n1,n2=bisection(zKuSJ,z13c[i+1]-10*dnv[i+1])
                    zetaS+=attKuSJ[n1]*10^dnv[i+1]/10^(z13c[i+1]*0.1*beta)*10^(z13obs[i+1]*0.1*beta)*dr
                end
                zeta1d[i+1]=zetaS
            end
            if i<=nodes[3] && i>nodes[2]
                if z13obs[i+1]>0
                    f=(nodes[3]-i)/(nodes[3]-nodes[2])
                    if z13c[i+1]+0.5<z13c[nodes[3]+1] && z13c[i+1]-0.5>z13c[nodes[2]+1]
                        f=(z13c[nodes[3]+1]-z13c[i+1])/(z13c[nodes[3]+1]-z13c[nodes[2]+1])
                    end
                    n1,n2=bisection2(10 .^zKuSJ,10 .^zKuBBJ,f,10^(z13c[i+1]-10*dnv[i+1]))
                    attKu=(f*attKuSJ[n1]+(1-f)*attKuBBJ[n1])*10^dnv[i+1]
                    zetaS+=attKu/10^(z13c[i+1]*0.1*beta)*10^(z13obs[i+1]*0.1*beta)*dr
                end
                zeta1d[i+1]=zetaS
            end
            if i<=nodes[4] && i>nodes[3]
                if z13obs[i+1]>0
                    f=(nodes[4]-i)/(nodes[4]-nodes[3])
                    n1,n2=bisection2(zKuBBJ,zKuJ,f,z13c[i+1]-10*dnv[i+1])
                    attKu=(f*attKuBBJ[n1]+(1-f)*attKuJ[n1])*10^dnv[i+1]
                    zetaS+=attKu/10^(z13c[i+1]*0.1*beta)*10^(z13obs[i+1]*0.1*beta)*dr
                end
                zeta1d[i+1]=zetaS
            end
            if i>nodes[4]
                if z13obs[i+1]>0
                    n1,n2=bisection(zKuJ,z13c[i+1]-10*dnv[i+1])
                    zetaS+=attKuJ[n1]*10^dnv[i+1]/10^(z13c[i+1]*0.1*beta)*10^(z13obs[i+1]*0.1*beta)*dr
                end
                zeta1d[i+1]=zetaS
            end
        end
        eps=min(1,(1-10^(-0.1*beta*piamax))/(q*beta*zetaS))
        for i=nodes[1]:nodes[5]
            z13c[i+1]=z13obs[i+1]-10/beta*log10(1-eps*q*beta*zeta1d[i+1])
        end
    end
    n1,n2=bisection(zKuJ,z13c[nodes[5]+1]-10*dnv[nodes[5]+1])
    println(z13c[nodes[5]+1]-z13obs[nodes[5]+1])
    return z13obs,z13c,eps,dmJ[n1],rrateJ[n1]*10^dnv[nodes[5]+1],
    -10/beta*log10(1-eps*q*beta*zeta1d[nodes[5]+1])
end

function hb_bb(zku,zka,bcf,bst,brs,bzd,bbPeak,pType,newTables,dr,dnv)
    zKuJ,attKuJ,pwcJ,rrateJ,dmJ,
    zKuBBJ,attKuBBJ,pwcBBJ,rrateBBJ,dmBBJ,
    zKuSJ,attKuSJ,pwcSJ,rrateSJ,dmSJ,nbinsj,
    zKaJ,attKaJ,zKaBBJ,attKaBBJ,zKaSJ,attKaSJ,dNw,dNwBB,dNwSJ=newTables
    nodes=[0,0,0,0,0]
    nodes[1]=Int32(trunc(bst[i1,j1]/2))
    nodes[3]=Int32(trunc((bzd[i1,j1])/2))
    if bbPeak[i1,j1]>0
        nodes[3]=Int32(trunc(bbPeak[i1,j1]/2))
    end
    nodes[2]=nodes[3]-3
    nodes[4]=nodes[3]+2
    nodes[5]=Int32(trunc(bcf[i1,j1]/2))
    if nodes[1]>nodes[2]
        nodes[1]=nodes[2]-3
    end
    z13obs=zeros(88)
    z35obs=zeros(88)
    dm=zeros(88)
    for k=1:88
        z13obs[k]=Cfloat(10.0*log10(1e-9+(10^(0.1*zku[i1,j1,2*k])+10^(0.1*zku[i1,j1,2*k-1]))/2))
        z35obs[k]=Cfloat(10.0*log10(1e-9+(10^(0.1*zka[i1,j1-12,2*k])+10^(0.1*zka[i1,j1-12,2*k-1]))/2))
    end
    #println(z35obs)
    z13c=copy(z13obs)
    iSurf=Int32(trunc(brs[i1,j1]/2))
    zeta1d=zeros(88)
    beta=0.76
    if z13obs[nodes[5]+1]<0
        z13obs[nodes[5]+1]=0
    end
    piamax=maximum(z13obs[nodes[1]+1:nodes[5]+1])+3-z13obs[nodes[5]+1]

    zetaS=0.0
    q=0.2*log(10)
    eps=1
    piaKa=0.0
    zKaSim=zeros(88).-99
    for it=1:2
        zetaS=0.0
        piaKa=0.0
        for i=nodes[1]:nodes[5]
            if i<=nodes[2]
                if z13obs[i+1]>0
                    n1,n2=bisection(zKuSJ,z13c[i+1]-10*dnv[i+1])
                    zetaS+=attKuSJ[n1]*10^dnv[i+1]/10^(z13c[i+1]*0.1*beta)*10^(z13obs[i+1]*0.1*beta)*dr
                end
                zeta1d[i+1]=zetaS
            end
            if i<=nodes[3] && i>nodes[2]
                if z13obs[i+1]>0
                    f=(nodes[3]-i)/(nodes[3]-nodes[2])
                    if z13c[i+1]+0.5<z13c[nodes[3]+1] && z13c[i+1]-0.5>z13c[nodes[2]+1]
                        f=(z13c[nodes[3]+1]-z13c[i+1])/(z13c[nodes[3]+1]-z13c[nodes[2]+1])
                    end
                    n1,n2=bisection2(10 .^zKuSJ,10 .^zKuBBJ,f,10^(z13c[i+1]-10*dnv[i+1]))
                    attKu=(f*attKuSJ[n1]+(1-f)*attKuBBJ[n1])*10^dnv[i+1]
                    zetaS+=attKu/10^(z13c[i+1]*0.1*beta)*10^(z13obs[i+1]*0.1*beta)*dr
                end
                zeta1d[i+1]=zetaS
            end
            if i<=nodes[4] && i>nodes[3]
                if z13obs[i+1]>0
                    f=(nodes[4]-i)/(nodes[4]-nodes[3])
                    n1,n2=bisection2(zKuBBJ,zKuJ,f,z13c[i+1]-10*dnv[i+1])
                    attKu=(f*attKuBBJ[n1]+(1-f)*attKuJ[n1])*10^dnv[i+1]
                    zetaS+=attKu/10^(z13c[i+1]*0.1*beta)*10^(z13obs[i+1]*0.1*beta)*dr
                end
                zeta1d[i+1]=zetaS
            end
            if i>nodes[4]
                if z13obs[i+1]>0
                    n1,n2=bisection(zKuJ,z13c[i+1]-10*dnv[i+1])
                    zetaS+=attKuJ[n1]*10^dnv[i+1]/10^(z13c[i+1]*0.1*beta)*10^(z13obs[i+1]*0.1*beta)*dr
                end
                zeta1d[i+1]=zetaS
            end
        end
        eps=min(1,(1-10^(-0.1*beta*piamax))/(q*zetaS))
        for i=nodes[1]:nodes[5]
            z13c[i+1]=z13obs[i+1]-10/beta*log10(1-eps*q*beta*zeta1d[i+1])
        end
    end
    retrate=zeros(88)
end

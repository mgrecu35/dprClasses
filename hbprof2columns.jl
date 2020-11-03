function prof(btop,bzd,bcf,bsfc,zKuL,zKaL,dr,pia,newTables)
    zKuJ,attKuJ,pwcJ,rrateJ,dmJ,zKuBBJ,attKuBBJ,pwcBBJ,rrateBBJ,dmBBJ,
    zKuSJ,attKuSJ,pwcSJ,rrateSJ,dmSJ,nbinsj,
    zKaJ,attKaJ,zKaBBJ,attKaBBJ,zKaSJ,attKaSJ,dNw,dNwBB,dNwSJ,
    zKuHJ,attKuHJ,pwcHJ,rrateHJ,dmHJ,
    zKaHJ,attKaHJ,pwcHJ2,rrateHJ2,dmHJ2=newTables
    npCoeff=[0.012061996352769808,-0.2876799692273608]
    piaB=pia
    dm1d=zeros(176)
    dn1d=zeros(176)
    rrate1d=zeros(176)
    rrate1ds=zeros(176)
    piaKu=0
    zKuC=copy(zKuL)
    zKuSim=zeros(176,4)
    piaKuV=zeros(4)
    dmCoeff=[0.027773772993636318,-0.6624076086959071]
    dns=zeros(176)
    zKaSim=zeros(176)
    dZKa=zeros(176)
    for it=1:2
        piaKu=0.0
        piaKa=0.0
        for k=btop:bzd-1
            if zKuC[k+1]>10
                if k<bzd
                    ztrue=zKuC[k+1]
                else
                    f=(bzd-k)/6.0
                    ztrue=10*log10(10^(0.1*zKuC[k+1])*f+1)
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
                zKuC[k+1]=zKuL[k+1]+piaKu
                zKaSim[k+1]=10*log10(10^(0.1*zKaSJ[n1]+dns[k+1])+10^(0.1*zKaHJ[n1H]))-piaKa
                if(zKaL[k+1]>10)
                    dZKa[k+1]=zKaSim[k+1]-zKaL[k+1]
                end
                piaKu+=(attKuS+attKuH)*dr
                piaKa+=(attKaS+attKaH)*dr
                dm1d[k+1]=0.5*(dmSJ[n1]+dmHJ[n1H])
                dn1d[k+1]=0.0
                #rrate1d[k+1]=rrateSJ[n1]*10^dns[k+1]+rrateHJ[n1H]
                rrate1ds[k+1]=rrateSJ[n1]*10^dns[k+1]+rrateHJ[n1H]
            end
            zKuSim[k+1,:].=zKuC[k+1]
        end
        if it==1
            dns[btop+1:bzd].=dns[btop+1:bzd].-0.6*(dZKa[btop+1:bzd])
            dns[dns.>0.75].=0.75
            dns[dns.<-1.5].=-1.5
        end
    end
    piaKuS=piaKu+0.0
    piaKaS=piaKa+0.0
    zeta1d=zeros(176)
    beta=0.76

    piamax=maximum(zKuL[bzd+1:bcf+1])-zKuL[bcf+1]+5

    q=0.2*log(10)
    zKuC[bzd+1:bcf+1].=zKuL[bzd+1:bcf+1].+piaKuS
    eps=1.0
    epst=eps
    attKu=0.0
    attKa=0.0
    for it=1:3
        zetaS=0.0
        piaKa=piaKaS+0.0
        for k=bzd:bcf
            if zKuC[k+1]>0
                ztrue=zKuC[k+1]
                if k<bzd+6
                    ftran=0.85+1/6.0*(k-bzd)
                else
                    ftran=1.0
                end
                ftran=1
                dm=ftran*1.05*exp(dmCoeff[1]*ztrue+dmCoeff[2])/epst^0.25
                if dm>3.0
                    dm=3.0
                end
                if dm<0.25
                    dm=0.25
                end
                n1,n2=bisection(dmJ,dm)
                dn=(ztrue-zKuJ[n1])/10.0
                attKu=attKuJ[n1]*10^dn
                attKa=attKaJ[n1]*10^dn
                zetaS+=attKu/10^(zKuC[k+1]*0.1*beta)*10^(zKuL[k+1]*0.1*beta)*dr
                rrate1d[k+1]=rrateJ[n1]*10^dn
                piaKa+=attKa*2*dr
            end
            zeta1d[k+1]=zetaS
        end
        eps=min(1,(1-10^(-0.1*beta*piamax))/(q*beta*zetaS))
        epst*=eps
        for k=bzd:bcf
            zKuC[k+1]=zKuL[k+1]+piaKuS-10/beta*log10(1-eps*q*beta*zeta1d[k+1])
        end
    end
    piaKu=(-10/beta*log10(1-eps*q*beta*zeta1d[bcf+1]))+attKu*(bsfc-bcf)*2*dr
    piaKa+=attKa*(bsfc-bcf)*2*dr
    println(eps,",$(epst), $(piamax),$(-10/beta*log10(1-eps*q*beta*zeta1d[bcf+1]))")
    println("$(piaKu)")
    return  dm1d, dn1d,piaKu, rrate1d,rrate1ds,zKuC,piaKu,piaKa

end

function profPIA(btop,bzd,bcf,bsfc,zKuL,zKaL,dr,pia,newTables)
    zKuJ,attKuJ,pwcJ,rrateJ,dmJ,zKuBBJ,attKuBBJ,pwcBBJ,rrateBBJ,dmBBJ,
    zKuSJ,attKuSJ,pwcSJ,rrateSJ,dmSJ,nbinsj,
    zKaJ,attKaJ,zKaBBJ,attKaBBJ,zKaSJ,attKaSJ,dNw,dNwBB,dNwSJ,
    zKuHJ,attKuHJ,pwcHJ,rrateHJ,dmHJ,
    zKaHJ,attKaHJ,pwcHJ2,rrateHJ2,dmHJ2=newTables
    npCoeff=[0.012061996352769808,-0.2876799692273608]
    piaB=pia
    dm1d=zeros(176)
    dn1d=zeros(176)
    rrate1d=zeros(176)
    rrate1ds=zeros(176)
    piaKu=0
    zKuC=copy(zKuL)
    zKuSim=zeros(176,4)
    piaKuV=zeros(4)
    dmCoeff=[0.027773772993636318,-0.6624076086959071]
    dns=zeros(176)
    zKaSim=zeros(176)
    dZKa=zeros(176)
    for it=1:2
        piaKu=0.0
        piaKa=0.0
        for k=btop:bzd-1
            if zKuC[k+1]>10
                if k<bzd
                    ztrue=zKuC[k+1]
                else
                    f=(bzd-k)/6.0
                    ztrue=10*log10(10^(0.1*zKuC[k+1])*f+1)
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
                zKuC[k+1]=zKuL[k+1]+piaKu
                zKaSim[k+1]=10*log10(10^(0.1*zKaSJ[n1]+dns[k+1])+10^(0.1*zKaHJ[n1H]))-piaKa
                if(zKaL[k+1]>10)
                    dZKa[k+1]=zKaSim[k+1]-zKaL[k+1]
                end
                piaKu+=(attKuS+attKuH)*dr
                piaKa+=(attKaS+attKaH)*dr
                dm1d[k+1]=0.5*(dmSJ[n1]+dmHJ[n1H])
                dn1d[k+1]=0.0
                #rrate1d[k+1]=rrateSJ[n1]*10^dns[k+1]+rrateHJ[n1H]
                rrate1ds[k+1]=rrateSJ[n1]*10^dns[k+1]+rrateHJ[n1H]
            end
            zKuSim[k+1,:].=zKuC[k+1]
        end
        if it==1
            dns[btop+1:bzd].=dns[btop+1:bzd].-0.6*(dZKa[btop+1:bzd])
            dns[dns.>0.75].=0.75
            dns[dns.<-1.5].=-1.5
        end
    end
    piaKuS=piaKu+0.0
    piaKaS=piaKa+0.0
    zeta1d=zeros(176)
    beta=0.76

    #piamax=maximum(zKuL[bzd+1:bcf+1])-zKuL[bcf+1]+5
    piamax=pia
    q=0.2*log(10)
    zKuC[bzd+1:bcf+1].=zKuL[bzd+1:bcf+1].+piaKuS
    eps=1.0
    epst=eps
    attKu=0.0
    attKa=0.0
    for it=1:3
        zetaS=0.0
        piaKa=piaKaS+0.0
        for k=bzd:bcf
            if zKuC[k+1]>0
                ztrue=zKuC[k+1]
                if k<bzd+6
                    ftran=0.85+1/6.0*(k-bzd)
                else
                    ftran=1.0
                end
                ftran=1
                dm=ftran*1.05*exp(dmCoeff[1]*ztrue+dmCoeff[2])/epst^0.25
                if dm>3.0
                    dm=3.0
                end
                if dm<0.25
                    dm=0.25
                end
                n1,n2=bisection(dmJ,dm)
                dn=(ztrue-zKuJ[n1])/10.0
                attKu=attKuJ[n1]*10^dn
                attKa=attKaJ[n1]*10^dn
                zetaS+=attKu/10^(zKuC[k+1]*0.1*beta)*10^(zKuL[k+1]*0.1*beta)*dr
                rrate1d[k+1]=rrateJ[n1]*10^dn
                piaKa+=attKa*2*dr
            end
            zeta1d[k+1]=zetaS
        end
        eps=(1-10^(-0.1*beta*piamax))/(q*beta*zetaS)
        epst*=eps
        for k=bzd:bcf
            zKuC[k+1]=zKuL[k+1]+piaKuS-10/beta*log10(1-eps*q*beta*zeta1d[k+1])
        end
    end
    piaKu=(-10/beta*log10(1-eps*q*beta*zeta1d[bcf+1]))+attKu*(bsfc-bcf)*2*dr
    piaKa+=attKa*(bsfc-bcf)*2*dr
    #println(eps,",$(epst), $(piamax),$(-10/beta*log10(1-eps*q*beta*zeta1d[bcf+1]))")
    #println("$(piaKu)")
    return  dm1d, dn1d,piaKu, rrate1d,rrate1ds,zKuC,piaKu,piaKa

end

function iter_prof(btop,bzd,bcf,bsfc,zKuL,zKaL,dr,pia,newTables)
    zKuJ,attKuJ,pwcJ,rrateJ,dmJ,zKuBBJ,attKuBBJ,pwcBBJ,rrateBBJ,dmBBJ,
    zKuSJ,attKuSJ,pwcSJ,rrateSJ,dmSJ,nbinsj,
    zKaJ,attKaJ,zKaBBJ,attKaBBJ,zKaSJ,attKaSJ,dNw,dNwBB,dNwSJ,
    zKuHJ,attKuHJ,pwcHJ,rrateHJ,dmHJ,
    zKaHJ,attKaHJ,pwcHJ2,rrateHJ2,dmHJ2=newTables
    npCoeff=[0.012061996352769808,-0.2876799692273608]
    piaB=pia
    dm1d=zeros(176)
    dn1d=zeros(176)
    rrate1d=zeros(176)
    piaKu=0
    zKuC=copy(zKuL)
    zKuSim=zeros(176,4)
    piaKuV=zeros(4)
    dmCoeff=[0.027773772993636318,-0.6624076086959071]
    dns=zeros(176)
    zKaSim=zeros(176)
    dZKa=zeros(176)
    for it=1:6
        piaKu=0.0
        piaKa=0.0
        for k=btop:bzd-1
            if zKuC[k+1]>10
                ztrue=zKuC[k+1]
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
                zKuC[k+1]=zKuL[k+1]+piaKu
                zKaSim[k+1]=10*log10(10^(0.1*zKaSJ[n1]+dns[k+1])+10^(0.1*zKaHJ[n1H]))-piaKa
                if(zKaL[k+1]>10)
                    dZKa[k+1]=zKaSim[k+1]-zKaL[k+1]
                end
                piaKu+=(attKuS+attKuH)*dr
                piaKa+=(attKaS+attKaH)*dr
                dm1d[k+1]=0.5*(dmSJ[n1]+dmHJ[n1H])
                dn1d[k+1]=0.0
                rrate1d[k+1]=rrateSJ[n1]*10^dns[k+1]+rrateHJ[n1H]
            end
            zKuSim[k+1,:].=zKuC[k+1]
        end
        if it==1
            dns[btop+1:bzd].=dns[btop+1:bzd].-0.6*(dZKa[btop+1:bzd])
            dns[dns.>1].=1
            dns[dns.<-1].=-1
        end
        for k=bzd:bcf
            if zKuC[k+1]>0
                ztrue=zKuC[k+1]
                dm=0.8*exp(dmCoeff[1]*ztrue+dmCoeff[2])
                if dm>3.0
                    dm=3.0
                end
                if dm<0.25
                    dm=0.25
                end
                n1,n2=bisection(dmJ,dm)
                dn=(ztrue-zKuJ[n1])/10.0
                attKu=attKuJ[n1]*10^dn
                piaKu+=attKu*dr
                zKuC[k+1]=zKuL[k+1]+piaKu
                piaKu+=attKu*dr
                dm1d[k+1]=dm
                dn1d[k+1]=dn
                rrate1d[k+1]=rrateJ[n1]*10^dn
            end
        end
    end

    return  dm1d, dn1d,piaKu, rrate1d,zKuC

end


function profw(btop,bzd,bcf,bsfc,zKuL,dr,pia,newTables)
    zKuJ,attKuJ,pwcJ,rrateJ,dmJ,zKuBBJ,attKuBBJ,pwcBBJ,rrateBBJ,dmBBJ,
    zKuSJ,attKuSJ,pwcSJ,rrateSJ,dmSJ,nbinsj,
    zKaJ,attKaJ,zKaBBJ,attKaBBJ,zKaSJ,attKaSJ,dNw,dNwBB,dNwSJ,
    zKuHJ,attKuHJ,pwcHJ,rrateHJ,dmHJ,
    zKaHJ,attKaHJ,pwcHJ2,rrateHJ2,dmHJ2=newTables
    npCoeff=[0.012061996352769808,-0.2876799692273608]
    dmCoeff=[0.027773772993636318,-0.6624076086959071]
    piaB=pia
    dm1d=zeros(176)
    dn1d=zeros(176)
    rrate1d=zeros(176)
    piaKu=0
    zKuC=copy(zKuL)
    zKuSim=zeros(176,4)
    piaKuV=zeros(4)
    vf=[2.,1,0.8,0.7]
    w=[1,1,2,2]
    dns=zeros(176)
    zKaSim=zeros(176)
    dZKa=zeros(176)
    eps=exp(np.random.randn()*0.5)
    #eps=1.0
    piaBZD=0.0
    piaBackw=zeros(4)
    dnsub=zeros(4,176)
    dms=zeros(4,176)

    piaKu=0.0
    piaKa=0.0
    zeta1d=zeros(176)
    beta=0.76

    piamax=maximum(zKuL[bzd+1:bcf+1])+3
    zetaS=0.0
    q=0.2*log(10)
    for it=1:2
        for k=bzd:bcf
            if zKuC[k+1]>0
                ztrue=zKuC[k+1]
                dm=0.99*exp(dmCoeff[1]*ztrue+dmCoeff[2])
                if dm>3.0
                    dm=3.0
                end
                if dm<0.25
                    dm=0.25
                end
                n1,n2=bisection(dmJ,dm)
                dn=(ztrue-zKuJ[n1])/10.0
                attKu=attKuJ[n1]*10^dn
                zetaS+=attKu/10^(zKuC[k+1]*0.1*beta)*10^(zKuL[k+1]*0.1*beta)*dr
                rrate1d[k+1]=rrateJ[n1]*10^dn
            end
            zeta1d[k+1]=zetaS
        end
        eps=min(1,(1-10^(-0.1*beta*piamax))/(q*beta*zetaS))
        for k=bzd:bcf
            zKuC[k+1]=zKuL[k+1]-10/beta*log10(1-eps*q*beta*zeta1d[k+1])
        end
    end
    return  dm1d, dn1d,piaKu, rrate1d,zKuC
end

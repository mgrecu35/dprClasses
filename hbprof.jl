function hb_bb(zku,zka,bcf,bst,brs,bzd,bbPeak,pType,newTables,dr,dnv)
    zKuJ,attKuJ,pwcJ,rrateJ,dmJ,
    zKuBBJ,attKuBBJ,pwcBBJ,rrateBBJ,dmBBJ,
    zKuSJ,attKuSJ,pwcSJ,rrateSJ,dmSJ,nbinsj,
    zKaJ,attKaJ,zKaBBJ,attKaBBJ,zKaSJ,attKaSJ,dNw,dNwBB,dNwSJ=newTables
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

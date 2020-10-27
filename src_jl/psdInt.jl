#using Dierckx
using Interpolations
using SpecialFunctions

function nw_lambd(swc,nc,mu,bscat,ext,scat,g,vfall,Deq,wl)
    rhow=1e6
    lambd=(nc*rhow*pi*gamma(4+mu)/gamma(1+mu)/6.0/swc)^(0.333)  # m-1
    n0=nc*lambd/gamma(1+mu) # m-4
    n0*=1e-3 # mm-1 m-3
    lambd*=1e-2 # cm-1
    #print(swc,nc,lambd)
    W=0
    dD=0.05
    rhow=1 #gcm-3
    #Dint=arange(160)*dD+dD/2.0
    #bscatInt=interp(Dint,Deq,bscat)
    #extInt=exp(interp(Dint,Deq,log(ext)))  #m^2
    #vfallInt=interp(Dint,Deq,vfall) Dint=0:159 Dint=Dint.*dD
    Dint=Dint.+dD/2.0
    #spl1=Spline1(Deq,log.(bscat)) #m^2
    spl1=interpolate((Deq,),log.(bscat),Gridded(Linear())) #m^2
    bscatInt=exp.(spl1(Dint))
    #spl2=Spline1D(Deq,log.(ext)) #m^2
    spl2=interpolate((Deq,),log.(ext),Gridded(Linear())) #m^2
    extInt=exp.(spl2(Dint))
    spl5=interpolate((Deq,),vfall,Gridded(Linear())) #m/s
    vfallInt=spl5(Dint)
    #spl3=Spline1D(Deq,log.(scat))  #m^2
    #scatInt=exp.(spl3(Dint))
    #spl6=Spline1D(Deq_r,vfall_r)  #m/s
    #vfallInt_r=spl6(Dint)
    #spl4=Spline1D(Deq,g)  #m^2
    #gInt=spl4(Dint)

    Z=0.0
    att=0.
    fact=1e6/pi^5/0.93*wl^4
    nc0=0
    Vol=0
    vdop=0
    m3=0
    m4=0
    for i=0:159
        d=dD*i+dD/2
        Nd=exp(-lambd*d*0.1)*(lambd*0.1*d)^mu*dD #(mm)
        W+=n0*Nd*(0.1*d)^3*pi/6*rhow #(g/m3)
        m4=n0*Nd*(0.1*d)^4
        m3=n0*Nd*(0.1*d)^3
        Z+=n0*Nd*bscatInt[i+1]
        vdop+=n0*Nd*bscatInt[i+1]*vfallInt[i+1]
        att+=n0*Nd*extInt[i+1]*1e3 #(/km)
        nc0+=n0*Nd
        Vol+=n0*Nd*(1e-3*d)^3*pi/6
    #print(nc0,nc)
    end
    return W, log10(Z*fact)*10., att, log10(n0/0.08e5), nc0, vdop/Z, m4/m3
end

function nw_lambd_1M(swc,dn,mu,bscat,ext,scat,g,vfall,Deq,wl)
    rhow=1e6
    #print(dn)
    #n0=nc*lambd/gammama(1+mu) # m-4
    n0=8e6*dn #m-1 m-3\n",
    lambd=(n0*rhow*pi*gamma(4+mu)/6.0/swc)^0.25
    n0*=1e-3 # mm-1 m-3
    lambd*=1e-2 # cm-1
    #print(swc,nc,lambd)
    W=0
    dD=0.05
    rhow=1 #gcm-3
    Dint=0:159
    Dint=Dint.*dD
    Dint=Dint.+dD/2.0
    spl1=interpolate((Deq,),log.(bscat),Gridded(Linear())) #m^2
    bscatInt=exp.(spl1(Dint))
    spl2=interpolate((Deq,),log.(ext),Gridded(Linear()))  #m^2
    extInt=exp.(spl2(Dint))
    spl3=interpolate((Deq,),log.(scat),Gridded(Linear())) #m^2
    scatInt=exp.(spl3(Dint))
    spl4=interpolate((Deq,),g,Gridded(Linear())) #m^2
    gInt=spl4(Dint)
    spl5=interpolate((Deq,),vfall,Gridded(Linear()))  #m/s
    vfallInt=spl5(Dint)

    #Dint=arange(160)*dD+dD/2.0
    #bscatInt=interp(Dint,Deq,bscat)
    #extInt=exp(interp(Dint,Deq,log(ext)))  #m^2
    #vfallInt=interp(Dint,Deq,vfall)
    Z=0.0
    att=0.
    scattPSD=0.
    gPSD=0.
    fact=1e6/pi^5/0.93*wl^4
    nc0=0
    Vol=0
    vdop=0
    #println(size(scatInt))
    #exit()
    m4=0
    m3=0
    rr=0.0
    for i=0:159
        d=dD*i+dD/2

        Nd=exp(-lambd*d*0.1)*(lambd*0.1*d)^mu*dD #(mm)
        W+=n0*Nd*(0.1*d)^3*pi/6*rhow #(g/m3)
        Z+=n0*Nd*bscatInt[i+1]
        vdop+=n0*Nd*bscatInt[i+1]*vfallInt[i+1]
        m4+=n0*Nd*(0.1*d)^4
        m3+=n0*Nd*(0.1*d)^3
        att+=n0*Nd*extInt[i+1]*1e3 #(/km)
        scattPSD+=n0*Nd*scatInt[i+1]*1e3 #(/km)
        gPSD+=n0*Nd*gInt[i+1]*scatInt[i+1]*1e3 #(/km)
        nc0+=n0*Nd
        Vol+=n0*Nd*(1e-3*d)^3*pi/6
        #print(nc0,nc)
        rr+=n0*Nd*(0.1*d)^3*pi/6*vfallInt[i+1]*rhow
    end
    if Z>0
        Zn=log10(Z*fact)*10.
    else
        Zn=-99
    end
    return W, Zn, att, log10(n0/0.08e8), nc0, vdop/Z, scattPSD, gPSD/scattPSD, m4/m3, rr
end

function nw_lambd_1Ms(swc,dn,mu,bscat,ext,scat,g,vfall,vfall_r,Deq,Deq_r,wl)
    rhow=1e6
    #print(dn)
    #n0=nc*lambd/gammama(1+mu) # m-4
    n0=8e6*dn #m-1 m-3\n",
    lambd=(n0*rhow*pi*gamma(4+mu)/6.0/swc)^0.25
    n0*=1e-3 # mm-1 m-3
    lambd*=1e-2 # cm-1
    #print(swc,nc,lambd)
    W=0
    dD=0.05
    rhow=1 #gcm-3
    Dint=0:159
    Dint=Dint.*dD
    Dint=Dint.+dD/2.0
    spl1=interpolate((Deq,),log.(bscat),Gridded(Linear())) #m^2
    bscatInt=exp.(spl1(Dint))
    spl2=interpolate((Deq,),log.(ext),Gridded(Linear()))  #m^2
    extInt=exp.(spl2(Dint))
    spl3=interpolate((Deq,),log.(scat),Gridded(Linear())) #m^2
    scatInt=exp.(spl3(Dint))
    spl4=interpolate((Deq,),g,Gridded(Linear())) #m^2
    gInt=spl4(Dint)
    spl5=interpolate((Deq,),vfall,Gridded(Linear()))  #m/s
    vfallInt=spl5(Dint)
    spl5_r=interpolate((Deq_r,),vfall_r,Gridded(Linear()))  #m/s
    vfallInt_r=spl5_r(Dint)

    #Dint=arange(160)*dD+dD/2.0
    #bscatInt=interp(Dint,Deq,bscat)
    #extInt=exp(interp(Dint,Deq,log(ext)))  #m^2
    #vfallInt=interp(Dint,Deq,vfall)
    Z=0.0
    att=0.
    scattPSD=0.
    gPSD=0.
    fact=1e6/pi^5/0.93*wl^4
    nc0=0
    Vol=0
    vdop=0
    #println(size(scatInt))
    #exit()
    m4=0
    m3=0
    rr=0
    for i=0:159
        d=dD*i+dD/2

        Nd=exp(-lambd*d*0.1)*(lambd*0.1*d)^mu*dD*vfallInt[i+1]/vfallInt_r[i+1] #(mm)
        W+=n0*Nd*(0.1*d)^3*pi/6*rhow #(g/m3)
        rr+=n0*Nd*(0.1*d)^3*pi/6*rhow*vfallInt[i+1] #(g/m3*m/s)
        Z+=n0*Nd*bscatInt[i+1]
        vdop+=n0*Nd*bscatInt[i+1]*vfallInt[i+1]
        m4+=n0*Nd*(0.1*d)^4
        m3+=n0*Nd*(0.1*d)^3
        att+=n0*Nd*extInt[i+1]*1e3 #(/km)
        scattPSD+=n0*Nd*scatInt[i+1]*1e3 #(/km)
        gPSD+=n0*Nd*gInt[i+1]*scatInt[i+1]*1e3 #(/km)
        nc0+=n0*Nd
        Vol+=n0*Nd*(1e-3*d)^3*pi/6
        #print(nc0,nc)
    end
    if Z>0
        Zn=log10(Z*fact)*10.
    else
        Zn=-99
    end
    return W, Zn, att, log10(n0/0.08e8), nc0, vdop/Z, scattPSD, gPSD/scattPSD, m4/m3, rr
end
function get_fZr(zObs,bscat_r,scat_r,ext_r,g_r,Deq_r,vfall_r,wl,dn,mu)
    rwc=0.1
    for it=1:6
        w,Z,att,dn1,nc,vdop,
        scatInt,gInt,dm,rr=nw_lambd_1M(rwc,dn,mu,bscat_r[9,:],ext_r[9,:],
                                        scat_r[9,:],g_r[9,:],
                                        vfall_r[:],Deq_r[:],wl)

        rwc1=rwc+0.005

        w1,Z1,att1,dn11,nc1,vdop1,
        scatInt,gInt,dm,rr=nw_lambd_1M(rwc1,dn,mu,bscat_r[9,:],
                                 ext_r[9,:],
                                 scat_r[9,:],g_r[9,:],
                                 vfall_r[:],Deq_r[:],wl)

        gradZ=(Z1-Z)/(log(rwc1/rwc))

        #println(rwc," Z= ",Z, " Z1= ",Z1," ",zObs, " dn=",dn)

        #println("$w $(w1)")
        rwc=rwc*exp(0.75*(zObs-Z)*gradZ/(gradZ*gradZ+0.1))
        if rwc<0.00001
            rwc=0.00001
        end

    end
    #exit()
    w,Z,att,dn1,nc,vdop,
    scatInt,gInt,dm,rr=nw_lambd_1M(rwc,dn,mu,bscat_r[9,:],ext_r[9,:],
                             scat_r[9,:],g_r[9,:],
                             vfall_r[:],Deq_r[:],wl)

    return rwc, Z, att,scatInt,gInt, vdop, dm, rr

end



function get_fZm(f,zObs,ns,bscat_r,scat_r,ext_r,g_r,Deq_r,vfall_r,
                bscat,scat,ext,g,Deq,vfall,wl,dn,mu)
    rwc=0.1
    for it=1:6
        wr,Zr,attr,dnr1,ncr,vdopr,
        scatIntr,gIntr,dm,rr=nw_lambd_1M(f*rwc,dn,mu,bscat_r[9,:],ext_r[9,:],
                                   scat_r[9,:],g_r[9,:],
                                   vfall_r[:],Deq_r[:],wl)

        ws,Zs,atts,dns1,ncs,vdops,
        scatInts,gInts,dm,rr=nw_lambd_1M((1-f)*rwc,dn,mu,bscat[end,ns,:],ext[end,ns,:],
                                   scat[end,ns,:],g[end,ns,:],
                                   vfall[ns,:],Deq[ns,:],wl)
        rwc1=rwc+0.01
        Z=10*log10(10^(0.1*Zr)+10^(0.1*Zs))

        wr1,Zr1,attr1,dnr11,ncr1,vdopr1,
        scatIntr,gIntr,dm,rr=nw_lambd_1M(f*rwc1,dn,mu,bscat_r[9,:],
                                 ext_r[9,:],
                                 scat_r[9,:],g_r[9,:],
                                   vfall_r[:],Deq_r[:],wl)
        ws1,Zs1,atts1,dns11,ncs1,vdops1,
        scatInts1,gInts1,dm,rr=nw_lambd_1M((1-f)*rwc1,dn,mu,bscat[end,ns,:],ext[end,ns,:],
                                   scat[end,ns,:],g[end,ns,:],
                                   vfall[ns,:],Deq[ns,:],wl)

        Z1=10*log10(10^(0.1*Zr1)+10^(0.1*Zs1))
        gradZ=(Z1-Z)/(log(rwc1/rwc))


        #println(rwc," ",Z, " ",zObs)
        rwc=rwc*exp((zObs-Z)*gradZ/(gradZ*gradZ+0.1))
        if rwc<0.00001
            rwc=0.00001
        end

    end

    wr,Zr,attr,dn1r,ncr,vdopr,
    scatIntr,gIntr,dm,rr=nw_lambd_1M(f*rwc,dn,mu,bscat_r[9,:],ext_r[9,:],
                             scat_r[9,:],g_r[9,:],
                               vfall_r[:],Deq_r[:],wl)

    ws,Zs,atts,dns1,ncs,vdops,
    scatInts,gInts,dm,rr=nw_lambd_1M((1-f)*rwc,dn,mu,bscat[end,ns,:],ext[end,ns,:],
                               scat[end,ns,:],g[end,ns,:],
                               vfall[ns,:],Deq[ns,:],wl)

    Z=10*log10(10^(0.1*Zr)+10^(0.1*Zs))
    att=attr+atts
    scatInt=scatIntr+scatInts
    gInt=(gIntr*scatIntr+gInts*scatInts)/scatInt
    return rwc, Z, att,scatInt,gInt

end


function get_fZs(zObs,ns,bscat,scat,ext,g,Deq,Deq_r,vfall,vfall_r,wl,dn,mu)
    rwc=0.1
    #println(ns)
    #println(wl)
    #println(dn)
    #println(mu)
    for it=1:6
        w,Z,att,dn1,nc,vdop,
        scatInt,gInt,dm,rr=nw_lambd_1Ms(rwc,dn,mu,bscat[end,ns,:],ext[end,ns,:],
                                        scat[end,ns,:],g[end,ns,:],
                                        vfall[ns,:],vfall_r,Deq[ns,:],Deq_r,wl)

        rwc1=rwc+0.01

        w1,Z1,att1,dn11,nc1,vdop1,
        scatInt,gInt,dm,rr=nw_lambd_1Ms(rwc1,dn,mu,bscat[end,ns,:],
                                              ext[end,ns,:],
                                              scat[end,ns,:],g[end,ns,:],
                                              vfall[ns,:],vfall_r,Deq[ns,:],Deq_r,wl)

        gradZ=(Z1-Z)/(log(rwc1/rwc))

        #println(rwc," Z= ",Z, " Z1= ",Z1," ",zObs, " dn=",dn)
        #println("snow ",rwc," ",Z, " ",zObs)
        rwc=rwc*exp((zObs-Z)*gradZ/(gradZ*gradZ+0.1))
        if rwc<0.00001
            rwc=0.00001
        end

    end

    w,Z,att,dn1,nc,vdop,
    scatInt,gInt,dm,rr=nw_lambd_1Ms(rwc,dn,mu,bscat[end,ns,:],
                             ext[end,ns,:],
                             scat[end,ns,:],g[end,ns,:],
                             vfall[ns,:],vfall_r,Deq[ns,:],Deq_r,wl)

    #println(Z)
    #println(att)
    return w, Z, att, scatInt,gInt,vdop,dm,rr

end



function pwcFromZHBv(piaR,zObs,dn1d,bscatKu,scatKu,extKu,gKu,DeqKu,
                     vfallKu,bscatKu_r,scatKu_r,extKu_r,gKu_r,
                     DeqKu_r,vfallKu_r,wlKu,hint,dr,
                     n1,n2,ns,mu)
    pia=0.0
    dn=-99.9
    dm=-99.9

    beta=0.72
    nz=size(zObs)[1]
    zeta1d=zeros(nz)

    kext1d=zeros(nz)
    scatt1d=zeros(nz)
    vdop1d=zeros(nz)
    g1d=zeros(nz)

    zT=copy(zObs)
    q=0.2*log(10)
    dn1=1.0
    eps=1.0
    pwc=zeros(nz)
    piaHB=0.0
    for it=1:3
        zeta=0.0
        for i=1:nz
            dn=dn1d[i]*dn1
            Z=0.0
            sca1s=0.0
            kext1s=0.0
            g1s=0.0
            sca1r=0.0
            kext1r=0.0
            g1r=0.0
            if hint[i]>4.0
                swc,Z,att,sca1s,g1s=get_Zs(zT[i],ns,bscatKu,scatKu,extKu,gKu,
                                 DeqKu,
                                 vfallKu,wlKu,dn,mu)
                zeta=zeta+att*4.343/10.0^(0.1*zT[i]*beta)*10.0^(0.1*zObs[i]*beta)*dr
                kext1s=att
                zeta1d[i]=zeta
                pwc[i]=swc
            else

                rwc,Z,att,sca1r,g1r,vdop=get_Zr(zT[i],bscatKu_r,scatKu_r,
                                 extKu_r,gKu_r,DeqKu_r,
                                 vfallKu_r,wlKu,dn,mu)
                zeta=zeta+att*4.343/10.0^(0.1*zT[i]*beta)*10.0^(0.1*zObs[i]*beta)*dr
                zeta1d[i]=zeta
                pwc[i]=rwc
                kext1r=att
                vdop1d[i]=vdop
            end
            kext1d[i]=kext1s+kext1r
            scatt1d[i]=sca1s+sca1r
            g1d[i]=(sca1s*g1s+sca1r*g1r)/(scatt1d[i]+1e-10)
            zT[i]=Z
            zeta1d[i]=zeta
        end
        eps=(1-10^(-0.1*piaR*beta))/(q*beta*zeta1d[nz])
        dn1=dn1*eps^(1.0/(1-beta))
        if q*beta*zeta1d[nz]!= q*beta*zeta1d[nz]
            println("dn1=",dn1)
            println(it, " ", piaR)
            println("eps=",eps)
            println("zeta=",zeta1d[:])
            exit(1)
        end
        #println(q*beta*zeta1d[nz])
        for i=1:nz
            zT[i]=zObs[i]-10/beta*log10(1-eps*q*beta*zeta1d[i])
        end
        piaHB=-10/beta*log10(1-eps*q*beta*zeta1d[nz])
    end

    return piaHB,zT,dn1,pwc, kext1d, scatt1d, g1d
end


function pwcFromZHBmixed(piaR,zObs,dn1d,bscatKu,scatKu,extKu,gKu,DeqKu,
                         vfallKu,bscatKu_r,scatKu_r,extKu_r,gKu_r,
                         DeqKu_r,vfallKu_r,wlKu,hint,dr,
                         ns,mu,nmTop,nmBot)
    pia=0.0
    dn=-99.9
    dm=-99.9

    beta=0.72
    nz=size(zObs)[1]
    zeta1d=zeros(nz)

    kext1d=zeros(nz)
    scatt1d=zeros(nz)
    g1d=zeros(nz)
    vdop1d=zeros(nz)

    zT=copy(zObs)
    q=0.2*log(10)
    dn1=1.0
    eps=1.0
    pwc=zeros(nz)
    piaHB=0.0
    for it=1:1
        zeta=0.0
        for i=1:nz
            if(zT[i]>-10 && zObs[i]>-10)
                dn=dn1d[i]*dn1
                #println("dn= ",dn)
                Z=0.0
                sca1s=0.0
                kext1s=0.0
                g1s=0.0
                sca1r=0.0
                kext1r=0.0
                g1r=0.0
                if i<=nmTop
                    swc,Z,att,sca1s,g1s=get_Zs(zT[i],ns,bscatKu,scatKu,extKu,gKu,
                                               DeqKu,
                                               vfallKu,wlKu,dn,mu)
                    if(swc>10)
                        println(zT[i], " ", zObs[i], " ", dn, " ", dn1," ",swc, "it=",it)
                        #exit()
                    end
                    zeta=zeta+att*4.343/10.0^(0.1*zT[i]*beta)*10.0^(0.1*zObs[i]*beta)*dr
                    kext1s=att
                    zeta1d[i]=zeta
                    pwc[i]=swc
                else
                    if i>=nmBot
                        rwc,Z,att,sca1r,g1r,vdopr=get_Zr(zT[i],bscatKu_r,scatKu_r,
                                                         extKu_r,gKu_r,DeqKu_r,
                                                         vfallKu_r,wlKu,dn,mu)
                        zeta=zeta+att*4.343/10.0^(0.1*zT[i]*beta)*10.0^(0.1*zObs[i]*beta)*dr
                        zeta1d[i]=zeta
                        pwc[i]=rwc
                        kext1r=att
                        vdop1d[i]=vdopr
                        if q*beta*zeta1d[i]<0.9999
                            atten=-10/beta*log10(1-q*beta*zeta1d[i])
                        else
                            atten=-10/beta*log10(1-0.9999)
                        end
                        println("Z($i)=$(zT[i]) dn=$(dn), $(zeta1d[i]) $(atten)")
                    else
                        f=(i-nmTop)/(nmBot-nmTop)
                        # get_Zm(f,zObs,ns,bscat_r,scat_r,ext_r,g_r,Deq_r,vfall_r,
                        #        bscat,scat,ext,g,Deq,vfall_r,wl,dn,mu)
                        rwc,Z,att,sca1r,g1r=get_Zm(f,zT[i],ns,bscatKu_r,scatKu_r,
                                                   extKu_r,gKu_r,DeqKu_r,
                                                   vfallKu_r,
                                                   bscatKu,scatKu,extKu,gKu,
                                                   DeqKu,
                                                   vfallKu,wlKu,dn,mu)
                        zeta=zeta+att*4.343/10.0^(0.1*zT[i]*beta)*10.0^(0.1*zObs[i]*beta)*dr
                        zeta1d[i]=zeta
                        pwc[i]=rwc
                        kext1r=att
                    end
                end
                kext1d[i]=kext1s+kext1r
                scatt1d[i]=sca1s+sca1r
                g1d[i]=(sca1s*g1s+sca1r*g1r)/(scatt1d[i]+1e-10)
                zT[i]=Z
                zeta1d[i]=zeta
            end
        end
        piaMax=25
        if q*beta*zeta1d[nz]>0.99999
            eps=(1-10^(-0.1*piaMax*beta))/(q*beta*zeta1d[nz])
            dn1=dn1*eps#^(1.0/(1-beta))
        end
        #eps=(1-10^(-0.1*piaR*beta))/(q*beta*zeta1d[nz])
        #dn1=dn1*eps^(1.0/(1-beta))
        if q*beta*zeta1d[nz]!= q*beta*zeta1d[nz]
            println("dn1=",dn1)
            println(it, " ", piaR)
            println("eps=",eps)
            println("zeta=",zeta1d[:])
            exit(1)
        end
        #println(q*beta*zeta1d[nz])
        for i=1:nz
            zT[i]=zObs[i]-10/beta*log10(1-eps*q*beta*zeta1d[i])
        end
        piaHB=-10/beta*log10(1-eps*q*beta*zeta1d[nz])
        println("$(piaHB) $(it)")
    end
    #println(dn1,"  ", eps)
    return piaHB,zT,dn1,pwc,vdop1d, kext1d, scatt1d, g1d
end



function profIterativeMixed(zObs,dn1d,bscatKu,scatKu,extKu,gKu,DeqKu,
                         vfallKu,bscatKu_r,scatKu_r,extKu_r,gKu_r,
                         DeqKu_r,vfallKu_r,wlKu,hint,dr,
                         ns,mu,nmTop,nmBot)
    pia=0.0
    dn=-99.9
    dm=-99.9

    beta=0.72
    nz=size(zObs)[1]
    zeta1d=zeros(nz)

    kext1d=zeros(nz)
    scatt1d=zeros(nz)
    g1d=zeros(nz)
    vdop1d=zeros(nz)

    zT=copy(zObs)
    q=0.2*log(10)
    dn1=1.0
    eps=1.0
    pwc=zeros(nz)
    piaHB=0.0
    cumPIA=zeros(nz)
    for it=1:2
        cumPIA=cumPIA.*0
        pia=0
        for i=1:nz
            if(zT[i]>-10 && zObs[i]>-10)
                dn=dn1d[i]*dn1
                #println("dn= ",dn)
                Z=0.0
                sca1s=0.0
                kext1s=0.0
                g1s=0.0
                sca1r=0.0
                kext1r=0.0
                g1r=0.0
                if i<=nmTop
                    swc,Z,att,sca1s,g1s=get_Zs(zT[i],ns,bscatKu,scatKu,extKu,gKu,
                                               DeqKu,
                                               vfallKu,wlKu,dn,mu)
                    if(swc>10)
                        println(zT[i], " ", zObs[i], " ", dn, " ", dn1," ",swc, "it=",it)
                        #exit()
                    end
                    kext1s=att
                    pwc[i]=swc
                else
                    if i>=nmBot
                        rwc,Z,att,sca1r,g1r,vdopr=get_Zr(zT[i],bscatKu_r,scatKu_r,
                                                         extKu_r,gKu_r,DeqKu_r,
                                                         vfallKu_r,wlKu,dn,mu)
                        pwc[i]=rwc
                        kext1r=att
                        vdop1d[i]=vdopr
                        println("$(zT[i]) $(zObs[i]) $(rwc) $(pia)")
                    else
                        f=(i-nmTop)/(nmBot-nmTop)
                        # get_Zm(f,zObs,ns,bscat_r,scat_r,ext_r,g_r,Deq_r,vfall_r,
                        #        bscat,scat,ext,g,Deq,vfall_r,wl,dn,mu)
                        rwc,Z,att,sca1r,g1r=get_Zm(f,zT[i],ns,bscatKu_r,scatKu_r,
                                                   extKu_r,gKu_r,DeqKu_r,
                                                   vfallKu_r,
                                                   bscatKu,scatKu,extKu,gKu,
                                                   DeqKu,
                                                   vfallKu,wlKu,dn,mu)
                        pwc[i]=rwc
                        kext1r=att
                    end
                end
                kext1d[i]=kext1s+kext1r
                scatt1d[i]=sca1s+sca1r
                g1d[i]=(sca1s*g1s+sca1r*g1r)/(scatt1d[i]+1e-10)
                pia=pia+4.343*kext1d[i]*dr
                zT[i]=zObs[i]+pia
                pia=pia+4.343*kext1d[i]*dr

            end
        end
        println(cumPIA[end], " ",it)
    end
    piaHB=pia
    return piaHB,zT,dn1,pwc,vdop1d, kext1d, scatt1d, g1d
end


function kZprofiling(dn,kzCoeffs,kzKuCoeffs,pzCoeffs,zObs,dr,nmTop,nmBot,vdopL,zL)
    pia=0.0
    beta=0.72
    nz=size(zObs)[1]
    zeta1d=zeros(nz)

    kext1d=zeros(nz)
    scatt1d=zeros(nz)
    g1d=zeros(nz)
    vdop1d=zeros(nz)

    zT=copy(zObs)
    zKuObs=copy(zObs)
    zKu=copy(zObs)
    q=0.2*log(10)
    eps=1.0
    pwc=zeros(nz)
    piaHB=0.0
    dn1=dn
    beta=kzCoeffs[11,1]*10.0
    piaKu=0.0
    for it=1:3
        zeta=0.0
        for i=1:nz
            if(zT[i]>-10 && zObs[i]>-10)
                if i<=nmTop
                    att=1^(1-10*kzCoeffs[1,1])*10^(kzCoeffs[1,1]*zT[i]+kzCoeffs[1,2])
                else
                    if i>=nmBot
                        att=dn1^(1-10*kzCoeffs[11,1])*10^(kzCoeffs[11,1]*zT[i]+kzCoeffs[11,2])
                    else
                        f=(i-nmTop)/(nmBot-nmTop)
                        i0=Int(trunc(f*10))+1
                        if i0>11
                            i0=11
                        end
                        dn11=(1-f)+f*dn1
                        att=dn11^(1-10*kzCoeffs[i0,1])*10^(kzCoeffs[i0,1]*zT[i]+kzCoeffs[i0,2])
                    end
                end
                zeta=zeta+att*4.343/10.0^(0.1*zT[i]*beta)*10.0^(0.1*zObs[i]*beta)*dr
                zeta1d[i]=zeta
            end
        end
        #println("$(dn1) $(q*beta*zeta1d[nz])")
        piaMax=25
        if q*beta*zeta1d[nz]>0.99999
            eps=(1-10^(-0.1*piaMax*beta))/(q*beta*zeta1d[nz])
            dn1=dn1*eps#^(1.0/(1-beta))
        end
        #println(dn1, " ", eps)
        #exit()
        #eps=(1-10^(-0.1*piaR*beta))/(q*beta*zeta1d[nz])
        #dn1=dn1*eps^(1.0/(1-beta))
        if q*beta*zeta1d[nz]!= q*beta*zeta1d[nz]
            println("dn1=",dn1)
            println(it, " ", piaR)
            println("eps=",eps)
            println("zeta=",zeta1d[:])
            exit(1)
        end
        #println(q*beta*zeta1d[nz])
        piaKu=0.0
        for i=1:nz
            zT[i]=zObs[i]-10/beta*log10(1-eps*q*beta*zeta1d[i])
            zTn=zT[i]-10*log10(dn1)
            zTn=max(-10,zTn)
            zTn=min(50,zTn)

            if(zTn>-10 && zTn<50)
                if i<=nmTop
                    zTn=zT[i]
                    zKu[i]=pzCoeffs[1,1]*zTn^2+pzCoeffs[1,2]*zTn+pzCoeffs[1,3]
                    attKu=1^(1-10*kzKuCoeffs[1,1])*10^(kzKuCoeffs[1,1]*zT[i]+kzKuCoeffs[1,2])
                else
                    if i>=nmBot
                        zTn=max(0,zTn)
                        zTn=min(zTn,45)
                        zKu[i]=pzCoeffs[11,1]*zTn^2+pzCoeffs[11,2]*zTn+pzCoeffs[11,3]+10*log10(dn1)
                        attKu=dn1^(1-10*kzKuCoeffs[11,1])*10^(kzKuCoeffs[11,1]*zT[i]+kzKuCoeffs[11,2])
                    else
                        f=(i-nmTop)/(nmBot-nmTop)
                        i0=Int(trunc(f*10))+1
                        if i0>11
                            i0=11
                        end
                        dn11=(1-f)+f*dn1
                        zTn=zT[i]-10*log(dn11)
                        zKu[i]=pzCoeffs[i0,1]*zTn^2+pzCoeffs[i0,2]*zTn+pzCoeffs[i0,3]+10*log10(dn11)
                        attKu=dn11^(1-10*kzKuCoeffs[i0,1])*10^(kzKuCoeffs[i0,1]*zT[i]+kzKuCoeffs[i0,2])
                    end
                end
                piaKu=piaKu+4.343*attKu*dr
                zKu[i]-=piaKu
                piaKu=piaKu+4.343*attKu*dr
            else
                #println(zTn," ",zT[i], " ")
                zKu[i]=-99
            end
        end
        piaHB=-10/beta*log10(1-eps*q*beta*zeta1d[nz])
        #println("$(piaHB) $(it)")
    end
    z1=zT[nz]-10*log(dn1)
    iz=Int(trunc(z1))+11
    iz=max(1,iz)
    if iz>40
        iz=40
    end
    #println(dn1,"  ", eps)
    return piaHB,piaKu,zT,zKu,dn1,vdopL[iz]
end

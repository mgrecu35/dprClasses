
push!(PyVector(pyimport("sys")["path"]), "./")
#using .scatTables
readS=pyimport("readScattProf")

function readvar2d(fh,g1,g2,vname,nst)
    var=fh.groups[g1].groups[g2].variables[vname]
    nt1,nr=var.shape
    var=var._get([nst,0],[nt1,nr],[1,1])
    return var
end

function readvar3d(fh,g1,g2,vname,nst)
    var=fh.groups[g1].groups[g2].variables[vname]
    nt1,nr,nz=var.shape
    var=var._get([nst,0,0],[nt1,nr,nz],[1,1,1])
    return var
end
function readvar4dj_cmb(fh,g1,vname,nst)
    var=fh.groups[g1].variables[vname]
    nt1,nr,nz,nf=var.shape
    var=var._get([nst,0,0,0],[nt1,nr,nz,nf],[1,1,1,1])
    return var
end
function readnetcdfj(fname)
    println(fname)
    netcdf=pyimport("netCDF4")
    fh=netcdf.Dataset.(fname,"r")
    #println(fh)
    zKu=readvar3d(fh,"NS","PRE","zFactorMeasured",0)
    zKa=readvar3d(fh,"MS","PRE","zFactorMeasured",0)
    zKuC=readvar3d(fh,"NS","SLV","zFactorCorrected",0)
    sfcPrecip=readvar2d(fh,"NS","SLV","precipRateNearSurface",0)
    precip3D=readvar3d(fh,"NS","SLV","precipRate",0)
    piaHB=readvar2d(fh,"NS","SLV","piaFinal",0)
    #exit(-1)
    #lon=fh['NS/Longitude'][a[0][b],:]
    #lat=fh['NS/Latitude'][a[0][b],:]
    zsfc=readvar2d(fh,"NS","PRE","elevation",0)
    hzero=readvar2d(fh,"NS","VER","heightZeroDeg",0)
    bcf=readvar2d(fh,"NS","PRE","binClutterFreeBottom",0)
    brs=readvar2d(fh,"NS","PRE","binRealSurface",0)
    bst=readvar2d(fh,"NS","PRE","binStormTop",0)
    bzd=readvar2d(fh,"NS","VER","binZeroDeg",0)
    bbPeak=readvar2d(fh,"NS","CSF","binBBPeak",0)
    pType=readvar2d(fh,"NS","CSF","typePrecip",0)
    pType=Int.(trunc.(pType./1e7))
    sfcType=readvar2d(fh,"NS","PRE","landSurfaceType",0)
    reliabF=readvar2d(fh,"NS","SRT","reliabFlag",0)
    piaKu=readvar2d(fh,"NS","SRT","PIAhybrid",0)
    locZAngle=readvar2d(fh,"NS","PRE","localZenithAngle",0)
    fh.close()
    GC.gc()
    return zKu,zKa,sfcPrecip,zsfc,hzero,bcf,bst,brs,bzd,bbPeak,pType,
    sfcType,reliabF,piaKu,locZAngle,piaHB,precip3D,zKuC
end

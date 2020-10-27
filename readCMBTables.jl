include("src_jl/procTables.jl")
function readCMBTables()
    nmu=Ref{Int32}(5)
    nmfreq=Ref{Int32}(8)

    ccall((:readtablesliang2_,"combAlg"),Cvoid,(Ref{Int32},Ref{Int32}),
      nmu,nmfreq)
      ccall((:cloud_init_,"./combAlg"),Cvoid,(Ref{Cint},),nmfreq)
      ccall((:__nbinmod_MOD_init_nbin,"./combAlg"),Cvoid,(),)
      sysdN=Cfloat(-0.25) # if stratiform -=0.1
      nmemb=50
      ccall((:initgeophys_,"./combAlg"),Cvoid,(Ref{Cint},Ref{Cint},Ref{Cfloat}),
            nmemb,nmfreq,sysdN)
    dnr=-0.3
    dnbb=-0.3
    dns=-0.5
    newTables=get_CMB_Tables(dnr,dnbb,dns)
    return newTables
end

import matplotlib.pyplot as plt

from netCDF4 import *

fh=Dataset('era5_news.nc')
def readvar(nc,name):
    flux_x=nc[name]
    s=flux_x.scale_factor
    a=flux_x.add_offset
    flux_x=flux_x[:,:,:]*s+a
    print(s,a)
    return flux_x

te=fh['p86.162'][:,:,:]
ke=fh['p82.162'][:,:,:]
divgp=fh['p85.162'][:,:,:]
divth=fh['p83.162'][:,:,:]
divqv=fh['p84.162'][:,:,:] #kg m**-2 s**-1"

slhf=fh['slhf'][:,:,:]
#slhf:standard_name = "surface_upward_latent_heat_flux" ;
ssr=fh['ssr'][:,:,:]
#ssr:standard_name = "surface_net_downward_shortwave_flux" ;
str=fh['str'][:,:,:]
#str:"surface_net_upward_longwave_flux" ;
sshf=fh['sshf'][:,:,:]
#sshf:standard_name = "surface_upward_sensible_heat_flux" ;
tsr=fh['tsr'][:,:,:]
#"toa_net_upward_shortwave_flux" ;
ttr=fh['ttr'][:,:,:]
#ttr:standard_name = "toa_outgoing_longwave_flux" ;

plt.figure()
plt.pcolormesh(divgp.mean(axis=0),cmap='RdBu',vmax=500,vmin=-500)


plt.figure()
plt.pcolormesh(divth.mean(axis=0),cmap='RdBu',vmax=500,vmin=-500)
#plt.show()

plt.figure()
plt.pcolormesh(divgp.mean(axis=0)+divth.mean(axis=0)+2.5e6*divqv.mean(axis=0),cmap='RdBu',vmax=500,vmin=-500)
#plt.show()

plt.figure()
plt.pcolormesh(te.mean(axis=0),cmap='RdBu',vmax=500,vmin=-500)
#plt.show()

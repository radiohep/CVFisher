#!/usr/bin/env python
import classy
import numpy as np
from castorina import castorinaBias, castorinaPn
from pritchard import pritchardTSpin
#redshift ranges
zlow=np.arange(0,6.01,0.1)
zhigh=np.arange(30,50.01,1)

# Cosmology definituin
# I could do this:
#
#   from astropy.cosmology import Planck15 as cosmo
#
# but it doesn't have amplitude, so might as well copy everything from a table.
#
# So I took numbers from 
# https://wiki.cosmos.esa.int/planckpla2015/index.php/Cosmological_Parameters
# file https://wiki.cosmos.esa.int/planckpla2015/images/f/f7/Baseline_params_table_2015_limit68.pdf
#
# Setion 2.7 Planck +EE + BAO, bestfit values
#
# Note that I set neutrino mass to zero, since presumably some codes have internal
# hubble's law and probably don't want to deal with distribution function for neutrinos.
#

params = {
    'output': 'mPk',
    'A_s': 2.234e-9,
    'n_s': 0.96708, 
    'h': 0.67610,
    'omega_b': 0.022319,
    'omega_cdm': 0.11910,
    'tau_reio': 0.0865,
    'P_k_max_h/Mpc': 20.,
    'YHe': 0.245370,
    'N_ncdm': 0,
    'z_pk':", ".join(map(str,np.concatenate((zlow,zhigh))))
}

# Create an instance of the CLASS wrapper
cldr = classy.Class()
# set parameters and go
cldr.set(params)
cldr.compute()

def Tspin(z,OmegaH, OmegaM):
## Formula from Tzu-Ching, in mK, arXiv:0709.3672
    return 0.3*(OmegaH/1e-3)*np.sqrt((1+z)/(2.5)*0.29/(OmegaM+(1.-OmegaM)/(1+z)**3))

OmegaH=1e-3
OmegaM=(params['omega_b']+params['omega_cdm'])/params['h']**2
h=params['h']
# first let's compute background values at our zs
for name,zs in [("zlow",zlow),("zhigh",zhigh)]:
    f=open("background_%s.dat"%(name),'w')
    f.write("# z distance(z)[Mpc] growth(z)[normalised at z=0] Tspin(z) [mK]  b(z) N(z) [Mpc^3] \n")
    for z in zs:
        dist=cldr.luminosity_distance(z)/(1+z) # this is what class has
        gf=np.sqrt(cldr.pk(0.1,z)/cldr.pk(0.1,0)) # evaluate at k=0.1 h/Mpc, good enough
        if z<10:
            Ts=Tspin(z,OmegaH,OmegaM)
            bias=castorinaBias(z)
            Pn=castorinaPn(z)/h**3 # converting from Mpc/h^3 to Mpc^3
        else:
            Ts=pritchardTSpin(z)
            bias=1.0
            Pn=0.0
        f.write ("%3.2f %g %g %g %g %g  \n"%(z,dist,gf,Ts,bias,Pn))
    f.close()
# next write power spectrum
f=open("pk0.dat","w")
f.write("# k [1/Mpc] Pk [Mpc^3] \n")
for k in np.logspace(-3,1,200):
    f.write("%g %g\n"%(k,cldr.pk(k,0)))
f.close()
cldr.struct_cleanup()
cldr.empty()



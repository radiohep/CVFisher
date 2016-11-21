#!/usr/bin/env python
import classy
import numpy as np
# I could do this, but it doesn't have amplitude,
# so might as well copy everything
# from astropy.cosmology import Planck15 as cosmo
#
# So I took numbers from 
# https://wiki.cosmos.esa.int/planckpla2015/index.php/Cosmological_Parameters
# file https://wiki.cosmos.esa.int/planckpla2015/images/f/f7/Baseline_params_table_2015_limit68.pdf
# Setion 2.7 Planck +EE + BAO:
#

zs=np.concatenate((np.arange(0,6,0.1),np.arange(30,50,1)))
params = {
    'output': 'mPk',
    'A_s': 2.234e-9,
    'n_s': 0.96708, 
    'h': 0.67610,
    'omega_b': 0.022319,
    'omega_cdm': 0.11910,
    'tau_reio': 0.0865,
    'P_k_max_h/Mpc': 10.,
    'YHe': 0.245370,
    'N_ncdm': 0,
    'z_pk':", ".join(map(str,zs))
}

# Create an instance of the CLASS wrapper
cldr = classy.Class()
# set parameters and go
cldr.set(params)
cldr.compute()

Pk = cldr.pk(0.1,0)

print Pk


cldr.struct_cleanup()
cldr.empty()



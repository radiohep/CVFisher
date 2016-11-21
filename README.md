# CVFisher

This repo is to co-ordinate running of Fisher codes for the 21-cm CV group.
The structure is as follows:

## assumptions

Three files in this directory:

 * `background_zlow.dat`: distance and growth factor for z=0..6
 * `background_zhigh.dat`: distance and growth factor for z=30-50
 * `pk0.dat`: linear power spectrum as z=0 -- to be assumed to just scale with growth factor^2

The assumed cosmology was Planck15 + BAO *with massless* neutrinos. The motivation here is that perhaps some codes do have
an internal Friedman equation and don't want to properly deal with neutrino distribution function.

Finally, the `assumptions/src/make_cosmo.py` is a script used to generate the above. In needs `classy` module that
interfaces CLASS code to python (imperfectly).

Cosmology used for the above is defined as LCDM(don't forget to update this is `make_cosmo` changes!!) with:
```
   'A_s': 2.234e-9,
   'n_s': 0.96708,
   'h': 0.67610,
   'omega_b': 0.022319,
   'omega_cdm': 0.11910,
   'tau_reio': 0.0865,
```





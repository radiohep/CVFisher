# S/N forecasts

The code in here will generate power spectrum S/N forecasts (currently only for
low-z experiment A).

## Requirements

- numpy
- [cora](http://github.com/radiocosmology/cora) - for a few helper routines

## Output

The output is HDF5 files containing the S/N on each power spectrum bin. By S/N I
mean the ratio of the signal in each PS bin to the expected error on it (the
inverse fractional error). This is essentially the square root of the scaled
Fisher matrix diagonal.

There are a few datasets:

- `sn_all` - the S/N on the 2D power spectrum combined across all redshifts.
  Packed as P[k_par, k_perp].
- `sn_band_X` - the S/N for a single frequency sub-band. The central redshift of
  the band is stored as an attribute on the dataset.
- `k_bin` - the binning used for the power spectrum, this array gives the lower
  and upper edges of each power spectrum bin.

As github isn't a great place to share the large binary files produced, for the
moment I've put them up on my departmental web space:

 - [Experiment A, T_recv=50 K](http://www.phas.ubc.ca/~richard/sn_lowz_expA_50K.h5)
 - [Experiment A, T_recv=10 K](http://www.phas.ubc.ca/~richard/sn_lowz_expA_10K.h5)

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from scipy.integrate import quad
import py_cosmo_mad as csm
import sys

NU_21=1420.405751786
CLIGHT=299.792458
FWHM2G=0.42466090014
H0=0.68
NSIDE=32
NDISH=NSIDE*NSIDE
DDISH=10.
TYR=5.
EFF=1.
THR=TYR*365*24*EFF
FSKY=0.5

if len(sys.argv)!=2 :
    print "Usage: 4cast.py T_inst [K]"
    exit(1)
TINST=float(sys.argv[1])
EXPNAME='A%d'%(int(TINST))
print TINST,EXPNAME

#Generate cosmology object
pcs=csm.PcsPar()
pcs.background_set(0.3,0.7,0.05,-1.,0.,H0,2.7255)

#Read power spectrum
ka,pka=np.loadtxt("../assumptions/pk0.dat",unpack=True)
pkf=interp1d(np.log(ka),pka,bounds_error=False,fill_value=0)

#Read background quantities
za,dza,gza,fza,tza,bza,nza=np.loadtxt("../assumptions/background_zlow.dat",unpack=True)
dzf=interp1d(za,dza,bounds_error=False,fill_value=0)
gzf=interp1d(za,gza,bounds_error=False,fill_value=0)
fzf=interp1d(za,fza,bounds_error=False,fill_value=0)
tzf=interp1d(za,tza,bounds_error=False,fill_value=0)
bzf=interp1d(za,bza,bounds_error=False,fill_value=0)
nzf=interp1d(za,nza,bounds_error=False,fill_value=0)

def get_grid_square(ns,dl) :
    #Gets antenna positions for a grid of ns x ns with a spacing dl
    pos_1=(np.arange(ns)+0.5)*dl
    pos=np.array([(pos_1[:,None]*np.ones([ns,ns])).flatten(),
                  (pos_1[None,:]*np.ones([ns,ns])).flatten()])

    return pos

def get_basedist(pos,nbins,dmin=None,dmax=None,fname=None,dec=90.,weigh_dist=False) :
    #Gets baseline distribution given antenna positions for observations at declination dec
    fac_dec=np.sin(np.pi*dec/180)
    dpos=pos[:,:,None]-pos[:,None,:];
    dpos_b=np.array([dpos[0].flatten(),dpos[1].flatten()])
    dposr=np.sqrt(dpos_b[0]**2+dpos_b[1]**2)*fac_dec;
    ind0=np.where(dposr!=0)[0]; dposr=dposr[ind0]
    if dmin==None :
        dmin=np.amin(dposr)
    if dmax==None :
        dmax=np.amax(dposr)
    counts,bins=np.histogram(dposr,bins=nbins,range=[dmin,dmax])
    rm=(bins[1:]+bins[:-1])*0.5
    if weigh_dist :
        hr,bins=np.histogram(dposr,bins=nbins,range=[dmin,dmax],weights=dposr)
        rm[np.where(counts>0)]=(hr/counts)[np.where(counts>0)]
    vol=np.pi*(bins[1:]**2-bins[:-1]**2)
    n_d=counts/vol/2.;
    #print np.amin(rm), np.amax(rm), np.amin(dposr), np.amax(dposr), np.sum(n_d*vol*2)-len(pos[0])*(len(pos[0])-1.)

    if fname!=None :
        np.savetxt(fname,np.transpose([rm,n_d]))
    return rm,n_d,counts

def get_basedist_grid(nside,d_min,fac=1,fname=None,decmax=80.,decstep=0.5,weigh_dist=False) :
    #Wrapper of the previous function for grid-like antenna positions
    pos=get_grid_square(nside,d_min)
    ndec=int(decmax/decstep)+1
    decarr=90.-decmax*np.arange(ndec)/(ndec-1.)
    nbins=int(fac*nside*0.75)

    #Here we average over declinations to smooth out the distribution a bit
    nbarr=np.zeros(nbins); darr=np.zeros(nbins); dcumul=np.zeros(nbins); cumul=np.zeros(nbins)
    for dec in decarr :
        d,n,c=get_basedist(pos,nbins,dmin=0,dmax=1.5*nside*d_min,fname=fname,dec=dec,weigh_dist=weigh_dist)
        darr+=d; nbarr+=n; cumul+=c; dcumul+=d*c
    darr/=ndec; darr[cumul>0]=(dcumul/cumul)[cumul>0]; nbarr/=ndec
    return darr,nbarr

#Generate baseline density for this experiment
darr,nbarr=get_basedist_grid(NSIDE,DDISH,decmax=30.)
ndistint=interp1d(darr,nbarr*darr*2*np.pi,bounds_error=False,fill_value=0.)
norm=0.5*NDISH*(NDISH-1.)/quad(ndistint,darr[0],darr[-1])[0];
nbarr*=norm; ndistf=interp1d(darr,nbarr,bounds_error=False,fill_value=0.)

def tsky(nu) :
    #Returns sky temperature
    return 2E6*(100./nu)**2.4+2.7E3

def get_noisepower(tsys,nu,ndish,fsky,t_total,d_dish,typ='both') :
    #Returns noise power spectrum (and the values of k_perp at which it's evaluated)
    z=NU_21/nu-1.
    sigma2_noise=((tsky(nu)+tsys*1000)/tzf(z))**2
    sigma2_noise*=4*np.pi*fsky/(3.6E9*t_total*NU_21)
    lam0=CLIGHT/nu

    beam_fwhm=lam0/d_dish
    l_arr=np.arange(10000)
    if (typ=='interferometer') or (typ=='both') :
        n_baselines=ndistf(l_arr*lam0/(2*np.pi))
        factor_beam_if=n_baselines*(lam0/beam_fwhm)**2
    else :
        factor_beam_if=1E-16*np.ones_like(l_arr)
    if (typ=='single_dish') or (typ=='both') :
        beam_rad=beam_fwhm*FWHM2G
        factor_beam_sd=ndish*np.exp(-(l_arr*beam_rad)**2)
    else :
        factor_beam_sd=1E-16*np.ones_like(l_arr)
    factor_beam=factor_beam_sd+factor_beam_if

    cl_noise=sigma2_noise/factor_beam
    chi=pcs.radial_comoving_distance(1./(1+z))
    pk_noise=cl_noise*(chi*(1+z))**2/pcs.hubble(1./(1+z))
    k_arr=l_arr/chi
    return k_arr,pk_noise

def pk2d_signal(kpar,kperp,z) :
    #Signal power spectrum (units of (Mpc/h)^3)
    k2=kpar*kpar+kperp*kperp
    mu2=kpar*kpar/k2
    return H0**3*(gzf(z)*(bzf(z)+fzf(z)*mu2))**2*pkf(0.5*np.log(k2))

def pk2d_noise(kpar,kperp,z,typ='both') :
    #Noise power spectrum (units of (Mpc/h)^3)
    k,nk=get_noisepower(TINST,NU_21/(1+z),NDISH,FSKY,THR,DDISH,typ=typ)
    nkf=interp1d(k,nk,bounds_error=False,fill_value=0)
    return nkf(kperp)

#Generate power spectrum arrays for this experiment
nz=11
nk=128
zarr=0.5*(np.arange(nz)+1)
lkarr=-3+3*np.arange(nk)/(nk-1.)
kpar_arr=10.**(lkarr[:,None]*(np.ones(nk))[None,:])
kperp_arr=10.**(lkarr[None,:]*(np.ones(nk))[:,None])
pks_arr=np.zeros([nz,nk,nk])
pkni_arr=np.zeros([nz,nk,nk])
pkns_arr=np.zeros([nz,nk,nk])
pknb_arr=np.zeros([nz,nk,nk])

for i,z in enumerate(zarr) :
    pks_arr[i,:,:]=pk2d_signal(kpar_arr,kperp_arr,z)
    pkni_arr[i,:,:]=pk2d_noise(kpar_arr,kperp_arr,z,typ='interferometer')
    pkns_arr[i,:,:]=pk2d_noise(kpar_arr,kperp_arr,z,typ='single_dish')
    pknb_arr[i,:,:]=pk2d_noise(kpar_arr,kperp_arr,z,typ='both')

np.savez("sn_pk_exp"+EXPNAME,
         t_inst=TINST,z_arr=zarr,kpar_arr=kpar_arr[:,0],kperp_arr=kperp_arr[0,:],
         pk_signal=pks_arr,
         pk_noise_interferometer=pkni_arr,
         pk_noise_singledish=pkns_arr,
         pk_noise_hybrid=pknb_arr)

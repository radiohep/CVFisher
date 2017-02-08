#!/usr/bin/env python
#
# module assumptions, effectively deals with
#
# quantities specifie in assumptions
#
# all units are Mpc rather than Mpc/h!
#

import numpy as np
from scipy.interpolate import interp1d

class Assumptions:
    def __init__ (self,ldir="./assumptions"):
        kk,pk=np.loadtxt(ldir+"/pk0.dat",skiprows=1,unpack=True)
        lk=np.log(kk)
        lpk=np.log(pk)
        self.lpki=interp1d(lk,lpk)
        z,d,g,Ts,b,N=np.loadtxt(ldir+"/background_zlow.dat",unpack=True)
        self.di=interp1d(z,d)
        self.gi=interp1d(z,g)
        dz=z[1]-z[0]
        dlngdz=(g[2:]-g[0:-2])/(2*dz)/g[1:-1]
        zc=z[1:-1]
        f=dlngdz*(-(1+zc)) ## d/dlna = d/dz * (1+z)
        self.fi=interp1d(zc,f)
        self.Tsi=interp1d(z,Ts)
        self.bi=interp1d(z,b)
        self.Ni=interp1d(z,N)

    def KBinModes(self,z1,z2,k1,k2, A=4*np.pi):
        """ 
        Gives the number of comoving modes between z1 and z2 and k1 and k2
        for a full sky experiment, unless A is specified
        """
        vol=A/3.*(self.di(z2)**3-self.di(z1)**3)
        volk=4*np.pi/3*(k2**3-k1**3)
        return vol*volk


    def KBinModes2Ddk(self,z1,z2,kperp, dk,  A=4*np.pi):
        """ 
        Gives the number of comoving modes between z1 and z2 and k1 and k2
        for a full sky experiment, unless A is specified
        """
        vol=A/3.*(self.di(z2)**3-self.di(z1)**3)
        volk=2*np.pi*kperp*dk*dk
        return vol*volk
    
    def PkisoDM(self,k,z):
        """ Returns isotropic power spectrum for DM at redshift z"""
        return self.lpki(np.log(k))*self.gi(z)
    
    def Pkmuz(self, k,mu,z, verbose=False):
        """ Returns anisotropic 21 cm power spectrum at k, mu and z
         in mK^2"""
        Pk=self.PkisoDM(k,z)
        f=self.fi(z)
        bias=self.bi(z)
        Tspin=self.Tsi(z)
        if verbose:
            print "At z=",z
            print "Tspin=",Tspin
            print "bias=",bias
            print "f=",f
            print "Pk=",Pk
        return Tspin**2*(bias + f*mu**2)**2 * Pk
    
if __name__=="__main__":
    A=Assumptions()
    tmp=A.Pkmuz(0.1,0.5,5,True)
    print "Power=",tmp
    

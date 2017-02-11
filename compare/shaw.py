#!/usr/bin/env python
import h5py
import numpy as np
class Shaw:
    def __init__ (self,A,fname="forecast_shaw/sn_lowz_expA_50K.h5"):
        sh=h5py.File(fname)
        self.ks=sh['k_bin'].value[:-1] ## we have bin centres
        self.dk=self.ks[1]
        self.ks+=self.dk/2 ## bin centres
        maxk=self.dk*(len(sh['k_bin'].value)-1)
        bi=0
        vals=[]
        zs=[]
        while True:
            try:
                ext=sh['sn_band_%i'%bi]
            except KeyError:
                break
            zs.append(ext.attrs['z'])
            vals.append(ext.value)
            bi+=1
        ## get band edges
        vals=vals[::-1]
        zs=zs[::-1]
        print "Have ",bi," z bins:",zs
        self.zs=np.array(zs)
        ztmp=list(0.5*(self.zs[:-1]+self.zs[1:]))
        self.zlow=np.array([0]+ztmp)
        self.zhigh=np.array(ztmp+[self.zs[-2]+(self.zs[-1]-self.zs[-2])*2])
        print self.zlow,self.zhigh
        # calculate SNR per mode
        kperp=np.outer(np.ones(len(self.ks)),self.ks)
        valspermode=[]
        for i,z in enumerate(self.zs):
            zl=self.zlow[i]
            zh=self.zhigh[i]
            Nm=A.KBinModes2Ddk(zl,zh,kperp,self.dk, self.dk)
            valspermode.append(vals[i]*np.sqrt(Nm))
        self.vals=vals
        self.valspermode=valspermode

        
    def getSNRperMode (self,kpar,kperp,z,A):
        #first find relevant indices
        #nearest grid point in k
        i=int(kperp/self.dk+0.5)
        j=int(kpar/self.dk+0.5)
        zi=0
        ## this will extrpolate beyon last point
        if (z>self.zs[-1]):
            zi=len(self.zs)-2
        else:
            for zi in range(len(self.zs)-1):
                if (z>=self.zs[zi]) and (z<self.zs[zi+1]): #ok, excessive i know, but clear
                    break
        DZ=self.zs[zi+1]-self.zs[zi]
        w=(z-self.zs[zi])/DZ
        return self.valspermode[zi][j,i]*(1-w)+self.valspermode[zi+1][j,i]*w

if __name__=="__main__":
    import assumptions as A
    A=A.Assumptions()
    s=Shaw(A)
    for z in np.arange(0.1,5,0.1):
        print z,s.getSNRperMode(0.1,0.1,z,A), s.getSNRperMode(0.2,0.2,z,A)
            
        

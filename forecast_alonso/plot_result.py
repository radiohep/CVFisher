import numpy as np
import matplotlib.pyplot as plt
import sys

if len(sys.argv)!=2 :
    print "Usage: plot_result.py fname_data"
    exit(1)
fname_data=sys.argv[1]

data=np.load(fname_data)
#Fields included in the files:
# t_inst <- instrument temperature in Kelvin
# z_arr <- array of redshifts
# kpar_arr <- array of k_parallel values
# kperp_arr <- array of k_perpendicular values
# pk_signal <- 3D array containing the values of the signal power spectrum.
#              data['pk_signal'][iz,ik1,ik2] contains the power spectrum at
#              z=data['z_arr'][iz], k_par=data['kpar_arr'][ik1] and k_perp=data['kperp_arr'][ik2]
# pk_noise_interferometer <- Noise power spectrum assuming interferometer mode
# pk_noise_single         <- Noise power spectrum assuming single-dish mode
# pk_noise_hybrid         <- Noise power spectrum assuming interferometer and single-dish observations

print "Instrument temperature is %.1lf K"%(data['t_inst'])
for i,z in enumerate(data['z_arr']) :
    plt.figure()
    plt.title("$z=%lf$"%z)
    plt.plot(data['kperp_arr'],data['pk_signal'][i,0,:]              ,'k-',label='Signal power spectrum')
    plt.plot(data['kperp_arr'],data['pk_noise_interferometer'][i,0,:],'r-',label='Noise, interferometer')
    plt.plot(data['kperp_arr'],data['pk_noise_singledish'][i,0,:]    ,'g-',label='Noise, single-dish')
    plt.plot(data['kperp_arr'],data['pk_noise_hybrid'][i,0,:]        ,'b-',label='Noise, both')
    plt.xlabel('$k_\\perp\\,[h\\,{\\rm Mpc}^{-1}]$',fontsize=16)
    plt.ylabel('$P(k_\\parallel=%.4lf,k_\\perp)\\,[{\\rm Mpc}/h)^3]$'%(data['kpar_arr'][0]),fontsize=16)
    plt.ylim([0.5*np.amin(data['pk_signal']),2*np.amax(data['pk_signal'])])
    plt.loglog()
    plt.legend(loc='center left',frameon=False,labelspacing=0.1)
plt.show()

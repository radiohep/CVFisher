#
# Bias and shot noise from Emanuele Castorina
#
# Evolution of HI bias and shot noise as a function of redshift 
# From equation 4 and 5 of https://arxiv.org/abs/1609.05157. In equation 7, alpha=1 and M_{min} = 3.5*10^10 Msun/h to fit DLA bias at z=2.3
# The only ASSUMPTION is that M_{min} does not depend on redshift, note also that bias and shot noise are independent from \Omega_{HI}
# | z | b_{HI}(z) | P_{SN}(z)
#
# AS: Added point at z=0 with extrapolated bias and Pn the same as first point
#
import numpy as np
from scipy.interpolate import interp1d
z_, bias_, Pn_=np.array(map(lambda x:map(float,x.split()),
"""0       1.42    528.1
0.8 	1.6526	528.1 
0.95	1.7082	430.63
1.1 	1.7671	356.19
1.25	1.8294	298.82
1.4 	1.895 	254.19
1.55	1.9638	219.13
1.7 	2.0357	191.36
1.85	2.1108	169.17
2.  	2.1891	151.33
2.15	2.2704	136.9 
2.3 	2.3549	125.18
2.45	2.4424	115.63
2.6 	2.533 	107.86
2.75	2.6268	101.54
2.9 	2.7235	96.437
3.05	2.8234	92.366
3.2 	2.9264	89.18 
3.35	3.0324	86.766
3.5 	3.1415	85.039
3.65	3.2536	83.932
3.8 	3.3688	83.397
3.95	3.4871	83.399
4.1 	3.6084	83.917
4.25	3.7328	84.938
4.4 	3.8603	86.46 
4.55	3.9908	88.49 
4.7 	4.1244	91.044
4.85	4.261 	94.145
5.  	4.4007	97.825
5.15	4.5434	102.13
5.3 	4.6892	107.1 
5.45	4.8381	112.81
5.6 	4.99  	119.32
5.75	5.145 	126.73
5.9 	5.303 	135.13
6.05	5.464 	144.64
6.2 	5.6282	155.4 """.split("\n"))).T
castorinaBias=interp1d(z_,bias_)
castorinaPn=interp1d(z_,Pn_)

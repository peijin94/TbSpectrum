### Peijin 2021 11 15


import astropy
import numpy as np
from astropy.io import ascii
from scipy.interpolate import interp2d

data=ascii.read('gauntff.dat',data_start=42,data_end=188,comment='\s@')
#dataErr=ascii.read('gauntff.dat',data_start=192,data_end=338,comment='\s@')
gauntArr=np.array([np.array(data['col'+str(i)]) for i in np.arange(81)+1 ])
#gauntErr=np.array([np.array(dataErr['col'+str(i)]) for i in np.arange(81)+1 ])

def gauntff(nu,tempr,gauntArr=gauntArr):

    z = 1
    hh = 6.62607551e-27
    kk = 1.38006504e-16
    ry = 2.17987e-11
    ngd = 81
    nud = 146
    gam20 = -6.
    u0 = -16.
    dex = 2e-1
    g = np.arange(ngd)*dex+gam20
    u = np.arange(nud)*dex+u0

    find_u = np.log10(hh/kk*nu/tempr)
    find_g = np.log10(z**2*ry/kk/tempr)
    
    gffFunc = interp2d(g,u,gauntArr.T)
    
    return(np.array([gffFunc(find_g[i],find_u[i]) 
                     for i in range(find_g.shape[0])] ))
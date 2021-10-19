import numpy as np
import math as mt
import matplotlib.pyplot as plt
import scipy.integrate as integrate
import scipy.optimize as opt
from locale import atof

## Open file
fp = open('dens.txt', "r")
line = fp.readline()

def column_value(line,col):
    a = line.split(' ')
    idx = 0
    for x in a:
        if x!='':
            idx=idx+1
        if idx==col:
            return atof(x)
    return 0

def enclosed_mass(r,dens,r_max=-1.):
    r_avg = (r[1:]+r[0:-1])/2.
    dr    = (r[1:]-r[0:-1])
    dens_avg = (dens[1:]+dens[0:-1])/2.
    if r_max<0. or r_avg[-1]<=r_max:
        return 4* np.pi* np.sum(r_avg**2*dr*dens_avg)
    else:
        id_max = np.where(r_avg>r_max)[0][0]
        r_avg = r_avg[0:id_max]
        dr    = dr[0:id_max]
        dens_avg = dens_avg[0:id_max]

        return 4* np.pi* np.sum(r_avg**2*dr*dens_avg)

center = 0.5
r=[]
dens=[]
while line:
    line = fp.readline()
    r.append(np.abs(column_value(line,4)-center)*3**0.5)
    dens.append(column_value(line,7)**0.5)
r = np.sort(r)
dens = -np.sort(-np.array(dens))

rm = 0.3

M = enclosed_mass(r,dens,r_max=rm)
print('enclosed mass:'+str(M))
print('circular velocity:'+str((M/rm)**0.5))
fp.close()
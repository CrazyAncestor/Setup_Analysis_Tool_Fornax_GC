import numpy as np
import math as mt
import matplotlib.pyplot as plt
import scipy.integrate as integrate

def King_dens(x,xt):
    C = (1+xt**2)**(-0.5)
    return ((1+x**2)**(-0.5) -C)

def density(rho_s,r0,rt,r,model_name):
    x  =  r/r0
    xt = rt/r0
    if model_name =="King":
        return rho_s*King_dens(x,xt)

def clustermass(rho_s,r0,rt,r,model_name):
    x = r/r0
    xt = rt/r0
    def massbase(x):
        if model_name =="King":
            return 4*np.pi*x*x*(r0**3) *King_dens(x,xt)
       
    f = integrate.nquad(massbase, [[0,x]])[0] *rho_s
    return f

def getrho(mass,r0,rt,model_name):
    rho_0 =1.0
    mass_init = clustermass(rho_0,r0,rt,rt,model_name)
    rho_0 = rho_0*mass/mass_init
    return rho_0

def create_table(table_filename,r,dens):
    f = open(table_filename, "w")
    f.write("GC Density Profile\n")

    for i in range(len(r)):
        if dens[i]!=0.:
            a = str(r[i])+" "+str(dens[i])+"\n"
            f.write(a)

    f.close()

if __name__== '__main__':

    #unit of length = 1kpc
    #unit of mass   = 1M Mass of sun
    GC_mass = 0.75 
    GC_rs = 0.01
    GC_tr = 0.059

    rho_GC = getrho(GC_mass,GC_rs,GC_tr,'King')
    
    r = np.logspace(np.log10(GC_rs/10),np.log10(GC_tr),1000)
    dens = np.array([density(rho_GC,GC_rs,GC_tr,r[i],'King') for i in range(len(r))])
    create_table('profile_GC.txt',r,dens)
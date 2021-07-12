import numpy as np
import scipy.integrate as integrate
import matplotlib.pyplot as plt

def NFW_dens(x):
    return (1/((1+x)*(1+x)*x))
def Burkert_dens(x):
    return ((1/((1+x)*(1+x*x))))
def Plummer_dens(x):
    return (1+x*x)**(-2.5)

def density(rho_s,r0,rt,r,model_name):
    x  =  r/r0
    if model_name =="Burkert":
        return rho_s*Burkert_dens(x)
    elif model_name =="NFW":
        return rho_s*NFW_dens(x)
    elif model_name =="Plummer":
        return rho_s*Plummer_dens(x)

def clustermass(rho_s,r0,rt,r,model_name):
    x = r/r0
    xt = rt/r0
    def massbase(x):
        if model_name =="Burkert":
            return 4*np.pi*x*x*(r0**3) *Burkert_dens(x)
        elif model_name =="NFW":
            return 4*np.pi*x*x*(r0**3) *NFW_dens(x)
        elif model_name =="Plummer":
            return 4*np.pi*x*x*(r0**3) *Plummer_dens(x)
    f = integrate.nquad(massbase, [[0,x]])[0] *rho_s
    return f

def getrho(mass,r0,rt,model_name):
    rho_0 =1.0
    mass_init = clustermass(rho_0,r0,rt,rt,model_name)
    rho_0 = rho_0*mass/mass_init
    return rho_0

def create_table(table_filename,r,dens):
    f = open(table_filename, "w")
    f.write("Fornax Density Profile\n")

    for i in range(len(r)):
        if dens[i]!=0.:
            a = str(r[i])+" "+str(dens[i])+"\n"
            f.write(a)

    f.close()

if __name__== '__main__':

    #unit of length = 1kpc
    #unit of mass   = 1M Mass of sun
    ri = 1.0
    rt = 6
    halo_rs = 0.25
    halo_mass = 318
    halo_rt = rt
    halo_type = "Burkert"

    stellar_rs = 0.668
    stellar_mass = 38.2
    stellar_rt = rt
    stellar_type = "Plummer"

    rho_halo= getrho(halo_mass,halo_rs,halo_rt,halo_type)
    rho_stellar= getrho(stellar_mass,stellar_rs,stellar_rt,stellar_type)

    print('scaling density:')
    print("DM Halo Density :"+str(rho_halo))
    print("Stellar Density :"+str(rho_stellar))

    enclosed_mass = clustermass(rho_halo,halo_rs,halo_rt,ri,halo_type) +clustermass(rho_stellar,stellar_rs,stellar_rt,ri,stellar_type)
    v = (enclosed_mass/ri)**0.5

    print('-------------------------------------------')
    print('intial velocity:')
    print("GC V :"+str(v/(2**0.5) *(stellar_mass+halo_mass)/(stellar_mass+halo_mass+0.75)))
    print("Halo V :"+str(v/(2**0.5) *(-0.75)/(stellar_mass+halo_mass+0.75)))

    r = np.logspace(np.log10(halo_rs/10),np.log10(rt),1000)
    dens = np.array([density(rho_halo,halo_rs,halo_rt,r[i],halo_type)+density(rho_stellar,stellar_rs,stellar_rt,r[i],stellar_type) for i in range(len(r))])
    create_table('profile_Fornax.txt',r,dens)
import numpy as np
import scipy.integrate as integrate
import matplotlib.pyplot as plt

def NFW_dens(x):
    return (1/((1+x)*(1+x)*x))
def Burkert_dens(x):
    return ((1/((1+x)*(1+x*x))))
def Plummer_dens(x):
    return (1+x*x)**(-2.5)
def LC_dens(x):
    gamma_0 = 0.07
    gamma_inf = 4.65
    eta = 3.7
    r0 = 1.4
    return x**(-gamma_0)*(1+x**eta)**((gamma_0-gamma_inf)/eta)
def King_dens(x):
    return ((1+x**2)**(-0.5))

def density(rho_s,r0,r,model_name):
    x  =  r/r0
    if model_name =="Burkert":
        return rho_s*Burkert_dens(x)
    elif model_name =="NFW":
        return rho_s*NFW_dens(x)
    elif model_name =="Plummer":
        return rho_s*Plummer_dens(x)
    elif model_name =="LC":
        return rho_s*LC_dens(x)
    elif model_name =="King":
        return rho_s*King_dens(x)
    else:
        return 0

def clustermass(rho_s,r0,r,model_name):
    
    x = r/r0
    def massbase(x):
        if model_name =="Burkert":
            return 4*np.pi*x*x*(r0**3) *Burkert_dens(x)
        elif model_name =="NFW":
            return 4*np.pi*x*x*(r0**3) *NFW_dens(x)
        elif model_name =="Plummer":
            return 4*np.pi*x*x*(r0**3) *Plummer_dens(x)
        elif model_name =="LC":
            return 4*np.pi*x*x*(r0**3) *LC_dens(x)
        elif model_name =="King":
            return 4*np.pi*x*x*(r0**3) *King_dens(x)
        else:
            return 0
    f = integrate.nquad(massbase, [[0,x]])[0] *rho_s
    return f

def getrho(mass,r0,rt,model_name):
    if model_name=="None":
        return 0
    rho_0 =1.0
    mass_init = clustermass(rho_0,r0,rt,model_name)
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

def SetCriteria(Halo, Stellar, GC, cell_size_base_level, max_level, ratio):

    halo_rs, halo_mass, halo_rt, halo_type = Halo
    stellar_rs, stellar_mass, stellar_rt, stellar_type = Stellar
    GC_rs, GC_mass, GC_rt, GC_type = GC

    ri = r_orbit
    
    rho_halo = getrho(halo_mass,halo_rs,halo_rt,halo_type)
    rho_stellar = getrho(stellar_mass,stellar_rs,stellar_rt,stellar_type)
    rho_GC = getrho(GC_mass,GC_rs,GC_rt,GC_type)

    cell_size = np.array([cell_size_base_level*2**(-i) for i in range(max_level+1)])
    dens = np.array([density(rho_halo,halo_rs,cell_size[i]*ratio,halo_type)+density(rho_stellar,stellar_rs,cell_size[i]*ratio,stellar_type) for i in range(max_level+1)])
    criteria = [dens[i]*(cell_size[i]**3) for i in range(max_level+1)]

    for i in range(len(criteria)):
        if(i>2):
            criteria[i] = rho_GC*(cell_size[i]**3) 
    print("criteria for each level:")
    print(criteria)


if __name__== '__main__':

    #unit of length = 1kpc
    #unit of mass   = 1M Mass of sun

    cell_size_base_level = 0.0390625*4
    max_level = 6
    ratio = 5

    Paper = 'Cole et al. 2012'

    if Paper == 'Cole et al. 2012':
        r_orbit = 0.4
        r_cutoff = 10.

        halo_rs = 1.4
        halo_mass = 412
        halo_rt = 1.8
        halo_type = "LC"
        Halo = halo_rs, halo_mass, halo_rt, halo_type

        stellar_rs = 0
        stellar_mass = 0
        stellar_rt = 0
        stellar_type = "None"
        Stellar = stellar_rs, stellar_mass, stellar_rt, stellar_type

        GC_mass = 0.13
        GC_rs = 0.005
        GC_rt = 0.08
        GC_type = "Plummer"
        GC = GC_rs, GC_mass, GC_rt, GC_type

        SetCriteria(Halo, Stellar, GC, cell_size_base_level, max_level, ratio)

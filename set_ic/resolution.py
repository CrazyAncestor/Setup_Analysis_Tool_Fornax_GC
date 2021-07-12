import numpy as np
import scipy.integrate as integrate
import matplotlib.pyplot as plt

def NFW_dens(x):
    return (1/((1+x)*(1+x)*x))
def Burkert_dens(x):
    return ((1/((1+x)*(1+x*x))))
def Plummer_dens(x):
    return (1+x*x)**(-2.5)
def King_dens(x,xt):
    C = (1+xt**2)**(-0.5)
    return ((1+x**2)**(-0.5) -C)


def density(rho_s,r0,rt,r,model_name):
    x  =  r/r0
    xt = rt/r0
    if model_name =="Burkert":
        return rho_s*Burkert_dens(x)
    elif model_name =="NFW":
        return rho_s*NFW_dens(x)
    elif model_name =="Plummer":
        return rho_s*Plummer_dens(x)
    if model_name =="King":
        return rho_s*King_dens(x,xt)

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
        elif model_name =="King":
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
    f.write("Fornax Density Profile\n")

    for i in range(len(r)):
        if dens[i]!=0.:
            a = str(r[i])+" "+str(dens[i])+"\n"
            f.write(a)

    f.close()

if __name__== '__main__':

    #unit of length = 1kpc
    #unit of mass   = 1M Mass of sun

    cell_size_base_level = 0.0390625*4
    max_level = 6
    ratio = 5

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

    GC_mass = 0.75 
    GC_rs = 0.01
    GC_tr = 0.059

    rho_GC = getrho(GC_mass,GC_rs,GC_tr,'King')
    rho_halo= getrho(halo_mass,halo_rs,halo_rt,halo_type)
    rho_stellar= getrho(stellar_mass,stellar_rs,stellar_rt,stellar_type)

    cell_size = np.array([cell_size_base_level*2**(-i) for i in range(max_level+1)])
    dens = np.array([density(rho_halo,halo_rs,halo_rt,cell_size[i]*ratio,halo_type)+density(rho_stellar,stellar_rs,stellar_rt,cell_size[i]*ratio,stellar_type) for i in range(max_level+1)])
    criteria = [dens[i]*(cell_size[i]**3) for i in range(max_level+1)]

    for i in range(len(criteria)):
        if(i>2):
            criteria[i] = rho_GC*(cell_size[i]**3) 
    print("criteria for each level:")
    print(criteria)

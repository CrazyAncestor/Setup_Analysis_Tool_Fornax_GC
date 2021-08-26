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
    #return x**(-gamma_0)*(1+x**eta)**((gamma_0-gamma_inf)/eta)
    if x<=1:
        return 1.#flat core
    else:
        return x**(-3)

def density(rho_s,r0,rt,r,model_name):
    x  =  r/r0
    if model_name =="Burkert":
        return rho_s*Burkert_dens(x)
    elif model_name =="NFW":
        return rho_s*NFW_dens(x)
    elif model_name =="Plummer":
        return rho_s*Plummer_dens(x)
    elif model_name =="LC":
        return rho_s*LC_dens(x)

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
        elif model_name =="LC":
            return 4*np.pi*x*x*(r0**3) *LC_dens(x)
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
    ri = 0.4
    rt = 1.8
    halo_rs = 1.4
    halo_mass = 412
    halo_rt = rt
    halo_type = "LC"

    """stellar_rs = 0.668
    stellar_mass = 38.2
    stellar_rt = rt
    #stellar_type = "Plummer
    """
    rho_halo= getrho(halo_mass,halo_rs,halo_rt,halo_type)
    #rho_stellar= getrho(stellar_mass,stellar_rs,stellar_rt,stellar_type)

    print('scaling density:')
    print("DM Halo Density :"+str(rho_halo))
    #print("Stellar Density :"+str(rho_stellar))

    enclosed_mass = clustermass(rho_halo,halo_rs,halo_rt,ri,halo_type)# +clustermass(rho_stellar,stellar_rs,stellar_rt,ri,stellar_type)
    v = (enclosed_mass/ri)**0.5
    theta = np.pi/4.
    rt = 10.

    print('-------------------------------------------')
    print("Enclosed Mass :"+str(enclosed_mass))
    print('intial x,y:')
    print("GC x :"+str(ri*np.cos(theta) ))
    print("GC y :"+str(ri*np.sin(theta) ))
    print('intial velocity:')
    M_GC = 131825/1e6 
    halo_mass =  clustermass(rho_halo,halo_rs,halo_rt,rt,halo_type)
    print("GC Vx :"+str(-v*np.sin(theta) *(halo_mass)/(halo_mass+M_GC)))
    print("GC Vy :"+str( v*np.cos(theta) *(halo_mass)/(halo_mass+M_GC)))
    print("Halo Vx :"+str(v*np.sin(theta) *(M_GC)/(halo_mass+M_GC)))
    print("Halo Vy :"+str(-v*np.cos(theta) *(M_GC)/(halo_mass+M_GC)))

    r = np.logspace(np.log10(halo_rs/100),np.log10(10),1000)
    dens = np.array([density(rho_halo,halo_rs,10,r[i],halo_type) for i in range(len(r))])
    create_table('profile_Fornax.txt',r,dens)
    print(halo_mass)

    plt.plot(r,dens)
    plt.xlabel('r (kpc)')
    plt.ylabel('dens (1e06 $M_{\odot}}/kpc^3$)')
    #plt.xscale('log')
    #plt.yscale('log')
    #plt.show()

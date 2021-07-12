import yt
import numpy as np
from yt.data_objects.particle_filters import add_particle_filter
from matplotlib import pyplot as plt

idx_start = 0
idx_end = 25
filein  = [ 'Data_%06d'%idx for idx in range(idx_start, idx_end+1) ] 
fileout = "fig__cluster__radius"
nbin    = 50
dpi     = 150

V_cluster = []
Time = []
for i in range(idx_start, idx_end+1):
   # load data
   ds = yt.load( filein[i] )
   # define the particle filter for the newly formed stars
   def halo_star( pfilter, data ):
      filter = (data[ "all", "particle_mass" ] >0)
      return filter

   add_particle_filter( "halo_star", function=halo_star, filtered_type="all", requires=["particle_mass"] )
   ds.add_particle_filter( "halo_star" )


   # get the mass and creation time of the new stars
   ad            = ds.all_data()

   velocity_x        = ad[ "halo_star", "particle_velocity_x" ].in_units( "kpc/Gyr" )
   velocity_y        = ad[ "halo_star", "particle_velocity_y" ].in_units( "kpc/Gyr" )
   velocity_z        = ad[ "halo_star", "particle_velocity_z" ].in_units( "kpc/Gyr" )

   velocity_tot =np.array([np.sum(velocity_x),np.sum(velocity_y),np.sum(velocity_z)])/len(velocity_x)
   velocity_mag =((np.sum(velocity_x**2)+np.sum(velocity_y**2)+np.sum(velocity_z**2))/len(velocity_x))**0.5
   v = np.sum(velocity_tot**2)**0.5
   print(str(v)+',')
   Time.append(i*0.47)
   V_cluster.append(v)#*3.08/3.16)
   #print(velocity_mag)

plt.plot(Time, V_cluster, color ='blue')
plt.xlabel('Time(Gyr)')
plt.ylabel('V(km/s)')

plt.savefig('V_T.png')

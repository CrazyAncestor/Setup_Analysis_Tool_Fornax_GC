import yt
import numpy as np
from yt.data_objects.particle_filters import add_particle_filter
from matplotlib import pyplot as plt


idx_start = 0
idx_end = 25
filein  = [ 'Data_%06d'%idx for idx in range(idx_start, idx_end+1) ] 

nbin    = 50
dpi     = 150

# define the particle filter for the newly formed stars


for i in range(idx_start, idx_end+1):
   # load data
   ds = yt.load( filein[i] )
   
   
   def halo_star( pfilter, data ):
      filter = (data[ "all", "particle_mass" ] >6.28524e35)
      return filter

   
   add_particle_filter( "halo_star", function=halo_star, filtered_type="all", requires=["particle_mass"] )
   
   ds.add_particle_filter( "halo_star" )


   

   ad            = ds.all_data()
   centre_v = np.array([np.average(ad[ "halo_star", "particle_velocity_x" ].in_units( "kpc/Gyr" ) ),
                      np.average(ad[ "halo_star", "particle_velocity_y" ].in_units( "kpc/Gyr" ) ),
                      np.average(ad[ "halo_star", "particle_velocity_z" ].in_units( "kpc/Gyr" ) )])
   
   
   deviation_v = (np.sum((np.array(centre_v))**2))**0.5

   #print(centre)
   print(deviation_v)
   
   

import yt
import numpy as np
from yt.data_objects.particle_filters import add_particle_filter
from matplotlib import pyplot as plt


idx_start = 100
idx_end = 100
filein  = [ 'Data_%06d'%idx for idx in range(idx_start, idx_end+1) ] 

nbin    = 50
dpi     = 150

# define the particle filter for the newly formed stars


for i in range(len(filein)):
   # load data
   ds = yt.load( filein[i] )
   df = yt.load( filein[i] )
   def gc_star( pfilter, data ):
      filter = (data[ "all", "particle_mass" ] <7.096752e35)
      return filter
   

   add_particle_filter( "gc_star", function=gc_star, filtered_type="all", requires=["particle_mass"] )
   
   ds.add_particle_filter( "gc_star" )
   


   # get the mass and creation time of the new stars
   ad            = ds.all_data()
   x = np.array(ad[ "gc_star", "particle_position_x" ].in_units( "kpc" ) )
   y = np.array(ad[ "gc_star", "particle_position_y" ].in_units( "kpc" ) )
   z = np.array(ad[ "gc_star", "particle_position_z" ].in_units( "kpc" ) )
   vx        = np.array(ad[ "gc_star", "particle_velocity_x" ].in_units( "kpc/Gyr" ))*0.47
   vy        = np.array(ad[ "gc_star", "particle_velocity_y" ].in_units( "kpc/Gyr" ))*0.47
   vz        = np.array(ad[ "gc_star", "particle_velocity_x" ].in_units( "kpc/Gyr" ))*0.47
   
   for i in range(len(x)):
       print(x[i])
       print(y[i])
       print(z[i])
       print(vx[i])
       print(vy[i])
       print(vz[i])

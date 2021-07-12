import yt
import numpy as np
from yt.data_objects.particle_filters import add_particle_filter
from matplotlib import pyplot as plt


idx_start = 0
idx_end = 100
filein  = [ 'Data_%06d'%idx for idx in range(idx_start, idx_end+1) ] 

nbin    = 50
dpi     = 150

# define the particle filter for the newly formed stars


for i in range(idx_start, idx_end+1):
   # load data
   ds = yt.load( filein[i] )
   df = yt.load( filein[i] )
   def gc_star( pfilter, data ):
      filter = (data[ "all", "particle_mass" ] >1e38)
      return filter
   

   add_particle_filter( "gc_star", function=gc_star, filtered_type="all", requires=["particle_mass"] )
   
   ds.add_particle_filter( "gc_star" )
   


   # get the mass and creation time of the new stars
   ad            = ds.all_data()
   
   centre = np.array([np.average(ad[ "gc_star", "particle_position_x" ].in_units( "kpc" ) ),
                      np.average(ad[ "gc_star", "particle_position_y" ].in_units( "kpc" ) ),
                      np.average(ad[ "gc_star", "particle_position_z" ].in_units( "kpc" ) )])
   
   print(str(centre[0])+",")
   
   

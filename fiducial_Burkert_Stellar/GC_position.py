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
   df = yt.load( filein[i] )
   def gc_star( pfilter, data ):
      filter = (data[ "all", "particle_mass" ] <6.28524e35)
      return filter
   def halo_star( pfilter, data ):
      filter = (data[ "all", "particle_mass" ] >6.28524e35)
      return filter

   add_particle_filter( "gc_star", function=gc_star, filtered_type="all", requires=["particle_mass"] )
   add_particle_filter( "halo_star", function=halo_star, filtered_type="all", requires=["particle_mass"] )
   ds.add_particle_filter( "gc_star" )
   df.add_particle_filter( "halo_star" )


   # get the mass and creation time of the new stars
   ad            = ds.all_data()
   mass          = ad[ "gc_star", "particle_mass"    ].in_units( "g" )

   ad2            = df.all_data()
   centre = np.array([np.average(ad2[ "halo_star", "particle_position_x" ].in_units( "kpc" ) ),
                      np.average(ad2[ "halo_star", "particle_position_y" ].in_units( "kpc" ) ),
                      np.average(ad2[ "halo_star", "particle_position_z" ].in_units( "kpc" ) )])
   
   position_x        = ad[ "gc_star", "particle_position_x" ].in_units( "kpc" )
   position_y        = ad[ "gc_star", "particle_position_y" ].in_units( "kpc" )
   position_z        = ad[ "gc_star", "particle_position_z" ].in_units( "kpc" )
   position_ave =np.array([np.sum(position_x),np.sum(position_y),np.sum(position_z)])/len(position_x)

   
   r = np.sum((position_ave -centre)**2)**0.5
   print(centre)
   print(position_ave)
   print(r)
   
   

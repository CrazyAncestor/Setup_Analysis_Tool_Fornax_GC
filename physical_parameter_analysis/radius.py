import yt
import numpy as np
from yt.data_objects.particle_filters import add_particle_filter
from matplotlib import pyplot as plt


idx_start = 68
idx_end = 68
filein  = [ 'Data_%06d'%idx for idx in range(idx_start, idx_end+1) ] 

nbin    = 50
dpi     = 150

# define the particle filter for the newly formed stars


for i in range(0, len(filein)):
   # load data
   ds = yt.load( filein[i] )
   df = yt.load( filein[i] )
   def gc_star( pfilter, data ):
      filter = (data[ "all", "particle_mass" ] >7.1e35)
      return filter

   def halo_star( pfilter, data ):
      filter = (data[ "all", "particle_mass" ] <7.1e35)
      return filter

   add_particle_filter( "gc_star", function=gc_star, filtered_type="all", requires=["particle_mass"] )
   add_particle_filter( "halo_star", function=halo_star, filtered_type="all", requires=["particle_mass"] )

   ds.add_particle_filter( "gc_star" )
   df.add_particle_filter( "halo_star" )

   # get the mass and creation time of the new stars
   ad            = ds.all_data()
   af            = df.all_data()
   x_gc = np.array(ad[ "gc_star", "particle_position_x" ].in_units( "kpc" ) )
   y_gc = np.array(ad[ "gc_star", "particle_position_y" ].in_units( "kpc" ) )
   z_gc = np.array(ad[ "gc_star", "particle_position_z" ].in_units( "kpc" ) )
   centre_gc = np.array([np.average(ad[ "gc_star", "particle_position_x" ].in_units( "kpc" ) ),
                      np.average(ad[ "gc_star", "particle_position_y" ].in_units( "kpc" ) ),
                      np.average(ad[ "gc_star", "particle_position_z" ].in_units( "kpc" ) )])
   x_halo = np.array(af[ "halo_star", "particle_position_x" ].in_units( "kpc" ) )
   y_halo = np.array(af[ "halo_star", "particle_position_y" ].in_units( "kpc" ) )
   z_halo = np.array(af[ "halo_star", "particle_position_z" ].in_units( "kpc" ) )
   centre_halo = np.array([np.average(af[ "halo_star", "particle_position_x" ].in_units( "kpc" ) ),
                      np.average(af[ "halo_star", "particle_position_y" ].in_units( "kpc" ) ),
                      np.average(af[ "halo_star", "particle_position_z" ].in_units( "kpc" ) )])
   
   
   radius = np.sum((centre_gc-centre_halo)**2)**0.5

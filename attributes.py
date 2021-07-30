import yt
import numpy as np
from yt.data_objects.particle_filters import add_particle_filter
from matplotlib import pyplot as plt
import argparse
import sys

# load the command-line parameters
parser = argparse.ArgumentParser( description='Plot the gas slices and projections' )

parser.add_argument( '-p', action='store', required=False, type=str, dest='prefix',
                     help='path prefix [%(default)s]', default='../' )
parser.add_argument( '-s', action='store', required=True,  type=int, dest='idx_start',
                     help='first data index' )
parser.add_argument( '-e', action='store', required=True,  type=int, dest='idx_end',
                     help='last data index' )
parser.add_argument( '-d', action='store', required=False, type=int, dest='didx',
                     help='delta data index [%(default)d]', default=1 )

args=parser.parse_args()

idx_start   = args.idx_start
idx_end     = args.idx_end

filein  = [ '../Data_%06d'%idx for idx in range(idx_start, idx_end+1) ] 

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

   # write into files
   f = open("Attributes_%06d" %i)
   
   for i in range(len(x)):
      f.write(x[i])
      f.write(y[i])
      f.write(z[i])
      f.write(vx[i])
      f.write(vy[i])
      f.write(vz[i])
   f.close()

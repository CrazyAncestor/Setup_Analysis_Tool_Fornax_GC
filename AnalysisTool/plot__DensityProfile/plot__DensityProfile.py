import argparse
import sys
import yt
import matplotlib
matplotlib.use('agg') 
import matplotlib.pyplot as plt
from yt.data_objects.particle_filters import add_particle_filter
import numpy as np

# load the command-line parameters
parser = argparse.ArgumentParser( description='Plot the gas slices and projections' )

parser.add_argument( '-p', action='store', required=False, type=str, dest='prefix',
                     help='path prefix [%(default)s]', default='./' )
parser.add_argument( '-s', action='store', required=True,  type=int, dest='idx_start',
                     help='first data index' )
parser.add_argument( '-e', action='store', required=True,  type=int, dest='idx_end',
                     help='last data index' )
parser.add_argument( '-d', action='store', required=False, type=int, dest='didx',
                     help='delta data index [%(default)d]', default=1 )

args=parser.parse_args()

# take note
print( '\nCommand-line arguments:' )
print( '-------------------------------------------------------------------' )
for t in range( len(sys.argv) ):
   print str(sys.argv[t]),
print( '' )
print( '-------------------------------------------------------------------\n' )


idx_start   = args.idx_start
idx_end     = args.idx_end
didx        = args.didx
prefix      = args.prefix

center_mode = 'c'
dpi         = 150

delta_t = 2.9375 # Unit: Myr
MassThreshold = 1.989e+34 # Unit: g

yt.enable_parallelism()

d0= yt.load(  prefix+'Data_000000' )
files = [ prefix+'Data_%06d'%idx for idx in range(idx_start, idx_end+1, didx) ]

for i in range(len(files)):

   f = files[i]
   ds  = yt.load(f)
   
   def SelectedParticle( pfilter, data ):
      filter = (data[ "all", "particle_mass" ] >MassThreshold)
      return filter

   add_particle_filter( "SelectedParticle", function=SelectedParticle, filtered_type="all", requires=["particle_mass"] )
   ds.add_particle_filter( "SelectedParticle" )
   ad            = ds.all_data()
   center = 10.
   x = np.array(ad[ "SelectedParticle", "particle_position_x" ].in_units( "kpc" ) )
   y = np.array(ad[ "SelectedParticle", "particle_position_y" ].in_units( "kpc" ) )
   z = np.array(ad[ "SelectedParticle", "particle_position_z" ].in_units( "kpc" ) )
   center_position = np.ones(len(x))*center

   x -= center_position
   y -= center_position
   z -= center_position

   Par_R = (x**2+y**2+z**2)**0.5
   r = np.linspace(0.01,10,100)
   counts, radius = np.histogram(Par_R, bins=r)

   radius = (r[0:-1] + r[1:])/2.
   dens =  counts/radius**2

   fig = plt.figure()
   plt.plot(radius,dens,label=str(i*delta_t)+'Myr')

   plt.xscale('log')
   plt.yscale('log')

   plt.xlabel('radius(kpc)')
   plt.ylabel('particle density(Msun/kpc^3)')
   plt.legend()

plt.savefig(f+'_DensityProfile.png')


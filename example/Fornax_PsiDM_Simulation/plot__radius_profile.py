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

filein  = [ 'Data_%06d'%idx for idx in range(idx_start, idx_end+1) ] 

delta_t = 0.0152 # Unit: Gyr
MassThreshold = 0 # Unit: g
center_pos = 7.5

def GetRadiusAndTime(filein):

    r = []
    t = []

    for i in range(0, len(filein)):
        # load data
        ds = yt.load( filein[i] )
        df = yt.load( filein[i] )
        def gc_star( pfilter, data ):
            filter = (data[ "all", "particle_mass" ] >MassThreshold)
            return filter

        add_particle_filter( "gc_star", function=gc_star, filtered_type="all", requires=["particle_mass"] )
        ds.add_particle_filter( "gc_star" )

        # get the mass and creation time of the new stars
        ad            = ds.all_data()
        af            = df.all_data()
        x_gc = np.array(ad[ "gc_star", "particle_position_x" ].in_units( "kpc" ) )
        y_gc = np.array(ad[ "gc_star", "particle_position_y" ].in_units( "kpc" ) )
        z_gc = np.array(ad[ "gc_star", "particle_position_z" ].in_units( "kpc" ) )
        centre_gc = np.array([np.average(ad[ "gc_star", "particle_position_x" ].in_units( "kpc" ) ),
                            np.average(ad[ "gc_star", "particle_position_y" ].in_units( "kpc" ) ),
                            np.average(ad[ "gc_star", "particle_position_z" ].in_units( "kpc" ) )])
        
        centre_halo = np.array([center_pos,center_pos,center_pos])
        radius = np.sum((centre_gc-centre_halo)**2)**0.5
        r.append(radius)
        t.append(delta_t*i)

    return t,r

t,r = GetRadiusAndTime(filein)

plt.plot(t,r)
plt.xlabel('t (Gyr)')
plt.ylabel('r (kpc)')
plt.legend(loc = 'best')

plt.savefig('r_t.png')

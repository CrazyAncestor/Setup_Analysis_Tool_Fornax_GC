import argparse
import sys
import yt
from yt.data_objects.particle_filters import add_particle_filter
import numpy as np

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
didx        = args.didx
prefix      = args.prefix

colormap    = 'arbre'
dpi         = 150
field       = 'particle_mass'


yt.enable_parallelism()

ts = yt.load( [ prefix+'/Data_%06d'%idx for idx in range(idx_start, idx_end+1, didx) ] )
p = []
for idx in range(idx_start,idx_end+1,1):
   ds = yt.load('Data_%06d'%idx)
   def gc_star( pfilter, data ):
      filter = ((data[ "all", "particle_mass" ] <6.0e+33)&(data[ "all", "particle_mass" ] >5.0e+33))
      return filter
   add_particle_filter( "gc_star", function=gc_star, filtered_type="all", requires=["particle_mass"] )
   ds.add_particle_filter( "gc_star" )
   p = yt.ParticlePlot(ds,  ("gc_star","particle_position_x") ,  ("gc_star","particle_position_y") ,("gc_star","particle_mass") , center='c',width=1.0)
   
   p.set_background_color( ("gc_star","particle_mass") )
   p.set_unit('particle_mass', 'code_mass')
   p.set_axes_unit( 'code_length' )
   p.annotate_timestamp( time_unit='code_time', corner='upper_right', text_args={'color':'k'} )
   p.save( mpl_kwargs={"dpi":dpi} )
   p.save()

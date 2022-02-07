from __future__ import print_function
import argparse
import sys
import yt


# hard-coded parameters (in code units)
field     = 'density'
proj_axis = 'z'
center    = [75,75,75]
width_x   = 150.
width_y   = 150.
width_z   = 150.
#colormap  = 'arbre'
colormap  = 'inferno'
dpi       = 150


# load the command-line parameters
parser = argparse.ArgumentParser( description='Plot the gas slices and projections' )

parser.add_argument( '-i', action='store', required=False, type=str, dest='prefix_in',
                     help='input path prefix [%(default)s]', default='./' )
parser.add_argument( '-o', action='store', required=False, type=str, dest='prefix_out',
                     help='output filename prefix [%(default)s]', default='fig__thin-slice-projection' )
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
   print( str(sys.argv[t]), end=' ' )
print( '' )
print( '-------------------------------------------------------------------\n' )

idx_start  = args.idx_start
idx_end    = args.idx_end
didx       = args.didx
prefix_in  = args.prefix_in
prefix_out = args.prefix_out


yt.enable_parallelism()
ts = yt.load( [ prefix_in+'/Data_%06d'%idx for idx in range(idx_start, idx_end+1, didx) ] )

for ds in ts.piter():

#  create a region object with the target geometry
   corner_L = ds.arr( [center[0]-0.5*width_x,
                       center[1]-0.5*width_y,
                       center[2]-0.5*width_z], 'code_length' )
   corner_R = ds.arr( [center[0]+0.5*width_x,
                       center[1]+0.5*width_y,
                       center[2]+0.5*width_z], 'code_length' )
   region   = ds.box( corner_L, corner_R )


#  plot
   plot = yt.ProjectionPlot( ds, 'z', field, center='c', weight_field='density' )
   #plot = yt.ProjectionPlot( ds, proj_axis, field, data_source=region,
   #                          center=center, width=[width_x,width_y] )
   plot.set_unit( field, 'Msun/kpc**3' )
   plot.set_zlim( field, 1.0e4, 1.0e8)
#  plot.set_log( field, False )
   plot.set_cmap( field, colormap )
   plot.annotate_timestamp( time_unit='code_time', corner='upper_right' )

#  time = ds.current_time
#  plt.suptitle( "t = %6.2f %s"%(time.d, time.units) )

   plot.annotate_particles( width_z, data_source=region )

#  save the image
   plot.save( prefix_out+'_'+ds.basename+'.png', mpl_kwargs={'dpi':dpi} )


import yt
import numpy as np
from yt.data_objects.particle_filters import add_particle_filter
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
from locale import atof


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

delta_t = 0.47/4. # Unit: Gyr
MassThreshold = 1.989e+34 # Unit: g

for i in range(0, len(filein)):

   # load data
   ds = yt.load( filein[i] )
   
   def SelectedParticle( pfilter, data ):
      filter = (data[ "all", "particle_mass" ] <MassThreshold)
      return filter

   add_particle_filter( "SelectedParticle", function=SelectedParticle, filtered_type="all", requires=["particle_mass"] )

   ds.add_particle_filter( "SelectedParticle" )

   # get the mass and creation time of the new stars
   ad            = ds.all_data()
   
   # position
   x_halo = np.array(ad[ "SelectedParticle", "particle_position_x" ].in_units( "code_length" ) )
   y_halo = np.array(ad[ "SelectedParticle", "particle_position_y" ].in_units( "code_length" ) )
   z_halo = np.array(ad[ "SelectedParticle", "particle_position_z" ].in_units( "code_length" ) )
  
   centre_halo = np.array([np.average(ad[ "SelectedParticle", "particle_position_x" ].in_units( "code_length" ) ),
                      np.average(ad[ "SelectedParticle", "particle_position_y" ].in_units( "code_length" ) ),
                      np.average(ad[ "SelectedParticle", "particle_position_z" ].in_units( "code_length" ) )])
   center = np.ones(len(x_halo))
   x_halo = x_halo - center*centre_halo[0]
   y_halo = y_halo - center*centre_halo[1]
   z_halo = z_halo - center*centre_halo[2]

   # velocity
   vx_halo = np.array(ad[ "SelectedParticle", "particle_velocity_x" ].in_units( "code_length/code_time" ) )
   vy_halo = np.array(ad[ "SelectedParticle", "particle_velocity_y" ].in_units( "code_length/code_time" ) )
   vz_halo = np.array(ad[ "SelectedParticle", "particle_velocity_z" ].in_units( "code_length/code_time" ) )

   bulk_vel = np.array([np.average(ad[ "SelectedParticle", "particle_velocity_x" ].in_units( "code_length/code_time" ) ),
                      np.average(ad[ "SelectedParticle", "particle_velocity_y" ].in_units( "code_length/code_time" ) ),
                      np.average(ad[ "SelectedParticle", "particle_velocity_z" ].in_units( "code_length/code_time" ) )])
   bulk = np.ones(len(x_halo))
   vx_halo = vx_halo - bulk*bulk_vel[0]
   vy_halo = vy_halo - bulk*bulk_vel[1]
   vz_halo = vz_halo - bulk*bulk_vel[2]

   v_radial = (x_halo*vx_halo + y_halo*vy_halo + z_halo*vz_halo)
   r   = (x_halo*x_halo+y_halo*y_halo+z_halo*z_halo)**0.5
   v_radial = v_radial/r

   v_2 = (vx_halo*vx_halo + vy_halo*vy_halo + vz_halo*vz_halo)
   v_transverse = (-y_halo*vx_halo + x_halo*vy_halo)
   rplane = (x_halo*x_halo+y_halo*y_halo)**0.5
   v_transverse = v_transverse/rplane

   plt.figure()
   plt.hist(v_radial, bins=1000,range=(-0.2,0.2),label='Radial Velocity Distribution')
   plt.hist(v_transverse, bins=1000,range=(-0.2,0.2),label='Transverse Velocity Distribution')

   plt.legend(loc = 'best')
   print('ratio:'+str(np.var(v_radial)/np.var(v_transverse)))
   plt.savefig(filein[i]+'__VelocityHistogram.png')

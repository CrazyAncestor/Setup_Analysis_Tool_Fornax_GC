#include<time.h>
#include<particle.h>
#include<file.h>
void init_pars(particle GC_particles[],int size,double pos[][3],double vel[][3],double pos_tot[3]){
    for (int i = 0; i < size; i++) {
        GC_particles[i].init(pos[i], vel[i]);
        for (int j = 0; j < 3; j++) {
            pos_tot[j] += GC_particles[i].get_position(j);
        }
    }
    for (int j = 0; j < 3; j++) {
        pos_tot[j] /= size;
    }
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < 3; j++) {
            pos[i][j] -= pos_tot[j];
        }
        GC_particles[i].init(pos[i], vel[i]);
    }
}
int main()
{ 
    int size = file_size()/6;
    double pos[size][3];
    double vel[size][3];
    read_file(size,pos,vel);
    cout<<"file read"<<endl;
    particle GC_particles[size];
    
    //initialize particle numbers
    double pos_tot[3];
    init_pars(GC_particles,size,pos,vel,pos_tot);

    //record time
    clock_t t;
    t = clock();
    
    //calculate potential
    cout<<"add potential"<<endl;    
#   pragma omp parallel num_threads( 16 )
{
#   pragma omp for collapse(1)
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
#           pragma omp critical
            if (i != j)GC_particles[i].add_potential(&GC_particles[j],true);
            
        }
        
    }
}
    t = clock()-t;
    cout<<"duration:"<<t/CLOCKS_PER_SEC<<"sec"<<endl;

    //remove unbounded particles
    int id = 0;
    int round = 0;
    
    while (true) {
        
        int idl = 0;
        double v_div[3] = { 0,0,0 };
        for (int i = 0; i < size; i++) {

            if (!GC_particles[i].bound() & !GC_particles[i].kicked) {
                double v_dis[3];
                for (int j = 0; j < 3; j++) {
                    v_dis[j] = -GC_particles[i].get_velocity(j)/ (size - id-idl-1);
                }
                for (int j = 0; j < size; j++) {
                    GC_particles[j].add_potential(&GC_particles[i], false);
                    GC_particles[j].change_velocity(v_dis);
                }
                GC_particles[i].kick();
                idl++;
            }
            if (!GC_particles[i].kicked) {
                for (int j = 0; j < 3; j++) {
                    v_div[j] += GC_particles[i].get_velocity(j);
                }
            }
        }//for (int i = 0; i < size; i++)
        if (idl == 0)break;
        id += idl;
        
        for (int j = 0; j < 3; j++) {
            v_div[j] /= -(size-id);
        }//for (int j = 0; j < 3; j++)

        round++;
        
    }//while (true)

    //print out mass loss rate
    cout <<"mass loss rate:"<< float((id+0.0)/(size+0.0)) << endl;

    //print GC position
    double r_all[3];
    for (int i = 0; i < size; i++) {
        GC_particles[i].init(pos[i], vel[i]);
        for (int j = 0; j < 3; j++) {
            if(!GC_particles[i].kicked)r_all[j] += GC_particles[i].get_position(j)-10.;
        }//for (int j = 0; j < 3; j++)
    }//for (int i = 0; i < size; i++)
    double r_GC = pow(r_all[0] * r_all[0] + r_all[1] * r_all[1] + r_all[2] * r_all[2],0.5)/(size-id) ;
    cout << r_GC << endl;
}//int main()

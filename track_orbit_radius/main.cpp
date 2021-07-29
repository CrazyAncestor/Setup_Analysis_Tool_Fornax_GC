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

    clock_t t;
    t = clock();

    int cell_num = 100;
    double cell_size=1e-3;
    double center_pos[3]={5,5,5};
    int ***par_num = new int**[cell_num];
    for(int i=0;i<cell_num;i++){
        par_num[i] = new int*[cell_num];
        for(int j=0;j<cell_num;j++){
            par_num[i][j] = new int[cell_num];
        }//for(int j=0;j<cell_num;j++)
    }//for(int i=0;i<cell_num;i++)

    //initialize particle numbers in cells

    for(int i=0;i<cell_num;i++){
        for(int j=0;j<cell_num;j++){
            for(int k=0;k<cell_num;k++){
                par_num[i][j][k]=0;
            }//for(int k=0;k<cell_num;k++)
        }//for(int j=0;j<cell_num;j++)      
    }//for(int i=0;i<cell_num;i++)



    //record particle numbers in cells
    for (int i = 0; i < size; i++) {
        int x[3];
        for(int j=0;j<3;j++){
            x[j] = int(GC_particles[i].get_position(j)/cell_size)+cell_num/2;
        }//for(j=0;j<3;j++)
        if(x[0]>0 && x[1]>0 && x[2]>0 && x[0]<cell_num && x[1]<cell_num && x[2]<cell_num){
            par_num[x[0]][x[1]][x[2]]++;
        }//if(x[0]>0 && x[1]>0 && x[2]>0 && x[0]<cell_num && x[1]<cell_num && x[2]<cell_num)
    }//for (int i = 0; i < size; i++)

    
    int max_num=0;
    int arg_max_index[3];
    int n_tot = 0;

    for(int i=0;i<cell_num;i++){
        for(int j=0;j<cell_num;j++){
            for(int k=0;k<cell_num;k++){
                n_tot += par_num[i][j][k];
                if(max_num<=par_num[i][j][k]){
                    max_num = par_num[i][j][k];
                    arg_max_index[0] = i;
                    arg_max_index[1] = j;
                    arg_max_index[2] = k;
                }//if(max_num<=par_num[i][j][k])
                //cout<<par_num[i][j][k]<<' ';
            }//for(int k=0;k<cell_num;k++)
        }//for(int j=0;j<cell_num;j++)    
        //cout<<endl;  
    }//for(int i=0;i<cell_num;i++)


    cout<<arg_max_index[0]<<','<<arg_max_index[1]<<','<<arg_max_index[2]<<endl;
    cout<<max_num<<endl;
    cout<<n_tot<<endl;
    cout<<pos_tot[0]<<','<<pos_tot[1]<<','<<pos_tot[2]<<endl;
    double pos_final[3];
    double radius=0;
    for(int j=0;j<3;j++){
        pos_final[j]=pos_tot[j]+(arg_max_index[j]-cell_num/2)*cell_size;
        cout<<(pos_final[j]-center_pos[j])<<endl;
        radius += (pos_final[j]-center_pos[j])*(pos_final[j]-center_pos[j]);
    }
    cout<<radius<<endl;
    t = clock()-t;
    cout<<"duration:"<<double((t+0.0)/CLOCKS_PER_SEC)<<"sec"<<endl;
/*   
    clock_t t;
    t = clock();
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
    cout <<"mass loss rate:"<< float((id+0.0)/(size+0.0)) << endl;

    double r_all[3];
    for (int i = 0; i < size; i++) {
        GC_particles[i].init(pos[i], vel[i]);
        for (int j = 0; j < 3; j++) {
            if(!GC_particles[i].kicked)r_all[j] += GC_particles[i].get_position(j)-10.;
        }//for (int j = 0; j < 3; j++)
    }//for (int i = 0; i < size; i++)
    double r_GC = pow(r_all[0] * r_all[0] + r_all[1] * r_all[1] + r_all[2] * r_all[2],0.5)/(size-id) ;
    cout << r_GC << endl;*/
}//int main()

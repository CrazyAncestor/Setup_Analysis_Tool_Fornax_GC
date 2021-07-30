#include<time.h>
#include<particle.h>
#include<file.h>

//#define DEBUG
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
double GC_orbit_radius(const char* filename){
    int size = file_size(filename)/6;
    double pos[size][3];
    double vel[size][3];
    read_file(filename,size,pos,vel);
    #ifdef DEBUG
    cout<<"file read"<<endl;
    #endif
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
            }//for(int k=0;k<cell_num;k++)
        }//for(int j=0;j<cell_num;j++)    
    }//for(int i=0;i<cell_num;i++)

    #ifdef DEBUG
    cout<<arg_max_index[0]<<','<<arg_max_index[1]<<','<<arg_max_index[2]<<endl;
    cout<<max_num<<endl;
    cout<<n_tot<<endl;
    cout<<pos_tot[0]<<','<<pos_tot[1]<<','<<pos_tot[2]<<endl;
    #endif
    double pos_final[3];
    double radius=0;
    for(int j=0;j<3;j++){
        pos_final[j]=pos_tot[j]+(arg_max_index[j]-cell_num/2)*cell_size;
        #ifdef DEBUG
        cout<<(pos_final[j]-center_pos[j])<<endl;
        #endif
        radius += (pos_final[j]-center_pos[j])*(pos_final[j]-center_pos[j]);
    }

    t = clock()-t;
    #ifdef DEBUG
    cout<<radius<<endl;
    cout<<"duration:"<<double((t+0.0)/CLOCKS_PER_SEC)<<"sec"<<endl;
    #endif
    return radius;
}
int main(int argc, char *argv[]) {
    printf("We have %d arguments:\n", argc);
    for (int i = 0; i < argc; ++i) {
        printf("[%d] %s\n", i, argv[i]);
    }

    fstream file;
    int idx_start=atoi(argv[1]);
    int idx_end  =atoi(argv[2]);
    file.open("r.txt",ios::out);
    for(int idx=idx_start;idx<idx_end+1;idx++){
	    char filename[100];
	    sprintf(filename,"../Attributes_%06d",idx);
        file<<GC_orbit_radius(filename)<<endl;
    }
    file.close();
}//int main()

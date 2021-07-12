//Statistics
#include<stdlib.h>
#include<iostream>
#include<fstream>
using namespace std;
#include<time.h>
#define PI 3.1415926535897
#include<math.h>
#include<omp.h>

//Physics Constant
double G = 1.0; // Gravitational Constant
double m = 3.6620343e-05; // particle mass
int file_size() {
    fstream file;
    file.open("Attributes", ios::in);        //�N�ɮ׶}�Ҭ���J���A

    if (!file)     //�ˬd�ɮ׬O�_���\�}��
    {
        cerr << "Can't open file!\n";
        exit(1);     //�b�����`���ΤU�A���_�{��������
    }
    string a;
    int id = 0;
    while (file >> a) {     //Ū���O���A�YŪ�����ɮ׵����h�Ǧ^0
        double s = atof(a.c_str());
        id++;
    }

    file.close();
    return id;
}
class particle
{
    
    double pos[3];
    double vel[3];
   
    double energy;
public:
    bool kicked;
    //double pos[3];
    //double vel[3];
    void init(double p[3], double v[3]) {
        energy = 0;
        for (int i = 0; i < 3; i++) {
            pos[i] = p[i];
            vel[i] = v[i];
            energy += m*(vel[i] * vel[i])/2.; //kinetic energy
        }
        
    }

    void display() {
        cout << "pos: (" << pos[0] << "," << pos[1] << "," << pos[2] << ")" << endl;
        cout << "vel: (" << vel[0] << "," << vel[1] << "," << vel[2] << ")" << endl;
        cout << "energy: " << energy  << endl;
    }
    
    void add_potential(particle* p,bool plus) {
        double r = 0;
        for (int i = 0; i < 3; i++) {
            r += pow((this->pos[i]-p->pos[i]),2);
        }
        r = pow(r, 0.5);
        if(plus)energy += -G * m * m / r;
        else energy -= -G * m * m / r;
        
    }
    double get_velocity(int i) {
        return vel[i];
    }
     void change_velocity(double dv[3]) {
         double potential = energy - 0.5 * m * (vel[0] * vel[0] + vel[1] * vel[1] + vel[2] * vel[2]);
         for (int i = 0; i < 3; i++) {
             vel[i] += dv[i];
         }
         energy = potential + 0.5 * m * (vel[0] * vel[0] + vel[1] * vel[1] + vel[2] * vel[2]);
    }
    bool bound() {
        if (this->energy < 0) return true;
        else  return false;
    }
    void kick() {

        kicked = true;
    }
    double get_position(int i) {
        return pos[i];
    }

    
};


int main()
{ 
    int size = file_size()/6;
    double pos[size][3];
    double vel[size][3];
    fstream file;

    file.open("Attributes", ios::in);        //�N�ɮ׶}�Ҭ���J���A

    if (!file)     //�ˬd�ɮ׬O�_���\�}��
    {
        cerr << "Can't open file!\n";
        exit(1);     //�b�����`���ΤU�A���_�{��������
    }
    string a;
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < 3; j++) {
            file >> a;
            pos[i][j] = atof(a.c_str());
        }
        for (int j = 0; j < 3; j++) {
            file >> a;
            vel[i][j] = atof(a.c_str());
        }
    }
    file.close();
    
    particle GC_particles[size];
    
    double v_tot[3];
    for (int i = 0; i < size; i++) {
        GC_particles[i].init(pos[i], vel[i]);
        for (int j = 0; j < 3; j++) {
            v_tot[j] += GC_particles[i].get_velocity(j);
        }
    }
    for (int j = 0; j < 3; j++) {
        v_tot[j] /= size;
    }
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < 3; j++) {
            vel[i][j] -= v_tot[j];
        }
        GC_particles[i].init(pos[i], vel[i]);
    }

    clock_t t;
    t = clock();
    
#   pragma omp parallel num_threads( 8 )
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
                //GC_particles[i].display();
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
        }
        if (idl == 0)break;
        id += idl;
        
        for (int j = 0; j < 3; j++) {
            v_div[j] /= -(size-id);
        }
        /*for (int i = 0; i < size; i++) {
            GC_particles[i].change_velocity(v_div);
        }*/
        //cout << "round" << round << endl;
       // cout << v_div[0] << endl;
       // cout << v_div[1] << endl;
       // cout << v_div[2] << endl;
        round++;
        
    }
    cout <<"mass loss rate:"<< float((id+0.0)/(size+0.0)) << endl;
    //cout << size << endl;

    double r_all[3];
    for (int i = 0; i < size; i++) {
        GC_particles[i].init(pos[i], vel[i]);
        for (int j = 0; j < 3; j++) {
            if(!GC_particles[i].kicked)r_all[j] += GC_particles[i].get_position(j)-10.;
        }
    }
    double r_GC = pow(r_all[0] * r_all[0] + r_all[1] * r_all[1] + r_all[2] * r_all[2],0.5)/(size-id) ;
    cout << r_GC << endl;
}

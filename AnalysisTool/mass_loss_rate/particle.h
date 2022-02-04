#ifndef PARTICLE_H
#define PARTICLE_H

#include<stdlib.h>
#include<iostream>
using namespace std;
#include<math.h>


//Physics Constant
#define PI 3.1415926535897
#define G 1.0 // Gravitational Constant
#define m 3.6620343e-05 // particle mass

#ifdef OMP_PARALLEL
#include<omp.h>
#endif

using namespace std;

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


#endif // #ifndef PARTICLE_H
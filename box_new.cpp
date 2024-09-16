#include <iostream>

#include <cmath>
#include <vector>
#include <fstream>
#include "particle.h"
#include "random.h"

using namespace std;

int main(){
    // Parameters
    double q = 0.5;
    double k = 15.0, Dt = 1.0, Dr = 10.0, l0 = 0.0, omega = 50.0;
    double Lx = 100.0; double Ly = 100.0;

    //Intializing particles
    particle p1, p2;
    p1.x = 50.0; p1.y = 50.0; p1.phi = 0.0;
    p2.x = 45.0; p2.y = 45.0;

    // Other variables
    double t = 0.0; double tmax = 3000000.0; double dt = 0.0001; int t_step = 10000; int bins = 25; double bin_size = Lx/bins;
    double dx, dy, theta, l, Fx, Fy, eta_x1, eta_y1, eta_phi, eta_x2, eta_y2;

    // Averaging variables
    int i=0; double xc; int count[25]; int j;
    for(int k = 0; k<25; k++){
        count[k] = 0;
    }

    // Opening files
    ofstream out1, out2, out3;
    out1.open("particle1.dat");
    out2.open("particle2.dat");
    out3.open("pdf_omega_large_50.dat");

    // Time loop
    seed_drand48();
    while(t<tmax){
        i++;
        t = t+dt;
        
        //theta = angle(p1,p2);
        l = length(p1, p2);

        Fx = -k*(l-l0)*(p1.x - p2.x)/l;
        Fy = -k*(l-l0)*(p1.y - p2.y)/l;

        eta_x1 = Gaussian();
        eta_y1 = Gaussian();
        eta_phi = Gaussian();
        eta_x2 = Gaussian();
        eta_y2 = Gaussian();

        p1.x += dt*(Fx + fs(p1.x, Lx)*cos(p1.phi)) + sqrt(2*Dt*dt)*eta_x1; p1.xa = p1.x;
        p1.y += dt*(Fy + fs(p1.x, Lx)*sin(p1.phi)) + sqrt(2*Dt*dt)*eta_y1; p1.ya = p1.y;
        p1.phi += omega*dt + sqrt(2*Dr*dt)*eta_phi;
        p2.x += -dt*Fx/q + sqrt(2*Dt*dt/q)*eta_x2; //p2.xa = p2.x;
        p2.y += -dt*Fy/q + sqrt(2*Dt*dt/q)*eta_y2; //p2.ya = p2.y;

        // Periodic Boundary Condition
        /*while(p1.xa<0) {p1.xa += Lx;}
        while(p1.xa>Lx) {p1.xa -= Lx;}
        while(p1.ya<0) {p1.ya += Ly;}
        while(p1.ya>Lx) {p1.ya -= Ly;}
        while(p2.xa<0) {p2.xa += Lx;}
        while(p2.xa>Lx) {p2.xa -= Lx;}
        while(p2.ya<0) {p2.ya += Ly;}
        while(p2.ya>Lx) {p2.ya -= Ly;}*/

        //out1<<p1.x<<"\t"<<p1.y<<endl;
        //out2<<p2.x<<"\t"<<p2.y<<endl;

        if(i%t_step == 0){
            xc = (p1.x + q*p2.x)/(1+q);
            while(xc<0) {xc += Lx;}
            while(xc>Lx) {xc -= Lx;}
            j = int(xc/4);
            count[j]++;
        }
    }

    for(int k=0; k<25; k++){
        out3<<4*k + 2<<"\t"<<double((count[k]*dt*t_step*Lx)/(tmax*bin_size))<<endl;
    }
}

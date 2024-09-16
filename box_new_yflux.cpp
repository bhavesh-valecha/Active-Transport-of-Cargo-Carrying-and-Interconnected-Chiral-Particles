#include <iostream>

#include <cmath>
#include <vector>
#include <fstream>
#include "particle.h"
#include "random.h"
#include "activity.h"

using namespace std;

int main(){
    // Parameters
    double q = 2;
    double k = 15.0, Dt = 1.0, Dr = 10.0, l0 = 0.0, omega = 12.0;
    double Lx = 100.0; double Ly = 100.0;

    //Intializing particles
    particle p1, p2;
    p1.x = 50.0; p1.y = 50.0; p1.phi = 0.0;
    p2.x = 45.0; p2.y = 45.0;

    // Other variables
    double t = 0.0; double tmax = 30000000.0; double dt = 0.0001; double jy_val = 0.0; int t_step = 1000; int bins = 25; double bin_size = Lx/bins; double rho_bulk = 0.0;
    double dx, dy, theta, l, Fx, Fy, eta_x1, eta_y1, eta_phi, eta_x2, eta_y2, eta_j;
    // Averaging variables
    int i=0; double xc; double xb; int count[bins]; int j;int jy; double countjy[bins];
    for(int k = 0; k<bins; k++){
	count[k] = 0;
	countjy[k] = 0;
    }

    // Opening files
    ofstream out1, out2, out3;
    out1.open("particle1.dat");
    out2.open("particle2.dat");
    out3.open("pdf_omega_12_yflux_25bins_largeact_onewave_final_long1.dat");

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

	p1.x += dt*(Fx + fs(p1.x, Lx)*cos(p1.phi)) + sqrt(2*Dt*dt)*eta_x1;
        p1.y += dt*(Fy + fs(p1.x, Lx)*sin(p1.phi)) + sqrt(2*Dt*dt)*eta_y1;
        p1.phi += omega*dt + sqrt(2*Dr*dt)*eta_phi;
        p2.x += -dt*Fx/q + sqrt(2*Dt*dt/q)*eta_x2;
        p2.y += -dt*Fy/q + sqrt(2*Dt*dt/q)*eta_y2;

        // Periodic Boundary Condition

        if(i%t_step==0){
            xc = (p1.x + q*p2.x)/(1+q);
            while(xc<0) {xc += Lx;}
            while(xc>Lx) {xc -= Lx;}
            j = int(xc/bin_size);
            count[j]++;
            countjy[j] += (1/(1+q))*fs(p1.x, Lx)*sin(p1.phi) + (1/(1+q))*sqrt(2*Dt*dt)*eta_y1 + (q/(1+q))*sqrt(2*Dt*dt/q)*eta_y2;
        }
    }

    for(int l=0; l<bins; l++){
        rho_bulk += count[l]*dt*t_step/tmax;
}

    for(int k=0; k<bins; k++){
        out3<<bin_size*k <<"\t"<<double((countjy[k]*dt*t_step*Lx)/(tmax*bin_size*rho_bulk))<<endl;
    }
}

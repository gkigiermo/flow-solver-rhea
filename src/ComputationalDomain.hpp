#ifndef _ComputationalDomain_
#define _ComputationalDomain_
#include <iostream>
#include "MacroParameters.hpp"
using namespace std;

class ComputationalDomain{
    public:
        ComputationalDomain(){};
        ComputationalDomain(double, double, double,double, double,double, double, double, double, int, int, int);
        
        void printDomain();
        void setBocos(int b[6]);
        int getGNx(){return gNx;};
        int getGNy(){return gNy;};
        int getGNz(){return gNz;};
        double getGlobx(int id){return globx[id];};
        double getGloby(int id){return globy[id];};
        double getGlobz(int id){return globz[id];};

        void set_x(int id,double val){ x[id]=val;};
        void set_y(int id,double val){ y[id]=val;};
        void set_z(int id,double val){ z[id]=val;};

        int getBoco(int bocoid){return bc[bocoid];};
        void calculateLocalGrid(int lNx,int lNy, int lNz);
        
        //Local positions
        double *x;
        double *y;
        double *z;


    protected:
        //Global Number of cells in each direction
        int gNx;
        int gNy;
        int gNz;

        //Cell dimensions
        double *gdx;
        double *gdy;
        double *gdz;

        //Local number of cells in each direction
        int Nx;
        int Ny;
        int Nz;

        //Local cell dimensions
        double *dx;
        double *dy;
        double *dz;

        //Boundary conditions

        int bc[6];

        //Positions globally
        double *globx;
        double *globy;
        double *globz;

        //Dimensions of the ComputationalDomain
        double L_x;
        double L_y;
        double L_z;

        //Origin 
        double x_0;
        double y_0;
        double z_0;

        // Stretching factors
        double A_x;
        double A_y;
        double A_z;

        void calculateGlobalGrid();
};


#endif

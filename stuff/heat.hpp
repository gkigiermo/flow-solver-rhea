#include <stdio.h>
#include <iostream>
#include <mpi.h>
#include "src/parameters.h"
#include "src/domain.h"
#include "src/comm_scheme.h"
#include "src/parvec.h"

int main(int argc, char** argv)
{

    MPI_Init(&argc, &argv);

    int _DEBUG_ = 0; /* This should be a parameter of the input file*/
    double Lx,Ly,Lz;
    Lx = 1.0;
    Ly = 1.0;
    Lz = 1.0;
    double origen_x,origen_y,origen_z;
    origen_x = 0.0;
    origen_y = 0.0;
    origen_z = 0.0;
   

    double kappa=0.1;
    double dt = 0.001;

    int bocos[6];
    bocos[_WEST_]  = _DIRICHLET_;
    bocos[_EAST_]  = _DIRICHLET_;
    bocos[_SOUTH_] = _NEUMANN_;
    bocos[_NORTH_] = _NEUMANN_;
    bocos[_BACK_]  = _NEUMANN_;
    bocos[_FRONT_] = _NEUMANN_;

    //This is just for testing, later would be read from a file
    domain dom(Lx, Ly, Lz, origen_x, origen_y, origen_z, RHEA_NX, RHEA_NY, RHEA_NZ);

    //To add the boundary conditions,  later would be read from a file
    dom.setBocos(bocos);

     //Number of procs in x,y,z
    int npx, npy, npz;
    npx=2;
    npy=1;
    npz=1;

    comm_scheme topo(&dom, npx, npy, npz);

    if(topo.getRank() == 0)
        dom.printDomain();

   for(int p = 0; p < npx*npy*npz; p++)
       topo.printCommSchemeToFile(p);

    parvec T(&topo);
    parvec Tnew(&topo);


    int _lNx_ = topo.getlNx();
    int _lNy_ = topo.getlNy();
    int _lNz_ = topo.getlNz();


    T = 0.0;
    Tnew = 0.0;

    int maxiter= 20000;


    for(int it=0; it < maxiter; it++)
    {

        for(int i = Tnew.ini_x; i <=  Tnew.fin_x; i++)
            for(int j = Tnew.ini_y; j <=  Tnew.fin_y; j++)
                for(int k = Tnew.ini_z; k <=  Tnew.fin_z; k++){
                    Tnew[I1D(i,j,k)] =  T[I1D(i,j,k)] + dt*kappa*( (2.0/(dom.x[i+1]-dom.x[i-1]))*( (T[I1D(i+1,j,k)] -T[I1D(i,j,k)])/(dom.x[i+1]-dom.x[i]) -  (T[I1D(i,j,k)] -T[I1D(i-1,j,k)])/(dom.x[i]-dom.x[i-1])) + 
                                                                   (2.0/(dom.y[j+1]-dom.y[j-1]))*( (T[I1D(i,j+1,k)] -T[I1D(i,j,k)])/(dom.y[j+1]-dom.y[j]) -  (T[I1D(i,j,k)] -T[I1D(i,j-1,k)])/(dom.y[j]-dom.y[j-1])) +        
                                                                   (2.0/(dom.z[k+1]-dom.z[k-1]))*( (T[I1D(i,j,k+1)] -T[I1D(i,j,k)])/(dom.z[k+1]-dom.z[k]) -  (T[I1D(i,j,k)] -T[I1D(i,j,k-1)])/(dom.z[k]-dom.z[k-1])) );    
                }
    
        Tnew.update();
        //WEST 
        for(int i = topo.iter_bound[_WEST_][_INIX_]; i <= topo.iter_bound[_WEST_][_ENDX_]; i++)
            for(int j = topo.iter_bound[_WEST_][_INIY_]; j <= topo.iter_bound[_WEST_][_ENDY_]; j++)
                for(int k = topo.iter_bound[_WEST_][_INIZ_]; k <=  topo.iter_bound[_WEST_][_ENDZ_]; k++)
                {
                    Tnew[I1D(i,j,k)] = 2*0.0 - Tnew[I1D(i+1,j,k)];
                }

        //EAST 
        for(int i = topo.iter_bound[_EAST_][_INIX_]; i <= topo.iter_bound[_EAST_][_ENDX_]; i++)
            for(int j = topo.iter_bound[_EAST_][_INIY_]; j <= topo.iter_bound[_EAST_][_ENDY_]; j++)
                for(int k = topo.iter_bound[_EAST_][_INIZ_]; k <=  topo.iter_bound[_EAST_][_ENDZ_]; k++)
                {
                    Tnew[I1D(i,j,k)] = 2*1.0 - Tnew[I1D(i-1,j,k)];
                }

        //SOUTH 
        for(int i = topo.iter_bound[_SOUTH_][_INIX_]; i <= topo.iter_bound[_SOUTH_][_ENDX_]; i++)
            for(int j = topo.iter_bound[_SOUTH_][_INIY_]; j <= topo.iter_bound[_SOUTH_][_ENDY_]; j++)
                for(int k = topo.iter_bound[_SOUTH_][_INIZ_]; k <=  topo.iter_bound[_SOUTH_][_ENDZ_]; k++)
                {
                    Tnew[I1D(i,j,k)] = Tnew[I1D(i,j+1,k)];
                }

        //NORTH 
        for(int i = topo.iter_bound[_NORTH_][_INIX_]; i <= topo.iter_bound[_NORTH_][_ENDX_]; i++)
            for(int j = topo.iter_bound[_NORTH_][_INIY_]; j <= topo.iter_bound[_NORTH_][_ENDY_]; j++)
                for(int k = topo.iter_bound[_NORTH_][_INIZ_]; k <=  topo.iter_bound[_NORTH_][_ENDZ_]; k++)
                {
                    Tnew[I1D(i,j,k)] = Tnew[I1D(i,(j-1),k)];
                }

        //BACK 
        for(int i = topo.iter_bound[_BACK_][_INIX_]; i <= topo.iter_bound[_BACK_][_ENDX_]; i++)
            for(int j = topo.iter_bound[_BACK_][_INIY_]; j <= topo.iter_bound[_BACK_][_ENDY_]; j++)
                for(int k = topo.iter_bound[_BACK_][_INIZ_]; k <=  topo.iter_bound[_BACK_][_ENDZ_]; k++)
                {
//                    Tnew[I1D(i,j,0)] = Tnew[I1D(i,j,_lNz_)];
                     Tnew[I1D(i,j,k)] = Tnew[I1D(i,j,k+1)];
                }

        //FRONT 
        for(int i = topo.iter_bound[_FRONT_][_INIX_]; i <= topo.iter_bound[_FRONT_][_ENDX_]; i++)
            for(int j = topo.iter_bound[_FRONT_][_INIY_]; j <= topo.iter_bound[_FRONT_][_ENDY_]; j++)
                for(int k = topo.iter_bound[_FRONT_][_INIZ_]; k <=  topo.iter_bound[_FRONT_][_ENDZ_]; k++)
                {
//                   Tnew[I1D(i,j,_lNz_+1)] = Tnew[I1D(i,j,1)];
                    Tnew[I1D(i,j,k)] = Tnew[I1D(i,j,k-1)];
               }



//        for(int i = topo.iter_common[_ALL_][_INIX_]; i <= topo.iter_common[_ALL_][_ENDX_]; i++)
//            for(int j =  topo.iter_common[_ALL_][_INIY_]; j <=  topo.iter_common[_ALL_][_ENDY_]; j++)
//                for(int k =  topo.iter_common[_ALL_][_INIZ_]; k <=  topo.iter_common[_ALL_][_ENDZ_]; k++){
//                    T[I1D(i,j,k)] = Tnew[I1D(i,j,k)];
//                }
 
        for(int l=0; l < topo.getSize(); l++)
        {
            T[l] = Tnew[l];
        }
 
    }



    if(topo.getRank() == 0 && _DEBUG_ ==0 )
        for(int l=0; l < topo.getSize(); l++)
            cout<<l<<":  "<<Tnew[l]<<endl;



    MPI_Finalize();
}

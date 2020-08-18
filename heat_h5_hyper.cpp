#include <stdio.h>
#include <iostream>
#include <mpi.h>
#include <hdf5.h>
//#include "h5Cpp.h"
#include "src/parameters.h"
#include "src/domain.h"
#include "src/comm_scheme.h"
#include "src/parvec.h"
#include "src/printer.h"

using namespace std;
//using namespace H5;

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
    double stretch_x,stretch_y,stretch_z;
    stretch_x = 0.0;
    stretch_y = 0.0;
    stretch_z = 0.0;   

    double kappa=0.1;
    double dt = 0.001;

    int bocos[6];
    bocos[_WEST_]  = _DIRICHLET_;
    bocos[_EAST_]  = _DIRICHLET_;
    bocos[_SOUTH_] = _NEUMANN_;
    bocos[_NORTH_] = _NEUMANN_;
    bocos[_BACK_]  = _PERIODIC_;
    bocos[_FRONT_] = _PERIODIC_;

    //This is just for testing, later would be read from a file
    domain dom(Lx, Ly, Lz, origen_x, origen_y, origen_z, stretch_x, stretch_y, stretch_z, RHEA_NX, RHEA_NY, RHEA_NZ);

    //To add the boundary conditions,  later would be read from a file
    dom.setBocos(bocos);

     //Number of procs in x,y,z
    int npx, npy, npz;
    npx=1;
    npy=1;
    npz=1;

    comm_scheme topo(&dom, npx, npy, npz);

    if(topo.getRank() == 0)
        dom.printDomain();

   for(int p = 0; p < npx*npy*npz; p++)
       topo.printCommSchemeToFile(p);


    /* Note that now the parallel vector also include a name as an input */
    /* That name is used to creating the output files */

    parvec T(&topo,"T");
    parvec Tnew(&topo,"Tnew");


    int _lNx_ = topo.getlNx();
    int _lNy_ = topo.getlNy();
    int _lNz_ = topo.getlNz();


    T = 0.0;
    Tnew = 0.0;
    cout<<"Before update "<<endl;
    Tnew.update();
    cout<<"After update "<<endl;
    int maxiter= 20000;


    parvec meshx(&topo,"X");

    for(int i = topo.iter_common[_ALL_][_INIX_]; i <= topo.iter_common[_ALL_][_ENDX_]; i++)
       for(int j =  topo.iter_common[_ALL_][_INIY_]; j <=  topo.iter_common[_ALL_][_ENDY_]; j++)
           for(int k =  topo.iter_common[_ALL_][_INIZ_]; k <=  topo.iter_common[_ALL_][_ENDZ_]; k++){
               meshx[I1D(i,j,k)] = dom.x[i];
           }
 
    /* Printer Class */
    /* It needs to be relabeled */
    
    /* Has input needs the topo and the name of the outputs*/
    /* The topo is used to obtain the info of the hyperslab */
    /* The name is used to create variable output names using the it number */
    printer output(&topo,"salida");

    /* By default it doesn't store anything*/
    /* Fields that you want to store need to be added explicitely only once */
    /* Just give a pointer of the field as a parameter */
    /* Those fields will be written or read acording to the function used */
    output.addField(&T);
    output.addField(&Tnew);
    output.addField(&meshx);


    /* This is an example on how it works the read */
    /* note that the output salida_12500.h5 should exist in order to work */
    /* if not the code will crash */
    /* This will only be needed when a restart is used (some if should be included in flowRhea*/
    //output.read(12500);
    
    /* Here I test printing the read values into a new file with name salida_9.h5 */
    /* You can verify that both files 12500 and 9 should be the same*/
    // output.write(9);


    for(int it=0; it < maxiter; it++)
    {

        for(int i = Tnew.ini_x; i <=  Tnew.fin_x; i++)
            for(int j = Tnew.ini_y; j <=  Tnew.fin_y; j++)
                for(int k = Tnew.ini_z; k <=  Tnew.fin_z; k++){
                    Tnew[I1D(i,j,k)] =  T[I1D(i,j,k)] + dt*kappa*( (2.0/(dom.x[i+1]-dom.x[i-1]))*( (T[I1D(i+1,j,k)] -T[I1D(i,j,k)])/(dom.x[i+1]-dom.x[i]) -  (T[I1D(i,j,k)] -T[I1D(i-1,j,k)])/(dom.x[i]-dom.x[i-1])) + 
                                                                   (2.0/(dom.y[j+1]-dom.y[j-1]))*( (T[I1D(i,j+1,k)] -T[I1D(i,j,k)])/(dom.y[j+1]-dom.y[j]) -  (T[I1D(i,j,k)] -T[I1D(i,j-1,k)])/(dom.y[j]-dom.y[j-1])) +        
                                                                   (2.0/(dom.z[k+1]-dom.z[k-1]))*( (T[I1D(i,j,k+1)] -T[I1D(i,j,k)])/(dom.z[k+1]-dom.z[k]) -  (T[I1D(i,j,k)] -T[I1D(i,j,k-1)])/(dom.z[k]-dom.z[k-1])) );    
                }
        
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

        Tnew.update();

 
        for(int l=0; l < topo.getSize(); l++)
        {
            T[l] = Tnew[l];
        }

        /* Suposing you want to save the fields every 2500 steps */
        /* Input of write is the iteration number than later is used to create the name */
        if(it%2500 ==0)
          output.write(it,dt);

 
    }

    if(topo.getRank() == 2 && _DEBUG_ ==0 )
        for(int l=0; l < topo.getSize(); l++)
        { 
            if(l%4 == 0){
                cout<<endl;
            }
            if(l%16 == 0)
                cout<<"------------------------"<<endl;
           cout<<" "<<Tnew[l];

        }
    
    MPI_Barrier(MPI_COMM_WORLD);

    MPI_Finalize();
}

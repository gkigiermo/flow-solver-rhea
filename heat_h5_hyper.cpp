#include <stdio.h>
#include <iostream>
#include <mpi.h>
#include <hdf5.h>
//#include "h5Cpp.h"
#include "src/parameters.h"
#include "src/domain.h"
#include "src/comm_scheme.h"
#include "src/parvec.h"

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
    domain dom(Lx, Ly, Lz, origen_x, origen_y, origen_z, RHEA_NX, RHEA_NY, RHEA_NZ);

    //To add the boundary conditions,  later would be read from a file
    dom.updateBocos(bocos);

     //Number of procs in x,y,z
    int npx, npy, npz;
    npx=2;
    npy=2;
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
    cout<<"Before update "<<endl;
    Tnew.update();
    cout<<"After update "<<endl;
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



   parvec meshx(&topo);

   for(int i = topo.iter_common[_ALL_][_INIX_]; i <= topo.iter_common[_ALL_][_ENDX_]; i++)
       for(int j =  topo.iter_common[_ALL_][_INIY_]; j <=  topo.iter_common[_ALL_][_ENDY_]; j++)
           for(int k =  topo.iter_common[_ALL_][_INIZ_]; k <=  topo.iter_common[_ALL_][_ENDZ_]; k++){
               meshx[I1D(i,j,k)] = dom.x[i];
           }
 
  /*  int lolo[10];

    hid_t plist_id;
    plist_id = H5Pcreate(H5P_FILE_ACCESS);

    MPI_Info info  = MPI_INFO_NULL;

    H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD, info);

    plist_id = H5Pcreate (H5P_DATASET_XFER);

    H5Pset_dxpl_mpio (plist_id, H5FD_MPIO_COLLECTIVE);

    status = H5Dwrite (dset_id, H5T_NATIVE_INT, memspace, filespace, plist_id, lolo );

*/

    char filename[100];

    int rank = topo.getRank();

  //sprintf(filename,"filetest-%d.h5",rank);
    MPI_Comm comm  = MPI_COMM_WORLD;
    MPI_Info info  = MPI_INFO_NULL;
    MPI_Info_create(&info);
    sprintf(filename,"hyperslab.h5");

    cout<<"writing  "<<filename<<endl;
    hid_t       file_id;   /* file identifier */
    herr_t      status;

    /* Create a new file using default properties. */
    hid_t fa_plist_id, fa_plist_id2, dx_plist_id;

    fa_plist_id = H5Pcreate(H5P_FILE_ACCESS);

    herr_t ret =  H5Pset_fapl_mpio(fa_plist_id, comm , info);

    cout<<" ret " <<ret<<" "<<-1<<endl;

    file_id = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, fa_plist_id);

    //Collective call

    fa_plist_id2 = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(fa_plist_id2, H5FD_MPIO_COLLECTIVE);

    int num_dims=3;
	hsize_t dim_gsize[3];               // dataset dimensions
	hsize_t dim_lsize[3];               // dataset dimensions
    hsize_t dim_offset[3];
	dim_gsize[0] = dom.getGNz() ; //  _lNz_;
    dim_gsize[1] = dom.getGNy() + 2;
    dim_gsize[2] = dom.getGNx() + 2;

	dim_lsize[0] = topo.lenslabz-2; //  _lNz_;
    dim_lsize[1] = topo.lenslaby;
    dim_lsize[2] = topo.lenslabx;

	dim_offset[0] = topo.offslab_z; //  _lNz_;
    dim_offset[1] = topo.offslab_y;
    dim_offset[2] = topo.offslab_x;


    cout<<"local offset "<<dim_offset[0]<<" "<<dim_offset[1]<<" "<<dim_offset[2]<<endl;
    cout<<"local sizes "<<dim_lsize[0]<<" "<<dim_lsize[1]<<" "<<dim_lsize[2]<<endl;


    hsize_t stride[3];
    hsize_t block[3];

    stride[0] = 1;
    stride[1] = 1;
    stride[2] = 1;

    block[0] = 1;
    block[1] = 1;
    block[2] = 1;

    hid_t filespace_id = H5Screate_simple(num_dims, dim_gsize, NULL);
    hid_t memspace_id  = H5Screate_simple(num_dims, dim_lsize, NULL);

    hid_t dataset_id1 = H5Dcreate2(file_id, "T1", H5T_NATIVE_DOUBLE, filespace_id, 
                          H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    status = H5Sselect_hyperslab(filespace_id,H5S_SELECT_SET,dim_offset,stride,dim_lsize,block);


    cout<<" status select hyperslab "<<status<<endl;

    double *slab_vec;
    slab_vec = new double[topo.lenslab];


    int l=0;
    for(int k = topo.iter_slab[_INIZ_]; k <= topo.iter_slab[_ENDZ_]; k++)
       for(int j =  topo.iter_slab[_INIY_]; j <=  topo.iter_slab[_ENDY_]; j++)
           for(int i =  topo.iter_slab[_INIX_]; i <=  topo.iter_slab[_ENDX_]; i++){
              slab_vec[l] =  T.vector[I1D(i,j,k)];
             cout<<" "<<slab_vec[l]<<endl;
               l=l+1;
           }
 


     status = H5Dwrite(dataset_id1, H5T_NATIVE_DOUBLE, memspace_id, filespace_id, fa_plist_id2,
                     slab_vec);

/*    file_id = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
*/
/*
    ret = H5Pclose(fa_plist_id);
    cout<<" ret2 " <<ret<<" "<<-1<<endl;


    dx_plist_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(dx_plist_id, H5FD_MPIO_COLLECTIVE);


    cout<<"after set dxpl"<<endl;
	// Create the data space for the dataset.
	hsize_t dim_size[1];               // dataset dimensions
	dim_size[0] = _lNx_*_lNy_*_lNz_;
    int num_dims=1;
*/

    // Create the data space for the dataset.
/*	hsize_t dim_size[3];               // dataset dimensions
	dim_size[0] = _lNz_;
    dim_size[1] = _lNy_;
    dim_size[2] = _lNx_;
    int num_dims=3;

    cout<<"dim 1: "<<dim_size[0]<<" 2: "<<dim_size[1]<<" 3: "<<dim_size[2]<<endl;

    hsize_t array_tonto[1];
    array_tonto[0] = 1;

    hid_t dataspace_id = H5Screate_simple(num_dims, dim_size, NULL);
*/
//	DataSpace dataspace(num_dims, dim_size);

 //   hid_t attspace_id = H5Screate_simple(num_dims, array_tonto, NULL);


	// Create the dataset.      
//	DataSet dataset = file.createDataSet("T1", H5T_NATIVE_DOUBLE, dataspace);
/*
    hid_t dataset_id1 = H5Dcreate2(file_id, "T1", H5T_NATIVE_DOUBLE, dataspace_id, 
                          H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

//    hid_t dataset_id1 = H5Dcreate(file_id, "T1", H5T_NATIVE_DOUBLE, dataspace_id, 
//                           H5P_DEFAULT);

 
    hid_t dataset_id2 = H5Dcreate2(file_id, "T2", H5T_NATIVE_DOUBLE, dataspace_id, 
                          H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);


//    hid_t attri1 = H5Acreate(file_id,"Time",H5T_NATIVE_DOUBLE, attspace_id,H5P_DEFAULT,H5P_DEFAULT);

    cout<<"Antes del write "<<endl;
*/
/*    status = H5Dwrite(dataset_id1, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, dx_plist_id,
                     T.vector);
    
    cout<<" status " <<status<<" -1 "<<endl;
*/
/*
    status = H5Dwrite(dataset_id1, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                     T.vector);


    status = H5Dwrite(dataset_id2, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                     meshx.vector);
*/
/*    status = H5Dwrite(dataset_id2, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, dx_plist_id,
                     meshx.vector);
*/
                     cout<<" status " <<status<<" -1 "<<endl;



    double tiempo=3.5;
    
  //  status = H5Awrite(attri1, H5T_NATIVE_DOUBLE, &tiempo);

   /* End access to the dataset and release resources used by it. */

  //  status = H5Aclose(attri1);

   status = H5Dclose(dataset_id1);

//   status = H5Dclose(dataset_id2);


   /* Terminate access to the data space. */ 
   status = H5Sclose(filespace_id);
   status = H5Sclose(memspace_id);


    /* Terminate access to the file. */
    status = H5Fclose(file_id); 


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

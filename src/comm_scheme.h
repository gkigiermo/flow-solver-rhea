#ifndef _comm_scheme_
#define _comm_scheme_
#include<mpi.h>
#include<iostream>
#include<fstream>


#include "domain.h"
#include "parameters.h"

using namespace std;

class comm_scheme{
    
    public:
        comm_scheme(domain*, int , int ,int );
        void printCommScheme(int);
        void printCommSchemeToFile(int);
        void update(double*);
        int getNB(int id){return neighb[id];};
        int getSize(){return len;};
        int getlNx(){return lNx;};
        int getlNy(){return lNy;};
        int getlNz(){return lNz;};
        int getRank(){return rank;};

        int getOffx(){return offx;};
        int getOffy(){return offy;};
        int getOffz(){return offz;};





        // Basic iterators to move on the domain

        int iter_common[2][6]; 
        int iter_bound[26][6];
        int iter_halo[26][6];
        int iter_toSend[26][6];
        int iter_toRecv[26][6];

        int iter_glob_ind[6];



    protected:
        int rank;
        int np;
        MPI_Comm RHEA_3DCOMM;
        int npx;
        int npy;
        int npz;

        // local number of inner cells
        int lNx;
        int lNy;
        int lNz;

        // local cells plus boundaries and ghosts /*not sure if both are needed*/
        int lncellsx;
        int lncellsy;
        int lncellsz;

        // offset in the global ids /* currently under the assumption equal partition */
        int offx; 
        int offy;
        int offz;

        // total size of a 1D array
        int len;

        // IDs of the neighbours
        int neighb[26];

        // IDs of the boundaries
        int lbounds[26];


        // 2nd level info

        int *info_2nd;
        int off_nb[6];


        // 3er level info

        int proc_z_start;
        int proc_z_end;

        // For the communication
        int len_xy;
        int len_xz;
        int len_yz;

        double *pack_send_w;
        double *pack_send_e;
        double *pack_send_s;
        double *pack_send_n;
        double *pack_send_b;
        double *pack_send_f;

        double *pack_recv_w;
        double *pack_recv_e;
        double *pack_recv_s;
        double *pack_recv_n;
        double *pack_recv_b;
        double *pack_recv_f;


        void create_common_iters();
        void create_basic_bound_iters();
        void create_halo_iters();
        void create_toSend_iters();
        void create_toRecv_iters();
        void create_global_iters();
        void create_comm_arrays();

        void exchange_2nd_level_neighbour_info();
        void find_extra_neighbours();
        void pack(double*);
        void halo_exchange();
        void unpack(double*);

};


#endif

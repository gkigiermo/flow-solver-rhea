#ifndef _ParallelTopology_
#define _ParallelTopology_

#include <mpi.h>
#include <iostream>
#include <fstream>


#include "ComputationalDomain.hpp"
#include "MacroParameters.hpp"

class ParallelTopology{
    
    public:
        ParallelTopology(ComputationalDomain*, int , int ,int );
        void printCommScheme(int);
        void printCommSchemeToFile(int);
        void update(double*);
        void update_simple(double*);
        int getNB(int id){return neighb[id];};
        int getSize(){return len;};
        int getlNx(){return lNx;};
        int getlNy(){return lNy;};
        int getlNz(){return lNz;};
        int getRank(){return rank;};

        int getOffx(){return offx;};
        int getOffy(){return offy;};
        int getOffz(){return offz;};

        
        ComputationalDomain* getMesh(){return mymesh; };


        // Basic iterators to move on the ComputationalDomain

        int iter_common[2][6]; 
        int iter_bound[26][6];
        int iter_halo[26][6];
        int iter_toSend[26][6];
        int iter_toRecv[26][6];

        int iter_glob_ind[6];

        int iter_slab[6];

        int offslab_x;
        int offslab_y;
        int offslab_z;

        int lenslabx;
        int lenslaby;
        int lenslabz;

        int lenslab;

    protected:
        ComputationalDomain* mymesh;
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

        // Tag Id
        int tagid_s[26];
        int tagid_r[26];


        // 2nd level info
        int is_inix;
        int is_endx;
        int is_iniy;
        int is_endy;
        int is_iniz;
        int is_endz;

        int *info_2nd;
        int off_nb[6];


        // 3er level info

        int proc_z_start;
        int proc_z_end;

        // For the communication
        int len_xy;
        int len_xz;
        int len_yz;

        int len_1Dx;
        int len_1Dy;
        int len_1Dz;

        int len_1pt;

        double *pack_send_w;
        double *pack_send_e;
        double *pack_send_s;
        double *pack_send_n;
        double *pack_send_b;
        double *pack_send_f;

        double *pack_send_ws;
        double *pack_send_wn;
        double *pack_send_wb;
        double *pack_send_wf;

        double *pack_send_es;
        double *pack_send_en;
        double *pack_send_eb;
        double *pack_send_ef;

        double *pack_send_sb;
        double *pack_send_sf;
        double *pack_send_nb;
        double *pack_send_nf;

        double *pack_send_wsb;
        double *pack_send_wnb;
        double *pack_send_wsf;
        double *pack_send_wnf;

        double *pack_send_esb;
        double *pack_send_enb;
        double *pack_send_esf;
        double *pack_send_enf;

        double *pack_recv_w;
        double *pack_recv_e;
        double *pack_recv_s;
        double *pack_recv_n;
        double *pack_recv_b;
        double *pack_recv_f;

        double *pack_recv_ws;
        double *pack_recv_wn;
        double *pack_recv_wb;
        double *pack_recv_wf;

        double *pack_recv_es;
        double *pack_recv_en;
        double *pack_recv_eb;
        double *pack_recv_ef;

        double *pack_recv_sb;
        double *pack_recv_sf;
        double *pack_recv_nb;
        double *pack_recv_nf;

        double *pack_recv_wsb;
        double *pack_recv_wnb;
        double *pack_recv_wsf;
        double *pack_recv_wnf;

        double *pack_recv_esb;
        double *pack_recv_enb;
        double *pack_recv_esf;
        double *pack_recv_enf;



        void create_common_iters();
        void create_basic_bound_iters();
        void create_halo_iters();
        void create_toSend_iters();
        void create_toRecv_iters();
        void create_global_iters();
        void create_comm_arrays();

        void exchange_2nd_level_neighbour_info();
        void find_extra_neighbours();
        void create_complex_bound_iters();
        void create_complex_halo_iters();
        void create_complex_toRecv_iters();
        void create_complex_toSend_iters();
        void create_complex_comm_arrays();

        void pack(double*);
        void pack_simple(double*);
        void halo_exchange();
        void halo_exchange_simple();
        void unpack(double*);
        void unpack_simple(double*);

        void calculate_tags();


};


#endif

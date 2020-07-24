//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                              //
// RHEA - an open-source Reproducible Hybrid-architecture flow solver Engineered for Academia                   //
//                                                                                                              //
// Rhea was the Titaness great Mother of the Gods, and goddess of female fertility, motherhood, and generation. //
// Her name means "flow" and "ease", representing the eternal flow of time and generations with ease.           //
//                                                                                                              //
//                                                                                                              //
// REHA is released under the MIT License:                                                                      //
//                                                                                                              //
// Copyright (c) 2020 Lluis Jofre Cruanyes & Guillermo Oyarzun Altamirano.                                      //
//                                                                                                              //
// Permission is hereby granted, free of charge, to any person obtaining a copy                                 //
// of this software and associated documentation files (the "Software"), to deal                                //
// in the Software without restriction, including without limitation the rights                                 //
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell                                    //
// copies of the Software, and to permit persons to whom the Software is                                        //
// furnished to do so, subject to the following conditions:                                                     //
//                                                                                                              //
// The above copyright notice and this permission notice shall be included in all                               //
// copies or substantial portions of the Software.                                                              //
//                                                                                                              //
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR                                   //
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,                                     //
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE                                  //
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER                                       //
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,                                //
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE                                //
// SOFTWARE.                                                                                                    //
//                                                                                                              //
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef _rhea_flow_solver_
#define _rhea_flow_solver_

#include <stdio.h>
#include <math.h>
#include <iostream>
#include <string>
#include <mpi.h>
#include "src/parameters.h"
#include "src/domain.h"
#include "src/comm_scheme.h"
#include "src/parvec.h"

using namespace std;

class rhea_flow_solver{
   
    ////////// VARIABLES DESCRIPTION //////////

    // Primitive variables:
    //   - Density rho
    //   - Velocities u, v, w
    //   - Specific total energy E
    //     ... E = e + ke
    //     ... is the sum of internal energy e
    //     ... and kinetic energy ke = (u*u + v*v + w*w)/2

    // Conserved variables:
    //   - Mass rho
    //   - Momentum rho*u, rho*v, rho*w
    //   - Total energy rho*E

    // Thermodynamic state:
    //   - Pressure P
    //   - Temperature T
    //   - Speed of sound sos

    // Thermophysical properties:
    //   - Specific gas constant R_specific
    //   - Ratio of heat capacities gamma
    //   - Dynamic viscosity mu
    //   - Thermal conductivity kappa
 
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
        int iter_bound[6][6];
        int iter_halo[6][6];
        int iter_toSend[6][6];
        int iter_toRecv[6][6];

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
        int neighb[6];

        // IDs of the boundaries
        int lbounds[6];

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

        void pack(double*);
        void halo_exchange();
        void unpack(double*);

};

#endif

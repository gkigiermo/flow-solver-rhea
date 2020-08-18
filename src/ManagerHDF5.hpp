#ifndef _printer_
#define _printer_

#include <iostream>
#include <hdf5.h>
#include <list>

#include "parameters.h"
#include "comm_scheme.h"
#include "parvec.h"
using namespace std;


class printer{

    public:
        printer(){};
        printer(comm_scheme*,char const*);
        void addField(parvec*);
        void printOnScreen();
        void write(int,double);
        void read(int);
    private:
        
        comm_scheme* myTopo;
        char outname[50];

        std::list<parvec*> fieldList;

        int num_dims;
        hsize_t dim_gsize[3];               // dataset dimensions
        hsize_t dim_lsize[3];               // dataset dimensions
        hsize_t dim_offset[3];
        hsize_t stride[3];
        hsize_t block[3];

        double *slab_vec;

        int _lNx_;
        int _lNy_;
        int _lNz_;

        // Variables to store the attribute

        hsize_t array_1D[1];
        int num_dim1D;



};


#endif

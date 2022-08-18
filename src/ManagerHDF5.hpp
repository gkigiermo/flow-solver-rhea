#ifndef _ManagerHDF5_
#define _ManagerHDF5_

#include <iostream>
#include <hdf5.h>
#include <list>
#include <string.h>
#include <map>
#include "MacroParameters.hpp"
#include "ParallelTopology.hpp"
#include "DistributedArray.hpp"

class ManagerHDF5{

    public:
        ManagerHDF5(){};
        ManagerHDF5(ParallelTopology*,char const*,bool);
        void addField(DistributedArray*);
        void printOnScreen();
        void write(int);
        void read(char const*);

        void addAttributeDouble(std::string str){ double dval=0.0; dattrib[str]=dval;};
        void addAttributeInt(std::string str){ int ival =0; iattrib[str] = ival; };

        void setAttribute(std::string str,double dval){ dattrib[str]=dval;};
        void setAttribute(std::string str,int ival){ iattrib[str] = ival; };

        double getAttributeDouble(std::string str){ return dattrib[str];};
        int    getAttributeInt(std::string str){ return iattrib[str];};


    private:
        
        ParallelTopology* myTopo;
        char outname[50];

        std::list<DistributedArray*> fieldList;

        bool gen_xdmf;

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


        //Attribute list

	std::map<std::string,double> dattrib;
	std::map<std::string,int> iattrib;


};


#endif

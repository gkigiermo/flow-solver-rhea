#ifndef _DistributedArray_
#define _DistributedArray_

#include <iostream>
#include "MacroParameters.hpp"
#include "ParallelTopology.hpp"

using namespace std;

class DistributedArray{
    public:
        DistributedArray(){};
        DistributedArray(ParallelTopology*,char const*);
        void update();
        void update_simple();
        void fillEdgeCornerBoundaries();

	#pragma acc routine
        double& operator[](int);

	#pragma acc routine
        void operator = (double);

        void setTopology(ParallelTopology* topo, char const* );
        char* printName(){return fieldName;};

        double* vector;
        int size;
        int ini_x;
        int ini_y;
        int ini_z;
        int fin_x;
        int fin_y;
        int fin_z;
    protected:

        ParallelTopology* mydomain;
        char  fieldName[30];
};

#endif


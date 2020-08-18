#ifndef _parvec_
#define _parvec_

#include <iostream>
#include "parameters.h"
#include "comm_scheme.h"

using namespace std;

class parvec{
    public:
        parvec(){};
        parvec(comm_scheme*,char const*);
        void update();
        void update_simple();
        void fillEdgeCornerBoundaries();
        double& operator[](int);
        void operator= (double);

        void setTopology(comm_scheme* topo, char const* );
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

        comm_scheme* mydomain;
        char  fieldName[30];
};

#endif


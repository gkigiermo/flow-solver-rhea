#ifndef _parvec_
#define _parvec_

#include <iostream>
#include "parameters.h"
#include "comm_scheme.h"

using namespace std;

class parvec{
    public:
        parvec(){};
        parvec(comm_scheme*);
        void update();
        double& operator[](int);
        void operator= (double);


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
};

#endif


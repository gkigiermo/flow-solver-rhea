#include "parvec.h"

parvec :: parvec(comm_scheme* topo)
{
    size = topo->getSize();
    vector = new double[size];
    ini_x = topo->iter_common[_INNER_][_INIX_] ;
    ini_y = topo->iter_common[_INNER_][_INIY_] ;
    ini_z = topo->iter_common[_INNER_][_INIZ_] ;
    fin_x = topo->iter_common[_INNER_][_ENDX_] ;
    fin_y = topo->iter_common[_INNER_][_ENDY_] ;
    fin_z = topo->iter_common[_INNER_][_ENDZ_] ;

    mydomain = topo;   
}

void parvec::update()
{
    mydomain->update(vector);
}

double& parvec::operator[](int idx)
{
    return vector[idx];
}

void parvec::operator=(double val)
{
    for(int l = 0; l < size; l++)
        vector[l] = val;
}

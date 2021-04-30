#include "DistributedArray.hpp"

DistributedArray::DistributedArray(ParallelTopology* topo, char const* auxName)
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

    sprintf(fieldName,"%s",auxName);

}

void DistributedArray::setTopology(ParallelTopology* topo, char const * auxName) {

    size = topo->getSize();
    vector = new double[size];
    ini_x = topo->iter_common[_INNER_][_INIX_] ;
    ini_y = topo->iter_common[_INNER_][_INIY_] ;
    ini_z = topo->iter_common[_INNER_][_INIZ_] ;
    fin_x = topo->iter_common[_INNER_][_ENDX_] ;
    fin_y = topo->iter_common[_INNER_][_ENDY_] ;
    fin_z = topo->iter_common[_INNER_][_ENDZ_] ;

    mydomain = topo;   

    sprintf(fieldName,"%s",auxName);


}

void DistributedArray::update()
{
    mydomain->update(vector);

    fillEdgeCornerBoundaries();
}

void DistributedArray::update_simple()
{
    mydomain->update_simple(vector);
    fillEdgeCornerBoundaries();
}


double& DistributedArray::operator[](int idx)
{
    return vector[idx];
}

void DistributedArray::operator=(double val)
{
    for(int l = 0; l < size; l++)
        vector[l] = val;
}






void DistributedArray::fillEdgeCornerBoundaries() {

    int _lNx_ = mydomain->getlNx();
    int _lNy_ = mydomain->getlNy();
    //int _lNz_ = mydomain->getlNz();


    /// West-South boundary points: rho, rhou, rhov, rhow and rhoE
    for(int i = mydomain->iter_bound[_WEST_S_][_INIX_]; i <= mydomain->iter_bound[_WEST_S_][_ENDX_]; i++) {
        for(int j = mydomain->iter_bound[_WEST_S_][_INIY_]; j <= mydomain->iter_bound[_WEST_S_][_ENDY_]; j++) {
            for(int k = mydomain->iter_bound[_WEST_S_][_INIZ_]; k <= mydomain->iter_bound[_WEST_S_][_ENDZ_]; k++) {
                vector[I1D(i,j,k)]  = ( 1.0/2.0 )*( vector[I1D(i+1,j,k)] + vector[I1D(i,j+1,k)] );
            }
        }
    }

    /// West-North boundary points: rho, rhou, rhov, rhow and rhoE
    for(int i = mydomain->iter_bound[_WEST_N_][_INIX_]; i <= mydomain->iter_bound[_WEST_N_][_ENDX_]; i++) {
        for(int j = mydomain->iter_bound[_WEST_N_][_INIY_]; j <= mydomain->iter_bound[_WEST_N_][_ENDY_]; j++) {
            for(int k = mydomain->iter_bound[_WEST_N_][_INIZ_]; k <= mydomain->iter_bound[_WEST_N_][_ENDZ_]; k++) {
                vector[I1D(i,j,k)]  = ( 1.0/2.0 )*( vector[I1D(i+1,j,k)] + vector[I1D(i,j-1,k)] );
            }
        }
    }

    /// West-Back boundary points: rho, rhou, rhov, rhow and rhoE
    for(int i = mydomain->iter_bound[_WEST_B_][_INIX_]; i <= mydomain->iter_bound[_WEST_B_][_ENDX_]; i++) {
        for(int j = mydomain->iter_bound[_WEST_B_][_INIY_]; j <= mydomain->iter_bound[_WEST_B_][_ENDY_]; j++) {
            for(int k = mydomain->iter_bound[_WEST_B_][_INIZ_]; k <= mydomain->iter_bound[_WEST_B_][_ENDZ_]; k++) {
                vector[I1D(i,j,k)]  = ( 1.0/2.0 )*( vector[I1D(i+1,j,k)] + vector[I1D(i,j,k+1)] );
            }
        }
    }

    /// West-Front boundary points: rho, rhou, rhov, rhow and rhoE
    for(int i = mydomain->iter_bound[_WEST_F_][_INIX_]; i <= mydomain->iter_bound[_WEST_F_][_ENDX_]; i++) {
        for(int j = mydomain->iter_bound[_WEST_F_][_INIY_]; j <= mydomain->iter_bound[_WEST_F_][_ENDY_]; j++) {
            for(int k = mydomain->iter_bound[_WEST_F_][_INIZ_]; k <= mydomain->iter_bound[_WEST_F_][_ENDZ_]; k++) {
                vector[I1D(i,j,k)]  = ( 1.0/2.0 )*( vector[I1D(i+1,j,k)] + vector[I1D(i,j,k-1)] );
            }
        }
    }

    /// East-South boundary points: rho, rhou, rhov, rhow and rhoE
    for(int i = mydomain->iter_bound[_EAST_S_][_INIX_]; i <= mydomain->iter_bound[_EAST_S_][_ENDX_]; i++) {
        for(int j = mydomain->iter_bound[_EAST_S_][_INIY_]; j <= mydomain->iter_bound[_EAST_S_][_ENDY_]; j++) {
            for(int k = mydomain->iter_bound[_EAST_S_][_INIZ_]; k <= mydomain->iter_bound[_EAST_S_][_ENDZ_]; k++) {
                vector[I1D(i,j,k)]  = ( 1.0/2.0 )*( vector[I1D(i-1,j,k)] + vector[I1D(i,j+1,k)] );
            }
        }
    }

    /// East-North boundary points: rho, rhou, rhov, rhow and rhoE
    for(int i = mydomain->iter_bound[_EAST_N_][_INIX_]; i <= mydomain->iter_bound[_EAST_N_][_ENDX_]; i++) {
        for(int j = mydomain->iter_bound[_EAST_N_][_INIY_]; j <= mydomain->iter_bound[_EAST_N_][_ENDY_]; j++) {
            for(int k = mydomain->iter_bound[_EAST_N_][_INIZ_]; k <= mydomain->iter_bound[_EAST_N_][_ENDZ_]; k++) {
                vector[I1D(i,j,k)]  = ( 1.0/2.0 )*( vector[I1D(i-1,j,k)] + vector[I1D(i,j-1,k)] );
            }
        }
    }

    /// East-Back boundary points: rho, rhou, rhov, rhow and rhoE
    for(int i = mydomain->iter_bound[_EAST_B_][_INIX_]; i <= mydomain->iter_bound[_EAST_B_][_ENDX_]; i++) {
        for(int j = mydomain->iter_bound[_EAST_B_][_INIY_]; j <= mydomain->iter_bound[_EAST_B_][_ENDY_]; j++) {
            for(int k = mydomain->iter_bound[_EAST_B_][_INIZ_]; k <= mydomain->iter_bound[_EAST_B_][_ENDZ_]; k++) {
                vector[I1D(i,j,k)]  = ( 1.0/2.0 )*( vector[I1D(i-1,j,k)] + vector[I1D(i,j,k+1)] );
                //cout<<" rank "<<mydomain->getRank()<<" Sets the variables as a corner "<<endl;
            }
        }
    }

    /// East-Front boundary points: rho, rhou, rhov, rhow and rhoE
    for(int i = mydomain->iter_bound[_EAST_F_][_INIX_]; i <= mydomain->iter_bound[_EAST_F_][_ENDX_]; i++) {
        for(int j = mydomain->iter_bound[_EAST_F_][_INIY_]; j <= mydomain->iter_bound[_EAST_F_][_ENDY_]; j++) {
            for(int k = mydomain->iter_bound[_EAST_F_][_INIZ_]; k <= mydomain->iter_bound[_EAST_F_][_ENDZ_]; k++) {
                vector[I1D(i,j,k)]  = ( 1.0/2.0 )*( vector[I1D(i-1,j,k)] + vector[I1D(i,j,k-1)] );
            }
        }
    }

    /// South-Back boundary points: rho, rhou, rhov, rhow and rhoE
    for(int i = mydomain->iter_bound[_SOUTH_B_][_INIX_]; i <= mydomain->iter_bound[_SOUTH_B_][_ENDX_]; i++) {
        for(int j = mydomain->iter_bound[_SOUTH_B_][_INIY_]; j <= mydomain->iter_bound[_SOUTH_B_][_ENDY_]; j++) {
            for(int k = mydomain->iter_bound[_SOUTH_B_][_INIZ_]; k <= mydomain->iter_bound[_SOUTH_B_][_ENDZ_]; k++) {
                vector[I1D(i,j,k)]  = ( 1.0/2.0 )*( vector[I1D(i,j+1,k)] + vector[I1D(i,j,k+1)] );
            }
        }
    }

    /// South-Front boundary points: rho, rhou, rhov, rhow and rhoE
    for(int i = mydomain->iter_bound[_SOUTH_F_][_INIX_]; i <= mydomain->iter_bound[_SOUTH_F_][_ENDX_]; i++) {
        for(int j = mydomain->iter_bound[_SOUTH_F_][_INIY_]; j <= mydomain->iter_bound[_SOUTH_F_][_ENDY_]; j++) {
            for(int k = mydomain->iter_bound[_SOUTH_F_][_INIZ_]; k <= mydomain->iter_bound[_SOUTH_F_][_ENDZ_]; k++) {
                vector[I1D(i,j,k)]  = ( 1.0/2.0 )*( vector[I1D(i,j+1,k)] + vector[I1D(i,j,k-1)] );
            }
        }
    }

    /// North-Back boundary points: rho, rhou, rhov, rhow and rhoE
    for(int i = mydomain->iter_bound[_NORTH_B_][_INIX_]; i <= mydomain->iter_bound[_NORTH_B_][_ENDX_]; i++) {
        for(int j = mydomain->iter_bound[_NORTH_B_][_INIY_]; j <= mydomain->iter_bound[_NORTH_B_][_ENDY_]; j++) {
            for(int k = mydomain->iter_bound[_NORTH_B_][_INIZ_]; k <= mydomain->iter_bound[_NORTH_B_][_ENDZ_]; k++) {
                vector[I1D(i,j,k)]  = ( 1.0/2.0 )*( vector[I1D(i,j-1,k)] + vector[I1D(i,j,k+1)] );
            }
        }
    }

    /// North-Front boundary points: rho, rhou, rhov, rhow and rhoE
    for(int i = mydomain->iter_bound[_NORTH_F_][_INIX_]; i <= mydomain->iter_bound[_NORTH_F_][_ENDX_]; i++) {
        for(int j = mydomain->iter_bound[_NORTH_F_][_INIY_]; j <= mydomain->iter_bound[_NORTH_F_][_ENDY_]; j++) {
            for(int k = mydomain->iter_bound[_NORTH_F_][_INIZ_]; k <= mydomain->iter_bound[_NORTH_F_][_ENDZ_]; k++) {
                vector[I1D(i,j,k)]  = ( 1.0/2.0 )*( vector[I1D(i,j-1,k)] + vector[I1D(i,j,k-1)] );
            }
        }
    }

    /// West-South-Back boundary point: rho, rhou, rhov, rhow and rhoE
    for(int i = mydomain->iter_bound[_WEST_S_B_][_INIX_]; i <= mydomain->iter_bound[_WEST_S_B_][_ENDX_]; i++) {
        for(int j = mydomain->iter_bound[_WEST_S_B_][_INIY_]; j <= mydomain->iter_bound[_WEST_S_B_][_ENDY_]; j++) {
            for(int k = mydomain->iter_bound[_WEST_S_B_][_INIZ_]; k <= mydomain->iter_bound[_WEST_S_B_][_ENDZ_]; k++) {
                vector[I1D(i,j,k)]  = ( 1.0/3.0 )*( vector[I1D(i+1,j,k)] + vector[I1D(i,j+1,k)] + vector[I1D(i,j,k+1)] );
            }
        }
    }

    /// West-North-Back boundary point: rho, rhou, rhov, rhow and rhoE
    for(int i = mydomain->iter_bound[_WEST_N_B_][_INIX_]; i <= mydomain->iter_bound[_WEST_N_B_][_ENDX_]; i++) {
        for(int j = mydomain->iter_bound[_WEST_N_B_][_INIY_]; j <= mydomain->iter_bound[_WEST_N_B_][_ENDY_]; j++) {
            for(int k = mydomain->iter_bound[_WEST_N_B_][_INIZ_]; k <= mydomain->iter_bound[_WEST_N_B_][_ENDZ_]; k++) {
                vector[I1D(i,j,k)]  = ( 1.0/3.0 )*( vector[I1D(i+1,j,k)] + vector[I1D(i,j-1,k)] + vector[I1D(i,j,k+1)] );
            }
        }
    }

    /// West-South-Front boundary point: rho, rhou, rhov, rhow and rhoE
    for(int i = mydomain->iter_bound[_WEST_S_F_][_INIX_]; i <= mydomain->iter_bound[_WEST_S_F_][_ENDX_]; i++) {
        for(int j = mydomain->iter_bound[_WEST_S_F_][_INIY_]; j <= mydomain->iter_bound[_WEST_S_F_][_ENDY_]; j++) {
            for(int k = mydomain->iter_bound[_WEST_S_F_][_INIZ_]; k <= mydomain->iter_bound[_WEST_S_F_][_ENDZ_]; k++) {
                vector[I1D(i,j,k)]  = ( 1.0/3.0 )*( vector[I1D(i+1,j,k)] + vector[I1D(i,j+1,k)] + vector[I1D(i,j,k-1)] );
            }
        }
    }

    /// West-North-Front boundary point: rho, rhou, rhov, rhow and rhoE
    for(int i = mydomain->iter_bound[_WEST_N_F_][_INIX_]; i <= mydomain->iter_bound[_WEST_N_F_][_ENDX_]; i++) {
        for(int j = mydomain->iter_bound[_WEST_N_F_][_INIY_]; j <= mydomain->iter_bound[_WEST_N_F_][_ENDY_]; j++) {
            for(int k = mydomain->iter_bound[_WEST_N_F_][_INIZ_]; k <= mydomain->iter_bound[_WEST_N_F_][_ENDZ_]; k++) {
                vector[I1D(i,j,k)]  = ( 1.0/3.0 )*( vector[I1D(i+1,j,k)] + vector[I1D(i,j-1,k)] + vector[I1D(i,j,k-1)] );
            }
        }
    }

    /// East-South-Back boundary point: rho, rhou, rhov, rhow and rhoE
    for(int i = mydomain->iter_bound[_EAST_S_B_][_INIX_]; i <= mydomain->iter_bound[_EAST_S_B_][_ENDX_]; i++) {
        for(int j = mydomain->iter_bound[_EAST_S_B_][_INIY_]; j <= mydomain->iter_bound[_EAST_S_B_][_ENDY_]; j++) {
            for(int k = mydomain->iter_bound[_EAST_S_B_][_INIZ_]; k <= mydomain->iter_bound[_EAST_S_B_][_ENDZ_]; k++) {
                vector[I1D(i,j,k)]  = ( 1.0/3.0 )*( vector[I1D(i-1,j,k)] + vector[I1D(i,j+1,k)] + vector[I1D(i,j,k+1)] );
            }
        }
    }

    /// East-North-Back boundary point: rho, rhou, rhov, rhow and rhoE
    for(int i = mydomain->iter_bound[_EAST_N_B_][_INIX_]; i <= mydomain->iter_bound[_EAST_N_B_][_ENDX_]; i++) {
        for(int j = mydomain->iter_bound[_EAST_N_B_][_INIY_]; j <= mydomain->iter_bound[_EAST_N_B_][_ENDY_]; j++) {
            for(int k = mydomain->iter_bound[_EAST_N_B_][_INIZ_]; k <= mydomain->iter_bound[_EAST_N_B_][_ENDZ_]; k++) {
                vector[I1D(i,j,k)]  = ( 1.0/3.0 )*( vector[I1D(i-1,j,k)] + vector[I1D(i,j-1,k)] + vector[I1D(i,j,k+1)] );
            }
        }
    }

    /// East-South-Front boundary point: rho, rhou, rhov, rhow and rhoE
    for(int i = mydomain->iter_bound[_EAST_S_F_][_INIX_]; i <= mydomain->iter_bound[_EAST_S_F_][_ENDX_]; i++) {
        for(int j = mydomain->iter_bound[_EAST_S_F_][_INIY_]; j <= mydomain->iter_bound[_EAST_S_F_][_ENDY_]; j++) {
            for(int k = mydomain->iter_bound[_EAST_S_F_][_INIZ_]; k <= mydomain->iter_bound[_EAST_S_F_][_ENDZ_]; k++) {
                vector[I1D(i,j,k)]  = ( 1.0/3.0 )*( vector[I1D(i-1,j,k)] + vector[I1D(i,j+1,k)] + vector[I1D(i,j,k-1)] );
            }
        }
    }

    /// East-North-Front boundary point: rho, rhou, rhov, rhow and rhoE
    for(int i = mydomain->iter_bound[_EAST_N_F_][_INIX_]; i <= mydomain->iter_bound[_EAST_N_F_][_ENDX_]; i++) {
        for(int j = mydomain->iter_bound[_EAST_N_F_][_INIY_]; j <= mydomain->iter_bound[_EAST_N_F_][_ENDY_]; j++) {
            for(int k = mydomain->iter_bound[_EAST_N_F_][_INIZ_]; k <= mydomain->iter_bound[_EAST_N_F_][_ENDZ_]; k++) {
                vector[I1D(i,j,k)]  = ( 1.0/3.0 )*( vector[I1D(i-1,j,k)] + vector[I1D(i,j-1,k)] + vector[I1D(i,j,k-1)] );
            }
        }
    }

};

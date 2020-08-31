#include "myRHEA.hpp"

using namespace std;


////////// myRHEA CLASS //////////

void myRHEA::setInitialConditions() {

    /// IMPORTANT: This method needs to be modified/overwritten according to the problem under consideration

    /// All (inner, boundary & halo) points: u, v, w, P and T
    for(int i = topo->iter_common[_ALL_][_INIX_]; i <= topo->iter_common[_ALL_][_ENDX_]; i++) {
        for(int j = topo->iter_common[_ALL_][_INIY_]; j <= topo->iter_common[_ALL_][_ENDY_]; j++) {
            for(int k = topo->iter_common[_ALL_][_INIZ_]; k <= topo->iter_common[_ALL_][_ENDZ_]; k++) {
		if( mesh->x[i] > 0.8 ) {
                    if( mesh->y[j] < 0.8 ) {
                        u_field[I1D(i,j,k)] = 0.0;
                        v_field[I1D(i,j,k)] = 1.206;
                        w_field[I1D(i,j,k)] = 0.0;
                        P_field[I1D(i,j,k)] = 0.3;
                        T_field[I1D(i,j,k)] = P_field[I1D(i,j,k)]/( 0.5323*thermodynamics->getSpecificGasConstant() );
		    } else {
                        u_field[I1D(i,j,k)] = 0.0;
                        v_field[I1D(i,j,k)] = 0.0;
                        w_field[I1D(i,j,k)] = 0.0;
                        P_field[I1D(i,j,k)] = 1.5;
                        T_field[I1D(i,j,k)] = P_field[I1D(i,j,k)]/( 1.5*thermodynamics->getSpecificGasConstant() );
		    }
		} else {
                    if( mesh->y[j] < 0.8 ) {
                        u_field[I1D(i,j,k)] = 1.206;
                        v_field[I1D(i,j,k)] = 1.206;
                        w_field[I1D(i,j,k)] = 0.0;
                        P_field[I1D(i,j,k)] = 0.029;
                        T_field[I1D(i,j,k)] = P_field[I1D(i,j,k)]/( 0.138*thermodynamics->getSpecificGasConstant() );
		    } else {
                        u_field[I1D(i,j,k)] = 1.206;
                        v_field[I1D(i,j,k)] = 0.0;
                        w_field[I1D(i,j,k)] = 0.0;
                        P_field[I1D(i,j,k)] = 0.3;
                        T_field[I1D(i,j,k)] = P_field[I1D(i,j,k)]/( 0.5323*thermodynamics->getSpecificGasConstant() );
		    }
		}		
            }
        }
    }

};

void myRHEA::calculateSourceTerms() {

    /// IMPORTANT: This method needs to be modified/overwritten according to the problem under consideration

    /// Inner points: f_rhou, f_rhov, f_rhow and f_rhoE
    for(int i = topo->iter_common[_INNER_][_INIX_]; i <= topo->iter_common[_INNER_][_ENDX_]; i++) {
        for(int j = topo->iter_common[_INNER_][_INIY_]; j <= topo->iter_common[_INNER_][_ENDY_]; j++) {
            for(int k = topo->iter_common[_INNER_][_INIZ_]; k <= topo->iter_common[_INNER_][_ENDZ_]; k++) {
                f_rhou_field[I1D(i,j,k)] = 0.0;
                f_rhov_field[I1D(i,j,k)] = 0.0;
                f_rhow_field[I1D(i,j,k)] = 0.0;
                f_rhoE_field[I1D(i,j,k)] = 0.0;
            }
        }
    }

};


////////// MAIN //////////
int main(int argc, char** argv) {

    /// Initialize MPI
    MPI_Init(&argc, &argv);

    /// Process command line arguments
    string configuration_file;
    if( argc >= 2 ) {
        configuration_file = argv[1];
    } else {
        cout << "Proper usage: RHEA.exe configuration_file.yaml" << endl;
        MPI_Abort( MPI_COMM_WORLD, 1 );
    }

    /// Construct my RHEA
    myRHEA my_RHEA( configuration_file );

    /// Execute my RHEA
    my_RHEA.execute();

    /// Destruct my RHEA ... destructor is called automatically

    /// Finalize MPI
    MPI_Finalize();

    /// Return exit code of program
    return 0;

}

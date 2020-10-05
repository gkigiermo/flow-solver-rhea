#include "myRHEA.hpp"

using namespace std;

/// pi number
const double pi = 2.0*asin( 1.0 );

/// PROBLEM PARAMETERS ///
const double gamma_0 = 1.4;					/// Reference ratio of heat capacities
//const double Re_0    = Pi;					/// Reynolds number
const double Ma_0    = ( 1.0e-2 )/sqrt( gamma_0 );		/// Mach number
const double rho_0   = 1.0;					/// Reference density	
const double U_0     = 1.0;					/// Reference velocity
//const double mu_0    = rho_0*U_0*Pi/Re_0;			/// Dynamic viscosity	
const double P_0     = ( rho_0/gamma_0 )*pow( U_0/Ma_0, 2.0 );	/// Reference pressure
const double L       = 2.0*pi;					/// Domain size			

////////// myRHEA CLASS //////////

void myRHEA::setInitialConditions() {

    /// IMPORTANT: This method needs to be modified/overwritten according to the problem under consideration

    /// All (inner, boundary & halo) points: u, v, w, P and T
    for(int i = topo->iter_common[_ALL_][_INIX_]; i <= topo->iter_common[_ALL_][_ENDX_]; i++) {
        for(int j = topo->iter_common[_ALL_][_INIY_]; j <= topo->iter_common[_ALL_][_ENDY_]; j++) {
            for(int k = topo->iter_common[_ALL_][_INIZ_]; k <= topo->iter_common[_ALL_][_ENDZ_]; k++) {
                u_field[I1D(i,j,k)] = U_0*( cos( 2.0*pi*mesh->x[i]/L )*sin( 2.0*pi*mesh->y[j]/L ) );
                v_field[I1D(i,j,k)] = ( -1.0 )*U_0*( sin( 2.0*pi*mesh->x[i]/L )*cos( 2.0*pi*mesh->y[j]/L ) );
                w_field[I1D(i,j,k)] = 0.0;
                P_field[I1D(i,j,k)] = P_0 + ( rho_0*U_0*U_0/4.0 )*( cos( 4.0*pi*mesh->x[i]/L ) + cos( 4.0*pi*mesh->y[j]/L ) );
                T_field[I1D(i,j,k)] = P_field[I1D(i,j,k)]/( rho_0*thermodynamics->getSpecificGasConstant() );
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

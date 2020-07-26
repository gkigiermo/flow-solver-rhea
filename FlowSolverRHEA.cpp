#include "FlowSolverRHEA.hpp"

using namespace std;

////////// FIXED PARAMETERS //////////
const double rk_order = 3;			// Time-integration Runge-Kutta order
const double epsilon  = 1.0e-15;		// Small epsilon number
const double Pi       = 2.0*asin(1.0);		// Pi number

////////// FlowSolverRHEA CLASS //////////
FlowSolverRHEA::FlowSolverRHEA() {};

FlowSolverRHEA::FlowSolverRHEA(const string name_configuration_file) : configuration_file(name_configuration_file) {

    /// Read configuration file
    // ... work in progress (YAML language)	

    // The lines below need to be introduced via configuration (input) file !!
    const int _DEBUG_    = 0;
    const double Re_tau  = 100.0;				// Friction Reynolds number [-]
    const double delta   = 1.0;					// Channel half-height [m]
    const double u_tau   = 1.0;					// Friction velocity [m/s]
    const double rho_ref = 1.0;					// Reference density [kg/m3]
    const double P_ref   = 101325.0;				// Reference pressure [Pa]

    R_specific        = 287.058;
    gamma             = 1.4;
    mu                = rho_ref*u_tau*delta/Re_tau;
    kappa             = 0.0;
    x_0               = 0.0;
    y_0               = 0.0;
    z_0               = 0.0;
    L_x               = 4.0*Pi*delta;
    L_y               = 2.0*delta;
    L_z               = 4.0*Pi*delta/3.0;
    initial_time      = 0.0;
    final_time        = 1.0e3;
    num_grid_x        = 64;
    num_grid_y        = 64;
    num_grid_z        = 64;
    CFL               = 0.9;
    max_num_time_iter = 1e6;
    output_iter       = 1e2;
    bocos[_WEST_]     = _PERIODIC_;
    bocos[_EAST_]     = _PERIODIC_;
    bocos[_SOUTH_]    = _PERIODIC_;
    bocos[_NORTH_]    = _PERIODIC_;
    bocos[_BACK_]     = _PERIODIC_;
    bocos[_FRONT_]    = _PERIODIC_;
    np_x              = 2;
    np_y              = 1;
    np_z              = 1;
    // The lines above need to be introduced via configuration (input) file !!

    /// Initialize (construct) computational domain
    dom = new domain(L_x, L_y, L_z, x_0, y_0, z_0, RHEA_NX, RHEA_NY, RHEA_NZ);

    /// Add boundary conditions to computational domain
    dom->updateBocos(bocos);

    /// Initialize (construct) communication scheme
    topo = new comm_scheme(dom, np_x, np_y, np_z);
    //if(topo->getRank() == 0) dom->printDomain();
    //for(int p = 0; p < np_x*np_y*np_z; p++) topo->printCommSchemeToFile(p);

    /// Allocate (parallel) primitive, conserved and thermodynamic variables	
    rho_field  = new parvec(topo);
    u_field    = new parvec(topo);
    v_field    = new parvec(topo);
    w_field    = new parvec(topo);
    E_field    = new parvec(topo);
    rhou_field = new parvec(topo);
    rhov_field = new parvec(topo);
    rhow_field = new parvec(topo);
    rhoE_field = new parvec(topo);
    P_field    = new parvec(topo);
    T_field    = new parvec(topo);
    sos_field  = new parvec(topo);

    /// Allocate (parallel) time-integration variables	
    rho_0_field  = new parvec(topo);
    rhou_0_field = new parvec(topo);
    rhov_0_field = new parvec(topo);
    rhow_0_field = new parvec(topo);
    rhoE_0_field = new parvec(topo);    

    /// Allocate (parallel) time-integration fluxes	
    rho_rk1_flux  = new parvec(topo);
    rho_rk2_flux  = new parvec(topo);
    rho_rk3_flux  = new parvec(topo);
    rhou_rk1_flux = new parvec(topo);
    rhou_rk2_flux = new parvec(topo);
    rhou_rk3_flux = new parvec(topo);    
    rhov_rk1_flux = new parvec(topo);
    rhov_rk2_flux = new parvec(topo);
    rhov_rk3_flux = new parvec(topo);    
    rhow_rk1_flux = new parvec(topo);
    rhow_rk2_flux = new parvec(topo);
    rhow_rk3_flux = new parvec(topo);
    rhoE_rk1_flux = new parvec(topo);
    rhoE_rk2_flux = new parvec(topo);
    rhoE_rk3_flux = new parvec(topo);

    /// Allocate (parallel) inviscid fluxes	
    rho_inv_flux  = new parvec(topo);
    rhou_inv_flux = new parvec(topo);
    rhov_inv_flux = new parvec(topo);
    rhow_inv_flux = new parvec(topo);
    rhoE_inv_flux = new parvec(topo);

    /// Allocate (parallel) viscous fluxes	
    rhou_vis_flux = new parvec(topo);
    rhov_vis_flux = new parvec(topo);
    rhow_vis_flux = new parvec(topo);
    rhoE_vis_flux = new parvec(topo);

    /// Allocate (parallel) volumetric & external forces
    f_rhou_field = new parvec(topo);
    f_rhov_field = new parvec(topo);
    f_rhow_field = new parvec(topo);
    f_rhoE_field = new parvec(topo);

};

FlowSolverRHEA::FlowSolverRHEA(const FlowSolverRHEA &in) {};

FlowSolverRHEA::~FlowSolverRHEA() {

    /// Free domain, boundary conditions & topology
    if(dom != NULL) free(dom);	
    if(bocos != NULL) free(bocos);	
    if(topo != NULL) free(topo);

    /// Free primitive, conserved and thermodynamic variables
    if(rho_field != NULL)  free(rho_field);	
    if(u_field != NULL)    free(u_field);	
    if(v_field != NULL)    free(v_field);	
    if(w_field != NULL)    free(w_field);	
    if(E_field != NULL)    free(E_field);	
    if(rhou_field != NULL) free(rhou_field);	
    if(rhov_field != NULL) free(rhov_field);	
    if(rhow_field != NULL) free(rhow_field);	
    if(rhoE_field != NULL) free(rhoE_field);	
    if(P_field != NULL)    free(P_field);	
    if(T_field != NULL)    free(T_field);	
    if(sos_field != NULL)  free(sos_field);	
        
    /// Free time-integration variables
    if(rho_0_field != NULL)  free(rho_0_field);	
    if(rhou_0_field != NULL) free(rhou_0_field);	
    if(rhov_0_field != NULL) free(rhov_0_field);	
    if(rhow_0_field != NULL) free(rhow_0_field);	
    if(rhoE_0_field != NULL) free(rhoE_0_field);	
        
    /// Free time-integration fluxes
    if(rho_rk1_flux != NULL)  free(rho_rk1_flux);	
    if(rho_rk2_flux != NULL)  free(rho_rk2_flux);	
    if(rho_rk3_flux != NULL)  free(rho_rk3_flux);
    if(rhou_rk1_flux != NULL) free(rhou_rk1_flux);	
    if(rhou_rk2_flux != NULL) free(rhou_rk2_flux);	
    if(rhou_rk3_flux != NULL) free(rhou_rk3_flux);    
    if(rhov_rk1_flux != NULL) free(rhov_rk1_flux);	
    if(rhov_rk2_flux != NULL) free(rhov_rk2_flux);	
    if(rhov_rk3_flux != NULL) free(rhov_rk3_flux);
    if(rhow_rk1_flux != NULL) free(rhow_rk1_flux);	
    if(rhow_rk2_flux != NULL) free(rhow_rk2_flux);	
    if(rhow_rk3_flux != NULL) free(rhow_rk3_flux);
    if(rhoE_rk1_flux != NULL) free(rhoE_rk1_flux);	
    if(rhoE_rk2_flux != NULL) free(rhoE_rk2_flux);	
    if(rhoE_rk3_flux != NULL) free(rhoE_rk3_flux);

    /// Free inviscid fluxes
    if(rho_inv_flux != NULL)  free(rho_inv_flux);
    if(rhou_inv_flux != NULL) free(rhou_inv_flux);
    if(rhov_inv_flux != NULL) free(rhov_inv_flux);
    if(rhow_inv_flux != NULL) free(rhow_inv_flux);
    if(rhoE_inv_flux != NULL) free(rhoE_inv_flux);

    /// Free viscous fluxes
    if(rhou_vis_flux != NULL) free(rhou_vis_flux);
    if(rhov_vis_flux != NULL) free(rhov_vis_flux);
    if(rhow_vis_flux != NULL) free(rhow_vis_flux);
    if(rhoE_vis_flux != NULL) free(rhoE_vis_flux);

    /// Free volumetric & external forces
    if(f_rhou_field != NULL) free(f_rhou_field);
    if(f_rhov_field != NULL) free(f_rhov_field);
    if(f_rhow_field != NULL) free(f_rhow_field);
    if(f_rhoE_field != NULL) free(f_rhoE_field);

};





/*
////////// FUNCTIONS //////////

// Initialize u, v, w, P and T variables
void initialize_uvwPT( u, v, w, P, T, grid ):

    # All points
    for i in range( 0, num_grid_x + 2 ):    
        for j in range( 0, num_grid_y + 2 ):    
            for k in range( 0, num_grid_z + 2 ):
                u[i][j][k] = u_tau*np.sin( grid[i][j][k][0] )
                v[i][j][k] = 0.0
                w[i][j][k] = u_tau*np.sin( grid[i][j][k][2] )
                P[i][j][k] = P_ref
                T[i][j][k] = ( 1.0/( rho_ref*R_specific ) )*P[i][j][k]
*/

////////// MAIN //////////
int main(int argc, char** argv) {

    /// Initialize MPI
    MPI_Init(&argc, &argv);

    /// Initialize (construct) flow solver RHEA
    FlowSolverRHEA flow_solver_RHEA("configuration_file.yaml");

/*    
    T = 0.0;
    Tnew = 0.0;

    int maxiter= 20000;


    for(int it=0; it < maxiter; it++)
    {

        for(int i = Tnew.ini_x; i <=  Tnew.fin_x; i++)
            for(int j = Tnew.ini_y; j <=  Tnew.fin_y; j++)
                for(int k = Tnew.ini_z; k <=  Tnew.fin_z; k++){
                    Tnew[I1D(i,j,k)] =  T[I1D(i,j,k)] + dt*kappa*( (2.0/(dom.x[i+1]-dom.x[i-1]))*( (T[I1D(i+1,j,k)] -T[I1D(i,j,k)])/(dom.x[i+1]-dom.x[i]) -  (T[I1D(i,j,k)] -T[I1D(i-1,j,k)])/(dom.x[i]-dom.x[i-1])) + 
                                                                   (2.0/(dom.y[j+1]-dom.y[j-1]))*( (T[I1D(i,j+1,k)] -T[I1D(i,j,k)])/(dom.y[j+1]-dom.y[j]) -  (T[I1D(i,j,k)] -T[I1D(i,j-1,k)])/(dom.y[j]-dom.y[j-1])) +        
                                                                   (2.0/(dom.z[k+1]-dom.z[k-1]))*( (T[I1D(i,j,k+1)] -T[I1D(i,j,k)])/(dom.z[k+1]-dom.z[k]) -  (T[I1D(i,j,k)] -T[I1D(i,j,k-1)])/(dom.z[k]-dom.z[k-1])) );    
                }
    
        Tnew.update();
        //WEST 
        for(int i = topo.iter_bound[_WEST_][_INIX_]; i <= topo.iter_bound[_WEST_][_ENDX_]; i++)
            for(int j = topo.iter_bound[_WEST_][_INIY_]; j <= topo.iter_bound[_WEST_][_ENDY_]; j++)
                for(int k = topo.iter_bound[_WEST_][_INIZ_]; k <=  topo.iter_bound[_WEST_][_ENDZ_]; k++)
                {
                    Tnew[I1D(i,j,k)] = 2*0.0 - Tnew[I1D(i+1,j,k)];
                }

        //EAST 
        for(int i = topo.iter_bound[_EAST_][_INIX_]; i <= topo.iter_bound[_EAST_][_ENDX_]; i++)
            for(int j = topo.iter_bound[_EAST_][_INIY_]; j <= topo.iter_bound[_EAST_][_ENDY_]; j++)
                for(int k = topo.iter_bound[_EAST_][_INIZ_]; k <=  topo.iter_bound[_EAST_][_ENDZ_]; k++)
                {
                    Tnew[I1D(i,j,k)] = 2*1.0 - Tnew[I1D(i-1,j,k)];
                }

        //SOUTH 
        for(int i = topo.iter_bound[_SOUTH_][_INIX_]; i <= topo.iter_bound[_SOUTH_][_ENDX_]; i++)
            for(int j = topo.iter_bound[_SOUTH_][_INIY_]; j <= topo.iter_bound[_SOUTH_][_ENDY_]; j++)
                for(int k = topo.iter_bound[_SOUTH_][_INIZ_]; k <=  topo.iter_bound[_SOUTH_][_ENDZ_]; k++)
                {
                    Tnew[I1D(i,j,k)] = Tnew[I1D(i,j+1,k)];
                }

        //NORTH 
        for(int i = topo.iter_bound[_NORTH_][_INIX_]; i <= topo.iter_bound[_NORTH_][_ENDX_]; i++)
            for(int j = topo.iter_bound[_NORTH_][_INIY_]; j <= topo.iter_bound[_NORTH_][_ENDY_]; j++)
                for(int k = topo.iter_bound[_NORTH_][_INIZ_]; k <=  topo.iter_bound[_NORTH_][_ENDZ_]; k++)
                {
                    Tnew[I1D(i,j,k)] = Tnew[I1D(i,(j-1),k)];
                }

        //BACK 
        for(int i = topo.iter_bound[_BACK_][_INIX_]; i <= topo.iter_bound[_BACK_][_ENDX_]; i++)
            for(int j = topo.iter_bound[_BACK_][_INIY_]; j <= topo.iter_bound[_BACK_][_ENDY_]; j++)
                for(int k = topo.iter_bound[_BACK_][_INIZ_]; k <=  topo.iter_bound[_BACK_][_ENDZ_]; k++)
                {
//                    Tnew[I1D(i,j,0)] = Tnew[I1D(i,j,_lNz_)];
                     Tnew[I1D(i,j,k)] = Tnew[I1D(i,j,k+1)];
                }

        //FRONT 
        for(int i = topo.iter_bound[_FRONT_][_INIX_]; i <= topo.iter_bound[_FRONT_][_ENDX_]; i++)
            for(int j = topo.iter_bound[_FRONT_][_INIY_]; j <= topo.iter_bound[_FRONT_][_ENDY_]; j++)
                for(int k = topo.iter_bound[_FRONT_][_INIZ_]; k <=  topo.iter_bound[_FRONT_][_ENDZ_]; k++)
                {
//                   Tnew[I1D(i,j,_lNz_+1)] = Tnew[I1D(i,j,1)];
                    Tnew[I1D(i,j,k)] = Tnew[I1D(i,j,k-1)];
               }



//        for(int i = topo.iter_common[_ALL_][_INIX_]; i <= topo.iter_common[_ALL_][_INIX_]; i++)
//            for(int j =  topo.iter_common[_ALL_][_INIY_]; j <=  topo.iter_common[_ALL_][_INIY_]; j++)
//                for(int k =  topo.iter_common[_ALL_][_INIZ_]; k <=  topo.iter_common[_ALL_][_INIZ_]; k++){
//                    T[I1D(i,j,k)] = Tnew[I1D(i,j,k)];
//                }
 
        for(int l=0; l < topo.getSize(); l++)
        {
            T[l] = Tnew[l];
        }
 
    }



    if(topo.getRank() == 0 && _DEBUG_ ==0 )
        for(int l=0; l < topo.getSize(); l++)
            cout<<l<<":  "<<Tnew[l]<<endl;
*/


    /// Finalize (destruct) flow solver RHEA ... destructor is called automatically

    /// Finalize MPI
    MPI_Finalize();

}

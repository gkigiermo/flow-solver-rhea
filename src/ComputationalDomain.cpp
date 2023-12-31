#include "ComputationalDomain.hpp"

using namespace std;

ComputationalDomain::ComputationalDomain(double sizex, double sizey, double sizez, double orig_x, double orig_y, double orig_z, double stretch_x, double stretch_y, double stretch_z, int ncellsx, int ncellsy, int ncellsz)
{
    gNx=ncellsx;
    gNy=ncellsy;
    gNz=ncellsz;

    L_x=sizex;
    L_y=sizey;
    L_z=sizez;

    x_0 = orig_x;
    y_0 = orig_y;
    z_0 = orig_z;

    A_x = stretch_x;
    A_y = stretch_y;
    A_z = stretch_z;

    // Position of the Global Cells
    globx = new double[gNx + 2];
    globy = new double[gNy + 2];
    globz = new double[gNz + 2];


    calculateGlobalGrid();
}

void ComputationalDomain::calculateGlobalGrid() {

    double eta_x, eta_y, eta_z;

    // Mesh in x
    for( int i = 0; i < (gNx + 2); i++ ) {
        eta_x = (i - 0.5)/gNx;
        globx[i] = x_0 + L_x*eta_x + A_x*( 0.5*L_x - L_x*eta_x )*( 1.0 - eta_x )*eta_x;

        if(i == 0) {
            eta_x = ( 1.0 - 0.5 )/gNx;
            globx[i] = x_0 - ( L_x*eta_x + A_x*( 0.5*L_x - L_x*eta_x )*( 1.0 - eta_x )*eta_x );
        }
        if( i == ( gNx + 1 ) ) {
            eta_x = ( gNx - 0.5 )/gNx;
            globx[i] = x_0 + 2.0*L_x - ( L_x*eta_x + A_x*( 0.5*L_x - L_x*eta_x )*( 1.0 - eta_x )*eta_x );
        }
    }

    // Mesh in y
    for( int j = 0; j < (gNy + 2); j++ ) {
        eta_y = (j - 0.5)/gNy;
        globy[j] = y_0 + L_y*eta_y + A_y*( 0.5*L_y - L_y*eta_y )*( 1.0 - eta_y )*eta_y;

        if( j == 0 ) {
            eta_y = ( 1.0 - 0.5 )/gNy;
            globy[j] = y_0 - ( L_y*eta_y + A_y*( 0.5*L_y - L_y*eta_y )*( 1.0 - eta_y )*eta_y );
        }
        if( j == (gNy + 1) ) {
            eta_y = ( gNy - 0.5 )/gNy;
            globy[j] = y_0 + 2.0*L_y - ( L_y*eta_y + A_y*( 0.5*L_y - L_y*eta_y )*( 1.0 - eta_y )*eta_y );
        }
    }

    // Mesh in z
    for( int k = 0; k < (gNz + 2); k++ ) {
        eta_z = (k - 0.5)/gNz;
        globz[k] = z_0 + L_z*eta_z + A_z*( 0.5*L_z - L_z*eta_z )*( 1.0 - eta_z )*eta_z;

        if( k == 0 ) {
            eta_z = ( 1.0 - 0.5 )/gNz;
            globz[k] = z_0 - ( L_z*eta_z + A_z*( 0.5*L_z - L_z*eta_z )*( 1.0 - eta_z )*eta_z );
        }
        if( k == (gNz + 1) ) {
            eta_z = ( gNz - 0.5 )/gNz;
            globz[k] = z_0 + 2.0*L_z - ( L_z*eta_z + A_z*( 0.5*L_z - L_z*eta_z )*( 1.0 - eta_z )*eta_z );
        }
    }

}

void ComputationalDomain::calculateLocalGrid(int lNx, int lNy, int lNz) {

    x = new double[lNx];
    y = new double[lNy];
    z = new double[lNz];

    Nx=lNx;
    Ny=lNy;
    Nz=lNz;
}

void ComputationalDomain::printDomain()
{

    cout<<"Global x: "<<endl;
    for(int i=0;i<gNx+2;i++)
        cout<<i<<" : "<<globx[i]<<" ; ";
    cout<<endl;

    cout<<"Global y: "<<endl;
    for(int i=0;i<gNy+2;i++)
        cout<<i<<" : "<<globy[i]<<" ; ";
    cout<<endl;

    cout<<"Global z: "<<endl;
    for(int i=0;i<gNz+2;i++)
        cout<<i<<" : "<<globz[i]<<" ; ";
    cout<<endl;

    cout<<"Local x: "<<endl;
    for(int i=0;i<Nx;i++)
        cout<<i<<" : "<<x[i]<<" ; ";
    cout<<endl;

    cout<<"Local y: "<<endl;
    for(int i=0;i<Ny;i++)
        cout<<i<<" : "<<y[i]<<" ; ";
    cout<<endl;

    cout<<"Local z: "<<endl;
    for(int i=0;i<Nz;i++)
        cout<<i<<" : "<<z[i]<<" ; ";
    cout<<endl;


}   

void ComputationalDomain::setBocos(int b[6])
{
    bc[_WEST_]  = b[_WEST_];
    bc[_EAST_]  = b[_EAST_];
    bc[_SOUTH_] = b[_SOUTH_];
    bc[_NORTH_] = b[_NORTH_];
    bc[_BACK_]  = b[_BACK_];
    bc[_FRONT_] = b[_FRONT_];
}


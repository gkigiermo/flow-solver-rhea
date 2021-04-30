/*ESTE FICHERO SERA CREADO POR UN PREPROCESS*/

#define R_INDX(i, j, k, N, M) (k*(N*M) + (j*N) + i)

//Common parameters

//Surface comms
#define _WEST_  0
#define _EAST_  1
#define _SOUTH_ 2
#define _NORTH_ 3
#define _BACK_  4
#define _FRONT_ 5

//Lineal comms
#define _WEST_S_  6
#define _WEST_N_  7
#define _WEST_B_  8
#define _WEST_F_  9

#define _EAST_S_  10
#define _EAST_N_  11
#define _EAST_B_  12
#define _EAST_F_  13

#define _SOUTH_B_  14
#define _SOUTH_F_  15
#define _NORTH_B_  16
#define _NORTH_F_  17

//Point comms
#define _WEST_S_B_  18
#define _WEST_N_B_  19
#define _WEST_S_F_  20
#define _WEST_N_F_  21

#define _EAST_S_B_  22
#define _EAST_N_B_  23
#define _EAST_S_F_  24
#define _EAST_N_F_  25

#define _INIX_ 0
#define _ENDX_ 1
#define _INIY_ 2
#define _ENDY_ 3
#define _INIZ_ 4
#define _ENDZ_ 5

#define _INNER_ 0
#define _ALL_   1

#define _NO_BOCO_   0
#define _DIRICHLET_ 1
#define _NEUMANN_   2
#define _PERIODIC_  3
#define _SUBSONIC_INFLOW_  4
#define _SUBSONIC_OUTFLOW_ 5

#define I1D(i, j, k) ((k)*(_lNx_*_lNy_) + ((j)*_lNx_) + (i))
//#define I1D(i, j, k) ((i)*(_lNz_*_lNy_) + ((j)*_lNz_) + (k))

#define _NO_NEIGHBOUR_ -1

/*ESTE FICHERO SERA CREADO POR UN PREPROCESS*/

#define RHEA_NX 10
#define RHEA_NY 1
#define RHEA_NZ 1

#define R_INDX(i, j, k, N, M) (k*(N*M) + (j*N) + i)

//Common parameters
#define _WEST_  0
#define _EAST_  1
#define _SOUTH_ 2
#define _NORTH_ 3
#define _BACK_  4
#define _FRONT_ 5

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

#define I1D(i, j, k) ((k)*(_lNx_*_lNy_) + ((j)*_lNx_) + (i))



#define _NO_NEIGHBOUR_ -1

#include "ParallelTopology.hpp"

ParallelTopology::ParallelTopology(ComputationalDomain* dom, int nprocsx, int nprocsy, int nprocsz)
{
    mymesh = dom;

    RHEA_3DCOMM = MPI_COMM_WORLD;    
    //MPI_Comm_dup(MPI_COMM_WORLD, &RHEA_3DCOMM);
    MPI_Comm_rank(RHEA_3DCOMM, &rank);
    MPI_Comm_size(RHEA_3DCOMM, &np);


    npx= nprocsx;
    npy= nprocsy;
    npz= nprocsz;

    if( np != npx*npy*npz ){

        cout<<"Mismatched in the ComputationalDomain partition and number of procs launched"<<endl;
        MPI_Abort(RHEA_3DCOMM,0);
    }

    //cout<<"testing 1"<<endl;
    MPI_Barrier(MPI_COMM_WORLD);
    int cellsx=dom->getGNx();
    int localNx= cellsx/npx;
    int divx = cellsx%npx;

    if(rank%npx < divx) 
        lNx=localNx+1;
    else
        lNx=localNx;


    lNx = lNx + 2;
    //cout<<"testing 2"<<endl;
    int cellsy=dom->getGNy();
    int localNy= cellsy/npy;
    int divy = cellsy%npy;

    int plane_rank= rank%(npx*npy);

    if(plane_rank/npx < divy) 
        lNy=localNy+1;
    else
        lNy=localNy;


    lNy = lNy + 2;
    //cout<<"testing 3"<<endl;
    int cellsz=dom->getGNz();
    int localNz= cellsz/npz;
    int divz = cellsz%npz;

    if(rank/(npx*npy) < divz) 
        lNz=localNz+1;
    else
        lNz=localNz;

    lNz = lNz + 2;


    len = lNx * lNy * lNz; 



    int mpirank;

    MPI_Comm_rank(MPI_COMM_WORLD, &mpirank);
    //cout<<" Rank Rhea "<<rank<<endl;
    //cout<<" Rank MPI  "<<mpirank<<endl;




    //Finding surface neighbours
    //By default there is no connection with other sub-ComputationalDomains
    for(int nb = 0 ; nb < 26 ; nb++)
        neighb[nb] = _NO_NEIGHBOUR_;

    is_inix=0;
    is_iniy=0;
    is_iniz=0;
    is_endx=0;
    is_endy=0;
    is_endz=0;

    int periodic_x;
    if(dom->getBoco(_WEST_) == _PERIODIC_)
        periodic_x = 1;
    else
        periodic_x = 0;

    //Direction X 
    //WEST
    if( rank%npx !=  0 ){
        neighb[_WEST_] = rank - 1;
    }
    else{
        if( periodic_x == 1 ) {
            neighb[_WEST_] = rank + npx - 1;
        }
        is_inix=1;
    }
    //EAST
    if(rank%npx != (npx - 1) ){
        neighb[_EAST_] = rank + 1;
    }
    else{
        if( periodic_x == 1 ){
            neighb[_EAST_] = rank - npx + 1;
            is_endx = 1;
        }
    }

   int periodic_y;
   if(dom->getBoco(_SOUTH_) == _PERIODIC_)
        periodic_y = 1;
   else
        periodic_y = 0;


    //SOUTH
    plane_rank = rank%(npx*npy);
    if(plane_rank/npx != 0)
    {
        neighb[_SOUTH_] = rank - npx;
    }
    else{
        if(periodic_y == 1) {
            neighb[_SOUTH_] = rank + (npy-1)*npx;
            is_iniy = 1;
        }
    }

    //NORTH
    if(plane_rank/npx != npy - 1){
        neighb[_NORTH_] = rank + npx;
    }
    else{
        if(periodic_y == 1) {
            neighb[_NORTH_] = rank - (npy-1)*npx;
            is_endy = 1;
        }
    }

    //BACK
    int periodic_z;
    if(dom->getBoco(_BACK_) == _PERIODIC_)
        periodic_z = 1;
    else
        periodic_z = 0;


    proc_z_start = 0;
    proc_z_end   = 0;


    plane_rank = rank/(npx*npy);
    if(plane_rank != 0){
        neighb[_BACK_] = rank - npx*npy;
    }
    else{
        if(periodic_z == 1 ){
            neighb[_BACK_] = rank + (npx*npy)*(npz-1);
            proc_z_start = 1;
            is_iniz = 1;
        }
    }

    //FRONT;
    if(plane_rank != npz-1){
        neighb[_FRONT_] = rank + npx*npy;
    }
    else{
        if(periodic_z == 1){
            neighb[_FRONT_] = rank - (npx*npy)*(npz-1);
            proc_z_end = 1;
            is_endz = 1;
        }
    }


    //Calculating Id offsets /* this should be updated to include the case with imbalanced partitions*/

    //offset x
    int factx = rank%npx;
    offx = factx*(lNx - 2 );

    //offset y
    int facty = (rank%(npx*npy));
    facty = (facty/npx);
    offy = facty*(lNy - 2 );

    //offset z
    int factz = (rank/(npx*npy));
    offz = factz*(lNz -2 );


    create_common_iters();
    create_halo_iters();
    create_basic_bound_iters();
    create_toRecv_iters();
    create_toSend_iters();
    create_global_iters();
    create_comm_arrays();

    dom->calculateLocalGrid(lNx,lNy,lNz);

    int l=0;
    for(int i=iter_glob_ind[_INIX_];i<=iter_glob_ind[_ENDX_];i++){
        dom->set_x( l, dom->getGlobx(i) );
        l++;
    }

    l=0;
    for(int j=iter_glob_ind[_INIY_];j<=iter_glob_ind[_ENDY_];j++){
        dom->set_y( l, dom->getGloby(j) );
        l++;
    }

    l=0;
    for(int k=iter_glob_ind[_INIZ_];k<=iter_glob_ind[_ENDZ_];k++){
        dom->set_z( l, dom->getGlobz(k) );
        l++;
    }

    // To create the extra connections

    exchange_2nd_level_neighbour_info();
    find_extra_neighbours();
    create_complex_bound_iters();
    create_complex_halo_iters();
    create_complex_toRecv_iters();
    create_complex_toSend_iters();
    create_complex_comm_arrays();
}

void ParallelTopology::exchange_2nd_level_neighbour_info()
{
    int pack[6];
    
    int count_nb=0; // number of 1st level neighbours

    MPI_Barrier(RHEA_3DCOMM);

    //Counting with whom am I connected
    for(int nb =_WEST_ ; nb <= _FRONT_ ; nb++){
        if( neighb[nb] != _NO_NEIGHBOUR_ ){
            off_nb[nb] = count_nb*6;
            count_nb++;
        }
    }
     
    //Packs the message that is sent to all the 1st level neighbours
    for(int l=_WEST_;l<=_FRONT_;l++)
        pack[l] = neighb[l];

    //int* recv_data;
    info_2nd = new int[count_nb*6];

    
    MPI_Request req_s[6];
    MPI_Request req_r[6];

    MPI_Status  stat_s[6];
    MPI_Status  stat_r[6];



    if(getNB(_WEST_) != _NO_NEIGHBOUR_ )
        MPI_Isend(&pack[0], 6, MPI_INT, getNB(_WEST_), 0, RHEA_3DCOMM, &req_s[_WEST_]);

    if(getNB(_EAST_) != _NO_NEIGHBOUR_ )
        MPI_Isend(&pack[0], 6, MPI_INT, getNB(_EAST_), 0, RHEA_3DCOMM, &req_s[_EAST_]);

    if(getNB(_SOUTH_) != _NO_NEIGHBOUR_ )
        MPI_Isend(&pack[0], 6, MPI_INT, getNB(_SOUTH_), 0, RHEA_3DCOMM, &req_s[_SOUTH_]);

    if(getNB(_NORTH_) != _NO_NEIGHBOUR_ )
        MPI_Isend(&pack[0], 6, MPI_INT, getNB(_NORTH_), 0, RHEA_3DCOMM, &req_s[_NORTH_]);

    if(getNB(_BACK_) != _NO_NEIGHBOUR_ )
        MPI_Isend(&pack[0], 6, MPI_INT, getNB(_BACK_), 0, RHEA_3DCOMM, &req_s[_BACK_]);

    if(getNB(_FRONT_) != _NO_NEIGHBOUR_ )
        MPI_Isend(&pack[0], 6, MPI_INT, getNB(_FRONT_), 0, RHEA_3DCOMM, &req_s[_FRONT_]);


    if(getNB(_WEST_) != _NO_NEIGHBOUR_ )
        MPI_Irecv(&info_2nd[off_nb[_WEST_]], 6, MPI_INT, getNB(_WEST_), 0, RHEA_3DCOMM, &req_r[_WEST_]);

    if(getNB(_EAST_) != _NO_NEIGHBOUR_ )
        MPI_Irecv(&info_2nd[off_nb[_EAST_]], 6, MPI_INT, getNB(_EAST_), 0, RHEA_3DCOMM, &req_r[_EAST_]);

    if(getNB(_SOUTH_) != _NO_NEIGHBOUR_ )
        MPI_Irecv(&info_2nd[off_nb[_SOUTH_]], 6, MPI_INT, getNB(_SOUTH_), 0, RHEA_3DCOMM, &req_r[_SOUTH_]);

    if(getNB(_NORTH_) != _NO_NEIGHBOUR_ )
        MPI_Irecv(&info_2nd[off_nb[_NORTH_]], 6, MPI_INT, getNB(_NORTH_), 0, RHEA_3DCOMM, &req_r[_NORTH_]);

    if(getNB(_BACK_) != _NO_NEIGHBOUR_ )
        MPI_Irecv(&info_2nd[off_nb[_BACK_]], 6, MPI_INT, getNB(_BACK_), 0, RHEA_3DCOMM, &req_r[_BACK_]);

    if(getNB(_FRONT_) != _NO_NEIGHBOUR_ )
        MPI_Irecv(&info_2nd[off_nb[_FRONT_]], 6, MPI_INT, getNB(_FRONT_), 0, RHEA_3DCOMM, &req_r[_FRONT_]);


    if(getNB(_WEST_) != _NO_NEIGHBOUR_ ){
        MPI_Wait(&req_r[_WEST_], &stat_r[_WEST_]);
        MPI_Wait(&req_s[_WEST_], &stat_s[_WEST_]);
    }
    if(getNB(_EAST_) != _NO_NEIGHBOUR_ ){
        MPI_Wait(&req_r[_EAST_], &stat_r[_EAST_]);
        MPI_Wait(&req_s[_EAST_], &stat_s[_EAST_]);
    }
    if(getNB(_SOUTH_) != _NO_NEIGHBOUR_ ){
        MPI_Wait(&req_r[_SOUTH_], &stat_r[_SOUTH_]);
        MPI_Wait(&req_s[_SOUTH_], &stat_s[_SOUTH_]);
    }
    if(getNB(_NORTH_) != _NO_NEIGHBOUR_ ){
        MPI_Wait(&req_r[_NORTH_], &stat_r[_NORTH_]);
        MPI_Wait(&req_s[_NORTH_], &stat_s[_NORTH_]);
    }
    if(getNB(_BACK_) != _NO_NEIGHBOUR_ ){
        MPI_Wait(&req_r[_BACK_], &stat_r[_BACK_]);
        MPI_Wait(&req_s[_BACK_], &stat_s[_BACK_]);
    }
    if(getNB(_FRONT_) != _NO_NEIGHBOUR_ ){
        MPI_Wait(&req_r[_FRONT_], &stat_r[_FRONT_]);
        MPI_Wait(&req_s[_FRONT_], &stat_s[_FRONT_]);
    }
}


void ParallelTopology::find_extra_neighbours()
{

// Lineal neighbours 

    //WEST SOUTH
    neighb[_WEST_S_] = neighb[_WEST_]; 
    if( neighb[_WEST_]  != _NO_NEIGHBOUR_ && neighb[_SOUTH_] != _NO_NEIGHBOUR_)
    {    
        neighb[_WEST_S_] =  info_2nd[ off_nb[_WEST_] + _SOUTH_ ]; 
    }
    else if ( neighb[_SOUTH_] != _NO_NEIGHBOUR_)
    {
        neighb[_WEST_S_] =  neighb[_SOUTH_]; 
    } 

    //WEST NORTH 
    neighb[_WEST_N_] = neighb[_WEST_]; 
    if( neighb[_WEST_]  != _NO_NEIGHBOUR_ && neighb[_NORTH_] != _NO_NEIGHBOUR_)
    {
        neighb[_WEST_N_] =  info_2nd[ off_nb[_WEST_] + _NORTH_ ]; 
    }
    else if ( neighb[_NORTH_] != _NO_NEIGHBOUR_)
    {
        neighb[_WEST_N_] =  neighb[_NORTH_]; 
    } 


    //WEST BACK 
    neighb[_WEST_B_] = neighb[_WEST_]; 
    if( neighb[_WEST_]  != _NO_NEIGHBOUR_ && neighb[_BACK_] != _NO_NEIGHBOUR_)
    {
        neighb[_WEST_B_] =  info_2nd[ off_nb[_WEST_] + _BACK_ ]; 
    }
    else if ( neighb[_BACK_] != _NO_NEIGHBOUR_)
    {
        neighb[_WEST_B_] =  neighb[_BACK_]; 
    } 


    //WEST FRONT 
    neighb[_WEST_F_] = neighb[_WEST_]; 
    if( neighb[_WEST_]  != _NO_NEIGHBOUR_ && neighb[_FRONT_] != _NO_NEIGHBOUR_)
    {
        neighb[_WEST_F_] =  info_2nd[ off_nb[_WEST_] + _FRONT_ ]; 
    }
    else if ( neighb[_FRONT_] != _NO_NEIGHBOUR_)
    {
        neighb[_WEST_F_] =  neighb[_FRONT_]; 
    } 



    //EAST SOUTH 
    neighb[_EAST_S_] = neighb[_EAST_]; 
    if( neighb[_EAST_]  != _NO_NEIGHBOUR_ && neighb[_SOUTH_] != _NO_NEIGHBOUR_)
    {
        neighb[_EAST_S_] = info_2nd[ off_nb[_EAST_] + _SOUTH_ ]; 
    }
    else if ( neighb[_SOUTH_] != _NO_NEIGHBOUR_)
    {
        neighb[_EAST_S_] =  neighb[_SOUTH_]; 
    } 


    //EAST NORTH 
    neighb[_EAST_N_] = neighb[_EAST_]; 
    if( neighb[_EAST_]  != _NO_NEIGHBOUR_ && neighb[_NORTH_] != _NO_NEIGHBOUR_)
    {
        neighb[_EAST_N_] = info_2nd[ off_nb[_EAST_] + _NORTH_ ]; 

    }
    else if ( neighb[_NORTH_] != _NO_NEIGHBOUR_)
    {
        neighb[_EAST_N_] =  neighb[_NORTH_]; 
    } 


    //EAST BACK 
    neighb[_EAST_B_] = neighb[_EAST_]; 
    if( neighb[_EAST_]  != _NO_NEIGHBOUR_ && neighb[_BACK_] != _NO_NEIGHBOUR_)
    {
        neighb[_EAST_B_] = info_2nd[ off_nb[_EAST_] + _BACK_ ]; 

    }
    else if ( neighb[_BACK_] != _NO_NEIGHBOUR_)
    {
        neighb[_EAST_B_] =  neighb[_BACK_]; 
    } 


    //EAST FRONT 
    neighb[_EAST_F_] = neighb[_EAST_]; 
    if( neighb[_EAST_]  != _NO_NEIGHBOUR_ && neighb[_FRONT_] != _NO_NEIGHBOUR_)
    {
        neighb[_EAST_F_] = info_2nd[ off_nb[_EAST_] + _FRONT_ ]; 

    }
    else if ( neighb[_FRONT_] != _NO_NEIGHBOUR_)
    {
        neighb[_EAST_F_] =  neighb[_FRONT_]; 
    } 



    //SOUTH BACK
    neighb[_SOUTH_B_] = neighb[_SOUTH_]; 
    if( neighb[_SOUTH_]  != _NO_NEIGHBOUR_ && neighb[_BACK_] != _NO_NEIGHBOUR_)
    {
        neighb[_SOUTH_B_] = info_2nd[ off_nb[_SOUTH_] + _BACK_ ];  
    }
    else if ( neighb[_BACK_] != _NO_NEIGHBOUR_)
    {
        neighb[_SOUTH_B_] =  neighb[_BACK_]; 
    } 


    //SOUTH FRONT
    neighb[_SOUTH_F_] = neighb[_SOUTH_]; 
    if( neighb[_SOUTH_]  != _NO_NEIGHBOUR_ && neighb[_FRONT_] != _NO_NEIGHBOUR_)
    {
        neighb[_SOUTH_F_] = info_2nd[ off_nb[_SOUTH_] + _FRONT_ ];  

    }
    else if ( neighb[_FRONT_] != _NO_NEIGHBOUR_)
    {
        neighb[_SOUTH_F_] =  neighb[_FRONT_]; 
    } 


    //NORTH BACK
    neighb[_NORTH_B_] = neighb[_NORTH_]; 
    if( neighb[_NORTH_]  != _NO_NEIGHBOUR_ && neighb[_BACK_] != _NO_NEIGHBOUR_)
    {
        neighb[_NORTH_B_] = info_2nd[ off_nb[_NORTH_] + _BACK_ ]; 
    }
    else if ( neighb[_BACK_] != _NO_NEIGHBOUR_)
    {
        neighb[_NORTH_B_] =  neighb[_BACK_]; 
    } 


    //NORTH FRONT
    neighb[_NORTH_F_] = neighb[_NORTH_]; 
    if( neighb[_NORTH_]  != _NO_NEIGHBOUR_ && neighb[_FRONT_] != _NO_NEIGHBOUR_)
    {
        neighb[_NORTH_F_] =  info_2nd[ off_nb[_NORTH_] + _FRONT_ ]; 
    }
    else if ( neighb[_FRONT_] != _NO_NEIGHBOUR_)
    {
        neighb[_NORTH_F_] =  neighb[_FRONT_]; 
    } 



    //3er level conexions
    int* lconx;
    int pack[6];
    //Packs the message that is sent to all the 1st level neighbours
    for(int l=_WEST_;l<=_FRONT_;l++)
        pack[l] = neighb[l];

   
    lconx = new int[6*npx*npy*npz]; // For each proc saves the 6 neighbours

    MPI_Allgather(pack, 6, MPI_INT, lconx, 6, MPI_INT, MPI_COMM_WORLD);


    //Point Neighbours

    //WEST SOUTH and NORTH BACK

    neighb[_WEST_S_B_] = neighb[_WEST_B_]; 
    neighb[_WEST_N_B_] = neighb[_WEST_B_]; 
    if( neighb[_WEST_B_] != _NO_NEIGHBOUR_)
    {
        if( lconx[neighb[_WEST_B_]*6 + _SOUTH_ ] != _NO_NEIGHBOUR_) 
            neighb[_WEST_S_B_] = lconx[neighb[_WEST_B_]*6 + _SOUTH_ ] ; 

        if( lconx[neighb[_WEST_B_]*6 + _NORTH_ ] != _NO_NEIGHBOUR_) 
            neighb[_WEST_N_B_] = lconx[neighb[_WEST_B_]*6 + _NORTH_ ] ; 
    }

    //WEST SOUTH and NORTH FRONT
    neighb[_WEST_S_F_] = neighb[_WEST_F_];
    neighb[_WEST_N_F_] = neighb[_WEST_F_];
    if( neighb[_WEST_F_] != _NO_NEIGHBOUR_)
    {
        if( lconx[neighb[_WEST_F_]*6 + _SOUTH_ ] != _NO_NEIGHBOUR_) 
            neighb[_WEST_S_F_] = lconx[neighb[_WEST_F_]*6 + _SOUTH_ ] ; 
 
        if( lconx[neighb[_WEST_F_]*6 + _NORTH_ ] != _NO_NEIGHBOUR_) 
            neighb[_WEST_N_F_] = lconx[neighb[_WEST_F_]*6 + _NORTH_ ] ; 
    }

    //EAST SOUTH and NORTH BACK
    neighb[_EAST_S_B_] = neighb[_EAST_B_];  
    neighb[_EAST_N_B_] = neighb[_EAST_B_];  
    if( neighb[_EAST_B_] != _NO_NEIGHBOUR_)
    {
        if( lconx[neighb[_EAST_B_]*6 + _SOUTH_ ] != _NO_NEIGHBOUR_) 
            neighb[_EAST_S_B_] = lconx[neighb[_EAST_B_]*6 + _SOUTH_ ] ; 

        if( lconx[neighb[_EAST_B_]*6 + _NORTH_ ] != _NO_NEIGHBOUR_) 
            neighb[_EAST_N_B_] = lconx[neighb[_EAST_B_]*6 + _NORTH_ ] ; 
    }

    //EAST SOUTH and NORTH FRONT 
    neighb[_EAST_S_F_] = neighb[_EAST_F_];  
    neighb[_EAST_N_F_] = neighb[_EAST_F_];  
    if( neighb[_EAST_F_] != _NO_NEIGHBOUR_)
    {
        if( lconx[neighb[_EAST_F_]*6 + _SOUTH_ ] != _NO_NEIGHBOUR_) 
           neighb[_EAST_S_F_] = lconx[neighb[_EAST_F_]*6 + _SOUTH_ ] ; 
 
        if( lconx[neighb[_EAST_F_]*6 + _NORTH_ ] != _NO_NEIGHBOUR_) 
            neighb[_EAST_N_F_] = lconx[neighb[_EAST_F_]*6 + _NORTH_ ] ; 
    }


    delete[] lconx; 

}


void ParallelTopology::printCommScheme(int proc)
{

    if(rank==proc){
        cout<<endl;
        cout<<"rank "<<rank<<endl;
        cout<<"lNx  "<<lNx<<" lNy "<<lNy<<" lNz "<<lNz<<endl;
        cout<<" X : WEST  "<<neighb[0]<<" EAST  "<<neighb[1]<<endl;
        cout<<" Y : SOUTH "<<neighb[2]<<" NORTH "<<neighb[3]<<endl;
        cout<<" Z : BACK  "<<neighb[4]<<" FRONT "<<neighb[5]<<endl;

     }
    MPI_Barrier(RHEA_3DCOMM);
}

void ParallelTopology::printCommSchemeToFile(int proc)
{
    char filename[100];
    sprintf(filename,"topo-%d.info",rank);
    if(rank==proc){
        ofstream myfile (filename);
        if(myfile.is_open()){
            myfile<<"--------"<<endl;
            myfile<<" Processor"<<endl;
            myfile<<" Rank "<<rank<<endl;
            myfile<<"--------"<<endl;
            myfile<<" Number of local cells"<<endl;
            myfile<<" lNx "<<lNx<<endl;
            myfile<<" lNy "<<lNy<<endl;
            myfile<<" lNz "<<lNz<<endl;
            myfile<<"--------"<<endl;
            myfile<<" Neighbour ids: "<<endl;
            myfile<<" X  : WEST  "<<neighb[_WEST_] <<" EAST  "<<neighb[_EAST_]<<endl;
            myfile<<" Y  : SOUTH "<<neighb[_SOUTH_]<<" NORTH "<<neighb[_NORTH_]<<endl;
            myfile<<" Z  : BACK  "<<neighb[_BACK_] <<" FRONT "<<neighb[_FRONT_]<<endl;
            myfile<<" 2nd level connections: "<<endl;
            myfile<<" XY : WEST  SOUTH  "<<neighb[_WEST_S_]  <<"  WEST NORTH "<<neighb[_WEST_N_]<<endl;
            myfile<<" XZ : WEST  BACK   "<<neighb[_WEST_B_]  <<"  WEST FRONT "<<neighb[_WEST_F_]<<endl;
            myfile<<" XY : EAST  SOUTH  "<<neighb[_EAST_S_]  <<"  EAST NORTH "<<neighb[_EAST_N_]<<endl;
            myfile<<" XZ : EAST  BACK   "<<neighb[_EAST_B_]  <<"  EAST FRONT "<<neighb[_EAST_F_]<<endl;
            myfile<<" YZ : SOUTH BACK   "<<neighb[_SOUTH_B_] <<" SOUTH FRONT "<<neighb[_SOUTH_F_]<<endl;
            myfile<<" YZ : NORTH BACK   "<<neighb[_NORTH_B_] <<" NORTH FRONT "<<neighb[_NORTH_F_]<<endl;
            myfile<<" 3er level connections: "<<endl;
            myfile<<" XYZ : WEST SOUTH BACK  "<<neighb[_WEST_S_B_]  <<"  WEST NORTH BACK  "<<neighb[_WEST_N_B_]<<endl;
            myfile<<" XYZ : WEST SOUTH FRONT "<<neighb[_WEST_S_F_]  <<"  WEST NORTH FRONT "<<neighb[_WEST_N_F_]<<endl;
            myfile<<" XYZ : EAST SOUTH BACK  "<<neighb[_EAST_S_B_]  <<"  EAST NORTH BACK  "<<neighb[_EAST_N_B_]<<endl;
            myfile<<" XYZ : EAST SOUTH FRONT "<<neighb[_EAST_S_F_]  <<"  EAST NORTH FRONT "<<neighb[_EAST_N_F_]<<endl;
            myfile<<"--------"<<endl;
            myfile<<" Internal iterators: "<<endl;
            myfile<<" INI_X  "<<iter_common[_INNER_][_INIX_]<<" END_X "<<iter_common[_INNER_][_ENDX_]<<endl;
            myfile<<" INI_Y  "<<iter_common[_INNER_][_INIY_]<<" END_Y "<<iter_common[_INNER_][_ENDY_]<<endl;
            myfile<<" INI_Z  "<<iter_common[_INNER_][_INIZ_]<<" END_Z "<<iter_common[_INNER_][_ENDZ_]<<endl;
            myfile<<" Total ComputationalDomain iterators: "<<endl;
            myfile<<" INI_X  "<<iter_common[_ALL_][_INIX_]<<" END_X "<<iter_common[_ALL_][_ENDX_]<<endl;
            myfile<<" INI_Y  "<<iter_common[_ALL_][_INIY_]<<" END_Y "<<iter_common[_ALL_][_ENDY_]<<endl;
            myfile<<" INI_Z  "<<iter_common[_ALL_][_INIZ_]<<" END_Z "<<iter_common[_ALL_][_ENDZ_]<<endl;
            myfile<<"--------"<<endl;
            myfile<<" Boundary iterators: "<<endl;
            myfile<<" WEST: "<<endl;
            myfile<<" INI_X  "<<iter_bound[_WEST_][_INIX_]<<" END_X "<<iter_bound[_WEST_][_ENDX_]<<endl;
            myfile<<" INI_Y  "<<iter_bound[_WEST_][_INIY_]<<" END_Y "<<iter_bound[_WEST_][_ENDY_]<<endl;
            myfile<<" INI_Z  "<<iter_bound[_WEST_][_INIZ_]<<" END_Z "<<iter_bound[_WEST_][_ENDZ_]<<endl;
            myfile<<endl;
            myfile<<" EAST: "<<endl;
            myfile<<" INI_X  "<<iter_bound[_EAST_][_INIX_]<<" END_X "<<iter_bound[_EAST_][_ENDX_]<<endl;
            myfile<<" INI_Y  "<<iter_bound[_EAST_][_INIY_]<<" END_Y "<<iter_bound[_EAST_][_ENDY_]<<endl;
            myfile<<" INI_Z  "<<iter_bound[_EAST_][_INIZ_]<<" END_Z "<<iter_bound[_EAST_][_ENDZ_]<<endl;
            myfile<<endl;
            myfile<<" SOUTH: "<<endl;
            myfile<<" INI_X  "<<iter_bound[_SOUTH_][_INIX_]<<" END_X "<<iter_bound[_SOUTH_][_ENDX_]<<endl;
            myfile<<" INI_Y  "<<iter_bound[_SOUTH_][_INIY_]<<" END_Y "<<iter_bound[_SOUTH_][_ENDY_]<<endl;
            myfile<<" INI_Z  "<<iter_bound[_SOUTH_][_INIZ_]<<" END_Z "<<iter_bound[_SOUTH_][_ENDZ_]<<endl;
            myfile<<endl;
            myfile<<" NORTH: "<<endl;
            myfile<<" INI_X  "<<iter_bound[_NORTH_][_INIX_]<<" END_X "<<iter_bound[_NORTH_][_ENDX_]<<endl;
            myfile<<" INI_Y  "<<iter_bound[_NORTH_][_INIY_]<<" END_Y "<<iter_bound[_NORTH_][_ENDY_]<<endl;
            myfile<<" INI_Z  "<<iter_bound[_NORTH_][_INIZ_]<<" END_Z "<<iter_bound[_NORTH_][_ENDZ_]<<endl;
            myfile<<endl;
            myfile<<" BACK: "<<endl;
            myfile<<" INI_X  "<<iter_bound[_BACK_][_INIX_]<<" END_X "<<iter_bound[_BACK_][_ENDX_]<<endl;
            myfile<<" INI_Y  "<<iter_bound[_BACK_][_INIY_]<<" END_Y "<<iter_bound[_BACK_][_ENDY_]<<endl;
            myfile<<" INI_Z  "<<iter_bound[_BACK_][_INIZ_]<<" END_Z "<<iter_bound[_BACK_][_ENDZ_]<<endl;
            myfile<<endl;
            myfile<<" FRONT: "<<endl;
            myfile<<" INI_X  "<<iter_bound[_FRONT_][_INIX_]<<" END_X "<<iter_bound[_FRONT_][_ENDX_]<<endl;
            myfile<<" INI_Y  "<<iter_bound[_FRONT_][_INIY_]<<" END_Y "<<iter_bound[_FRONT_][_ENDY_]<<endl;
            myfile<<" INI_Z  "<<iter_bound[_FRONT_][_INIZ_]<<" END_Z "<<iter_bound[_FRONT_][_ENDZ_]<<endl;
            myfile<<"--------"<<endl;
            myfile<<" Send iterators: "<<endl;
            myfile<<" WEST: "<<endl;
            myfile<<" INI_X  "<<iter_toSend[_WEST_][_INIX_]<<" END_X "<<iter_toSend[_WEST_][_ENDX_]<<endl;
            myfile<<" INI_Y  "<<iter_toSend[_WEST_][_INIY_]<<" END_Y "<<iter_toSend[_WEST_][_ENDY_]<<endl;
            myfile<<" INI_Z  "<<iter_toSend[_WEST_][_INIZ_]<<" END_Z "<<iter_toSend[_WEST_][_ENDZ_]<<endl;
            myfile<<endl;
            myfile<<" EAST: "<<endl;
            myfile<<" INI_X  "<<iter_toSend[_EAST_][_INIX_]<<" END_X "<<iter_toSend[_EAST_][_ENDX_]<<endl;
            myfile<<" INI_Y  "<<iter_toSend[_EAST_][_INIY_]<<" END_Y "<<iter_toSend[_EAST_][_ENDY_]<<endl;
            myfile<<" INI_Z  "<<iter_toSend[_EAST_][_INIZ_]<<" END_Z "<<iter_toSend[_EAST_][_ENDZ_]<<endl;
            myfile<<endl;
            myfile<<" SOUTH: "<<endl;
            myfile<<" INI_X  "<<iter_toSend[_SOUTH_][_INIX_]<<" END_X "<<iter_toSend[_SOUTH_][_ENDX_]<<endl;
            myfile<<" INI_Y  "<<iter_toSend[_SOUTH_][_INIY_]<<" END_Y "<<iter_toSend[_SOUTH_][_ENDY_]<<endl;
            myfile<<" INI_Z  "<<iter_toSend[_SOUTH_][_INIZ_]<<" END_Z "<<iter_toSend[_SOUTH_][_ENDZ_]<<endl;
            myfile<<endl;
            myfile<<" NORTH: "<<endl;
            myfile<<" INI_X  "<<iter_toSend[_NORTH_][_INIX_]<<" END_X "<<iter_toSend[_NORTH_][_ENDX_]<<endl;
            myfile<<" INI_Y  "<<iter_toSend[_NORTH_][_INIY_]<<" END_Y "<<iter_toSend[_NORTH_][_ENDY_]<<endl;
            myfile<<" INI_Z  "<<iter_toSend[_NORTH_][_INIZ_]<<" END_Z "<<iter_toSend[_NORTH_][_ENDZ_]<<endl;
            myfile<<endl;
            myfile<<" BACK: "<<endl;
            myfile<<" INI_X  "<<iter_toSend[_BACK_][_INIX_]<<" END_X "<<iter_toSend[_BACK_][_ENDX_]<<endl;
            myfile<<" INI_Y  "<<iter_toSend[_BACK_][_INIY_]<<" END_Y "<<iter_toSend[_BACK_][_ENDY_]<<endl;
            myfile<<" INI_Z  "<<iter_toSend[_BACK_][_INIZ_]<<" END_Z "<<iter_toSend[_BACK_][_ENDZ_]<<endl;
            myfile<<endl;
            myfile<<" FRONT: "<<endl;
            myfile<<" INI_X  "<<iter_toSend[_FRONT_][_INIX_]<<" END_X "<<iter_toSend[_FRONT_][_ENDX_]<<endl;
            myfile<<" INI_Y  "<<iter_toSend[_FRONT_][_INIY_]<<" END_Y "<<iter_toSend[_FRONT_][_ENDY_]<<endl;
            myfile<<" INI_Z  "<<iter_toSend[_FRONT_][_INIZ_]<<" END_Z "<<iter_toSend[_FRONT_][_ENDZ_]<<endl;
            myfile<<"--------"<<endl;
            myfile<<" Recv iterators: "<<endl;
            myfile<<" WEST: "<<endl;
            myfile<<" INI_X  "<<iter_toRecv[_WEST_][_INIX_]<<" END_X "<<iter_toRecv[_WEST_][_ENDX_]<<endl;
            myfile<<" INI_Y  "<<iter_toRecv[_WEST_][_INIY_]<<" END_Y "<<iter_toRecv[_WEST_][_ENDY_]<<endl;
            myfile<<" INI_Z  "<<iter_toRecv[_WEST_][_INIZ_]<<" END_Z "<<iter_toRecv[_WEST_][_ENDZ_]<<endl;
            myfile<<endl;
            myfile<<" EAST: "<<endl;
            myfile<<" INI_X  "<<iter_toRecv[_EAST_][_INIX_]<<" END_X "<<iter_toRecv[_EAST_][_ENDX_]<<endl;
            myfile<<" INI_Y  "<<iter_toRecv[_EAST_][_INIY_]<<" END_Y "<<iter_toRecv[_EAST_][_ENDY_]<<endl;
            myfile<<" INI_Z  "<<iter_toRecv[_EAST_][_INIZ_]<<" END_Z "<<iter_toRecv[_EAST_][_ENDZ_]<<endl;
            myfile<<endl;
            myfile<<" SOUTH: "<<endl;
            myfile<<" INI_X  "<<iter_toRecv[_SOUTH_][_INIX_]<<" END_X "<<iter_toRecv[_SOUTH_][_ENDX_]<<endl;
            myfile<<" INI_Y  "<<iter_toRecv[_SOUTH_][_INIY_]<<" END_Y "<<iter_toRecv[_SOUTH_][_ENDY_]<<endl;
            myfile<<" INI_Z  "<<iter_toRecv[_SOUTH_][_INIZ_]<<" END_Z "<<iter_toRecv[_SOUTH_][_ENDZ_]<<endl;
            myfile<<endl;
            myfile<<" NORTH: "<<endl;
            myfile<<" INI_X  "<<iter_toRecv[_NORTH_][_INIX_]<<" END_X "<<iter_toRecv[_NORTH_][_ENDX_]<<endl;
            myfile<<" INI_Y  "<<iter_toRecv[_NORTH_][_INIY_]<<" END_Y "<<iter_toRecv[_NORTH_][_ENDY_]<<endl;
            myfile<<" INI_Z  "<<iter_toRecv[_NORTH_][_INIZ_]<<" END_Z "<<iter_toRecv[_NORTH_][_ENDZ_]<<endl;
            myfile<<endl;
            myfile<<" BACK: "<<endl;
            myfile<<" INI_X  "<<iter_toRecv[_BACK_][_INIX_]<<" END_X "<<iter_toRecv[_BACK_][_ENDX_]<<endl;
            myfile<<" INI_Y  "<<iter_toRecv[_BACK_][_INIY_]<<" END_Y "<<iter_toRecv[_BACK_][_ENDY_]<<endl;
            myfile<<" INI_Z  "<<iter_toRecv[_BACK_][_INIZ_]<<" END_Z "<<iter_toRecv[_BACK_][_ENDZ_]<<endl;
            myfile<<endl;
            myfile<<" FRONT: "<<endl;
            myfile<<" INI_X  "<<iter_toRecv[_FRONT_][_INIX_]<<" END_X "<<iter_toRecv[_FRONT_][_ENDX_]<<endl;
            myfile<<" INI_Y  "<<iter_toRecv[_FRONT_][_INIY_]<<" END_Y "<<iter_toRecv[_FRONT_][_ENDY_]<<endl;
            myfile<<" INI_Z  "<<iter_toRecv[_FRONT_][_INIZ_]<<" END_Z "<<iter_toRecv[_FRONT_][_ENDZ_]<<endl;


            myfile.close();
        }
        else{
            cout<<"Unable to create file: "<<filename<<endl;
        }

       
     }
    MPI_Barrier(RHEA_3DCOMM);
}


/* This part assumes halo size 1*/
void ParallelTopology::create_common_iters()
{
    // By default non halos are created
    for(int bd=_INNER_; bd<= _ALL_; bd++)
    {
        iter_common[bd][_INIX_] =   0;
        iter_common[bd][_ENDX_] =  -1;
        iter_common[bd][_INIY_] =   0;
        iter_common[bd][_ENDY_] =  -1;
        iter_common[bd][_INIZ_] =   0;
        iter_common[bd][_ENDZ_] =  -1;
    }


    //INNER
    iter_common[_INNER_][_INIX_] =   1;
    iter_common[_INNER_][_ENDX_] = lNx-2;
    iter_common[_INNER_][_INIY_] =   1;
    iter_common[_INNER_][_ENDY_] = lNy-2;
    iter_common[_INNER_][_INIZ_] =   1;
    iter_common[_INNER_][_ENDZ_] = lNz-2;

    //ALL
    iter_common[_ALL_][_INIX_] =   0;
    iter_common[_ALL_][_ENDX_] = lNx-1;
    iter_common[_ALL_][_INIY_] =   0;
    iter_common[_ALL_][_ENDY_] = lNy-1;
    iter_common[_ALL_][_INIZ_] =   0;
    iter_common[_ALL_][_ENDZ_] = lNz-1;

}


/* This part assumes halo size 1*/
void ParallelTopology::create_basic_bound_iters()
{
    // By default non halos are created
    for(int bd=_WEST_; bd<= _EAST_N_F_; bd++)
    {
        iter_bound[bd][_INIX_] =   0;
        iter_bound[bd][_ENDX_] =  -1;
        iter_bound[bd][_INIY_] =   0;
        iter_bound[bd][_ENDY_] =  -1;
        iter_bound[bd][_INIZ_] =   0;
        iter_bound[bd][_ENDZ_] =  -1;
    }


    //WEST
    if(neighb[_WEST_]== _NO_NEIGHBOUR_) {
        iter_bound[_WEST_][_INIX_] =   0;
        iter_bound[_WEST_][_ENDX_] =   0;
        iter_bound[_WEST_][_INIY_] =   1;
        iter_bound[_WEST_][_ENDY_] = lNy-2;
        iter_bound[_WEST_][_INIZ_] =   1;
        iter_bound[_WEST_][_ENDZ_] = lNz-2;
    }

    //EAST
    if(neighb[_EAST_]== _NO_NEIGHBOUR_) {
        iter_bound[_EAST_][_INIX_] = lNx-1;
        iter_bound[_EAST_][_ENDX_] = lNx-1;
        iter_bound[_EAST_][_INIY_] =   1;
        iter_bound[_EAST_][_ENDY_] = lNy-2;
        iter_bound[_EAST_][_INIZ_] =   1;
        iter_bound[_EAST_][_ENDZ_] = lNz-2;
    }

    //SOUTH
    if(neighb[_SOUTH_]== _NO_NEIGHBOUR_) {
        iter_bound[_SOUTH_][_INIX_] =   1;
        iter_bound[_SOUTH_][_ENDX_] = lNx-2;
        iter_bound[_SOUTH_][_INIY_] =   0;
        iter_bound[_SOUTH_][_ENDY_] =   0;
        iter_bound[_SOUTH_][_INIZ_] =   1;
        iter_bound[_SOUTH_][_ENDZ_] = lNz-2;
    }

    //NORTH
    if(neighb[_NORTH_]== _NO_NEIGHBOUR_) {
        iter_bound[_NORTH_][_INIX_] =   1;
        iter_bound[_NORTH_][_ENDX_] =   lNx-2;
        iter_bound[_NORTH_][_INIY_] =   lNy-1;
        iter_bound[_NORTH_][_ENDY_] =   lNy-1;
        iter_bound[_NORTH_][_INIZ_] =   1;
        iter_bound[_NORTH_][_ENDZ_] = lNz-2;
    }

    //BACK
    if(neighb[_BACK_]== _NO_NEIGHBOUR_) {
        iter_bound[_BACK_][_INIX_] =   1;
        iter_bound[_BACK_][_ENDX_] =   lNx-2;
        iter_bound[_BACK_][_INIY_] =   1;
        iter_bound[_BACK_][_ENDY_] =   lNy-2;
        iter_bound[_BACK_][_INIZ_] =   0;
        iter_bound[_BACK_][_ENDZ_] =   0;
    }

    //FRONT
    if(neighb[_FRONT_]== _NO_NEIGHBOUR_) {
        iter_bound[_FRONT_][_INIX_] =   1;
        iter_bound[_FRONT_][_ENDX_] =   lNx-2;
        iter_bound[_FRONT_][_INIY_] =   1;
        iter_bound[_FRONT_][_ENDY_] =   lNy-2;
        iter_bound[_FRONT_][_INIZ_] =   lNz-1;
        iter_bound[_FRONT_][_ENDZ_] =   lNz-1;
    }
}
    
void ParallelTopology::create_complex_bound_iters(){

    // 2nd Level (lines)
    //WEST SOUTH
    if(neighb[_WEST_S_]== _NO_NEIGHBOUR_) {
        iter_bound[_WEST_S_][_INIX_] =   0;
        iter_bound[_WEST_S_][_ENDX_] =   0;
        iter_bound[_WEST_S_][_INIY_] =   0;
        iter_bound[_WEST_S_][_ENDY_] =   0;
        iter_bound[_WEST_S_][_INIZ_] =   1;
        iter_bound[_WEST_S_][_ENDZ_] = lNz-2;
    }

    //WEST NORTH
    if(neighb[_WEST_N_]== _NO_NEIGHBOUR_) {
        iter_bound[_WEST_N_][_INIX_] =   0;
        iter_bound[_WEST_N_][_ENDX_] =   0;
        iter_bound[_WEST_N_][_INIY_] = lNy-1;
        iter_bound[_WEST_N_][_ENDY_] = lNy-1;
        iter_bound[_WEST_N_][_INIZ_] =   1;
        iter_bound[_WEST_N_][_ENDZ_] = lNz-2;
    }

    //WEST BACK
    if(neighb[_WEST_B_]== _NO_NEIGHBOUR_) {
        iter_bound[_WEST_B_][_INIX_] =   0;
        iter_bound[_WEST_B_][_ENDX_] =   0;
        iter_bound[_WEST_B_][_INIY_] =   1;
        iter_bound[_WEST_B_][_ENDY_] = lNy-2;
        iter_bound[_WEST_B_][_INIZ_] =   0;
        iter_bound[_WEST_B_][_ENDZ_] =   0;
    }

    //WEST FRONT
    if(neighb[_WEST_F_]== _NO_NEIGHBOUR_) {
        iter_bound[_WEST_F_][_INIX_] =   0;
        iter_bound[_WEST_F_][_ENDX_] =   0;
        iter_bound[_WEST_F_][_INIY_] =   1;
        iter_bound[_WEST_F_][_ENDY_] = lNy-2;
        iter_bound[_WEST_F_][_INIZ_] = lNz-1;
        iter_bound[_WEST_F_][_ENDZ_] = lNz-1;
    }


    //EAST SOUTH
    if(neighb[_EAST_S_]== _NO_NEIGHBOUR_) {
        iter_bound[_EAST_S_][_INIX_] = lNx-1;
        iter_bound[_EAST_S_][_ENDX_] = lNx-1;
        iter_bound[_EAST_S_][_INIY_] =   0;
        iter_bound[_EAST_S_][_ENDY_] =   0;
        iter_bound[_EAST_S_][_INIZ_] =   1;
        iter_bound[_EAST_S_][_ENDZ_] = lNz-2;
    }

    //EAST NORTH
    if(neighb[_EAST_N_]== _NO_NEIGHBOUR_) {
        iter_bound[_EAST_N_][_INIX_] = lNx-1;
        iter_bound[_EAST_N_][_ENDX_] = lNx-1;
        iter_bound[_EAST_N_][_INIY_] = lNy-1;
        iter_bound[_EAST_N_][_ENDY_] = lNy-1;
        iter_bound[_EAST_N_][_INIZ_] =   1;
        iter_bound[_EAST_N_][_ENDZ_] = lNz-2;
    }

    //EAST BACK
    if(neighb[_EAST_B_]== _NO_NEIGHBOUR_) {
        iter_bound[_EAST_B_][_INIX_] = lNx-1;
        iter_bound[_EAST_B_][_ENDX_] = lNx-1;
        iter_bound[_EAST_B_][_INIY_] =   1;
        iter_bound[_EAST_B_][_ENDY_] = lNy-2;
        iter_bound[_EAST_B_][_INIZ_] =   0;
        iter_bound[_EAST_B_][_ENDZ_] =   0;
    }

    //EAST FRONT
    if(neighb[_EAST_F_]== _NO_NEIGHBOUR_) {
        iter_bound[_EAST_F_][_INIX_] = lNx-1;
        iter_bound[_EAST_F_][_ENDX_] = lNx-1;
        iter_bound[_EAST_F_][_INIY_] =   1;
        iter_bound[_EAST_F_][_ENDY_] = lNy-2;
        iter_bound[_EAST_F_][_INIZ_] = lNz-1;
        iter_bound[_EAST_F_][_ENDZ_] = lNz-1;
    }

    //SOUTH BACK
    if(neighb[_SOUTH_B_]== _NO_NEIGHBOUR_) {
        iter_bound[_SOUTH_B_][_INIX_] =   1;
        iter_bound[_SOUTH_B_][_ENDX_] = lNx-2;
        iter_bound[_SOUTH_B_][_INIY_] =   0;
        iter_bound[_SOUTH_B_][_ENDY_] =   0;
        iter_bound[_SOUTH_B_][_INIZ_] =   0;
        iter_bound[_SOUTH_B_][_ENDZ_] =   0;
    }

    //SOUTH FRONT
    if(neighb[_SOUTH_F_]== _NO_NEIGHBOUR_) {
        iter_bound[_SOUTH_F_][_INIX_] =   1;
        iter_bound[_SOUTH_F_][_ENDX_] = lNx-2;
        iter_bound[_SOUTH_F_][_INIY_] =   0;
        iter_bound[_SOUTH_F_][_ENDY_] =   0;
        iter_bound[_SOUTH_F_][_INIZ_] = lNz-1;
        iter_bound[_SOUTH_F_][_ENDZ_] = lNz-1;
    }


    //NORTH BACK
    if(neighb[_NORTH_B_]== _NO_NEIGHBOUR_) {
        iter_bound[_NORTH_B_][_INIX_] =   1;
        iter_bound[_NORTH_B_][_ENDX_] = lNx-2;
        iter_bound[_NORTH_B_][_INIY_] = lNy-1;
        iter_bound[_NORTH_B_][_ENDY_] = lNy-1;
        iter_bound[_NORTH_B_][_INIZ_] =   0;
        iter_bound[_NORTH_B_][_ENDZ_] =   0;
    }

    //NORTH FRONT
    if(neighb[_NORTH_F_]== _NO_NEIGHBOUR_) {
        iter_bound[_NORTH_F_][_INIX_] =   1;
        iter_bound[_NORTH_F_][_ENDX_] = lNx-2;
        iter_bound[_NORTH_F_][_INIY_] = lNy-1;
        iter_bound[_NORTH_F_][_ENDY_] = lNy-1;
        iter_bound[_NORTH_F_][_INIZ_] = lNz-1;
        iter_bound[_NORTH_F_][_ENDZ_] = lNz-1;
    }

// 3er Level (corners)

    //WEST SOUTH BACK
    if(neighb[_WEST_S_B_]== _NO_NEIGHBOUR_) {
        iter_bound[_WEST_S_B_][_INIX_] =   0;
        iter_bound[_WEST_S_B_][_ENDX_] =   0;
        iter_bound[_WEST_S_B_][_INIY_] =   0;
        iter_bound[_WEST_S_B_][_ENDY_] =   0;
        iter_bound[_WEST_S_B_][_INIZ_] =   0;
        iter_bound[_WEST_S_B_][_ENDZ_] =   0;
    }

    //WEST NORTH BACK
    if(neighb[_WEST_N_B_]== _NO_NEIGHBOUR_) {
        iter_bound[_WEST_N_B_][_INIX_] =   0;
        iter_bound[_WEST_N_B_][_ENDX_] =   0;
        iter_bound[_WEST_N_B_][_INIY_] = lNy-1;
        iter_bound[_WEST_N_B_][_ENDY_] = lNy-1;
        iter_bound[_WEST_N_B_][_INIZ_] =   0;
        iter_bound[_WEST_N_B_][_ENDZ_] =   0;
    }

    //WEST SOUTH FRONT
    if(neighb[_WEST_S_F_]== _NO_NEIGHBOUR_) {
        iter_bound[_WEST_S_F_][_INIX_] =   0;
        iter_bound[_WEST_S_F_][_ENDX_] =   0;
        iter_bound[_WEST_S_F_][_INIY_] =   0;
        iter_bound[_WEST_S_F_][_ENDY_] =   0;
        iter_bound[_WEST_S_F_][_INIZ_] = lNz-1;
        iter_bound[_WEST_S_F_][_ENDZ_] = lNz-1;
    }

    //WEST NORTH FRONT
    if(neighb[_WEST_N_F_]== _NO_NEIGHBOUR_) {
        iter_bound[_WEST_N_F_][_INIX_] =   0;
        iter_bound[_WEST_N_F_][_ENDX_] =   0;
        iter_bound[_WEST_N_F_][_INIY_] = lNy-1;
        iter_bound[_WEST_N_F_][_ENDY_] = lNy-1;
        iter_bound[_WEST_N_F_][_INIZ_] = lNz-1;
        iter_bound[_WEST_N_F_][_ENDZ_] = lNz-1;
    }

 
    //EAST SOUTH BACK
    if(neighb[_EAST_S_B_]== _NO_NEIGHBOUR_) {
        iter_bound[_EAST_S_B_][_INIX_] = lNx-1;
        iter_bound[_EAST_S_B_][_ENDX_] = lNx-1;
        iter_bound[_EAST_S_B_][_INIY_] =   0;
        iter_bound[_EAST_S_B_][_ENDY_] =   0;
        iter_bound[_EAST_S_B_][_INIZ_] =   0;
        iter_bound[_EAST_S_B_][_ENDZ_] =   0;
    }

    //EAST NORTH BACK
    if(neighb[_EAST_N_B_]== _NO_NEIGHBOUR_) {
        iter_bound[_EAST_N_B_][_INIX_] = lNx-1;
        iter_bound[_EAST_N_B_][_ENDX_] = lNx-1;
        iter_bound[_EAST_N_B_][_INIY_] = lNy-1;
        iter_bound[_EAST_N_B_][_ENDY_] = lNy-1;
        iter_bound[_EAST_N_B_][_INIZ_] =   0;
        iter_bound[_EAST_N_B_][_ENDZ_] =   0;
    }


    //EAST SOUTH FRONT
    if(neighb[_EAST_S_F_]== _NO_NEIGHBOUR_) {
        iter_bound[_EAST_S_F_][_INIX_] = lNx-1;
        iter_bound[_EAST_S_F_][_ENDX_] = lNx-1;
        iter_bound[_EAST_S_F_][_INIY_] =   0;
        iter_bound[_EAST_S_F_][_ENDY_] =   0;
        iter_bound[_EAST_S_F_][_INIZ_] = lNz-1;
        iter_bound[_EAST_S_F_][_ENDZ_] = lNz-1;
    }

    //EAST NORTH FRONT
    if(neighb[_EAST_N_F_]== _NO_NEIGHBOUR_) {
        iter_bound[_EAST_N_F_][_INIX_] = lNx-1;
        iter_bound[_EAST_N_F_][_ENDX_] = lNx-1;
        iter_bound[_EAST_N_F_][_INIY_] = lNy-1;
        iter_bound[_EAST_N_F_][_ENDY_] = lNy-1;
        iter_bound[_EAST_N_F_][_INIZ_] = lNz-1;
        iter_bound[_EAST_N_F_][_ENDZ_] = lNz-1;
    }


}


/* This part assumes halo size 1*/
/* This rarely will be used */
void ParallelTopology::create_halo_iters()
{
    // By default non halos are created
    for(int bd=_WEST_; bd<= _EAST_N_F_; bd++)
    {
        iter_halo[bd][_INIX_] =   0;
        iter_halo[bd][_ENDX_] =  -1;
        iter_halo[bd][_INIY_] =   0;
        iter_halo[bd][_ENDY_] =  -1;
        iter_halo[bd][_INIZ_] =   0;
        iter_halo[bd][_ENDZ_] =  -1;
    }


    //WEST
    if(neighb[_WEST_] != _NO_NEIGHBOUR_ ) {
        iter_halo[_WEST_][_INIX_] =   0;
        iter_halo[_WEST_][_ENDX_] =   0;
        iter_halo[_WEST_][_INIY_] =   1;
        iter_halo[_WEST_][_ENDY_] = lNy-2;
        iter_halo[_WEST_][_INIZ_] =   1;
        iter_halo[_WEST_][_ENDZ_] = lNz-2;
    }

    //EAST
    if(neighb[_EAST_] != _NO_NEIGHBOUR_ ) {
        iter_halo[_EAST_][_INIX_] = lNx-1;
        iter_halo[_EAST_][_ENDX_] = lNx-1;
        iter_halo[_EAST_][_INIY_] =   1;
        iter_halo[_EAST_][_ENDY_] = lNy-2;
        iter_halo[_EAST_][_INIZ_] =   1;
        iter_halo[_EAST_][_ENDZ_] = lNz-2;
    }

    //SOUTH
    if(neighb[_SOUTH_] != _NO_NEIGHBOUR_ ) {
        iter_halo[_SOUTH_][_INIX_] =   1;
        iter_halo[_SOUTH_][_ENDX_] = lNx-2;
        iter_halo[_SOUTH_][_INIY_] =   0;
        iter_halo[_SOUTH_][_ENDY_] =   0;
        iter_halo[_SOUTH_][_INIZ_] =   1;
        iter_halo[_SOUTH_][_ENDZ_] = lNz-2;
    }

    //NORTH
    if(neighb[_NORTH_] != _NO_NEIGHBOUR_ ) {
        iter_halo[_NORTH_][_INIX_] =   1;
        iter_halo[_NORTH_][_ENDX_] =   lNx-2;
        iter_halo[_NORTH_][_INIY_] =   lNy-1;
        iter_halo[_NORTH_][_ENDY_] =   lNy-1;
        iter_halo[_NORTH_][_INIZ_] =   1;
        iter_halo[_NORTH_][_ENDZ_] = lNz-2;
    }

    //BACK
    if(neighb[_BACK_] != _NO_NEIGHBOUR_ ) {
        iter_halo[_BACK_][_INIX_] =   1;
        iter_halo[_BACK_][_ENDX_] =   lNx-2;
        iter_halo[_BACK_][_INIY_] =   1;
        iter_halo[_BACK_][_ENDY_] =   lNy-2;
        iter_halo[_BACK_][_INIZ_] =   0;
        iter_halo[_BACK_][_ENDZ_] =   0;
    }

    //FRONT
    if(neighb[_FRONT_] != _NO_NEIGHBOUR_ ) {
        iter_halo[_FRONT_][_INIX_] =   1;
        iter_halo[_FRONT_][_ENDX_] =   lNx-2;
        iter_halo[_FRONT_][_INIY_] =   1;
        iter_halo[_FRONT_][_ENDY_] =   lNy-2;
        iter_halo[_FRONT_][_INIZ_] =   lNz-1;
        iter_halo[_FRONT_][_ENDZ_] =   lNz-1;
    }
}

void ParallelTopology::create_complex_halo_iters()
{

    //2nd Level Halos (lines)
    //WEST SOUTH
    if(neighb[_WEST_S_] != _NO_NEIGHBOUR_) {
        iter_halo[_WEST_S_][_INIX_] =   0;
        iter_halo[_WEST_S_][_ENDX_] =   0;
        iter_halo[_WEST_S_][_INIY_] =   0;
        iter_halo[_WEST_S_][_ENDY_] =   0;
        iter_halo[_WEST_S_][_INIZ_] =   1;
        iter_halo[_WEST_S_][_ENDZ_] = lNz-2;
    }

    //WEST NORTH
    if(neighb[_WEST_N_] != _NO_NEIGHBOUR_) {
        iter_halo[_WEST_N_][_INIX_] =   0;
        iter_halo[_WEST_N_][_ENDX_] =   0;
        iter_halo[_WEST_N_][_INIY_] = lNy-1;
        iter_halo[_WEST_N_][_ENDY_] = lNy-1;
        iter_halo[_WEST_N_][_INIZ_] =   1;
        iter_halo[_WEST_N_][_ENDZ_] = lNz-2;
    }

    //WEST BACK
    if(neighb[_WEST_B_] != _NO_NEIGHBOUR_) {
        iter_halo[_WEST_B_][_INIX_] =   0;
        iter_halo[_WEST_B_][_ENDX_] =   0;
        iter_halo[_WEST_B_][_INIY_] =   1;
        iter_halo[_WEST_B_][_ENDY_] = lNy-2;
        iter_halo[_WEST_B_][_INIZ_] =   0;
        iter_halo[_WEST_B_][_ENDZ_] =   0;
    }

    //WEST FRONT
    if(neighb[_WEST_F_] != _NO_NEIGHBOUR_) {
        iter_halo[_WEST_F_][_INIX_] =   0;
        iter_halo[_WEST_F_][_ENDX_] =   0;
        iter_halo[_WEST_F_][_INIY_] =   1;
        iter_halo[_WEST_F_][_ENDY_] = lNy-2;
        iter_halo[_WEST_F_][_INIZ_] = lNz-1;
        iter_halo[_WEST_F_][_ENDZ_] = lNz-1;
    }


    //EAST SOUTH
    if(neighb[_EAST_S_] != _NO_NEIGHBOUR_) {
        iter_halo[_EAST_S_][_INIX_] = lNx-1;
        iter_halo[_EAST_S_][_ENDX_] = lNx-1;
        iter_halo[_EAST_S_][_INIY_] =   0;
        iter_halo[_EAST_S_][_ENDY_] =   0;
        iter_halo[_EAST_S_][_INIZ_] =   1;
        iter_halo[_EAST_S_][_ENDZ_] = lNz-2;
    }

    //EAST NORTH
    if(neighb[_EAST_N_] != _NO_NEIGHBOUR_) {
        iter_halo[_EAST_N_][_INIX_] = lNx-1;
        iter_halo[_EAST_N_][_ENDX_] = lNx-1;
        iter_halo[_EAST_N_][_INIY_] = lNy-1;
        iter_halo[_EAST_N_][_ENDY_] = lNy-1;
        iter_halo[_EAST_N_][_INIZ_] =   1;
        iter_halo[_EAST_N_][_ENDZ_] = lNz-2;
    }

    //EAST BACK
    if(neighb[_EAST_B_] != _NO_NEIGHBOUR_) {
        iter_halo[_EAST_B_][_INIX_] = lNx-1;
        iter_halo[_EAST_B_][_ENDX_] = lNx-1;
        iter_halo[_EAST_B_][_INIY_] =   1;
        iter_halo[_EAST_B_][_ENDY_] = lNy-2;
        iter_halo[_EAST_B_][_INIZ_] =   0;
        iter_halo[_EAST_B_][_ENDZ_] =   0;
    }

    //EAST FRONT
    if(neighb[_EAST_F_] != _NO_NEIGHBOUR_) {
        iter_halo[_EAST_F_][_INIX_] = lNx-1;
        iter_halo[_EAST_F_][_ENDX_] = lNx-1;
        iter_halo[_EAST_F_][_INIY_] =   1;
        iter_halo[_EAST_F_][_ENDY_] = lNy-2;
        iter_halo[_EAST_F_][_INIZ_] = lNz-1;
        iter_halo[_EAST_F_][_ENDZ_] = lNz-1;
    }

    //SOUTH BACK
    if(neighb[_SOUTH_B_] != _NO_NEIGHBOUR_) {
        iter_halo[_SOUTH_B_][_INIX_] =   1;
        iter_halo[_SOUTH_B_][_ENDX_] = lNx-2;
        iter_halo[_SOUTH_B_][_INIY_] =   0;
        iter_halo[_SOUTH_B_][_ENDY_] =   0;
        iter_halo[_SOUTH_B_][_INIZ_] =   0;
        iter_halo[_SOUTH_B_][_ENDZ_] =   0;
    }

    //SOUTH FRONT
    if(neighb[_SOUTH_F_] != _NO_NEIGHBOUR_) {
        iter_halo[_SOUTH_F_][_INIX_] =   1;
        iter_halo[_SOUTH_F_][_ENDX_] = lNx-2;
        iter_halo[_SOUTH_F_][_INIY_] =   0;
        iter_halo[_SOUTH_F_][_ENDY_] =   0;
        iter_halo[_SOUTH_F_][_INIZ_] = lNz-1;
        iter_halo[_SOUTH_F_][_ENDZ_] = lNz-1;
    }


    //NORTH BACK
    if(neighb[_NORTH_B_] != _NO_NEIGHBOUR_) {
        iter_halo[_NORTH_B_][_INIX_] =   1;
        iter_halo[_NORTH_B_][_ENDX_] = lNx-2;
        iter_halo[_NORTH_B_][_INIY_] = lNy-1;
        iter_halo[_NORTH_B_][_ENDY_] = lNy-1;
        iter_halo[_NORTH_B_][_INIZ_] =   0;
        iter_halo[_NORTH_B_][_ENDZ_] =   0;
    }

    //NORTH FRONT
    if(neighb[_NORTH_F_] != _NO_NEIGHBOUR_) {
        iter_halo[_NORTH_F_][_INIX_] =   1;
        iter_halo[_NORTH_F_][_ENDX_] = lNx-2;
        iter_halo[_NORTH_F_][_INIY_] = lNy-1;
        iter_halo[_NORTH_F_][_ENDY_] = lNy-1;
        iter_halo[_NORTH_F_][_INIZ_] = lNz-1;
        iter_halo[_NORTH_F_][_ENDZ_] = lNz-1;
    }

// 3er Level (corners)

    //WEST SOUTH BACK
    if(neighb[_WEST_S_B_] != _NO_NEIGHBOUR_) {
        iter_halo[_WEST_S_B_][_INIX_] =   0;
        iter_halo[_WEST_S_B_][_ENDX_] =   0;
        iter_halo[_WEST_S_B_][_INIY_] =   0;
        iter_halo[_WEST_S_B_][_ENDY_] =   0;
        iter_halo[_WEST_S_B_][_INIZ_] =   0;
        iter_halo[_WEST_S_B_][_ENDZ_] =   0;
    }

    //WEST NORTH BACK
    if(neighb[_WEST_N_B_] != _NO_NEIGHBOUR_) {
        iter_halo[_WEST_N_B_][_INIX_] =   0;
        iter_halo[_WEST_N_B_][_ENDX_] =   0;
        iter_halo[_WEST_N_B_][_INIY_] = lNy-1;
        iter_halo[_WEST_N_B_][_ENDY_] = lNy-1;
        iter_halo[_WEST_N_B_][_INIZ_] =   0;
        iter_halo[_WEST_N_B_][_ENDZ_] =   0;
    }

    //WEST SOUTH FRONT
    if(neighb[_WEST_S_F_] != _NO_NEIGHBOUR_) {
        iter_halo[_WEST_S_F_][_INIX_] =   0;
        iter_halo[_WEST_S_F_][_ENDX_] =   0;
        iter_halo[_WEST_S_F_][_INIY_] =   0;
        iter_halo[_WEST_S_F_][_ENDY_] =   0;
        iter_halo[_WEST_S_F_][_INIZ_] = lNz-1;
        iter_halo[_WEST_S_F_][_ENDZ_] = lNz-1;
    }

    //WEST NORTH FRONT
    if(neighb[_WEST_N_F_] != _NO_NEIGHBOUR_) {
        iter_halo[_WEST_N_F_][_INIX_] =   0;
        iter_halo[_WEST_N_F_][_ENDX_] =   0;
        iter_halo[_WEST_N_F_][_INIY_] = lNy-1;
        iter_halo[_WEST_N_F_][_ENDY_] = lNy-1;
        iter_halo[_WEST_N_F_][_INIZ_] = lNz-1;
        iter_halo[_WEST_N_F_][_ENDZ_] = lNz-1;
    }

 
    //EAST SOUTH BACK
    if(neighb[_EAST_S_B_] != _NO_NEIGHBOUR_) {
        iter_halo[_EAST_S_B_][_INIX_] = lNx-1;
        iter_halo[_EAST_S_B_][_ENDX_] = lNx-1;
        iter_halo[_EAST_S_B_][_INIY_] =   0;
        iter_halo[_EAST_S_B_][_ENDY_] =   0;
        iter_halo[_EAST_S_B_][_INIZ_] =   0;
        iter_halo[_EAST_S_B_][_ENDZ_] =   0;
    }

    //EAST NORTH BACK
    if(neighb[_EAST_N_B_] != _NO_NEIGHBOUR_) {
        iter_halo[_EAST_N_B_][_INIX_] = lNx-1;
        iter_halo[_EAST_N_B_][_ENDX_] = lNx-1;
        iter_halo[_EAST_N_B_][_INIY_] = lNy-1;
        iter_halo[_EAST_N_B_][_ENDY_] = lNy-1;
        iter_halo[_EAST_N_B_][_INIZ_] =   0;
        iter_halo[_EAST_N_B_][_ENDZ_] =   0;
    }


    //EAST SOUTH FRONT
    if(neighb[_EAST_S_F_] != _NO_NEIGHBOUR_) {
        iter_halo[_EAST_S_F_][_INIX_] = lNx-1;
        iter_halo[_EAST_S_F_][_ENDX_] = lNx-1;
        iter_halo[_EAST_S_F_][_INIY_] =   0;
        iter_halo[_EAST_S_F_][_ENDY_] =   0;
        iter_halo[_EAST_S_F_][_INIZ_] = lNz-1;
        iter_halo[_EAST_S_F_][_ENDZ_] = lNz-1;
    }

    //EAST NORTH FRONT
    if(neighb[_EAST_N_F_] != _NO_NEIGHBOUR_) {
        iter_halo[_EAST_N_F_][_INIX_] = lNx-1;
        iter_halo[_EAST_N_F_][_ENDX_] = lNx-1;
        iter_halo[_EAST_N_F_][_INIY_] = lNy-1;
        iter_halo[_EAST_N_F_][_ENDY_] = lNy-1;
        iter_halo[_EAST_N_F_][_INIZ_] = lNz-1;
        iter_halo[_EAST_N_F_][_ENDZ_] = lNz-1;
    }

}


/* This part assumes halo size 1*/
void ParallelTopology::create_toRecv_iters()
{

    for(int bd=_WEST_; bd<= _EAST_N_F_; bd++)
    {
        iter_toRecv[bd][_INIX_] =   0;
        iter_toRecv[bd][_ENDX_] =  -1;
        iter_toRecv[bd][_INIY_] =   0;
        iter_toRecv[bd][_ENDY_] =  -1;
        iter_toRecv[bd][_INIZ_] =   0;
        iter_toRecv[bd][_ENDZ_] =  -1;
    }

    //WEST
    if(neighb[_WEST_] != _NO_NEIGHBOUR_ ) {
        iter_toRecv[_WEST_][_INIX_] =   0;
        iter_toRecv[_WEST_][_ENDX_] =   0;
        iter_toRecv[_WEST_][_INIY_] =   1;
        iter_toRecv[_WEST_][_ENDY_] = lNy-2;
        iter_toRecv[_WEST_][_INIZ_] =   1;
        iter_toRecv[_WEST_][_ENDZ_] = lNz-2;
    }

    //EAST
    if(neighb[_EAST_] != _NO_NEIGHBOUR_ ) {
        iter_toRecv[_EAST_][_INIX_] = lNx-1;
        iter_toRecv[_EAST_][_ENDX_] = lNx-1;
        iter_toRecv[_EAST_][_INIY_] =   1;
        iter_toRecv[_EAST_][_ENDY_] = lNy-2;
        iter_toRecv[_EAST_][_INIZ_] =   1;
        iter_toRecv[_EAST_][_ENDZ_] = lNz-2;
    }

    //SOUTH
    if(neighb[_SOUTH_] != _NO_NEIGHBOUR_ ) {
        iter_toRecv[_SOUTH_][_INIX_] =   1;
        iter_toRecv[_SOUTH_][_ENDX_] = lNx-2;
        iter_toRecv[_SOUTH_][_INIY_] =   0;
        iter_toRecv[_SOUTH_][_ENDY_] =   0;
        iter_toRecv[_SOUTH_][_INIZ_] =   1;
        iter_toRecv[_SOUTH_][_ENDZ_] = lNz-2;
    }

    //NORTH
    if(neighb[_NORTH_] != _NO_NEIGHBOUR_ ) {
        iter_toRecv[_NORTH_][_INIX_] =   1;
        iter_toRecv[_NORTH_][_ENDX_] =   lNx-2;
        iter_toRecv[_NORTH_][_INIY_] =   lNy-1;
        iter_toRecv[_NORTH_][_ENDY_] =   lNy-1;
        iter_toRecv[_NORTH_][_INIZ_] =   1;
        iter_toRecv[_NORTH_][_ENDZ_] = lNz-2;
    }

    //BACK
    if(neighb[_BACK_] != _NO_NEIGHBOUR_ ) {
        iter_toRecv[_BACK_][_INIX_] =   1;
        iter_toRecv[_BACK_][_ENDX_] =   lNx-2;
        iter_toRecv[_BACK_][_INIY_] =   1;
        iter_toRecv[_BACK_][_ENDY_] =   lNy-2;
        iter_toRecv[_BACK_][_INIZ_] =   0;
        iter_toRecv[_BACK_][_ENDZ_] =   0;
    }

    //FRONT
    if(neighb[_FRONT_] != _NO_NEIGHBOUR_ ) {
        iter_toRecv[_FRONT_][_INIX_] =   1;
        iter_toRecv[_FRONT_][_ENDX_] =   lNx-2;
        iter_toRecv[_FRONT_][_INIY_] =   1;
        iter_toRecv[_FRONT_][_ENDY_] =   lNy-2;
        iter_toRecv[_FRONT_][_INIZ_] =   lNz-1;
        iter_toRecv[_FRONT_][_ENDZ_] =   lNz-1;
    }

}


void ParallelTopology::create_complex_toRecv_iters()
{
    //2nd Level Halos (lines)
    //WEST SOUTH
    if(neighb[_WEST_S_] != _NO_NEIGHBOUR_) {
        iter_toRecv[_WEST_S_][_INIX_] =   0;
        iter_toRecv[_WEST_S_][_ENDX_] =   0;
        iter_toRecv[_WEST_S_][_INIY_] =   0;
        iter_toRecv[_WEST_S_][_ENDY_] =   0;
        iter_toRecv[_WEST_S_][_INIZ_] =   1;
        iter_toRecv[_WEST_S_][_ENDZ_] = lNz-2;
    }

    //WEST NORTH
    if(neighb[_WEST_N_] != _NO_NEIGHBOUR_) {
        iter_toRecv[_WEST_N_][_INIX_] =   0;
        iter_toRecv[_WEST_N_][_ENDX_] =   0;
        iter_toRecv[_WEST_N_][_INIY_] = lNy-1;
        iter_toRecv[_WEST_N_][_ENDY_] = lNy-1;
        iter_toRecv[_WEST_N_][_INIZ_] =   1;
        iter_toRecv[_WEST_N_][_ENDZ_] = lNz-2;
    }

    //WEST BACK
    if(neighb[_WEST_B_] != _NO_NEIGHBOUR_) {
        iter_toRecv[_WEST_B_][_INIX_] =   0;
        iter_toRecv[_WEST_B_][_ENDX_] =   0;
        iter_toRecv[_WEST_B_][_INIY_] =   1;
        iter_toRecv[_WEST_B_][_ENDY_] = lNy-2;
        iter_toRecv[_WEST_B_][_INIZ_] =   0;
        iter_toRecv[_WEST_B_][_ENDZ_] =   0;
    }

    //WEST FRONT
    if(neighb[_WEST_F_] != _NO_NEIGHBOUR_) {
        iter_toRecv[_WEST_F_][_INIX_] =   0;
        iter_toRecv[_WEST_F_][_ENDX_] =   0;
        iter_toRecv[_WEST_F_][_INIY_] =   1;
        iter_toRecv[_WEST_F_][_ENDY_] = lNy-2;
        iter_toRecv[_WEST_F_][_INIZ_] = lNz-1;
        iter_toRecv[_WEST_F_][_ENDZ_] = lNz-1;
    }


    //EAST SOUTH
    if(neighb[_EAST_S_] != _NO_NEIGHBOUR_) {
        iter_toRecv[_EAST_S_][_INIX_] = lNx-1;
        iter_toRecv[_EAST_S_][_ENDX_] = lNx-1;
        iter_toRecv[_EAST_S_][_INIY_] =   0;
        iter_toRecv[_EAST_S_][_ENDY_] =   0;
        iter_toRecv[_EAST_S_][_INIZ_] =   1;
        iter_toRecv[_EAST_S_][_ENDZ_] = lNz-2;
    }

    //EAST NORTH
    if(neighb[_EAST_N_] != _NO_NEIGHBOUR_) {
        iter_toRecv[_EAST_N_][_INIX_] = lNx-1;
        iter_toRecv[_EAST_N_][_ENDX_] = lNx-1;
        iter_toRecv[_EAST_N_][_INIY_] = lNy-1;
        iter_toRecv[_EAST_N_][_ENDY_] = lNy-1;
        iter_toRecv[_EAST_N_][_INIZ_] =   1;
        iter_toRecv[_EAST_N_][_ENDZ_] = lNz-2;
    }

    //EAST BACK
    if(neighb[_EAST_B_] != _NO_NEIGHBOUR_) {
        iter_toRecv[_EAST_B_][_INIX_] = lNx-1;
        iter_toRecv[_EAST_B_][_ENDX_] = lNx-1;
        iter_toRecv[_EAST_B_][_INIY_] =   1;
        iter_toRecv[_EAST_B_][_ENDY_] = lNy-2;
        iter_toRecv[_EAST_B_][_INIZ_] =   0;
        iter_toRecv[_EAST_B_][_ENDZ_] =   0;
    }

    //EAST FRONT
    if(neighb[_EAST_F_] != _NO_NEIGHBOUR_) {
        iter_toRecv[_EAST_F_][_INIX_] = lNx-1;
        iter_toRecv[_EAST_F_][_ENDX_] = lNx-1;
        iter_toRecv[_EAST_F_][_INIY_] =   1;
        iter_toRecv[_EAST_F_][_ENDY_] = lNy-2;
        iter_toRecv[_EAST_F_][_INIZ_] = lNz-1;
        iter_toRecv[_EAST_F_][_ENDZ_] = lNz-1;
    }

    //SOUTH BACK
    if(neighb[_SOUTH_B_] != _NO_NEIGHBOUR_) {
        iter_toRecv[_SOUTH_B_][_INIX_] =   1;
        iter_toRecv[_SOUTH_B_][_ENDX_] = lNx-2;
        iter_toRecv[_SOUTH_B_][_INIY_] =   0;
        iter_toRecv[_SOUTH_B_][_ENDY_] =   0;
        iter_toRecv[_SOUTH_B_][_INIZ_] =   0;
        iter_toRecv[_SOUTH_B_][_ENDZ_] =   0;
    }

    //SOUTH FRONT
    if(neighb[_SOUTH_F_] != _NO_NEIGHBOUR_) {
        iter_toRecv[_SOUTH_F_][_INIX_] =   1;
        iter_toRecv[_SOUTH_F_][_ENDX_] = lNx-2;
        iter_toRecv[_SOUTH_F_][_INIY_] =   0;
        iter_toRecv[_SOUTH_F_][_ENDY_] =   0;
        iter_toRecv[_SOUTH_F_][_INIZ_] = lNz-1;
        iter_toRecv[_SOUTH_F_][_ENDZ_] = lNz-1;
    }


    //NORTH BACK
    if(neighb[_NORTH_B_] != _NO_NEIGHBOUR_) {
        iter_toRecv[_NORTH_B_][_INIX_] =   1;
        iter_toRecv[_NORTH_B_][_ENDX_] = lNx-2;
        iter_toRecv[_NORTH_B_][_INIY_] = lNy-1;
        iter_toRecv[_NORTH_B_][_ENDY_] = lNy-1;
        iter_toRecv[_NORTH_B_][_INIZ_] =   0;
        iter_toRecv[_NORTH_B_][_ENDZ_] =   0;
    }

    //NORTH FRONT
    if(neighb[_NORTH_F_] != _NO_NEIGHBOUR_) {
        iter_toRecv[_NORTH_F_][_INIX_] =   1;
        iter_toRecv[_NORTH_F_][_ENDX_] = lNx-2;
        iter_toRecv[_NORTH_F_][_INIY_] = lNy-1;
        iter_toRecv[_NORTH_F_][_ENDY_] = lNy-1;
        iter_toRecv[_NORTH_F_][_INIZ_] = lNz-1;
        iter_toRecv[_NORTH_F_][_ENDZ_] = lNz-1;
    }

// 3er Level (corners)

    //WEST SOUTH BACK
    if(neighb[_WEST_S_B_] != _NO_NEIGHBOUR_) {
        iter_toRecv[_WEST_S_B_][_INIX_] =   0;
        iter_toRecv[_WEST_S_B_][_ENDX_] =   0;
        iter_toRecv[_WEST_S_B_][_INIY_] =   0;
        iter_toRecv[_WEST_S_B_][_ENDY_] =   0;
        iter_toRecv[_WEST_S_B_][_INIZ_] =   0;
        iter_toRecv[_WEST_S_B_][_ENDZ_] =   0;
    }

    //WEST NORTH BACK
    if(neighb[_WEST_N_B_] != _NO_NEIGHBOUR_) {
        iter_toRecv[_WEST_N_B_][_INIX_] =   0;
        iter_toRecv[_WEST_N_B_][_ENDX_] =   0;
        iter_toRecv[_WEST_N_B_][_INIY_] = lNy-1;
        iter_toRecv[_WEST_N_B_][_ENDY_] = lNy-1;
        iter_toRecv[_WEST_N_B_][_INIZ_] =   0;
        iter_toRecv[_WEST_N_B_][_ENDZ_] =   0;
    }

    //WEST SOUTH FRONT
    if(neighb[_WEST_S_F_] != _NO_NEIGHBOUR_) {
        iter_toRecv[_WEST_S_F_][_INIX_] =   0;
        iter_toRecv[_WEST_S_F_][_ENDX_] =   0;
        iter_toRecv[_WEST_S_F_][_INIY_] =   0;
        iter_toRecv[_WEST_S_F_][_ENDY_] =   0;
        iter_toRecv[_WEST_S_F_][_INIZ_] = lNz-1;
        iter_toRecv[_WEST_S_F_][_ENDZ_] = lNz-1;
    }

    //WEST NORTH FRONT
    if(neighb[_WEST_N_F_] != _NO_NEIGHBOUR_) {
        iter_toRecv[_WEST_N_F_][_INIX_] =   0;
        iter_toRecv[_WEST_N_F_][_ENDX_] =   0;
        iter_toRecv[_WEST_N_F_][_INIY_] = lNy-1;
        iter_toRecv[_WEST_N_F_][_ENDY_] = lNy-1;
        iter_toRecv[_WEST_N_F_][_INIZ_] = lNz-1;
        iter_toRecv[_WEST_N_F_][_ENDZ_] = lNz-1;
    }

 
    //EAST SOUTH BACK
    if(neighb[_EAST_S_B_] != _NO_NEIGHBOUR_) {
        iter_toRecv[_EAST_S_B_][_INIX_] = lNx-1;
        iter_toRecv[_EAST_S_B_][_ENDX_] = lNx-1;
        iter_toRecv[_EAST_S_B_][_INIY_] =   0;
        iter_toRecv[_EAST_S_B_][_ENDY_] =   0;
        iter_toRecv[_EAST_S_B_][_INIZ_] =   0;
        iter_toRecv[_EAST_S_B_][_ENDZ_] =   0;
    }

    //EAST NORTH BACK
    if(neighb[_EAST_N_B_] != _NO_NEIGHBOUR_) {
        iter_toRecv[_EAST_N_B_][_INIX_] = lNx-1;
        iter_toRecv[_EAST_N_B_][_ENDX_] = lNx-1;
        iter_toRecv[_EAST_N_B_][_INIY_] = lNy-1;
        iter_toRecv[_EAST_N_B_][_ENDY_] = lNy-1;
        iter_toRecv[_EAST_N_B_][_INIZ_] =   0;
        iter_toRecv[_EAST_N_B_][_ENDZ_] =   0;
    }


    //EAST SOUTH FRONT
    if(neighb[_EAST_S_F_] != _NO_NEIGHBOUR_) {
        iter_toRecv[_EAST_S_F_][_INIX_] = lNx-1;
        iter_toRecv[_EAST_S_F_][_ENDX_] = lNx-1;
        iter_toRecv[_EAST_S_F_][_INIY_] =   0;
        iter_toRecv[_EAST_S_F_][_ENDY_] =   0;
        iter_toRecv[_EAST_S_F_][_INIZ_] = lNz-1;
        iter_toRecv[_EAST_S_F_][_ENDZ_] = lNz-1;
    }

    //EAST NORTH FRONT
    if(neighb[_EAST_N_F_] != _NO_NEIGHBOUR_) {
        iter_toRecv[_EAST_N_F_][_INIX_] = lNx-1;
        iter_toRecv[_EAST_N_F_][_ENDX_] = lNx-1;
        iter_toRecv[_EAST_N_F_][_INIY_] = lNy-1;
        iter_toRecv[_EAST_N_F_][_ENDY_] = lNy-1;
        iter_toRecv[_EAST_N_F_][_INIZ_] = lNz-1;
        iter_toRecv[_EAST_N_F_][_ENDZ_] = lNz-1;
    }



}

/* This part assumes halo size 1*/
void ParallelTopology::create_toSend_iters()
{

    for(int bd=_WEST_; bd<= _EAST_N_F_; bd++)
    {
        iter_toSend[bd][_INIX_] =   0;
        iter_toSend[bd][_ENDX_] =  -1;
        iter_toSend[bd][_INIY_] =   0;
        iter_toSend[bd][_ENDY_] =  -1;
        iter_toSend[bd][_INIZ_] =   0;
        iter_toSend[bd][_ENDZ_] =  -1;
    }

    //WEST
    if(neighb[_WEST_] != _NO_NEIGHBOUR_ ) {
        iter_toSend[_WEST_][_INIX_] =   1;// <= this
        iter_toSend[_WEST_][_ENDX_] =   1;
        iter_toSend[_WEST_][_INIY_] =   1;
        iter_toSend[_WEST_][_ENDY_] = lNy-2;
        iter_toSend[_WEST_][_INIZ_] =   1;
        iter_toSend[_WEST_][_ENDZ_] = lNz-2;
    }

    //EAST
    if(neighb[_EAST_] != _NO_NEIGHBOUR_ ) {
        iter_toSend[_EAST_][_INIX_] = lNx-2; // <=this
        iter_toSend[_EAST_][_ENDX_] = lNx-2;
        iter_toSend[_EAST_][_INIY_] =   1;
        iter_toSend[_EAST_][_ENDY_] = lNy-2;
        iter_toSend[_EAST_][_INIZ_] =   1;
        iter_toSend[_EAST_][_ENDZ_] = lNz-2;
    }

    //SOUTH
    if(neighb[_SOUTH_] != _NO_NEIGHBOUR_ ) {
        iter_toSend[_SOUTH_][_INIX_] =   1;
        iter_toSend[_SOUTH_][_ENDX_] = lNx-2;
        iter_toSend[_SOUTH_][_INIY_] =   1; // <=this
        iter_toSend[_SOUTH_][_ENDY_] =   1;
        iter_toSend[_SOUTH_][_INIZ_] =   1;
        iter_toSend[_SOUTH_][_ENDZ_] = lNz-2;
    }

    //NORTH
    if(neighb[_NORTH_] != _NO_NEIGHBOUR_ ) {
        iter_toSend[_NORTH_][_INIX_] =   1;
        iter_toSend[_NORTH_][_ENDX_] = lNx-2;
        iter_toSend[_NORTH_][_INIY_] = lNy-2; //<=this
        iter_toSend[_NORTH_][_ENDY_] = lNy-2;
        iter_toSend[_NORTH_][_INIZ_] =   1;
        iter_toSend[_NORTH_][_ENDZ_] = lNz-2;
    }

    //BACK
    if(neighb[_BACK_] != _NO_NEIGHBOUR_ ) {
        iter_toSend[_BACK_][_INIX_] =   1;
        iter_toSend[_BACK_][_ENDX_] = lNx-2;
        iter_toSend[_BACK_][_INIY_] =   1;
        iter_toSend[_BACK_][_ENDY_] = lNy-2;
        iter_toSend[_BACK_][_INIZ_] =   1; //<= this
        iter_toSend[_BACK_][_ENDZ_] =   1;
    }

    //FRONT
    if(neighb[_FRONT_] != _NO_NEIGHBOUR_) {
        iter_toSend[_FRONT_][_INIX_] =   1;
        iter_toSend[_FRONT_][_ENDX_] = lNx-2;
        iter_toSend[_FRONT_][_INIY_] =   1;
        iter_toSend[_FRONT_][_ENDY_] = lNy-2;
        iter_toSend[_FRONT_][_INIZ_] = lNz-2; //<= this
        iter_toSend[_FRONT_][_ENDZ_] = lNz-2;
    }

}

void ParallelTopology::create_complex_toSend_iters()
{


// 2nd Level Sends

    //WEST SOUTH
    if(neighb[_WEST_S_] != _NO_NEIGHBOUR_ ) {
        if(neighb[_WEST_] != _NO_NEIGHBOUR_ && neighb[_SOUTH_] != _NO_NEIGHBOUR_){  
            iter_toSend[_WEST_S_][_INIX_] =   1;// <= this
            iter_toSend[_WEST_S_][_ENDX_] =   1;
            iter_toSend[_WEST_S_][_INIY_] =   1;
            iter_toSend[_WEST_S_][_ENDY_] =   1;
            iter_toSend[_WEST_S_][_INIZ_] =   1;
            iter_toSend[_WEST_S_][_ENDZ_] = lNz-2;
        }
        else if(neighb[_WEST_] != _NO_NEIGHBOUR_ )
         {

                iter_toSend[_WEST_S_][_INIX_] =   1;// <= this
                iter_toSend[_WEST_S_][_ENDX_] =   1;
                iter_toSend[_WEST_S_][_INIY_] =   0;
                iter_toSend[_WEST_S_][_ENDY_] =   0;
                iter_toSend[_WEST_S_][_INIZ_] =   1;
                iter_toSend[_WEST_S_][_ENDZ_] = lNz-2;

         }
        else if( neighb[_SOUTH_] != _NO_NEIGHBOUR_)
         {
                 iter_toSend[_WEST_S_][_INIX_] =   0;// <= this
                 iter_toSend[_WEST_S_][_ENDX_] =   0;
                 iter_toSend[_WEST_S_][_INIY_] = lNy-2;
                 iter_toSend[_WEST_S_][_ENDY_] = lNy-2;
                 iter_toSend[_WEST_S_][_INIZ_] =   1;
                 iter_toSend[_WEST_S_][_ENDZ_] = lNz-2;
         }
    }
    //WEST NORTH
    if(neighb[_WEST_N_] != _NO_NEIGHBOUR_ ) {
        if(neighb[_WEST_] != _NO_NEIGHBOUR_ && neighb[_NORTH_] != _NO_NEIGHBOUR_){  

            iter_toSend[_WEST_N_][_INIX_] =   1;// <= this
            iter_toSend[_WEST_N_][_ENDX_] =   1;
            iter_toSend[_WEST_N_][_INIY_] = lNy-2;
            iter_toSend[_WEST_N_][_ENDY_] = lNy-2;
            iter_toSend[_WEST_N_][_INIZ_] =   1;
            iter_toSend[_WEST_N_][_ENDZ_] = lNz-2;
        }
        else if(neighb[_WEST_]  != _NO_NEIGHBOUR_){  

            iter_toSend[_WEST_N_][_INIX_] =   1;// <= this
            iter_toSend[_WEST_N_][_ENDX_] =   1;
            iter_toSend[_WEST_N_][_INIY_] = lNy-1;
            iter_toSend[_WEST_N_][_ENDY_] = lNy-1;
            iter_toSend[_WEST_N_][_INIZ_] =   1;
            iter_toSend[_WEST_N_][_ENDZ_] = lNz-2;
        }
        else if(neighb[_NORTH_] != _NO_NEIGHBOUR_){  

            iter_toSend[_WEST_N_][_INIX_] =   0;// <= this
            iter_toSend[_WEST_N_][_ENDX_] =   0;
            iter_toSend[_WEST_N_][_INIY_] = lNy-2;
            iter_toSend[_WEST_N_][_ENDY_] = lNy-2;
            iter_toSend[_WEST_N_][_INIZ_] =   1;
            iter_toSend[_WEST_N_][_ENDZ_] = lNz-2;

        }
    }



    //WEST BACK
    if(neighb[_WEST_B_] != _NO_NEIGHBOUR_ ) {
        if(neighb[_WEST_] != _NO_NEIGHBOUR_ && neighb[_BACK_] != _NO_NEIGHBOUR_){  
            iter_toSend[_WEST_B_][_INIX_] =   1;// <= this
            iter_toSend[_WEST_B_][_ENDX_] =   1;    
            iter_toSend[_WEST_B_][_INIY_] =   1;
            iter_toSend[_WEST_B_][_ENDY_] = lNy-2;
            iter_toSend[_WEST_B_][_INIZ_] =   1;
            iter_toSend[_WEST_B_][_ENDZ_] =   1;
        }
        else if(neighb[_WEST_] != _NO_NEIGHBOUR_){  
            iter_toSend[_WEST_B_][_INIX_] =   1;// <= this
            iter_toSend[_WEST_B_][_ENDX_] =   1;    
            iter_toSend[_WEST_B_][_INIY_] =   1;
            iter_toSend[_WEST_B_][_ENDY_] = lNy-2;
            iter_toSend[_WEST_B_][_INIZ_] =   0;
            iter_toSend[_WEST_B_][_ENDZ_] =   0;

        }
        else if(neighb[_BACK_] != _NO_NEIGHBOUR_){  
            iter_toSend[_WEST_B_][_INIX_] =   0;// <= this
            iter_toSend[_WEST_B_][_ENDX_] =   0;    
            iter_toSend[_WEST_B_][_INIY_] =   1;
            iter_toSend[_WEST_B_][_ENDY_] = lNy-2;
            iter_toSend[_WEST_B_][_INIZ_] =   1;
            iter_toSend[_WEST_B_][_ENDZ_] =   1;

        } 
    }

    //WEST FRONT
    if(neighb[_WEST_F_] != _NO_NEIGHBOUR_ ) {
        if(neighb[_WEST_] != _NO_NEIGHBOUR_ && neighb[_FRONT_] != _NO_NEIGHBOUR_){  
            iter_toSend[_WEST_F_][_INIX_] =   1;// <= this
            iter_toSend[_WEST_F_][_ENDX_] =   1;
            iter_toSend[_WEST_F_][_INIY_] =   1;
            iter_toSend[_WEST_F_][_ENDY_] = lNy-2;
            iter_toSend[_WEST_F_][_INIZ_] = lNz-2;
            iter_toSend[_WEST_F_][_ENDZ_] = lNz-2;
        }
        else if(neighb[_WEST_] != _NO_NEIGHBOUR_ ){  
            iter_toSend[_WEST_F_][_INIX_] =   1;// <= this
            iter_toSend[_WEST_F_][_ENDX_] =   1;
            iter_toSend[_WEST_F_][_INIY_] =   1;
            iter_toSend[_WEST_F_][_ENDY_] = lNy-2;
            iter_toSend[_WEST_F_][_INIZ_] = lNz-1;
            iter_toSend[_WEST_F_][_ENDZ_] = lNz-1;

        }
        else if(neighb[_FRONT_] != _NO_NEIGHBOUR_){  
            iter_toSend[_WEST_F_][_INIX_] =   0;// <= this
            iter_toSend[_WEST_F_][_ENDX_] =   0;
            iter_toSend[_WEST_F_][_INIY_] =   1;
            iter_toSend[_WEST_F_][_ENDY_] = lNy-2;
            iter_toSend[_WEST_F_][_INIZ_] = lNz-2;
            iter_toSend[_WEST_F_][_ENDZ_] = lNz-2;

        }
 
    }

    //EAST SOUTH
    if(neighb[_EAST_S_] != _NO_NEIGHBOUR_ ) {
        if(neighb[_EAST_] != _NO_NEIGHBOUR_ && neighb[_SOUTH_] != _NO_NEIGHBOUR_){  
            iter_toSend[_EAST_S_][_INIX_] = lNx-2; // <=this
            iter_toSend[_EAST_S_][_ENDX_] = lNx-2;
            iter_toSend[_EAST_S_][_INIY_] =   1;
            iter_toSend[_EAST_S_][_ENDY_] =   1;
            iter_toSend[_EAST_S_][_INIZ_] =   1;
            iter_toSend[_EAST_S_][_ENDZ_] = lNz-2;
        }
        else if(neighb[_EAST_] != _NO_NEIGHBOUR_){ 
            iter_toSend[_EAST_S_][_INIX_] = lNx-2; // <=this
            iter_toSend[_EAST_S_][_ENDX_] = lNx-2;
            iter_toSend[_EAST_S_][_INIY_] =   0;
            iter_toSend[_EAST_S_][_ENDY_] =   0;
            iter_toSend[_EAST_S_][_INIZ_] =   1;
            iter_toSend[_EAST_S_][_ENDZ_] = lNz-2;
       }
       else if(neighb[_SOUTH_] != _NO_NEIGHBOUR_){ 
            iter_toSend[_EAST_S_][_INIX_] = lNx-1; // <=this
            iter_toSend[_EAST_S_][_ENDX_] = lNx-1;
            iter_toSend[_EAST_S_][_INIY_] =   1;
            iter_toSend[_EAST_S_][_ENDY_] =   1;
            iter_toSend[_EAST_S_][_INIZ_] =   1;
            iter_toSend[_EAST_S_][_ENDZ_] = lNz-2;
       } 
    }

    //EAST NORTH
    if(neighb[_EAST_N_] != _NO_NEIGHBOUR_ ) {
        if(neighb[_EAST_] != _NO_NEIGHBOUR_ && neighb[_NORTH_] != _NO_NEIGHBOUR_){  
            iter_toSend[_EAST_N_][_INIX_] = lNx-2; // <=this
            iter_toSend[_EAST_N_][_ENDX_] = lNx-2;
            iter_toSend[_EAST_N_][_INIY_] = lNy-2;
            iter_toSend[_EAST_N_][_ENDY_] = lNy-2;
            iter_toSend[_EAST_N_][_INIZ_] =   1;
            iter_toSend[_EAST_N_][_ENDZ_] = lNz-2;
        }
        else if(neighb[_EAST_] != _NO_NEIGHBOUR_){  
            iter_toSend[_EAST_N_][_INIX_] = lNx-2; // <=this
            iter_toSend[_EAST_N_][_ENDX_] = lNx-2;
            iter_toSend[_EAST_N_][_INIY_] = lNy-1;
            iter_toSend[_EAST_N_][_ENDY_] = lNy-1;
            iter_toSend[_EAST_N_][_INIZ_] =   1;
            iter_toSend[_EAST_N_][_ENDZ_] = lNz-2;  
        }
        else if(neighb[_NORTH_] != _NO_NEIGHBOUR_){  
            iter_toSend[_EAST_N_][_INIX_] = lNx-1; // <=this
            iter_toSend[_EAST_N_][_ENDX_] = lNx-1;
            iter_toSend[_EAST_N_][_INIY_] = lNy-2;
            iter_toSend[_EAST_N_][_ENDY_] = lNy-2;
            iter_toSend[_EAST_N_][_INIZ_] =   1;
            iter_toSend[_EAST_N_][_ENDZ_] = lNz-2;      
        }
    }

    //EAST BACK
    if(neighb[_EAST_B_] != _NO_NEIGHBOUR_ ) {
        if(neighb[_EAST_] != _NO_NEIGHBOUR_ && neighb[_BACK_] != _NO_NEIGHBOUR_){  
            iter_toSend[_EAST_B_][_INIX_] = lNx-2; // <=this
            iter_toSend[_EAST_B_][_ENDX_] = lNx-2;
            iter_toSend[_EAST_B_][_INIY_] =   1;
            iter_toSend[_EAST_B_][_ENDY_] = lNy-2;
            iter_toSend[_EAST_B_][_INIZ_] =   1;
            iter_toSend[_EAST_B_][_ENDZ_] =   1;
        }
        else if ( neighb[_EAST_] != _NO_NEIGHBOUR_ )
        {
            iter_toSend[_EAST_B_][_INIX_] = lNx-2; // <=this
            iter_toSend[_EAST_B_][_ENDX_] = lNx-2;
            iter_toSend[_EAST_B_][_INIY_] =   1;
            iter_toSend[_EAST_B_][_ENDY_] = lNy-2;
            iter_toSend[_EAST_B_][_INIZ_] = lNz-1;
            iter_toSend[_EAST_B_][_ENDZ_] = lNz-1;


        }
        else if( neighb[_BACK_] != _NO_NEIGHBOUR_ ){
            iter_toSend[_EAST_B_][_INIX_] = lNx-1; // <=this
            iter_toSend[_EAST_B_][_ENDX_] = lNx-1;
            iter_toSend[_EAST_B_][_INIY_] =   1;
            iter_toSend[_EAST_B_][_ENDY_] = lNy-2;
            iter_toSend[_EAST_B_][_INIZ_] =   1;
            iter_toSend[_EAST_B_][_ENDZ_] =   1;
        }

    }

    //EAST FRONT
    if(neighb[_EAST_F_] != _NO_NEIGHBOUR_ ) {
        if(neighb[_EAST_] != _NO_NEIGHBOUR_ && neighb[_FRONT_] != _NO_NEIGHBOUR_){  
            iter_toSend[_EAST_F_][_INIX_] = lNx-2; // <=this
            iter_toSend[_EAST_F_][_ENDX_] = lNx-2;
            iter_toSend[_EAST_F_][_INIY_] =   1;
            iter_toSend[_EAST_F_][_ENDY_] = lNy-2;
            iter_toSend[_EAST_F_][_INIZ_] = lNz-2;
            iter_toSend[_EAST_F_][_ENDZ_] = lNz-2;
        }
        else if(neighb[_EAST_] != _NO_NEIGHBOUR_){ 
                iter_toSend[_EAST_F_][_INIX_] = lNx-2; // <=this
                iter_toSend[_EAST_F_][_ENDX_] = lNx-2;
                iter_toSend[_EAST_F_][_INIY_] =   1;
                iter_toSend[_EAST_F_][_ENDY_] = lNy-2;
                iter_toSend[_EAST_F_][_INIZ_] = lNz-1;
                iter_toSend[_EAST_F_][_ENDZ_] = lNz-1;       
            } 
        else if(neighb[_FRONT_] != _NO_NEIGHBOUR_){  
                iter_toSend[_EAST_F_][_INIX_] = lNx-1; // <=this
                iter_toSend[_EAST_F_][_ENDX_] = lNx-1;
                iter_toSend[_EAST_F_][_INIY_] =   1;
                iter_toSend[_EAST_F_][_ENDY_] = lNy-2;
                iter_toSend[_EAST_F_][_INIZ_] = lNz-2;
                iter_toSend[_EAST_F_][_ENDZ_] = lNz-2;      
            }
        }


    //SOUTH BACK
    if(neighb[_SOUTH_B_] != _NO_NEIGHBOUR_ ) {
        if(neighb[_SOUTH_] != _NO_NEIGHBOUR_ && neighb[_BACK_] != _NO_NEIGHBOUR_){  
            iter_toSend[_SOUTH_B_][_INIX_] =   1;
            iter_toSend[_SOUTH_B_][_ENDX_] = lNx-2;
            iter_toSend[_SOUTH_B_][_INIY_] =   1; // <=this
            iter_toSend[_SOUTH_B_][_ENDY_] =   1;
            iter_toSend[_SOUTH_B_][_INIZ_] =   1;
            iter_toSend[_SOUTH_B_][_ENDZ_] =   1;
        }
        else if(neighb[_SOUTH_] != _NO_NEIGHBOUR_){  
            iter_toSend[_SOUTH_B_][_INIX_] =   1;
            iter_toSend[_SOUTH_B_][_ENDX_] = lNx-2;
            iter_toSend[_SOUTH_B_][_INIY_] =   1; // <=this
            iter_toSend[_SOUTH_B_][_ENDY_] =   1;
            iter_toSend[_SOUTH_B_][_INIZ_] =   0;
            iter_toSend[_SOUTH_B_][_ENDZ_] =   0;

        }
        else if(neighb[_BACK_] != _NO_NEIGHBOUR_){  
            iter_toSend[_SOUTH_B_][_INIX_] =   1;
            iter_toSend[_SOUTH_B_][_ENDX_] = lNx-2;
            iter_toSend[_SOUTH_B_][_INIY_] =   0; // <=this
            iter_toSend[_SOUTH_B_][_ENDY_] =   0;
            iter_toSend[_SOUTH_B_][_INIZ_] =   1;
            iter_toSend[_SOUTH_B_][_ENDZ_] =   1;

        } 
    }

    //SOUTH FRONT
    if(neighb[_SOUTH_F_] != _NO_NEIGHBOUR_ ) {
        if(neighb[_SOUTH_] != _NO_NEIGHBOUR_ && neighb[_FRONT_] != _NO_NEIGHBOUR_){  
            iter_toSend[_SOUTH_F_][_INIX_] =   1;
            iter_toSend[_SOUTH_F_][_ENDX_] = lNx-2;
            iter_toSend[_SOUTH_F_][_INIY_] =   1; // <=this
            iter_toSend[_SOUTH_F_][_ENDY_] =   1;
            iter_toSend[_SOUTH_F_][_INIZ_] = lNz-2;
            iter_toSend[_SOUTH_F_][_ENDZ_] = lNz-2;
        }
        else if(neighb[_SOUTH_] != _NO_NEIGHBOUR_){  
            iter_toSend[_SOUTH_F_][_INIX_] =   1;
            iter_toSend[_SOUTH_F_][_ENDX_] = lNx-2;
            iter_toSend[_SOUTH_F_][_INIY_] =   1; // <=this
            iter_toSend[_SOUTH_F_][_ENDY_] =   1;
            iter_toSend[_SOUTH_F_][_INIZ_] = lNz-1;
            iter_toSend[_SOUTH_F_][_ENDZ_] = lNz-1;
        }
        else if(neighb[_FRONT_] != _NO_NEIGHBOUR_){  
            iter_toSend[_SOUTH_F_][_INIX_] =   1;
            iter_toSend[_SOUTH_F_][_ENDX_] = lNx-2;
            iter_toSend[_SOUTH_F_][_INIY_] =   0; // <=this
            iter_toSend[_SOUTH_F_][_ENDY_] =   0;
            iter_toSend[_SOUTH_F_][_INIZ_] = lNz-2;
            iter_toSend[_SOUTH_F_][_ENDZ_] = lNz-2;
        }
    }


    //NORTH BACK
    if(neighb[_NORTH_B_] != _NO_NEIGHBOUR_ ) {
        if(neighb[_NORTH_] != _NO_NEIGHBOUR_ && neighb[_BACK_] != _NO_NEIGHBOUR_){  
            iter_toSend[_NORTH_B_][_INIX_] =   1;
            iter_toSend[_NORTH_B_][_ENDX_] = lNx-2;
            iter_toSend[_NORTH_B_][_INIY_] = lNy-2; //<=this
            iter_toSend[_NORTH_B_][_ENDY_] = lNy-2;
            iter_toSend[_NORTH_B_][_INIZ_] =   1;
            iter_toSend[_NORTH_B_][_ENDZ_] =   1;
        }
        else if(neighb[_NORTH_] != _NO_NEIGHBOUR_){  
            iter_toSend[_NORTH_B_][_INIX_] =   1;
            iter_toSend[_NORTH_B_][_ENDX_] = lNx-2;
            iter_toSend[_NORTH_B_][_INIY_] = lNy-2; //<=this
            iter_toSend[_NORTH_B_][_ENDY_] = lNy-2;
            iter_toSend[_NORTH_B_][_INIZ_] =   0;
            iter_toSend[_NORTH_B_][_ENDZ_] =   0;
        }
        else if(neighb[_BACK_] != _NO_NEIGHBOUR_){  
            iter_toSend[_NORTH_B_][_INIX_] =   1;
            iter_toSend[_NORTH_B_][_ENDX_] = lNx-2;
            iter_toSend[_NORTH_B_][_INIY_] = lNy-1; //<=this
            iter_toSend[_NORTH_B_][_ENDY_] = lNy-1;
            iter_toSend[_NORTH_B_][_INIZ_] =   1;
            iter_toSend[_NORTH_B_][_ENDZ_] =   1;
        }
    }

    //NORTH FRONT
    if(neighb[_NORTH_F_] != _NO_NEIGHBOUR_ ) {
        if(neighb[_NORTH_] != _NO_NEIGHBOUR_ && neighb[_FRONT_] != _NO_NEIGHBOUR_){  
            iter_toSend[_NORTH_F_][_INIX_] =   1;
            iter_toSend[_NORTH_F_][_ENDX_] = lNx-2;
            iter_toSend[_NORTH_F_][_INIY_] = lNy-2; //<=this
            iter_toSend[_NORTH_F_][_ENDY_] = lNy-2;
            iter_toSend[_NORTH_F_][_INIZ_] = lNz-2;
            iter_toSend[_NORTH_F_][_ENDZ_] = lNz-2;
        }
        else if(neighb[_NORTH_] != _NO_NEIGHBOUR_){  
            iter_toSend[_NORTH_F_][_INIX_] =   1;
            iter_toSend[_NORTH_F_][_ENDX_] = lNx-2;
            iter_toSend[_NORTH_F_][_INIY_] = lNy-2; //<=this
            iter_toSend[_NORTH_F_][_ENDY_] = lNy-2;
            iter_toSend[_NORTH_F_][_INIZ_] = lNz-1;
            iter_toSend[_NORTH_F_][_ENDZ_] = lNz-1;
        }
        else if(neighb[_FRONT_] != _NO_NEIGHBOUR_){  
            iter_toSend[_NORTH_F_][_INIX_] =   1;
            iter_toSend[_NORTH_F_][_ENDX_] = lNx-2;
            iter_toSend[_NORTH_F_][_INIY_] = lNy-1; //<=this
            iter_toSend[_NORTH_F_][_ENDY_] = lNy-1;
            iter_toSend[_NORTH_F_][_INIZ_] = lNz-2;
            iter_toSend[_NORTH_F_][_ENDZ_] = lNz-2;
        }   
    }


// 3er level

//WEST SOUTH BACK
    if(neighb[_WEST_S_B_] != _NO_NEIGHBOUR_ ) {
        if(neighb[_WEST_] != _NO_NEIGHBOUR_ && neighb[_SOUTH_] != _NO_NEIGHBOUR_ &&  neighb[_BACK_] != _NO_NEIGHBOUR_ ){  

            iter_toSend[_WEST_S_B_][_INIX_] =   1;// <= this
            iter_toSend[_WEST_S_B_][_ENDX_] =   1;
            iter_toSend[_WEST_S_B_][_INIY_] =   1;
            iter_toSend[_WEST_S_B_][_ENDY_] =   1;
            iter_toSend[_WEST_S_B_][_INIZ_] =   1;
            iter_toSend[_WEST_S_B_][_ENDZ_] =   1;
        }
        else
        {
            if(neighb[_WEST_S_] != _NO_NEIGHBOUR_ && neighb[_BACK_] == _NO_NEIGHBOUR_ ){ 

                iter_toSend[_WEST_S_B_][_INIX_] =  iter_toSend[_WEST_S_][_INIX_];// <= this
                iter_toSend[_WEST_S_B_][_ENDX_] =  iter_toSend[_WEST_S_][_ENDX_];
                iter_toSend[_WEST_S_B_][_INIY_] =  iter_toSend[_WEST_S_][_INIY_];
                iter_toSend[_WEST_S_B_][_ENDY_] =  iter_toSend[_WEST_S_][_ENDY_];
                iter_toSend[_WEST_S_B_][_INIZ_] =  0;
                iter_toSend[_WEST_S_B_][_ENDZ_] =  0;
            }
            else if(neighb[_WEST_S_] != _NO_NEIGHBOUR_ && neighb[_BACK_] != _NO_NEIGHBOUR_ ){
                
                iter_toSend[_WEST_S_B_][_INIX_] =  iter_toSend[_WEST_S_][_INIX_];// <= this
                iter_toSend[_WEST_S_B_][_ENDX_] =  iter_toSend[_WEST_S_][_ENDX_];
                iter_toSend[_WEST_S_B_][_INIY_] =  iter_toSend[_WEST_S_][_INIY_];
                iter_toSend[_WEST_S_B_][_ENDY_] =  iter_toSend[_WEST_S_][_ENDY_];
                iter_toSend[_WEST_S_B_][_INIZ_] =  1;
                iter_toSend[_WEST_S_B_][_ENDZ_] =  1;
            }
             if(neighb[_WEST_B_] != _NO_NEIGHBOUR_ && neighb[_SOUTH_] == _NO_NEIGHBOUR_ ){ 

                iter_toSend[_WEST_S_B_][_INIX_] =  iter_toSend[_WEST_B_][_INIX_];// <= this
                iter_toSend[_WEST_S_B_][_ENDX_] =  iter_toSend[_WEST_B_][_ENDX_];
                iter_toSend[_WEST_S_B_][_INIY_] =  0;
                iter_toSend[_WEST_S_B_][_ENDY_] =  0;
                iter_toSend[_WEST_S_B_][_INIZ_] =  iter_toSend[_WEST_B_][_INIZ_];
                iter_toSend[_WEST_S_B_][_ENDZ_] =  iter_toSend[_WEST_B_][_ENDZ_];
            }
            else if(neighb[_WEST_B_] != _NO_NEIGHBOUR_ && neighb[_SOUTH_] != _NO_NEIGHBOUR_ ){
                
                iter_toSend[_WEST_S_B_][_INIX_] =  iter_toSend[_WEST_B_][_INIX_];// <= this
                iter_toSend[_WEST_S_B_][_ENDX_] =  iter_toSend[_WEST_B_][_ENDX_];
                iter_toSend[_WEST_S_B_][_INIY_] =  1;
                iter_toSend[_WEST_S_B_][_ENDY_] =  1;
                iter_toSend[_WEST_S_B_][_INIZ_] =  iter_toSend[_WEST_B_][_INIZ_];
                iter_toSend[_WEST_S_B_][_ENDZ_] =  iter_toSend[_WEST_B_][_ENDZ_];
            }
             if(neighb[_SOUTH_B_] != _NO_NEIGHBOUR_ && neighb[_WEST_] == _NO_NEIGHBOUR_ ){ 

                iter_toSend[_WEST_S_B_][_INIX_] =  0;// <= this
                iter_toSend[_WEST_S_B_][_ENDX_] =  0;
                iter_toSend[_WEST_S_B_][_INIY_] =  iter_toSend[_SOUTH_B_][_INIY_];
                iter_toSend[_WEST_S_B_][_ENDY_] =  iter_toSend[_SOUTH_B_][_ENDY_];
                iter_toSend[_WEST_S_B_][_INIZ_] =  iter_toSend[_SOUTH_B_][_INIZ_];
                iter_toSend[_WEST_S_B_][_ENDZ_] =  iter_toSend[_SOUTH_B_][_ENDZ_];
            }
            else if(neighb[_SOUTH_B_] != _NO_NEIGHBOUR_ && neighb[_WEST_] != _NO_NEIGHBOUR_ ){
                
                iter_toSend[_WEST_S_B_][_INIX_] =  1;// <= this
                iter_toSend[_WEST_S_B_][_ENDX_] =  1;
                iter_toSend[_WEST_S_B_][_INIY_] =  iter_toSend[_SOUTH_B_][_INIY_];
                iter_toSend[_WEST_S_B_][_ENDY_] =  iter_toSend[_SOUTH_B_][_ENDY_];
                iter_toSend[_WEST_S_B_][_INIZ_] =  iter_toSend[_SOUTH_B_][_INIZ_];
                iter_toSend[_WEST_S_B_][_ENDZ_] =  iter_toSend[_SOUTH_B_][_ENDZ_];
            }

        }
    }

    //WEST NORTH BACK
    if(neighb[_WEST_N_B_] != _NO_NEIGHBOUR_ ) {
         if(neighb[_WEST_] != _NO_NEIGHBOUR_ && neighb[_NORTH_] != _NO_NEIGHBOUR_ &&  neighb[_BACK_] != _NO_NEIGHBOUR_ ){  

       iter_toSend[_WEST_N_B_][_INIX_] =   1;// <= this
        iter_toSend[_WEST_N_B_][_ENDX_] =   1;
        iter_toSend[_WEST_N_B_][_INIY_] = lNy-2;
        iter_toSend[_WEST_N_B_][_ENDY_] = lNy-2;
        iter_toSend[_WEST_N_B_][_INIZ_] =   1;
        iter_toSend[_WEST_N_B_][_ENDZ_] =   1;
         }
        else
        {
            if(neighb[_WEST_N_] != _NO_NEIGHBOUR_ && neighb[_BACK_] == _NO_NEIGHBOUR_ ){ 

                iter_toSend[_WEST_N_B_][_INIX_] =  iter_toSend[_WEST_N_][_INIX_];// <= this
                iter_toSend[_WEST_N_B_][_ENDX_] =  iter_toSend[_WEST_N_][_ENDX_];
                iter_toSend[_WEST_N_B_][_INIY_] =  iter_toSend[_WEST_N_][_INIY_];
                iter_toSend[_WEST_N_B_][_ENDY_] =  iter_toSend[_WEST_N_][_ENDY_];
                iter_toSend[_WEST_N_B_][_INIZ_] =  0;
                iter_toSend[_WEST_N_B_][_ENDZ_] =  0;
            }
            else if(neighb[_WEST_N_] != _NO_NEIGHBOUR_ && neighb[_BACK_] != _NO_NEIGHBOUR_ ){
                
                iter_toSend[_WEST_N_B_][_INIX_] =  iter_toSend[_WEST_N_][_INIX_];// <= this
                iter_toSend[_WEST_N_B_][_ENDX_] =  iter_toSend[_WEST_N_][_ENDX_];
                iter_toSend[_WEST_N_B_][_INIY_] =  iter_toSend[_WEST_N_][_INIY_];
                iter_toSend[_WEST_N_B_][_ENDY_] =  iter_toSend[_WEST_N_][_ENDY_];
                iter_toSend[_WEST_N_B_][_INIZ_] =  1;
                iter_toSend[_WEST_N_B_][_ENDZ_] =  1;
            }
             if(neighb[_WEST_B_] != _NO_NEIGHBOUR_ && neighb[_NORTH_] == _NO_NEIGHBOUR_ ){ 

                iter_toSend[_WEST_N_B_][_INIX_] =  iter_toSend[_WEST_B_][_INIX_];// <= this
                iter_toSend[_WEST_N_B_][_ENDX_] =  iter_toSend[_WEST_B_][_ENDX_];
                iter_toSend[_WEST_N_B_][_INIY_] =  lNy-1;
                iter_toSend[_WEST_N_B_][_ENDY_] =  lNy-1;
                iter_toSend[_WEST_N_B_][_INIZ_] =  iter_toSend[_WEST_B_][_INIZ_];
                iter_toSend[_WEST_N_B_][_ENDZ_] =  iter_toSend[_WEST_B_][_ENDZ_];
            }
            else if(neighb[_WEST_B_] != _NO_NEIGHBOUR_ && neighb[_NORTH_] != _NO_NEIGHBOUR_ ){
                
                iter_toSend[_WEST_N_B_][_INIX_] =  iter_toSend[_WEST_B_][_INIX_];// <= this
                iter_toSend[_WEST_N_B_][_ENDX_] =  iter_toSend[_WEST_B_][_ENDX_];
                iter_toSend[_WEST_N_B_][_INIY_] =  lNy-2;
                iter_toSend[_WEST_N_B_][_ENDY_] =  lNy-2;
                iter_toSend[_WEST_N_B_][_INIZ_] =  iter_toSend[_WEST_B_][_INIZ_];
                iter_toSend[_WEST_N_B_][_ENDZ_] =  iter_toSend[_WEST_B_][_ENDZ_];
            }
             if(neighb[_NORTH_B_] != _NO_NEIGHBOUR_ && neighb[_WEST_] == _NO_NEIGHBOUR_ ){ 

                iter_toSend[_WEST_N_B_][_INIX_] =  0;// <= this
                iter_toSend[_WEST_N_B_][_ENDX_] =  0;
                iter_toSend[_WEST_N_B_][_INIY_] =  iter_toSend[_NORTH_B_][_INIY_];
                iter_toSend[_WEST_N_B_][_ENDY_] =  iter_toSend[_NORTH_B_][_ENDY_];
                iter_toSend[_WEST_N_B_][_INIZ_] =  iter_toSend[_NORTH_B_][_INIZ_];
                iter_toSend[_WEST_N_B_][_ENDZ_] =  iter_toSend[_NORTH_B_][_ENDZ_];
            }
            else if(neighb[_NORTH_B_] != _NO_NEIGHBOUR_ && neighb[_WEST_] != _NO_NEIGHBOUR_ ){
                
                iter_toSend[_WEST_N_B_][_INIX_] =  1;// <= this
                iter_toSend[_WEST_N_B_][_ENDX_] =  1;
                iter_toSend[_WEST_N_B_][_INIY_] =  iter_toSend[_NORTH_B_][_INIY_];
                iter_toSend[_WEST_N_B_][_ENDY_] =  iter_toSend[_NORTH_B_][_ENDY_];
                iter_toSend[_WEST_N_B_][_INIZ_] =  iter_toSend[_NORTH_B_][_INIZ_];
                iter_toSend[_WEST_N_B_][_ENDZ_] =  iter_toSend[_NORTH_B_][_ENDZ_];
            }

        }

    }

    //WEST SOUTH FRONT
    if(neighb[_WEST_S_F_] != _NO_NEIGHBOUR_ ) {
        if(neighb[_WEST_] != _NO_NEIGHBOUR_ && neighb[_SOUTH_] != _NO_NEIGHBOUR_ &&  neighb[_FRONT_] != _NO_NEIGHBOUR_ ){  

            iter_toSend[_WEST_S_F_][_INIX_] =   1;// <= this
            iter_toSend[_WEST_S_F_][_ENDX_] =   1;
            iter_toSend[_WEST_S_F_][_INIY_] =   1;
            iter_toSend[_WEST_S_F_][_ENDY_] =   1;
            iter_toSend[_WEST_S_F_][_INIZ_] = lNz-2;
            iter_toSend[_WEST_S_F_][_ENDZ_] = lNz-2;
        }
        else
        {
            if(neighb[_WEST_S_] != _NO_NEIGHBOUR_ && neighb[_FRONT_] == _NO_NEIGHBOUR_ ){ 

                iter_toSend[_WEST_S_F_][_INIX_] =  iter_toSend[_WEST_S_][_INIX_];// <= this
                iter_toSend[_WEST_S_F_][_ENDX_] =  iter_toSend[_WEST_S_][_ENDX_];
                iter_toSend[_WEST_S_F_][_INIY_] =  iter_toSend[_WEST_S_][_INIY_];
                iter_toSend[_WEST_S_F_][_ENDY_] =  iter_toSend[_WEST_S_][_ENDY_];
                iter_toSend[_WEST_S_F_][_INIZ_] =  lNz-1;
                iter_toSend[_WEST_S_F_][_ENDZ_] =  lNz-1;
            }
            else if(neighb[_WEST_S_] != _NO_NEIGHBOUR_ && neighb[_FRONT_] != _NO_NEIGHBOUR_ ){

                iter_toSend[_WEST_S_F_][_INIX_] =  iter_toSend[_WEST_S_][_INIX_];// <= this
                iter_toSend[_WEST_S_F_][_ENDX_] =  iter_toSend[_WEST_S_][_ENDX_];
                iter_toSend[_WEST_S_F_][_INIY_] =  iter_toSend[_WEST_S_][_INIY_];
                iter_toSend[_WEST_S_F_][_ENDY_] =  iter_toSend[_WEST_S_][_ENDY_];
                iter_toSend[_WEST_S_F_][_INIZ_] =  lNz-2;
                iter_toSend[_WEST_S_F_][_ENDZ_] =  lNz-2;
            }
            if(neighb[_WEST_F_] != _NO_NEIGHBOUR_ && neighb[_SOUTH_] == _NO_NEIGHBOUR_ ){ 

                iter_toSend[_WEST_S_F_][_INIX_] =  iter_toSend[_WEST_F_][_INIX_];// <= this
                iter_toSend[_WEST_S_F_][_ENDX_] =  iter_toSend[_WEST_F_][_ENDX_];
                iter_toSend[_WEST_S_F_][_INIY_] =  0;
                iter_toSend[_WEST_S_F_][_ENDY_] =  0;
                iter_toSend[_WEST_S_F_][_INIZ_] =  iter_toSend[_WEST_F_][_INIZ_];
                iter_toSend[_WEST_S_F_][_ENDZ_] =  iter_toSend[_WEST_F_][_ENDZ_];
            }
            else if(neighb[_WEST_F_] != _NO_NEIGHBOUR_ && neighb[_SOUTH_] != _NO_NEIGHBOUR_ ){

                iter_toSend[_WEST_S_F_][_INIX_] =  iter_toSend[_WEST_F_][_INIX_];// <= this
                iter_toSend[_WEST_S_F_][_ENDX_] =  iter_toSend[_WEST_F_][_ENDX_];
                iter_toSend[_WEST_S_F_][_INIY_] =  1;
                iter_toSend[_WEST_S_F_][_ENDY_] =  1;
                iter_toSend[_WEST_S_F_][_INIZ_] =  iter_toSend[_WEST_F_][_INIZ_];
                iter_toSend[_WEST_S_F_][_ENDZ_] =  iter_toSend[_WEST_F_][_ENDZ_];
            }
            if(neighb[_SOUTH_F_] != _NO_NEIGHBOUR_ && neighb[_WEST_] == _NO_NEIGHBOUR_ ){ 

                iter_toSend[_WEST_S_F_][_INIX_] =  0;// <= this
                iter_toSend[_WEST_S_F_][_ENDX_] =  0;
                iter_toSend[_WEST_S_F_][_INIY_] =  iter_toSend[_SOUTH_F_][_INIY_];
                iter_toSend[_WEST_S_F_][_ENDY_] =  iter_toSend[_SOUTH_F_][_ENDY_];
                iter_toSend[_WEST_S_F_][_INIZ_] =  iter_toSend[_SOUTH_F_][_INIZ_];
                iter_toSend[_WEST_S_F_][_ENDZ_] =  iter_toSend[_SOUTH_F_][_ENDZ_];
            }
            else if(neighb[_SOUTH_F_] != _NO_NEIGHBOUR_ && neighb[_WEST_] != _NO_NEIGHBOUR_ ){

                iter_toSend[_WEST_S_F_][_INIX_] =  1;// <= this
                iter_toSend[_WEST_S_F_][_ENDX_] =  1;
                iter_toSend[_WEST_S_F_][_INIY_] =  iter_toSend[_SOUTH_F_][_INIY_];
                iter_toSend[_WEST_S_F_][_ENDY_] =  iter_toSend[_SOUTH_F_][_ENDY_];
                iter_toSend[_WEST_S_F_][_INIZ_] =  iter_toSend[_SOUTH_F_][_INIZ_];
                iter_toSend[_WEST_S_F_][_ENDZ_] =  iter_toSend[_SOUTH_F_][_ENDZ_];
            }

        }

    }

    //WEST NORTH FRONT
    if(neighb[_WEST_N_F_] != _NO_NEIGHBOUR_ ) {
         if(neighb[_WEST_] != _NO_NEIGHBOUR_ && neighb[_NORTH_] != _NO_NEIGHBOUR_ &&  neighb[_FRONT_] != _NO_NEIGHBOUR_ ){  

       iter_toSend[_WEST_N_F_][_INIX_] =   1;// <= this
        iter_toSend[_WEST_N_F_][_ENDX_] =   1;
        iter_toSend[_WEST_N_F_][_INIY_] = lNy-2;
        iter_toSend[_WEST_N_F_][_ENDY_] = lNy-2;
        iter_toSend[_WEST_N_F_][_INIZ_] = lNz-2;
        iter_toSend[_WEST_N_F_][_ENDZ_] = lNz-2;
       }
        else
        {
            if(neighb[_WEST_N_] != _NO_NEIGHBOUR_ && neighb[_FRONT_] == _NO_NEIGHBOUR_ ){ 

                iter_toSend[_WEST_N_F_][_INIX_] =  iter_toSend[_WEST_N_][_INIX_];// <= this
                iter_toSend[_WEST_N_F_][_ENDX_] =  iter_toSend[_WEST_N_][_ENDX_];
                iter_toSend[_WEST_N_F_][_INIY_] =  iter_toSend[_WEST_N_][_INIY_];
                iter_toSend[_WEST_N_F_][_ENDY_] =  iter_toSend[_WEST_N_][_ENDY_];
                iter_toSend[_WEST_N_F_][_INIZ_] =  lNz-1;
                iter_toSend[_WEST_N_F_][_ENDZ_] =  lNz-1;
            }
            else if(neighb[_WEST_N_] != _NO_NEIGHBOUR_ && neighb[_FRONT_] != _NO_NEIGHBOUR_ ){

                iter_toSend[_WEST_N_F_][_INIX_] =  iter_toSend[_WEST_N_][_INIX_];// <= this
                iter_toSend[_WEST_N_F_][_ENDX_] =  iter_toSend[_WEST_N_][_ENDX_];
                iter_toSend[_WEST_N_F_][_INIY_] =  iter_toSend[_WEST_N_][_INIY_];
                iter_toSend[_WEST_N_F_][_ENDY_] =  iter_toSend[_WEST_N_][_ENDY_];
                iter_toSend[_WEST_N_F_][_INIZ_] =  lNz-2;
                iter_toSend[_WEST_N_F_][_ENDZ_] =  lNz-2;
            }
            if(neighb[_WEST_F_] != _NO_NEIGHBOUR_ && neighb[_NORTH_] == _NO_NEIGHBOUR_ ){ 

                iter_toSend[_WEST_N_F_][_INIX_] =  iter_toSend[_WEST_F_][_INIX_];// <= this
                iter_toSend[_WEST_N_F_][_ENDX_] =  iter_toSend[_WEST_F_][_ENDX_];
                iter_toSend[_WEST_N_F_][_INIY_] =  lNy-1;
                iter_toSend[_WEST_N_F_][_ENDY_] =  lNy-1;
                iter_toSend[_WEST_N_F_][_INIZ_] =  iter_toSend[_WEST_F_][_INIZ_];
                iter_toSend[_WEST_N_F_][_ENDZ_] =  iter_toSend[_WEST_F_][_ENDZ_];
            }
            else if(neighb[_WEST_F_] != _NO_NEIGHBOUR_ && neighb[_NORTH_] != _NO_NEIGHBOUR_ ){

                iter_toSend[_WEST_N_F_][_INIX_] =  iter_toSend[_WEST_F_][_INIX_];// <= this
                iter_toSend[_WEST_N_F_][_ENDX_] =  iter_toSend[_WEST_F_][_ENDX_];
                iter_toSend[_WEST_N_F_][_INIY_] =  lNy-2;
                iter_toSend[_WEST_N_F_][_ENDY_] =  lNy-2;
                iter_toSend[_WEST_N_F_][_INIZ_] =  iter_toSend[_WEST_F_][_INIZ_];
                iter_toSend[_WEST_N_F_][_ENDZ_] =  iter_toSend[_WEST_F_][_ENDZ_];
            }
            if(neighb[_NORTH_F_] != _NO_NEIGHBOUR_ && neighb[_WEST_] == _NO_NEIGHBOUR_ ){ 

                iter_toSend[_WEST_N_F_][_INIX_] =  0;// <= this
                iter_toSend[_WEST_N_F_][_ENDX_] =  0;
                iter_toSend[_WEST_N_F_][_INIY_] =  iter_toSend[_NORTH_F_][_INIY_];
                iter_toSend[_WEST_N_F_][_ENDY_] =  iter_toSend[_NORTH_F_][_ENDY_];
                iter_toSend[_WEST_N_F_][_INIZ_] =  iter_toSend[_NORTH_F_][_INIZ_];
                iter_toSend[_WEST_N_F_][_ENDZ_] =  iter_toSend[_NORTH_F_][_ENDZ_];
            }
            else if(neighb[_NORTH_F_] != _NO_NEIGHBOUR_ && neighb[_WEST_] != _NO_NEIGHBOUR_ ){

                iter_toSend[_WEST_N_F_][_INIX_] =  1;// <= this
                iter_toSend[_WEST_N_F_][_ENDX_] =  1;
                iter_toSend[_WEST_N_F_][_INIY_] =  iter_toSend[_NORTH_F_][_INIY_];
                iter_toSend[_WEST_N_F_][_ENDY_] =  iter_toSend[_NORTH_F_][_ENDY_];
                iter_toSend[_WEST_N_F_][_INIZ_] =  iter_toSend[_NORTH_F_][_INIZ_];
                iter_toSend[_WEST_N_F_][_ENDZ_] =  iter_toSend[_NORTH_F_][_ENDZ_];
            }

        }


    }

    //EAST SOUTH BACK
    if(neighb[_EAST_S_B_] != _NO_NEIGHBOUR_ ) {
         if(neighb[_EAST_] != _NO_NEIGHBOUR_ && neighb[_SOUTH_] != _NO_NEIGHBOUR_ &&  neighb[_BACK_] != _NO_NEIGHBOUR_ ){  

       iter_toSend[_EAST_S_B_][_INIX_] = lNx-2; // <=this
        iter_toSend[_EAST_S_B_][_ENDX_] = lNx-2;
        iter_toSend[_EAST_S_B_][_INIY_] =   1;
        iter_toSend[_EAST_S_B_][_ENDY_] =   1;
        iter_toSend[_EAST_S_B_][_INIZ_] =   1;
        iter_toSend[_EAST_S_B_][_ENDZ_] =   1;
         }
        else
        {
            if(neighb[_EAST_S_] != _NO_NEIGHBOUR_ && neighb[_BACK_] == _NO_NEIGHBOUR_ ){ 

                iter_toSend[_EAST_S_B_][_INIX_] =  iter_toSend[_EAST_S_][_INIX_];// <= this
                iter_toSend[_EAST_S_B_][_ENDX_] =  iter_toSend[_EAST_S_][_ENDX_];
                iter_toSend[_EAST_S_B_][_INIY_] =  iter_toSend[_EAST_S_][_INIY_];
                iter_toSend[_EAST_S_B_][_ENDY_] =  iter_toSend[_EAST_S_][_ENDY_];
                iter_toSend[_EAST_S_B_][_INIZ_] =  0;
                iter_toSend[_EAST_S_B_][_ENDZ_] =  0;
            }
            else if(neighb[_EAST_S_] != _NO_NEIGHBOUR_ && neighb[_BACK_] != _NO_NEIGHBOUR_ ){
                
                iter_toSend[_EAST_S_B_][_INIX_] =  iter_toSend[_EAST_S_][_INIX_];// <= this
                iter_toSend[_EAST_S_B_][_ENDX_] =  iter_toSend[_EAST_S_][_ENDX_];
                iter_toSend[_EAST_S_B_][_INIY_] =  iter_toSend[_EAST_S_][_INIY_];
                iter_toSend[_EAST_S_B_][_ENDY_] =  iter_toSend[_EAST_S_][_ENDY_];
                iter_toSend[_EAST_S_B_][_INIZ_] =  1;
                iter_toSend[_EAST_S_B_][_ENDZ_] =  1;
            }
             if(neighb[_EAST_B_] != _NO_NEIGHBOUR_ && neighb[_SOUTH_] == _NO_NEIGHBOUR_ ){ 

                iter_toSend[_EAST_S_B_][_INIX_] =  iter_toSend[_EAST_B_][_INIX_];// <= this
                iter_toSend[_EAST_S_B_][_ENDX_] =  iter_toSend[_EAST_B_][_ENDX_];
                iter_toSend[_EAST_S_B_][_INIY_] =  0;
                iter_toSend[_EAST_S_B_][_ENDY_] =  0;
                iter_toSend[_EAST_S_B_][_INIZ_] =  iter_toSend[_EAST_B_][_INIZ_];
                iter_toSend[_EAST_S_B_][_ENDZ_] =  iter_toSend[_EAST_B_][_ENDZ_];
            }
            else if(neighb[_EAST_B_] != _NO_NEIGHBOUR_ && neighb[_SOUTH_] != _NO_NEIGHBOUR_ ){
                
                iter_toSend[_EAST_S_B_][_INIX_] =  iter_toSend[_EAST_B_][_INIX_];// <= this
                iter_toSend[_EAST_S_B_][_ENDX_] =  iter_toSend[_EAST_B_][_ENDX_];
                iter_toSend[_EAST_S_B_][_INIY_] =  1;
                iter_toSend[_EAST_S_B_][_ENDY_] =  1;
                iter_toSend[_EAST_S_B_][_INIZ_] =  iter_toSend[_EAST_B_][_INIZ_];
                iter_toSend[_EAST_S_B_][_ENDZ_] =  iter_toSend[_EAST_B_][_ENDZ_];
            }
             if(neighb[_SOUTH_B_] != _NO_NEIGHBOUR_ && neighb[_EAST_] == _NO_NEIGHBOUR_ ){ 

                iter_toSend[_EAST_S_B_][_INIX_] =  lNx-1;// <= this
                iter_toSend[_EAST_S_B_][_ENDX_] =  lNx-1;
                iter_toSend[_EAST_S_B_][_INIY_] =  iter_toSend[_SOUTH_B_][_INIY_];
                iter_toSend[_EAST_S_B_][_ENDY_] =  iter_toSend[_SOUTH_B_][_ENDY_];
                iter_toSend[_EAST_S_B_][_INIZ_] =  iter_toSend[_SOUTH_B_][_INIZ_];
                iter_toSend[_EAST_S_B_][_ENDZ_] =  iter_toSend[_SOUTH_B_][_ENDZ_];
            }
            else if(neighb[_SOUTH_B_] != _NO_NEIGHBOUR_ && neighb[_EAST_] != _NO_NEIGHBOUR_ ){
                
                iter_toSend[_EAST_S_B_][_INIX_] =  lNx-2;// <= this
                iter_toSend[_EAST_S_B_][_ENDX_] =  lNx-2;
                iter_toSend[_EAST_S_B_][_INIY_] =  iter_toSend[_SOUTH_B_][_INIY_];
                iter_toSend[_EAST_S_B_][_ENDY_] =  iter_toSend[_SOUTH_B_][_ENDY_];
                iter_toSend[_EAST_S_B_][_INIZ_] =  iter_toSend[_SOUTH_B_][_INIZ_];
                iter_toSend[_EAST_S_B_][_ENDZ_] =  iter_toSend[_SOUTH_B_][_ENDZ_];
            }

        }

    }

    //EAST NORTH BACK
    if(neighb[_EAST_N_B_] != _NO_NEIGHBOUR_ ) {
          if(neighb[_EAST_] != _NO_NEIGHBOUR_ && neighb[_NORTH_] != _NO_NEIGHBOUR_ &&  neighb[_BACK_] != _NO_NEIGHBOUR_ ){  

       iter_toSend[_EAST_N_B_][_INIX_] = lNx-2; // <=this
        iter_toSend[_EAST_N_B_][_ENDX_] = lNx-2;
        iter_toSend[_EAST_N_B_][_INIY_] = lNy-2;
        iter_toSend[_EAST_N_B_][_ENDY_] = lNy-2;
        iter_toSend[_EAST_N_B_][_INIZ_] =   1;
        iter_toSend[_EAST_N_B_][_ENDZ_] =   1;
          }
        else
        {
            if(neighb[_EAST_N_] != _NO_NEIGHBOUR_ && neighb[_BACK_] == _NO_NEIGHBOUR_ ){ 

                iter_toSend[_EAST_N_B_][_INIX_] =  iter_toSend[_EAST_N_][_INIX_];// <= this
                iter_toSend[_EAST_N_B_][_ENDX_] =  iter_toSend[_EAST_N_][_ENDX_];
                iter_toSend[_EAST_N_B_][_INIY_] =  iter_toSend[_EAST_N_][_INIY_];
                iter_toSend[_EAST_N_B_][_ENDY_] =  iter_toSend[_EAST_N_][_ENDY_];
                iter_toSend[_EAST_N_B_][_INIZ_] =  0;
                iter_toSend[_EAST_N_B_][_ENDZ_] =  0;
            }
            else if(neighb[_EAST_N_] != _NO_NEIGHBOUR_ && neighb[_BACK_] != _NO_NEIGHBOUR_ ){
                
                iter_toSend[_EAST_N_B_][_INIX_] =  iter_toSend[_EAST_N_][_INIX_];// <= this
                iter_toSend[_EAST_N_B_][_ENDX_] =  iter_toSend[_EAST_N_][_ENDX_];
                iter_toSend[_EAST_N_B_][_INIY_] =  iter_toSend[_EAST_N_][_INIY_];
                iter_toSend[_EAST_N_B_][_ENDY_] =  iter_toSend[_EAST_N_][_ENDY_];
                iter_toSend[_EAST_N_B_][_INIZ_] =  1;
                iter_toSend[_EAST_N_B_][_ENDZ_] =  1;
            }
             if(neighb[_EAST_B_] != _NO_NEIGHBOUR_ && neighb[_NORTH_] == _NO_NEIGHBOUR_ ){ 

                iter_toSend[_EAST_N_B_][_INIX_] =  iter_toSend[_EAST_B_][_INIX_];// <= this
                iter_toSend[_EAST_N_B_][_ENDX_] =  iter_toSend[_EAST_B_][_ENDX_];
                iter_toSend[_EAST_N_B_][_INIY_] =  lNy-1;
                iter_toSend[_EAST_N_B_][_ENDY_] =  lNy-1;
                iter_toSend[_EAST_N_B_][_INIZ_] =  iter_toSend[_EAST_B_][_INIZ_];
                iter_toSend[_EAST_N_B_][_ENDZ_] =  iter_toSend[_EAST_B_][_ENDZ_];
            }
            else if(neighb[_EAST_B_] != _NO_NEIGHBOUR_ && neighb[_NORTH_] != _NO_NEIGHBOUR_ ){
                
                iter_toSend[_EAST_N_B_][_INIX_] =  iter_toSend[_EAST_B_][_INIX_];// <= this
                iter_toSend[_EAST_N_B_][_ENDX_] =  iter_toSend[_EAST_B_][_ENDX_];
                iter_toSend[_EAST_N_B_][_INIY_] =  lNy-2;
                iter_toSend[_EAST_N_B_][_ENDY_] =  lNy-2;
                iter_toSend[_EAST_N_B_][_INIZ_] =  iter_toSend[_EAST_B_][_INIZ_];
                iter_toSend[_EAST_N_B_][_ENDZ_] =  iter_toSend[_EAST_B_][_ENDZ_];
            }
             if(neighb[_NORTH_B_] != _NO_NEIGHBOUR_ && neighb[_EAST_] == _NO_NEIGHBOUR_ ){ 

                iter_toSend[_EAST_N_B_][_INIX_] =  lNx-1;// <= this
                iter_toSend[_EAST_N_B_][_ENDX_] =  lNx-1;
                iter_toSend[_EAST_N_B_][_INIY_] =  iter_toSend[_NORTH_B_][_INIY_];
                iter_toSend[_EAST_N_B_][_ENDY_] =  iter_toSend[_NORTH_B_][_ENDY_];
                iter_toSend[_EAST_N_B_][_INIZ_] =  iter_toSend[_NORTH_B_][_INIZ_];
                iter_toSend[_EAST_N_B_][_ENDZ_] =  iter_toSend[_NORTH_B_][_ENDZ_];
            }
            else if(neighb[_NORTH_B_] != _NO_NEIGHBOUR_ && neighb[_EAST_] != _NO_NEIGHBOUR_ ){
                
                iter_toSend[_EAST_N_B_][_INIX_] =  lNx-2;// <= this
                iter_toSend[_EAST_N_B_][_ENDX_] =  lNx-2;
                iter_toSend[_EAST_N_B_][_INIY_] =  iter_toSend[_NORTH_B_][_INIY_];
                iter_toSend[_EAST_N_B_][_ENDY_] =  iter_toSend[_NORTH_B_][_ENDY_];
                iter_toSend[_EAST_N_B_][_INIZ_] =  iter_toSend[_NORTH_B_][_INIZ_];
                iter_toSend[_EAST_N_B_][_ENDZ_] =  iter_toSend[_NORTH_B_][_ENDZ_];
            }

        }

    }

    //EAST SOUTH FRONT
    if(neighb[_EAST_S_F_] != _NO_NEIGHBOUR_ ) {
         if(neighb[_EAST_] != _NO_NEIGHBOUR_ && neighb[_SOUTH_] != _NO_NEIGHBOUR_ &&  neighb[_FRONT_] != _NO_NEIGHBOUR_ ){  

       iter_toSend[_EAST_S_F_][_INIX_] = lNx-2; // <=this
        iter_toSend[_EAST_S_F_][_ENDX_] = lNx-2;
        iter_toSend[_EAST_S_F_][_INIY_] =   1;
        iter_toSend[_EAST_S_F_][_ENDY_] =   1;
        iter_toSend[_EAST_S_F_][_INIZ_] = lNz-2;
        iter_toSend[_EAST_S_F_][_ENDZ_] = lNz-2;
         }
        else
        {
            if(neighb[_EAST_S_] != _NO_NEIGHBOUR_ && neighb[_FRONT_] == _NO_NEIGHBOUR_ ){ 

                iter_toSend[_EAST_S_F_][_INIX_] =  iter_toSend[_EAST_S_][_INIX_];// <= this
                iter_toSend[_EAST_S_F_][_ENDX_] =  iter_toSend[_EAST_S_][_ENDX_];
                iter_toSend[_EAST_S_F_][_INIY_] =  iter_toSend[_EAST_S_][_INIY_];
                iter_toSend[_EAST_S_F_][_ENDY_] =  iter_toSend[_EAST_S_][_ENDY_];
                iter_toSend[_EAST_S_F_][_INIZ_] =  lNz-1;
                iter_toSend[_EAST_S_F_][_ENDZ_] =  lNz-1;
            }
            else if(neighb[_EAST_S_] != _NO_NEIGHBOUR_ && neighb[_FRONT_] != _NO_NEIGHBOUR_ ){

                iter_toSend[_EAST_S_F_][_INIX_] =  iter_toSend[_EAST_S_][_INIX_];// <= this
                iter_toSend[_EAST_S_F_][_ENDX_] =  iter_toSend[_EAST_S_][_ENDX_];
                iter_toSend[_EAST_S_F_][_INIY_] =  iter_toSend[_EAST_S_][_INIY_];
                iter_toSend[_EAST_S_F_][_ENDY_] =  iter_toSend[_EAST_S_][_ENDY_];
                iter_toSend[_EAST_S_F_][_INIZ_] =  lNz-2;
                iter_toSend[_EAST_S_F_][_ENDZ_] =  lNz-2;
            }
            if(neighb[_EAST_F_] != _NO_NEIGHBOUR_ && neighb[_SOUTH_] == _NO_NEIGHBOUR_ ){ 

                iter_toSend[_EAST_S_F_][_INIX_] =  iter_toSend[_EAST_F_][_INIX_];// <= this
                iter_toSend[_EAST_S_F_][_ENDX_] =  iter_toSend[_EAST_F_][_ENDX_];
                iter_toSend[_EAST_S_F_][_INIY_] =  0;
                iter_toSend[_EAST_S_F_][_ENDY_] =  0;
                iter_toSend[_EAST_S_F_][_INIZ_] =  iter_toSend[_EAST_F_][_INIZ_];
                iter_toSend[_EAST_S_F_][_ENDZ_] =  iter_toSend[_EAST_F_][_ENDZ_];
            }
            else if(neighb[_EAST_F_] != _NO_NEIGHBOUR_ && neighb[_SOUTH_] != _NO_NEIGHBOUR_ ){

                iter_toSend[_EAST_S_F_][_INIX_] =  iter_toSend[_EAST_F_][_INIX_];// <= this
                iter_toSend[_EAST_S_F_][_ENDX_] =  iter_toSend[_EAST_F_][_ENDX_];
                iter_toSend[_EAST_S_F_][_INIY_] =  1;
                iter_toSend[_EAST_S_F_][_ENDY_] =  1;
                iter_toSend[_EAST_S_F_][_INIZ_] =  iter_toSend[_EAST_F_][_INIZ_];
                iter_toSend[_EAST_S_F_][_ENDZ_] =  iter_toSend[_EAST_F_][_ENDZ_];
            }
            if(neighb[_SOUTH_F_] != _NO_NEIGHBOUR_ && neighb[_EAST_] == _NO_NEIGHBOUR_ ){ 

                iter_toSend[_EAST_S_F_][_INIX_] =  lNx-1;// <= this
                iter_toSend[_EAST_S_F_][_ENDX_] =  lNx-1;
                iter_toSend[_EAST_S_F_][_INIY_] =  iter_toSend[_SOUTH_F_][_INIY_];
                iter_toSend[_EAST_S_F_][_ENDY_] =  iter_toSend[_SOUTH_F_][_ENDY_];
                iter_toSend[_EAST_S_F_][_INIZ_] =  iter_toSend[_SOUTH_F_][_INIZ_];
                iter_toSend[_EAST_S_F_][_ENDZ_] =  iter_toSend[_SOUTH_F_][_ENDZ_];
            }
            else if(neighb[_SOUTH_F_] != _NO_NEIGHBOUR_ && neighb[_EAST_] != _NO_NEIGHBOUR_ ){

                iter_toSend[_EAST_S_F_][_INIX_] =  lNx-2;// <= this
                iter_toSend[_EAST_S_F_][_ENDX_] =  lNx-2;
                iter_toSend[_EAST_S_F_][_INIY_] =  iter_toSend[_SOUTH_F_][_INIY_];
                iter_toSend[_EAST_S_F_][_ENDY_] =  iter_toSend[_SOUTH_F_][_ENDY_];
                iter_toSend[_EAST_S_F_][_INIZ_] =  iter_toSend[_SOUTH_F_][_INIZ_];
                iter_toSend[_EAST_S_F_][_ENDZ_] =  iter_toSend[_SOUTH_F_][_ENDZ_];
            }

        }

    }

    //EAST NORTH FRONT
    if(neighb[_EAST_N_F_] != _NO_NEIGHBOUR_ ) {
          if(neighb[_WEST_] != _NO_NEIGHBOUR_ && neighb[_NORTH_] != _NO_NEIGHBOUR_ &&  neighb[_FRONT_] != _NO_NEIGHBOUR_ ){  

       iter_toSend[_EAST_N_F_][_INIX_] = lNx-2; // <=this
        iter_toSend[_EAST_N_F_][_ENDX_] = lNx-2;
        iter_toSend[_EAST_N_F_][_INIY_] = lNy-2;
        iter_toSend[_EAST_N_F_][_ENDY_] = lNy-2;
        iter_toSend[_EAST_N_F_][_INIZ_] = lNz-2;
        iter_toSend[_EAST_N_F_][_ENDZ_] = lNz-2;
          }

        else
        {
            if(neighb[_EAST_N_] != _NO_NEIGHBOUR_ && neighb[_FRONT_] == _NO_NEIGHBOUR_ ){ 

                iter_toSend[_EAST_N_F_][_INIX_] =  iter_toSend[_EAST_N_][_INIX_];// <= this
                iter_toSend[_EAST_N_F_][_ENDX_] =  iter_toSend[_EAST_N_][_ENDX_];
                iter_toSend[_EAST_N_F_][_INIY_] =  iter_toSend[_EAST_N_][_INIY_];
                iter_toSend[_EAST_N_F_][_ENDY_] =  iter_toSend[_EAST_N_][_ENDY_];
                iter_toSend[_EAST_N_F_][_INIZ_] =  lNz-1;
                iter_toSend[_EAST_N_F_][_ENDZ_] =  lNz-1;
            }
            else if(neighb[_EAST_N_] != _NO_NEIGHBOUR_ && neighb[_FRONT_] != _NO_NEIGHBOUR_ ){

                iter_toSend[_EAST_N_F_][_INIX_] =  iter_toSend[_EAST_N_][_INIX_];// <= this
                iter_toSend[_EAST_N_F_][_ENDX_] =  iter_toSend[_EAST_N_][_ENDX_];
                iter_toSend[_EAST_N_F_][_INIY_] =  iter_toSend[_EAST_N_][_INIY_];
                iter_toSend[_EAST_N_F_][_ENDY_] =  iter_toSend[_EAST_N_][_ENDY_];
                iter_toSend[_EAST_N_F_][_INIZ_] =  lNz-2;
                iter_toSend[_EAST_N_F_][_ENDZ_] =  lNz-2;
            }
            if(neighb[_EAST_F_] != _NO_NEIGHBOUR_ && neighb[_NORTH_] == _NO_NEIGHBOUR_ ){ 

                iter_toSend[_EAST_N_F_][_INIX_] =  iter_toSend[_EAST_F_][_INIX_];// <= this
                iter_toSend[_EAST_N_F_][_ENDX_] =  iter_toSend[_EAST_F_][_ENDX_];
                iter_toSend[_EAST_N_F_][_INIY_] =  lNy-1;
                iter_toSend[_EAST_N_F_][_ENDY_] =  lNy-1;
                iter_toSend[_EAST_N_F_][_INIZ_] =  iter_toSend[_EAST_F_][_INIZ_];
                iter_toSend[_EAST_N_F_][_ENDZ_] =  iter_toSend[_EAST_F_][_ENDZ_];
            }
            else if(neighb[_EAST_F_] != _NO_NEIGHBOUR_ && neighb[_NORTH_] != _NO_NEIGHBOUR_ ){

                iter_toSend[_EAST_N_F_][_INIX_] =  iter_toSend[_EAST_F_][_INIX_];// <= this
                iter_toSend[_EAST_N_F_][_ENDX_] =  iter_toSend[_EAST_F_][_ENDX_];
                iter_toSend[_EAST_N_F_][_INIY_] =  lNy-2;
                iter_toSend[_EAST_N_F_][_ENDY_] =  lNy-2;
                iter_toSend[_EAST_N_F_][_INIZ_] =  iter_toSend[_EAST_F_][_INIZ_];
                iter_toSend[_EAST_N_F_][_ENDZ_] =  iter_toSend[_EAST_F_][_ENDZ_];
            }
            if(neighb[_NORTH_F_] != _NO_NEIGHBOUR_ && neighb[_EAST_] == _NO_NEIGHBOUR_ ){ 

                iter_toSend[_EAST_N_F_][_INIX_] =  lNx-1;// <= this
                iter_toSend[_EAST_N_F_][_ENDX_] =  lNx-1;
                iter_toSend[_EAST_N_F_][_INIY_] =  iter_toSend[_NORTH_F_][_INIY_];
                iter_toSend[_EAST_N_F_][_ENDY_] =  iter_toSend[_NORTH_F_][_ENDY_];
                iter_toSend[_EAST_N_F_][_INIZ_] =  iter_toSend[_NORTH_F_][_INIZ_];
                iter_toSend[_EAST_N_F_][_ENDZ_] =  iter_toSend[_NORTH_F_][_ENDZ_];
            }
            else if(neighb[_NORTH_F_] != _NO_NEIGHBOUR_ && neighb[_EAST_] != _NO_NEIGHBOUR_ ){

                iter_toSend[_EAST_N_F_][_INIX_] =  lNx-2;// <= this
                iter_toSend[_EAST_N_F_][_ENDX_] =  lNx-2;
                iter_toSend[_EAST_N_F_][_INIY_] =  iter_toSend[_NORTH_F_][_INIY_];
                iter_toSend[_EAST_N_F_][_ENDY_] =  iter_toSend[_NORTH_F_][_ENDY_];
                iter_toSend[_EAST_N_F_][_INIZ_] =  iter_toSend[_NORTH_F_][_INIZ_];
                iter_toSend[_EAST_N_F_][_ENDZ_] =  iter_toSend[_NORTH_F_][_ENDZ_];
            }

        }


    }

}


void ParallelTopology::create_global_iters()
{
    
    iter_glob_ind[_INIX_] = iter_common[_INNER_][_INIX_] + offx - 1;
    iter_glob_ind[_ENDX_] = iter_common[_INNER_][_ENDX_] + offx + 1;
    iter_glob_ind[_INIY_] = iter_common[_INNER_][_INIY_] + offy - 1;
    iter_glob_ind[_ENDY_] = iter_common[_INNER_][_ENDY_] + offy + 1;
    iter_glob_ind[_INIZ_] = iter_common[_INNER_][_INIZ_] + offz - 1;
    iter_glob_ind[_ENDZ_] = iter_common[_INNER_][_ENDZ_] + offz + 1;



    if(iter_glob_ind[_INIX_] == 0)
        iter_slab[_INIX_] = 0;
    else
        iter_slab[_INIX_] = 1;


    if(iter_glob_ind[_ENDX_] ==  (mymesh->getGNx() + 1) )
      iter_slab[_ENDX_] = lNx-1;
    else
      iter_slab[_ENDX_] = lNx-2;
 
    if(iter_glob_ind[_INIY_] == 0)
        iter_slab[_INIY_] = 0;
    else 
        iter_slab[_INIY_] = 1;
 
    if(iter_glob_ind[_ENDY_] ==  (mymesh->getGNy() + 1) )
        iter_slab[_ENDY_] = lNy-1;
    else
        iter_slab[_ENDY_] = lNy-2;


    if(iter_glob_ind[_INIZ_] == 0)
        iter_slab[_INIZ_] = 0;
    else
        iter_slab[_INIZ_] = 1;

    if(iter_glob_ind[_ENDZ_] ==  (mymesh->getGNz() + 1) )
        iter_slab[_ENDZ_] = lNz-1;
    else
        iter_slab[_ENDZ_] = lNz-2;

    //
    if(iter_glob_ind[_INIX_] == 0)
        offslab_x = offx;
    else
        offslab_x = offx + 1;

    if(iter_glob_ind[_INIY_] == 0)
        offslab_y = offy;
    else
        offslab_y = offy + 1;


    if(iter_glob_ind[_INIZ_] == 0)
        offslab_z = offz;
    else
        offslab_z = offz + 1;


    lenslabx = iter_slab[_ENDX_] -  iter_slab[_INIX_] + 1; 
    lenslaby = iter_slab[_ENDY_] -  iter_slab[_INIY_] + 1; 
    lenslabz = iter_slab[_ENDZ_] -  iter_slab[_INIZ_] + 1; 

    lenslab = lenslabx*lenslaby*lenslabz;




    //cout<<" Globals "<<rank<<" offz "<<offz<<endl;
    //cout<<" "<<iter_glob_ind[_INIX_] <<endl;
    //cout<<" "<<iter_glob_ind[_ENDX_] <<endl;
    //cout<<" "<<iter_glob_ind[_INIY_] <<endl;
    //cout<<" "<<iter_glob_ind[_ENDY_] <<endl;
    //cout<<" "<<iter_glob_ind[_INIZ_] <<endl;
    //cout<<" "<<iter_glob_ind[_ENDZ_] <<endl;





}

void ParallelTopology::create_comm_arrays()
{
    //Size of the faces

    len_xy = (lNx-2)*(lNy-2);
    len_xz = (lNx-2)*(lNz-2);
    len_yz = (lNy-2)*(lNz-2);

    if( len_yz != 0 ) {
        pack_send_w = new double [len_yz];    
        pack_send_e = new double [len_yz];    
    }
    if( len_xz != 0 ) {
        pack_send_s = new double [len_xz];    
        pack_send_n = new double [len_xz];    
    }
    if( len_xy != 0 ) {
        pack_send_b = new double [len_xy];    
        pack_send_f = new double [len_xy];    
    }

    pack_recv_w = new double [len_yz];    
    pack_recv_e = new double [len_yz];    
    pack_recv_s = new double [len_xz];    
    pack_recv_n = new double [len_xz];    
    pack_recv_b = new double [len_xy];    
    pack_recv_f = new double [len_xy];    
 
}

void ParallelTopology::create_complex_comm_arrays()
{

    len_1Dx = lNx-2;
    len_1Dy = lNy-2;
    len_1Dz = lNz-2;

    len_1pt = 1;

// pack send

    pack_send_ws = new double [len_1Dz];
    pack_send_wn = new double [len_1Dz];
    pack_send_wb = new double [len_1Dy];
    pack_send_wf = new double [len_1Dy];

    pack_send_es = new double [len_1Dz];
    pack_send_en = new double [len_1Dz];
    pack_send_eb = new double [len_1Dy];
    pack_send_ef = new double [len_1Dy];

    pack_send_sb = new double [len_1Dx];
    pack_send_sf = new double [len_1Dx];
    pack_send_nb = new double [len_1Dx];
    pack_send_nf = new double [len_1Dx];

    pack_send_wsb = new double [len_1pt];
    pack_send_wnb = new double [len_1pt];
    pack_send_wsf = new double [len_1pt];
    pack_send_wnf = new double [len_1pt];

    pack_send_esb = new double [len_1pt];
    pack_send_enb = new double [len_1pt];
    pack_send_esf = new double [len_1pt];
    pack_send_enf = new double [len_1pt];


// pack recv

    pack_recv_ws = new double [len_1Dz];
    pack_recv_wn = new double [len_1Dz];
    pack_recv_wb = new double [len_1Dy];
    pack_recv_wf = new double [len_1Dy];

    pack_recv_es = new double [len_1Dz];
    pack_recv_en = new double [len_1Dz];
    pack_recv_eb = new double [len_1Dy];
    pack_recv_ef = new double [len_1Dy];

    pack_recv_sb = new double [len_1Dx];
    pack_recv_sf = new double [len_1Dx];
    pack_recv_nb = new double [len_1Dx];
    pack_recv_nf = new double [len_1Dx];

    pack_recv_wsb = new double [len_1pt];
    pack_recv_wnb = new double [len_1pt];
    pack_recv_wsf = new double [len_1pt];
    pack_recv_wnf = new double [len_1pt];

    pack_recv_esb = new double [len_1pt];
    pack_recv_enb = new double [len_1pt];
    pack_recv_esf = new double [len_1pt];
    pack_recv_enf = new double [len_1pt];



}
void ParallelTopology::update(double *vec)
{
    pack(vec);
    halo_exchange();
    unpack(vec);
}

void ParallelTopology::update_simple(double *vec)
{
    pack_simple(vec);
    halo_exchange_simple();
    unpack_simple(vec);
}


void ParallelTopology::pack_simple(double *vec)
{

    //WEST
    int l=0;
    for(int i=iter_toSend[_WEST_][_INIX_]; i<= iter_toSend[_WEST_][_ENDX_]; i++)
        for(int j=iter_toSend[_WEST_][_INIY_]; j<= iter_toSend[_WEST_][_ENDY_]; j++)
            for(int k=iter_toSend[_WEST_][_INIZ_]; k<= iter_toSend[_WEST_][_ENDZ_]; k++)
            {
                pack_send_w[l] = vec[ R_INDX( i , j, k, lNx, lNy)];
                l++; 
            }


    //EAST
    l=0;
    for(int i=iter_toSend[_EAST_][_INIX_]; i<= iter_toSend[_EAST_][_ENDX_]; i++)
        for(int j=iter_toSend[_EAST_][_INIY_]; j<= iter_toSend[_EAST_][_ENDY_]; j++)
            for(int k=iter_toSend[_EAST_][_INIZ_]; k<= iter_toSend[_EAST_][_ENDZ_]; k++)
            {
                pack_send_e[l] = vec[ R_INDX( i , j, k, lNx, lNy)];
                l++;
            }


    //SOUTH
    l=0;
    for(int i=iter_toSend[_SOUTH_][_INIX_]; i<= iter_toSend[_SOUTH_][_ENDX_]; i++)
        for(int j=iter_toSend[_SOUTH_][_INIY_]; j<= iter_toSend[_SOUTH_][_ENDY_]; j++)
            for(int k=iter_toSend[_SOUTH_][_INIZ_]; k<= iter_toSend[_SOUTH_][_ENDZ_]; k++)
            {
                pack_send_s[l] = vec[ R_INDX( i , j, k, lNx, lNy)];
                l++;
            }


    //NORTH
    l=0;
    for(int i=iter_toSend[_NORTH_][_INIX_]; i<= iter_toSend[_NORTH_][_ENDX_]; i++)
        for(int j=iter_toSend[_NORTH_][_INIY_]; j<= iter_toSend[_NORTH_][_ENDY_]; j++)
            for(int k=iter_toSend[_NORTH_][_INIZ_]; k<= iter_toSend[_NORTH_][_ENDZ_]; k++)
            {
                pack_send_n[l] = vec[ R_INDX( i , j, k, lNx, lNy)];
                l++;
            }

    //BACK
    l=0;
    for(int i=iter_toSend[_BACK_][_INIX_]; i<= iter_toSend[_BACK_][_ENDX_]; i++)
        for(int j=iter_toSend[_BACK_][_INIY_]; j<= iter_toSend[_BACK_][_ENDY_]; j++)
            for(int k=iter_toSend[_BACK_][_INIZ_]; k<= iter_toSend[_BACK_][_ENDZ_]; k++)
            {
                pack_send_b[l] = vec[ R_INDX( i , j, k, lNx, lNy)];
                l++;
            }


    //FRONT 
    l=0;
    for(int i=iter_toSend[_FRONT_][_INIX_]; i<= iter_toSend[_FRONT_][_ENDX_]; i++)
        for(int j=iter_toSend[_FRONT_][_INIY_]; j<= iter_toSend[_FRONT_][_ENDY_]; j++)
            for(int k=iter_toSend[_FRONT_][_INIZ_]; k<= iter_toSend[_FRONT_][_ENDZ_]; k++)
            {
                pack_send_f[l] = vec[ R_INDX( i , j, k, lNx, lNy)];
                l++;
            }

}


void ParallelTopology::pack(double *vec)
{

    //WEST
    int l=0;
    for(int i=iter_toSend[_WEST_][_INIX_]; i<= iter_toSend[_WEST_][_ENDX_]; i++)
        for(int j=iter_toSend[_WEST_][_INIY_]; j<= iter_toSend[_WEST_][_ENDY_]; j++)
            for(int k=iter_toSend[_WEST_][_INIZ_]; k<= iter_toSend[_WEST_][_ENDZ_]; k++)
            {
                pack_send_w[l] = vec[ R_INDX( i , j, k, lNx, lNy)];
                l++; 
            }


    //EAST
    l=0;
    for(int i=iter_toSend[_EAST_][_INIX_]; i<= iter_toSend[_EAST_][_ENDX_]; i++)
        for(int j=iter_toSend[_EAST_][_INIY_]; j<= iter_toSend[_EAST_][_ENDY_]; j++)
            for(int k=iter_toSend[_EAST_][_INIZ_]; k<= iter_toSend[_EAST_][_ENDZ_]; k++)
            {
                pack_send_e[l] = vec[ R_INDX( i , j, k, lNx, lNy)];
                l++;
            }


    //SOUTH
    l=0;
    for(int i=iter_toSend[_SOUTH_][_INIX_]; i<= iter_toSend[_SOUTH_][_ENDX_]; i++)
        for(int j=iter_toSend[_SOUTH_][_INIY_]; j<= iter_toSend[_SOUTH_][_ENDY_]; j++)
            for(int k=iter_toSend[_SOUTH_][_INIZ_]; k<= iter_toSend[_SOUTH_][_ENDZ_]; k++)
            {
                pack_send_s[l] = vec[ R_INDX( i , j, k, lNx, lNy)];
                l++;
            }


    //NORTH
    l=0;
    for(int i=iter_toSend[_NORTH_][_INIX_]; i<= iter_toSend[_NORTH_][_ENDX_]; i++)
        for(int j=iter_toSend[_NORTH_][_INIY_]; j<= iter_toSend[_NORTH_][_ENDY_]; j++)
            for(int k=iter_toSend[_NORTH_][_INIZ_]; k<= iter_toSend[_NORTH_][_ENDZ_]; k++)
            {
                pack_send_n[l] = vec[ R_INDX( i , j, k, lNx, lNy)];
                l++;
            }

    //BACK
    l=0;
    for(int i=iter_toSend[_BACK_][_INIX_]; i<= iter_toSend[_BACK_][_ENDX_]; i++)
        for(int j=iter_toSend[_BACK_][_INIY_]; j<= iter_toSend[_BACK_][_ENDY_]; j++)
            for(int k=iter_toSend[_BACK_][_INIZ_]; k<= iter_toSend[_BACK_][_ENDZ_]; k++)
            {
                pack_send_b[l] = vec[ R_INDX( i , j, k, lNx, lNy)];
                l++;
            }


    //FRONT 
    l=0;
    for(int i=iter_toSend[_FRONT_][_INIX_]; i<= iter_toSend[_FRONT_][_ENDX_]; i++)
        for(int j=iter_toSend[_FRONT_][_INIY_]; j<= iter_toSend[_FRONT_][_ENDY_]; j++)
            for(int k=iter_toSend[_FRONT_][_INIZ_]; k<= iter_toSend[_FRONT_][_ENDZ_]; k++)
            {
                pack_send_f[l] = vec[ R_INDX( i , j, k, lNx, lNy)];
                l++;
            }

 
    //WEST SOUTH 
    l=0;
    for(int i=iter_toSend[_WEST_S_][_INIX_]; i<= iter_toSend[_WEST_S_][_ENDX_]; i++)
        for(int j=iter_toSend[_WEST_S_][_INIY_]; j<= iter_toSend[_WEST_S_][_ENDY_]; j++)
            for(int k=iter_toSend[_WEST_S_][_INIZ_]; k<= iter_toSend[_WEST_S_][_ENDZ_]; k++)
            {
                pack_send_ws[l] = vec[ R_INDX( i , j, k, lNx, lNy)];
                l++; 
            }

    //WEST NORTH
    l=0;
    for(int i=iter_toSend[_WEST_N_][_INIX_]; i<= iter_toSend[_WEST_N_][_ENDX_]; i++)
        for(int j=iter_toSend[_WEST_N_][_INIY_]; j<= iter_toSend[_WEST_N_][_ENDY_]; j++)
            for(int k=iter_toSend[_WEST_N_][_INIZ_]; k<= iter_toSend[_WEST_N_][_ENDZ_]; k++)
            {
                pack_send_wn[l] = vec[ R_INDX( i , j, k, lNx, lNy)];
                l++; 
            }

    //WEST BACK 
    l=0;
    for(int i=iter_toSend[_WEST_B_][_INIX_]; i<= iter_toSend[_WEST_B_][_ENDX_]; i++)
        for(int j=iter_toSend[_WEST_B_][_INIY_]; j<= iter_toSend[_WEST_B_][_ENDY_]; j++)
            for(int k=iter_toSend[_WEST_B_][_INIZ_]; k<= iter_toSend[_WEST_B_][_ENDZ_]; k++)
            {
                pack_send_wb[l] = vec[ R_INDX( i , j, k, lNx, lNy)];
                l++; 
            }

    //WEST FRONT 
    l=0;
    for(int i=iter_toSend[_WEST_F_][_INIX_]; i<= iter_toSend[_WEST_F_][_ENDX_]; i++)
        for(int j=iter_toSend[_WEST_F_][_INIY_]; j<= iter_toSend[_WEST_F_][_ENDY_]; j++)
            for(int k=iter_toSend[_WEST_F_][_INIZ_]; k<= iter_toSend[_WEST_F_][_ENDZ_]; k++)
            {
                pack_send_wf[l] = vec[ R_INDX( i , j, k, lNx, lNy)];
                l++; 
            }

    //EAST SOUTH 
    l=0;
    for(int i=iter_toSend[_EAST_S_][_INIX_]; i<= iter_toSend[_EAST_S_][_ENDX_]; i++)
        for(int j=iter_toSend[_EAST_S_][_INIY_]; j<= iter_toSend[_EAST_S_][_ENDY_]; j++)
            for(int k=iter_toSend[_EAST_S_][_INIZ_]; k<= iter_toSend[_EAST_S_][_ENDZ_]; k++)
            {
                pack_send_es[l] = vec[ R_INDX( i , j, k, lNx, lNy)];
                l++; 
            }

    //EAST NORTH
    l=0;
    for(int i=iter_toSend[_EAST_N_][_INIX_]; i<= iter_toSend[_EAST_N_][_ENDX_]; i++)
        for(int j=iter_toSend[_EAST_N_][_INIY_]; j<= iter_toSend[_EAST_N_][_ENDY_]; j++)
            for(int k=iter_toSend[_EAST_N_][_INIZ_]; k<= iter_toSend[_EAST_N_][_ENDZ_]; k++)
            {
                pack_send_en[l] = vec[ R_INDX( i , j, k, lNx, lNy)];
                l++; 
            }

    //EAST BACK 
    l=0;
    for(int i=iter_toSend[_EAST_B_][_INIX_]; i<= iter_toSend[_EAST_B_][_ENDX_]; i++)
        for(int j=iter_toSend[_EAST_B_][_INIY_]; j<= iter_toSend[_EAST_B_][_ENDY_]; j++)
            for(int k=iter_toSend[_EAST_B_][_INIZ_]; k<= iter_toSend[_EAST_B_][_ENDZ_]; k++)
            {
                pack_send_eb[l] = vec[ R_INDX( i , j, k, lNx, lNy)];
                l++; 
            }

    //EAST FRONT 
    l=0;
    for(int i=iter_toSend[_EAST_F_][_INIX_]; i<= iter_toSend[_EAST_F_][_ENDX_]; i++)
        for(int j=iter_toSend[_EAST_F_][_INIY_]; j<= iter_toSend[_EAST_F_][_ENDY_]; j++)
            for(int k=iter_toSend[_EAST_F_][_INIZ_]; k<= iter_toSend[_EAST_F_][_ENDZ_]; k++)
            {
                pack_send_ef[l] = vec[ R_INDX( i , j, k, lNx, lNy)];
                l++; 
            }

    //SOUTH BACK 
    l=0;
    for(int i=iter_toSend[_SOUTH_B_][_INIX_]; i<= iter_toSend[_SOUTH_B_][_ENDX_]; i++)
        for(int j=iter_toSend[_SOUTH_B_][_INIY_]; j<= iter_toSend[_SOUTH_B_][_ENDY_]; j++)
            for(int k=iter_toSend[_SOUTH_B_][_INIZ_]; k<= iter_toSend[_SOUTH_B_][_ENDZ_]; k++)
            {
                pack_send_sb[l] = vec[ R_INDX( i , j, k, lNx, lNy)];
                l++; 
            }

    //SOUTH FRONT 
    l=0;
    for(int i=iter_toSend[_SOUTH_F_][_INIX_]; i<= iter_toSend[_SOUTH_F_][_ENDX_]; i++)
        for(int j=iter_toSend[_SOUTH_F_][_INIY_]; j<= iter_toSend[_SOUTH_F_][_ENDY_]; j++)
            for(int k=iter_toSend[_SOUTH_F_][_INIZ_]; k<= iter_toSend[_SOUTH_F_][_ENDZ_]; k++)
            {
                pack_send_sf[l] = vec[ R_INDX( i , j, k, lNx, lNy)];
                l++; 
            }

    //NORTH BACK 
    l=0;
    for(int i=iter_toSend[_NORTH_B_][_INIX_]; i<= iter_toSend[_NORTH_B_][_ENDX_]; i++)
        for(int j=iter_toSend[_NORTH_B_][_INIY_]; j<= iter_toSend[_NORTH_B_][_ENDY_]; j++)
            for(int k=iter_toSend[_NORTH_B_][_INIZ_]; k<= iter_toSend[_NORTH_B_][_ENDZ_]; k++)
            {
                pack_send_nb[l] = vec[ R_INDX( i , j, k, lNx, lNy)];
                l++; 
            }

    //NORTH FRONT 
    l=0;
    for(int i=iter_toSend[_NORTH_F_][_INIX_]; i<= iter_toSend[_NORTH_F_][_ENDX_]; i++)
        for(int j=iter_toSend[_NORTH_F_][_INIY_]; j<= iter_toSend[_NORTH_F_][_ENDY_]; j++)
            for(int k=iter_toSend[_NORTH_F_][_INIZ_]; k<= iter_toSend[_NORTH_F_][_ENDZ_]; k++)
            {
                pack_send_nf[l] = vec[ R_INDX( i , j, k, lNx, lNy)];
                l++; 
            }
 
    //WEST SOUTH BACK 
    l=0;
    for(int i=iter_toSend[_WEST_S_B_][_INIX_]; i<= iter_toSend[_WEST_S_B_][_ENDX_]; i++)
        for(int j=iter_toSend[_WEST_S_B_][_INIY_]; j<= iter_toSend[_WEST_S_B_][_ENDY_]; j++)
            for(int k=iter_toSend[_WEST_S_B_][_INIZ_]; k<= iter_toSend[_WEST_S_B_][_ENDZ_]; k++)
            {
                pack_send_wsb[l] = vec[ R_INDX( i , j, k, lNx, lNy)];
                l++; 
            }

    //WEST NORTH BACK
    l=0;
    for(int i=iter_toSend[_WEST_N_B_][_INIX_]; i<= iter_toSend[_WEST_N_B_][_ENDX_]; i++)
        for(int j=iter_toSend[_WEST_N_B_][_INIY_]; j<= iter_toSend[_WEST_N_B_][_ENDY_]; j++)
            for(int k=iter_toSend[_WEST_N_B_][_INIZ_]; k<= iter_toSend[_WEST_N_B_][_ENDZ_]; k++)
            {
                pack_send_wnb[l] = vec[ R_INDX( i , j, k, lNx, lNy)];
                l++; 
            }

    //WEST SOUTH FRONT
    l=0;
    for(int i=iter_toSend[_WEST_S_F_][_INIX_]; i<= iter_toSend[_WEST_S_F_][_ENDX_]; i++)
        for(int j=iter_toSend[_WEST_S_F_][_INIY_]; j<= iter_toSend[_WEST_S_F_][_ENDY_]; j++)
            for(int k=iter_toSend[_WEST_S_F_][_INIZ_]; k<= iter_toSend[_WEST_S_F_][_ENDZ_]; k++)
            {
                pack_send_wsf[l] = vec[ R_INDX( i , j, k, lNx, lNy)];
                l++; 
            }

    //WEST NORTH FRONT
    l=0;
    for(int i=iter_toSend[_WEST_N_F_][_INIX_]; i<= iter_toSend[_WEST_N_F_][_ENDX_]; i++)
        for(int j=iter_toSend[_WEST_N_F_][_INIY_]; j<= iter_toSend[_WEST_N_F_][_ENDY_]; j++)
            for(int k=iter_toSend[_WEST_N_F_][_INIZ_]; k<= iter_toSend[_WEST_N_F_][_ENDZ_]; k++)
            {
                pack_send_wnf[l] = vec[ R_INDX( i , j, k, lNx, lNy)];
                l++; 
            }

    //EAST SOUTH BACK 
    l=0;
    for(int i=iter_toSend[_EAST_S_B_][_INIX_]; i<= iter_toSend[_EAST_S_B_][_ENDX_]; i++)
        for(int j=iter_toSend[_EAST_S_B_][_INIY_]; j<= iter_toSend[_EAST_S_B_][_ENDY_]; j++)
            for(int k=iter_toSend[_EAST_S_B_][_INIZ_]; k<= iter_toSend[_EAST_S_B_][_ENDZ_]; k++)
            {
                pack_send_esb[l] = vec[ R_INDX( i , j, k, lNx, lNy)];
                l++; 
            }

    //EAST NORTH BACK
    l=0;
    for(int i=iter_toSend[_EAST_N_B_][_INIX_]; i<= iter_toSend[_EAST_N_B_][_ENDX_]; i++)
        for(int j=iter_toSend[_EAST_N_B_][_INIY_]; j<= iter_toSend[_EAST_N_B_][_ENDY_]; j++)
            for(int k=iter_toSend[_EAST_N_B_][_INIZ_]; k<= iter_toSend[_EAST_N_B_][_ENDZ_]; k++)
            {
                pack_send_enb[l] = vec[ R_INDX( i , j, k, lNx, lNy)];
                l++; 
            }
 
    //EAST SOUTH FRONT
    l=0;
    for(int i=iter_toSend[_EAST_S_F_][_INIX_]; i<= iter_toSend[_EAST_S_F_][_ENDX_]; i++)
        for(int j=iter_toSend[_EAST_S_F_][_INIY_]; j<= iter_toSend[_EAST_S_F_][_ENDY_]; j++)
            for(int k=iter_toSend[_EAST_S_F_][_INIZ_]; k<= iter_toSend[_EAST_S_F_][_ENDZ_]; k++)
            {
                pack_send_esf[l] = vec[ R_INDX( i , j, k, lNx, lNy)];
                l++; 
            }

    //EAST NORTH FRONT
    l=0;
    for(int i=iter_toSend[_EAST_N_F_][_INIX_]; i<= iter_toSend[_EAST_N_F_][_ENDX_]; i++)
        for(int j=iter_toSend[_EAST_N_F_][_INIY_]; j<= iter_toSend[_EAST_N_F_][_ENDY_]; j++)
            for(int k=iter_toSend[_EAST_N_F_][_INIZ_]; k<= iter_toSend[_EAST_N_F_][_ENDZ_]; k++)
            {
                pack_send_enf[l] = vec[ R_INDX( i , j, k, lNx, lNy)];
                l++; 
            }

}


void ParallelTopology::unpack_simple(double *vec)
{
    //WEST
    int l=0;
    for(int i=iter_toRecv[_WEST_][_INIX_]; i<= iter_toRecv[_WEST_][_ENDX_]; i++)
        for(int j=iter_toRecv[_WEST_][_INIY_]; j<= iter_toRecv[_WEST_][_ENDY_]; j++)
            for(int k=iter_toRecv[_WEST_][_INIZ_]; k<= iter_toRecv[_WEST_][_ENDZ_]; k++)
            {
                vec[ R_INDX( i , j, k, lNx, lNy)] = pack_recv_w[l];
                l++;
            }

    //EAST
    l=0;
    for(int i=iter_toRecv[_EAST_][_INIX_]; i<= iter_toRecv[_EAST_][_ENDX_]; i++)
        for(int j=iter_toRecv[_EAST_][_INIY_]; j<= iter_toRecv[_EAST_][_ENDY_]; j++)
            for(int k=iter_toRecv[_EAST_][_INIZ_]; k<= iter_toRecv[_EAST_][_ENDZ_]; k++)
            {
                vec[ R_INDX( i , j, k, lNx, lNy)] = pack_recv_e[l];
                l++;
            }

    //SOUTH
    l=0;
    for(int i=iter_toRecv[_SOUTH_][_INIX_]; i<= iter_toRecv[_SOUTH_][_ENDX_]; i++)
        for(int j=iter_toRecv[_SOUTH_][_INIY_]; j<= iter_toRecv[_SOUTH_][_ENDY_]; j++)
            for(int k=iter_toRecv[_SOUTH_][_INIZ_]; k<= iter_toRecv[_SOUTH_][_ENDZ_]; k++)
            {
                vec[ R_INDX( i , j, k, lNx, lNy)] = pack_recv_s[l];
                l++;
            }

    //NORTH
    l=0;
    for(int i=iter_toRecv[_NORTH_][_INIX_]; i<= iter_toRecv[_NORTH_][_ENDX_]; i++)
        for(int j=iter_toRecv[_NORTH_][_INIY_]; j<= iter_toRecv[_NORTH_][_ENDY_]; j++)
            for(int k=iter_toRecv[_NORTH_][_INIZ_]; k<= iter_toRecv[_NORTH_][_ENDZ_]; k++)
            {
                vec[ R_INDX( i , j, k, lNx, lNy)] = pack_recv_n[l];
                l++;
            }

    //BACK
    l=0;
    for(int i=iter_toRecv[_BACK_][_INIX_]; i<= iter_toRecv[_BACK_][_ENDX_]; i++)
        for(int j=iter_toRecv[_BACK_][_INIY_]; j<= iter_toRecv[_BACK_][_ENDY_]; j++)
            for(int k=iter_toRecv[_BACK_][_INIZ_]; k<= iter_toRecv[_BACK_][_ENDZ_]; k++)
            {
                vec[ R_INDX( i , j, k, lNx, lNy)] = pack_recv_b[l];
                l++;
            }

    //FRONT 
    l=0;
    for(int i=iter_toRecv[_FRONT_][_INIX_]; i<= iter_toRecv[_FRONT_][_ENDX_]; i++)
        for(int j=iter_toRecv[_FRONT_][_INIY_]; j<= iter_toRecv[_FRONT_][_ENDY_]; j++)
            for(int k=iter_toRecv[_FRONT_][_INIZ_]; k<= iter_toRecv[_FRONT_][_ENDZ_]; k++)
            {
                vec[ R_INDX( i , j, k, lNx, lNy)] = pack_recv_f[l];
                l++;
            }
}

void ParallelTopology::unpack(double *vec)
{
    //WEST
    int l=0;
    for(int i=iter_toRecv[_WEST_][_INIX_]; i<= iter_toRecv[_WEST_][_ENDX_]; i++)
        for(int j=iter_toRecv[_WEST_][_INIY_]; j<= iter_toRecv[_WEST_][_ENDY_]; j++)
            for(int k=iter_toRecv[_WEST_][_INIZ_]; k<= iter_toRecv[_WEST_][_ENDZ_]; k++)
            {
                vec[ R_INDX( i , j, k, lNx, lNy)] = pack_recv_w[l];
                l++;
            }

    //EAST
    l=0;
    for(int i=iter_toRecv[_EAST_][_INIX_]; i<= iter_toRecv[_EAST_][_ENDX_]; i++)
        for(int j=iter_toRecv[_EAST_][_INIY_]; j<= iter_toRecv[_EAST_][_ENDY_]; j++)
            for(int k=iter_toRecv[_EAST_][_INIZ_]; k<= iter_toRecv[_EAST_][_ENDZ_]; k++)
            {
                vec[ R_INDX( i , j, k, lNx, lNy)] = pack_recv_e[l];
                l++;
            }

    //SOUTH
    l=0;
    for(int i=iter_toRecv[_SOUTH_][_INIX_]; i<= iter_toRecv[_SOUTH_][_ENDX_]; i++)
        for(int j=iter_toRecv[_SOUTH_][_INIY_]; j<= iter_toRecv[_SOUTH_][_ENDY_]; j++)
            for(int k=iter_toRecv[_SOUTH_][_INIZ_]; k<= iter_toRecv[_SOUTH_][_ENDZ_]; k++)
            {
                vec[ R_INDX( i , j, k, lNx, lNy)] = pack_recv_s[l];
                l++;
            }

    //NORTH
    l=0;
    for(int i=iter_toRecv[_NORTH_][_INIX_]; i<= iter_toRecv[_NORTH_][_ENDX_]; i++)
        for(int j=iter_toRecv[_NORTH_][_INIY_]; j<= iter_toRecv[_NORTH_][_ENDY_]; j++)
            for(int k=iter_toRecv[_NORTH_][_INIZ_]; k<= iter_toRecv[_NORTH_][_ENDZ_]; k++)
            {
                vec[ R_INDX( i , j, k, lNx, lNy)] = pack_recv_n[l];
                l++;
            }

    //BACK
    l=0;
    for(int i=iter_toRecv[_BACK_][_INIX_]; i<= iter_toRecv[_BACK_][_ENDX_]; i++)
        for(int j=iter_toRecv[_BACK_][_INIY_]; j<= iter_toRecv[_BACK_][_ENDY_]; j++)
            for(int k=iter_toRecv[_BACK_][_INIZ_]; k<= iter_toRecv[_BACK_][_ENDZ_]; k++)
            {
                vec[ R_INDX( i , j, k, lNx, lNy)] = pack_recv_b[l];
                l++;
            }

    //FRONT 
    l=0;
    for(int i=iter_toRecv[_FRONT_][_INIX_]; i<= iter_toRecv[_FRONT_][_ENDX_]; i++)
        for(int j=iter_toRecv[_FRONT_][_INIY_]; j<= iter_toRecv[_FRONT_][_ENDY_]; j++)
            for(int k=iter_toRecv[_FRONT_][_INIZ_]; k<= iter_toRecv[_FRONT_][_ENDZ_]; k++)
            {
                vec[ R_INDX( i , j, k, lNx, lNy)] = pack_recv_f[l];
                l++;
            }
    

    //2level comms
    //WEST SOUTH
    l=0;
    for(int i=iter_toRecv[_WEST_S_][_INIX_]; i<= iter_toRecv[_WEST_S_][_ENDX_]; i++)
        for(int j=iter_toRecv[_WEST_S_][_INIY_]; j<= iter_toRecv[_WEST_S_][_ENDY_]; j++)
            for(int k=iter_toRecv[_WEST_S_][_INIZ_]; k<= iter_toRecv[_WEST_S_][_ENDZ_]; k++)
            {
                vec[ R_INDX( i , j, k, lNx, lNy)] = pack_recv_ws[l];
                l++;
            }

    //WEST NORTH
    l=0;
    for(int i=iter_toRecv[_WEST_N_][_INIX_]; i<= iter_toRecv[_WEST_N_][_ENDX_]; i++)
        for(int j=iter_toRecv[_WEST_N_][_INIY_]; j<= iter_toRecv[_WEST_N_][_ENDY_]; j++)
            for(int k=iter_toRecv[_WEST_N_][_INIZ_]; k<= iter_toRecv[_WEST_N_][_ENDZ_]; k++)
            {
                vec[ R_INDX( i , j, k, lNx, lNy)] = pack_recv_wn[l];
                l++;
            }

    //WEST BACK
    l=0;
    for(int i=iter_toRecv[_WEST_B_][_INIX_]; i<= iter_toRecv[_WEST_B_][_ENDX_]; i++)
        for(int j=iter_toRecv[_WEST_B_][_INIY_]; j<= iter_toRecv[_WEST_B_][_ENDY_]; j++)
            for(int k=iter_toRecv[_WEST_B_][_INIZ_]; k<= iter_toRecv[_WEST_B_][_ENDZ_]; k++)
            {
                vec[ R_INDX( i , j, k, lNx, lNy)] = pack_recv_wb[l];
                l++;
            }

    //WEST FRONT
    l=0;
    for(int i=iter_toRecv[_WEST_F_][_INIX_]; i<= iter_toRecv[_WEST_F_][_ENDX_]; i++)
        for(int j=iter_toRecv[_WEST_F_][_INIY_]; j<= iter_toRecv[_WEST_F_][_ENDY_]; j++)
            for(int k=iter_toRecv[_WEST_F_][_INIZ_]; k<= iter_toRecv[_WEST_F_][_ENDZ_]; k++)
            {
                vec[ R_INDX( i , j, k, lNx, lNy)] = pack_recv_wf[l];
                l++;
            }

    //EAST SOUTH
    l=0;
    for(int i=iter_toRecv[_EAST_S_][_INIX_]; i<= iter_toRecv[_EAST_S_][_ENDX_]; i++)
        for(int j=iter_toRecv[_EAST_S_][_INIY_]; j<= iter_toRecv[_EAST_S_][_ENDY_]; j++)
            for(int k=iter_toRecv[_EAST_S_][_INIZ_]; k<= iter_toRecv[_EAST_S_][_ENDZ_]; k++)
            {
                vec[ R_INDX( i , j, k, lNx, lNy)] = pack_recv_es[l];
                l++;
            }

    //EAST NORTH
    l=0;
    for(int i=iter_toRecv[_EAST_N_][_INIX_]; i<= iter_toRecv[_EAST_N_][_ENDX_]; i++)
        for(int j=iter_toRecv[_EAST_N_][_INIY_]; j<= iter_toRecv[_EAST_N_][_ENDY_]; j++)
            for(int k=iter_toRecv[_EAST_N_][_INIZ_]; k<= iter_toRecv[_EAST_N_][_ENDZ_]; k++)
            {
                vec[ R_INDX( i , j, k, lNx, lNy)] = pack_recv_en[l];
                l++;
            }

    //EAST BACK
    l=0;
    for(int i=iter_toRecv[_EAST_B_][_INIX_]; i<= iter_toRecv[_EAST_B_][_ENDX_]; i++)
        for(int j=iter_toRecv[_EAST_B_][_INIY_]; j<= iter_toRecv[_EAST_B_][_ENDY_]; j++)
            for(int k=iter_toRecv[_EAST_B_][_INIZ_]; k<= iter_toRecv[_EAST_B_][_ENDZ_]; k++)
            {
                vec[ R_INDX( i , j, k, lNx, lNy)] = pack_recv_eb[l];
                l++;
            }

    //EAST FRONT
    l=0;
    for(int i=iter_toRecv[_EAST_F_][_INIX_]; i<= iter_toRecv[_EAST_F_][_ENDX_]; i++)
        for(int j=iter_toRecv[_EAST_F_][_INIY_]; j<= iter_toRecv[_EAST_F_][_ENDY_]; j++)
            for(int k=iter_toRecv[_EAST_F_][_INIZ_]; k<= iter_toRecv[_EAST_F_][_ENDZ_]; k++)
            {
                vec[ R_INDX( i , j, k, lNx, lNy)] = pack_recv_ef[l];
                l++;
            }

    //SOUTH BACK
    l=0;
    for(int i=iter_toRecv[_SOUTH_B_][_INIX_]; i<= iter_toRecv[_SOUTH_B_][_ENDX_]; i++)
        for(int j=iter_toRecv[_SOUTH_B_][_INIY_]; j<= iter_toRecv[_SOUTH_B_][_ENDY_]; j++)
            for(int k=iter_toRecv[_SOUTH_B_][_INIZ_]; k<= iter_toRecv[_SOUTH_B_][_ENDZ_]; k++)
            {
                vec[ R_INDX( i , j, k, lNx, lNy)] = pack_recv_sb[l];
                l++;
            }

    //SOUTH FRONT
    l=0;
    for(int i=iter_toRecv[_SOUTH_F_][_INIX_]; i<= iter_toRecv[_SOUTH_F_][_ENDX_]; i++)
        for(int j=iter_toRecv[_SOUTH_F_][_INIY_]; j<= iter_toRecv[_SOUTH_F_][_ENDY_]; j++)
            for(int k=iter_toRecv[_SOUTH_F_][_INIZ_]; k<= iter_toRecv[_SOUTH_F_][_ENDZ_]; k++)
            {
                vec[ R_INDX( i , j, k, lNx, lNy)] = pack_recv_sf[l];
                l++;
            }

    //NORTH BACK
    l=0;
    for(int i=iter_toRecv[_NORTH_B_][_INIX_]; i<= iter_toRecv[_NORTH_B_][_ENDX_]; i++)
        for(int j=iter_toRecv[_NORTH_B_][_INIY_]; j<= iter_toRecv[_NORTH_B_][_ENDY_]; j++)
            for(int k=iter_toRecv[_NORTH_B_][_INIZ_]; k<= iter_toRecv[_NORTH_B_][_ENDZ_]; k++)
            {
                vec[ R_INDX( i , j, k, lNx, lNy)] = pack_recv_nb[l];
                l++;
            }

    //NORTH FRONT
    l=0;
    for(int i=iter_toRecv[_NORTH_F_][_INIX_]; i<= iter_toRecv[_NORTH_F_][_ENDX_]; i++)
        for(int j=iter_toRecv[_NORTH_F_][_INIY_]; j<= iter_toRecv[_NORTH_F_][_ENDY_]; j++)
            for(int k=iter_toRecv[_NORTH_F_][_INIZ_]; k<= iter_toRecv[_NORTH_F_][_ENDZ_]; k++)
            {
                vec[ R_INDX( i , j, k, lNx, lNy)] = pack_recv_nf[l];
                l++;
            }

    //WEST SOUTH BACK 
    l=0;
    for(int i=iter_toRecv[_WEST_S_B_][_INIX_]; i<= iter_toRecv[_WEST_S_B_][_ENDX_]; i++)
        for(int j=iter_toRecv[_WEST_S_B_][_INIY_]; j<= iter_toRecv[_WEST_S_B_][_ENDY_]; j++)
            for(int k=iter_toRecv[_WEST_S_B_][_INIZ_]; k<= iter_toRecv[_WEST_S_B_][_ENDZ_]; k++)
            {
                vec[ R_INDX( i , j, k, lNx, lNy)] = pack_recv_wsb[l];
                l++;
            }

    //WEST NORTH BACK
    l=0;
    for(int i=iter_toRecv[_WEST_N_B_][_INIX_]; i<= iter_toRecv[_WEST_N_B_][_ENDX_]; i++)
        for(int j=iter_toRecv[_WEST_N_B_][_INIY_]; j<= iter_toRecv[_WEST_N_B_][_ENDY_]; j++)
            for(int k=iter_toRecv[_WEST_N_B_][_INIZ_]; k<= iter_toRecv[_WEST_N_B_][_ENDZ_]; k++)
            {
                vec[ R_INDX( i , j, k, lNx, lNy)] = pack_recv_wnb[l];
                l++;
            }

    //WEST SOUTH FRONT
    l=0;
    for(int i=iter_toRecv[_WEST_S_F_][_INIX_]; i<= iter_toRecv[_WEST_S_F_][_ENDX_]; i++)
        for(int j=iter_toRecv[_WEST_S_F_][_INIY_]; j<= iter_toRecv[_WEST_S_F_][_ENDY_]; j++)
            for(int k=iter_toRecv[_WEST_S_F_][_INIZ_]; k<= iter_toRecv[_WEST_S_F_][_ENDZ_]; k++)
            {
                vec[ R_INDX( i , j, k, lNx, lNy)] = pack_recv_wsf[l];
                l++;
            }

    //WEST NORTH FRONT
    l=0;
    for(int i=iter_toRecv[_WEST_N_F_][_INIX_]; i<= iter_toRecv[_WEST_N_F_][_ENDX_]; i++)
        for(int j=iter_toRecv[_WEST_N_F_][_INIY_]; j<= iter_toRecv[_WEST_N_F_][_ENDY_]; j++)
            for(int k=iter_toRecv[_WEST_N_F_][_INIZ_]; k<= iter_toRecv[_WEST_N_F_][_ENDZ_]; k++)
            {
                vec[ R_INDX( i , j, k, lNx, lNy)] = pack_recv_wnf[l];
                l++;
            }

    //EAST SOUTH BACK 
    l=0;
    for(int i=iter_toRecv[_EAST_S_B_][_INIX_]; i<= iter_toRecv[_EAST_S_B_][_ENDX_]; i++)
        for(int j=iter_toRecv[_EAST_S_B_][_INIY_]; j<= iter_toRecv[_EAST_S_B_][_ENDY_]; j++)
            for(int k=iter_toRecv[_EAST_S_B_][_INIZ_]; k<= iter_toRecv[_EAST_S_B_][_ENDZ_]; k++)
            {
                vec[ R_INDX( i , j, k, lNx, lNy)] = pack_recv_esb[l];
                l++;
            }

    //EAST NORTH BACK
    l=0;
    for(int i=iter_toRecv[_EAST_N_B_][_INIX_]; i<= iter_toRecv[_EAST_N_B_][_ENDX_]; i++)
        for(int j=iter_toRecv[_EAST_N_B_][_INIY_]; j<= iter_toRecv[_EAST_N_B_][_ENDY_]; j++)
            for(int k=iter_toRecv[_EAST_N_B_][_INIZ_]; k<= iter_toRecv[_EAST_N_B_][_ENDZ_]; k++)
            {
                vec[ R_INDX( i , j, k, lNx, lNy)] = pack_recv_enb[l];
                l++;
            }

    //EAST SOUTH FRONT
    l=0;
    for(int i=iter_toRecv[_EAST_S_F_][_INIX_]; i<= iter_toRecv[_EAST_S_F_][_ENDX_]; i++)
        for(int j=iter_toRecv[_EAST_S_F_][_INIY_]; j<= iter_toRecv[_EAST_S_F_][_ENDY_]; j++)
            for(int k=iter_toRecv[_EAST_S_F_][_INIZ_]; k<= iter_toRecv[_EAST_S_F_][_ENDZ_]; k++)
            {
                vec[ R_INDX( i , j, k, lNx, lNy)] = pack_recv_esf[l];
                l++;
            }

    //EAST NORTH FRONT
    l=0;
    for(int i=iter_toRecv[_EAST_N_F_][_INIX_]; i<= iter_toRecv[_EAST_N_F_][_ENDX_]; i++)
        for(int j=iter_toRecv[_EAST_N_F_][_INIY_]; j<= iter_toRecv[_EAST_N_F_][_ENDY_]; j++)
            for(int k=iter_toRecv[_EAST_N_F_][_INIZ_]; k<= iter_toRecv[_EAST_N_F_][_ENDZ_]; k++)
            {
                vec[ R_INDX( i , j, k, lNx, lNy)] = pack_recv_enf[l];
                l++;
            }

}


void ParallelTopology :: halo_exchange_simple()
{

    MPI_Request req_s[6];
    MPI_Request req_r[6];

    MPI_Status  stat_s[6];
    MPI_Status  stat_r[6];



    if(getNB(_WEST_) != _NO_NEIGHBOUR_ )
        MPI_Isend(pack_send_w, len_yz, MPI_DOUBLE, getNB(_WEST_), 0, RHEA_3DCOMM, &req_s[_WEST_]);

    if(getNB(_EAST_) != _NO_NEIGHBOUR_ )
        MPI_Isend(pack_send_e, len_yz, MPI_DOUBLE, getNB(_EAST_), 0, RHEA_3DCOMM, &req_s[_EAST_]);

    if(getNB(_SOUTH_) != _NO_NEIGHBOUR_ )
        MPI_Isend(pack_send_s, len_xz, MPI_DOUBLE, getNB(_SOUTH_), 0, RHEA_3DCOMM, &req_s[_SOUTH_]);

    if(getNB(_NORTH_) != _NO_NEIGHBOUR_ )
        MPI_Isend(pack_send_n, len_xz, MPI_DOUBLE, getNB(_NORTH_), 0, RHEA_3DCOMM, &req_s[_NORTH_]);

    if(getNB(_BACK_) != _NO_NEIGHBOUR_ )
        MPI_Isend(pack_send_b, len_xy, MPI_DOUBLE, getNB(_BACK_), 0, RHEA_3DCOMM, &req_s[_BACK_]);

    if(getNB(_FRONT_) != _NO_NEIGHBOUR_ )
        MPI_Isend(pack_send_f, len_xy, MPI_DOUBLE, getNB(_FRONT_), 0, RHEA_3DCOMM, &req_s[_FRONT_]);


    if(getNB(_WEST_) != _NO_NEIGHBOUR_ )
        MPI_Irecv(pack_recv_w, len_yz, MPI_DOUBLE, getNB(_WEST_), 0, RHEA_3DCOMM, &req_r[_WEST_]);

    if(getNB(_EAST_) != _NO_NEIGHBOUR_ )
        MPI_Irecv(pack_recv_e, len_yz, MPI_DOUBLE, getNB(_EAST_), 0, RHEA_3DCOMM, &req_r[_EAST_]);

    if(getNB(_SOUTH_) != _NO_NEIGHBOUR_ )
        MPI_Irecv(pack_recv_s, len_xz, MPI_DOUBLE, getNB(_SOUTH_), 0, RHEA_3DCOMM, &req_r[_SOUTH_]);

    if(getNB(_NORTH_) != _NO_NEIGHBOUR_ )
        MPI_Irecv(pack_recv_n, len_xz, MPI_DOUBLE, getNB(_NORTH_), 0, RHEA_3DCOMM, &req_r[_NORTH_]);

    if(getNB(_BACK_) != _NO_NEIGHBOUR_ )
        MPI_Irecv(pack_recv_b, len_xy, MPI_DOUBLE, getNB(_BACK_), 0, RHEA_3DCOMM, &req_r[_BACK_]);

    if(getNB(_FRONT_) != _NO_NEIGHBOUR_ )
        MPI_Irecv(pack_recv_f, len_xy, MPI_DOUBLE, getNB(_FRONT_), 0, RHEA_3DCOMM, &req_r[_FRONT_]);


    if(getNB(_WEST_) != _NO_NEIGHBOUR_ ) {
        MPI_Wait(&req_r[_WEST_], &stat_r[_WEST_]);
        MPI_Wait(&req_s[_WEST_], &stat_s[_WEST_]);
    }
    if(getNB(_EAST_) != _NO_NEIGHBOUR_ ) {
        MPI_Wait(&req_r[_EAST_], &stat_r[_EAST_]);
        MPI_Wait(&req_s[_EAST_], &stat_s[_EAST_]);
    }
    if(getNB(_SOUTH_) != _NO_NEIGHBOUR_ ) {
        MPI_Wait(&req_r[_SOUTH_], &stat_r[_SOUTH_]);
        MPI_Wait(&req_s[_SOUTH_], &stat_s[_SOUTH_]);
    }

    if(getNB(_NORTH_) != _NO_NEIGHBOUR_ ) {
        MPI_Wait(&req_r[_NORTH_], &stat_r[_NORTH_]);
        MPI_Wait(&req_s[_NORTH_], &stat_s[_NORTH_]);
    }
    if(getNB(_BACK_) != _NO_NEIGHBOUR_ ) {
        MPI_Wait(&req_r[_BACK_], &stat_r[_BACK_]);
        MPI_Wait(&req_s[_BACK_], &stat_s[_BACK_]);
    }
    if(getNB(_FRONT_) != _NO_NEIGHBOUR_ ){
        MPI_Wait(&req_r[_FRONT_], &stat_r[_FRONT_]);
        MPI_Wait(&req_s[_FRONT_], &stat_s[_FRONT_]);
    }


}


void ParallelTopology :: halo_exchange()
{

    MPI_Request req_s[26];
    MPI_Request req_r[26];

    MPI_Status  stat_s[26];
    MPI_Status  stat_r[26];



    if(getNB(_WEST_) != _NO_NEIGHBOUR_ )
        MPI_Isend(pack_send_w, len_yz, MPI_DOUBLE, getNB(_WEST_), 0, RHEA_3DCOMM, &req_s[_WEST_]);

    if(getNB(_EAST_) != _NO_NEIGHBOUR_ )
        MPI_Isend(pack_send_e, len_yz, MPI_DOUBLE, getNB(_EAST_), 0, RHEA_3DCOMM, &req_s[_EAST_]);

    if(getNB(_SOUTH_) != _NO_NEIGHBOUR_ )
        MPI_Isend(pack_send_s, len_xz, MPI_DOUBLE, getNB(_SOUTH_), 0, RHEA_3DCOMM, &req_s[_SOUTH_]);

    if(getNB(_NORTH_) != _NO_NEIGHBOUR_ )
        MPI_Isend(pack_send_n, len_xz, MPI_DOUBLE, getNB(_NORTH_), 0, RHEA_3DCOMM, &req_s[_NORTH_]);

    if(getNB(_BACK_) != _NO_NEIGHBOUR_ )
        MPI_Isend(pack_send_b, len_xy, MPI_DOUBLE, getNB(_BACK_), 0, RHEA_3DCOMM, &req_s[_BACK_]);

    if(getNB(_FRONT_) != _NO_NEIGHBOUR_ )
        MPI_Isend(pack_send_f, len_xy, MPI_DOUBLE, getNB(_FRONT_), 0, RHEA_3DCOMM, &req_s[_FRONT_]);

// 2nd level

    if(getNB(_WEST_S_) != _NO_NEIGHBOUR_ )
        MPI_Isend(pack_send_ws, len_1Dz, MPI_DOUBLE, getNB(_WEST_S_), 0, RHEA_3DCOMM, &req_s[_WEST_S_]);

    if(getNB(_WEST_N_) != _NO_NEIGHBOUR_ )
        MPI_Isend(pack_send_wn, len_1Dz, MPI_DOUBLE, getNB(_WEST_N_), 0, RHEA_3DCOMM, &req_s[_WEST_N_]);

    if(getNB(_WEST_B_) != _NO_NEIGHBOUR_ )
        MPI_Isend(pack_send_wb, len_1Dy, MPI_DOUBLE, getNB(_WEST_B_), 0, RHEA_3DCOMM, &req_s[_WEST_B_]);

    if(getNB(_WEST_F_) != _NO_NEIGHBOUR_ )
        MPI_Isend(pack_send_wf, len_1Dy, MPI_DOUBLE, getNB(_WEST_F_), 0, RHEA_3DCOMM, &req_s[_WEST_F_]);


    if(getNB(_EAST_S_) != _NO_NEIGHBOUR_ )
        MPI_Isend(pack_send_es, len_1Dz, MPI_DOUBLE, getNB(_EAST_S_), 0, RHEA_3DCOMM, &req_s[_EAST_S_]);

    if(getNB(_EAST_N_) != _NO_NEIGHBOUR_ )
        MPI_Isend(pack_send_en, len_1Dz, MPI_DOUBLE, getNB(_EAST_N_), 0, RHEA_3DCOMM, &req_s[_EAST_N_]);

    if(getNB(_EAST_B_) != _NO_NEIGHBOUR_ )
        MPI_Isend(pack_send_eb, len_1Dy, MPI_DOUBLE, getNB(_EAST_B_), 0, RHEA_3DCOMM, &req_s[_EAST_B_]);

    if(getNB(_EAST_F_) != _NO_NEIGHBOUR_ )
        MPI_Isend(pack_send_ef, len_1Dy, MPI_DOUBLE, getNB(_EAST_F_), 0, RHEA_3DCOMM, &req_s[_EAST_F_]);


    if(getNB(_SOUTH_B_) != _NO_NEIGHBOUR_ )
        MPI_Isend(pack_send_sb, len_1Dx, MPI_DOUBLE, getNB(_SOUTH_B_), 0, RHEA_3DCOMM, &req_s[_SOUTH_B_]);

    if(getNB(_SOUTH_F_) != _NO_NEIGHBOUR_ )
        MPI_Isend(pack_send_sf, len_1Dx, MPI_DOUBLE, getNB(_SOUTH_F_), 0, RHEA_3DCOMM, &req_s[_SOUTH_F_]);

    if(getNB(_NORTH_B_) != _NO_NEIGHBOUR_ )
        MPI_Isend(pack_send_nb, len_1Dx, MPI_DOUBLE, getNB(_NORTH_B_), 0, RHEA_3DCOMM, &req_s[_NORTH_B_]);

    if(getNB(_NORTH_F_) != _NO_NEIGHBOUR_ )
        MPI_Isend(pack_send_nf, len_1Dx, MPI_DOUBLE, getNB(_NORTH_F_), 0, RHEA_3DCOMM, &req_s[_NORTH_F_]);


// 3er level

    if(getNB(_WEST_S_B_) != _NO_NEIGHBOUR_ )
        MPI_Isend(pack_send_wsb, len_1pt, MPI_DOUBLE, getNB(_WEST_S_B_), 0, RHEA_3DCOMM, &req_s[_WEST_S_B_]);

    if(getNB(_WEST_N_B_) != _NO_NEIGHBOUR_ )
        MPI_Isend(pack_send_wnb, len_1pt, MPI_DOUBLE, getNB(_WEST_N_B_), 0, RHEA_3DCOMM, &req_s[_WEST_N_B_]);

    if(getNB(_WEST_S_F_) != _NO_NEIGHBOUR_ )
        MPI_Isend(pack_send_wsf, len_1pt, MPI_DOUBLE, getNB(_WEST_S_F_), 0, RHEA_3DCOMM, &req_s[_WEST_S_F_]);

    if(getNB(_WEST_N_F_) != _NO_NEIGHBOUR_ )
        MPI_Isend(pack_send_wnf, len_1pt, MPI_DOUBLE, getNB(_WEST_N_F_), 0, RHEA_3DCOMM, &req_s[_WEST_N_F_]);



    if(getNB(_EAST_S_B_) != _NO_NEIGHBOUR_ )
        MPI_Isend(pack_send_esb, len_1pt, MPI_DOUBLE, getNB(_EAST_S_B_), 0, RHEA_3DCOMM, &req_s[_EAST_S_B_]);

    if(getNB(_EAST_N_B_) != _NO_NEIGHBOUR_ )
        MPI_Isend(pack_send_enb, len_1pt, MPI_DOUBLE, getNB(_EAST_N_B_), 0, RHEA_3DCOMM, &req_s[_EAST_N_B_]);

    if(getNB(_EAST_S_F_) != _NO_NEIGHBOUR_ )
        MPI_Isend(pack_send_esf, len_1pt, MPI_DOUBLE, getNB(_EAST_S_F_), 0, RHEA_3DCOMM, &req_s[_EAST_S_F_]);

    if(getNB(_EAST_N_F_) != _NO_NEIGHBOUR_ )
        MPI_Isend(pack_send_enf, len_1pt, MPI_DOUBLE, getNB(_EAST_N_F_), 0, RHEA_3DCOMM, &req_s[_EAST_N_F_]);



///////

    if(getNB(_WEST_) != _NO_NEIGHBOUR_ )
        MPI_Irecv(pack_recv_w, len_yz, MPI_DOUBLE, getNB(_WEST_), 0, RHEA_3DCOMM, &req_r[_WEST_]);

    if(getNB(_EAST_) != _NO_NEIGHBOUR_ )
        MPI_Irecv(pack_recv_e, len_yz, MPI_DOUBLE, getNB(_EAST_), 0, RHEA_3DCOMM, &req_r[_EAST_]);

    if(getNB(_SOUTH_) != _NO_NEIGHBOUR_ )
        MPI_Irecv(pack_recv_s, len_xz, MPI_DOUBLE, getNB(_SOUTH_), 0, RHEA_3DCOMM, &req_r[_SOUTH_]);

    if(getNB(_NORTH_) != _NO_NEIGHBOUR_ )
        MPI_Irecv(pack_recv_n, len_xz, MPI_DOUBLE, getNB(_NORTH_), 0, RHEA_3DCOMM, &req_r[_NORTH_]);

    if(getNB(_BACK_) != _NO_NEIGHBOUR_ )
        MPI_Irecv(pack_recv_b, len_xy, MPI_DOUBLE, getNB(_BACK_), 0, RHEA_3DCOMM, &req_r[_BACK_]);

    if(getNB(_FRONT_) != _NO_NEIGHBOUR_ )
        MPI_Irecv(pack_recv_f, len_xy, MPI_DOUBLE, getNB(_FRONT_), 0, RHEA_3DCOMM, &req_r[_FRONT_]);

//2nd Level

    if(getNB(_WEST_S_) != _NO_NEIGHBOUR_ )
        MPI_Irecv(pack_recv_ws, len_1Dz, MPI_DOUBLE, getNB(_WEST_S_), 0, RHEA_3DCOMM, &req_r[_WEST_S_]);
    if(getNB(_WEST_N_) != _NO_NEIGHBOUR_ )
        MPI_Irecv(pack_recv_wn, len_1Dz, MPI_DOUBLE, getNB(_WEST_N_), 0, RHEA_3DCOMM, &req_r[_WEST_N_]);
    if(getNB(_WEST_B_) != _NO_NEIGHBOUR_ )
        MPI_Irecv(pack_recv_wb, len_1Dy, MPI_DOUBLE, getNB(_WEST_B_), 0, RHEA_3DCOMM, &req_r[_WEST_B_]);
    if(getNB(_WEST_F_) != _NO_NEIGHBOUR_ )
        MPI_Irecv(pack_recv_wf, len_1Dy, MPI_DOUBLE, getNB(_WEST_F_), 0, RHEA_3DCOMM, &req_r[_WEST_F_]);

    if(getNB(_EAST_S_) != _NO_NEIGHBOUR_ )
        MPI_Irecv(pack_recv_es, len_1Dz, MPI_DOUBLE, getNB(_EAST_S_), 0, RHEA_3DCOMM, &req_r[_EAST_S_]);
    if(getNB(_EAST_N_) != _NO_NEIGHBOUR_ )
        MPI_Irecv(pack_recv_en, len_1Dz, MPI_DOUBLE, getNB(_EAST_N_), 0, RHEA_3DCOMM, &req_r[_EAST_N_]);
    if(getNB(_EAST_B_) != _NO_NEIGHBOUR_ )
        MPI_Irecv(pack_recv_eb, len_1Dy, MPI_DOUBLE, getNB(_EAST_B_), 0, RHEA_3DCOMM, &req_r[_EAST_B_]);
    if(getNB(_EAST_F_) != _NO_NEIGHBOUR_ )
        MPI_Irecv(pack_recv_ef, len_1Dy, MPI_DOUBLE, getNB(_EAST_F_), 0, RHEA_3DCOMM, &req_r[_EAST_F_]);

    if(getNB(_SOUTH_B_) != _NO_NEIGHBOUR_ )
        MPI_Irecv(pack_recv_sb, len_1Dx, MPI_DOUBLE, getNB(_SOUTH_B_), 0, RHEA_3DCOMM, &req_r[_SOUTH_B_]);
    if(getNB(_SOUTH_F_) != _NO_NEIGHBOUR_ )
        MPI_Irecv(pack_recv_sf, len_1Dx, MPI_DOUBLE, getNB(_SOUTH_F_), 0, RHEA_3DCOMM, &req_r[_SOUTH_F_]);
    if(getNB(_NORTH_B_) != _NO_NEIGHBOUR_ )
        MPI_Irecv(pack_recv_nb, len_1Dx, MPI_DOUBLE, getNB(_NORTH_B_), 0, RHEA_3DCOMM, &req_r[_NORTH_B_]);
    if(getNB(_NORTH_F_) != _NO_NEIGHBOUR_ )
        MPI_Irecv(pack_recv_nf, len_1Dx, MPI_DOUBLE, getNB(_NORTH_F_), 0, RHEA_3DCOMM, &req_r[_NORTH_F_]);


// 3er level

    if(getNB(_WEST_S_B_) != _NO_NEIGHBOUR_ )
        MPI_Irecv(pack_recv_wsb, len_1pt, MPI_DOUBLE, getNB(_WEST_S_B_), 0, RHEA_3DCOMM, &req_r[_WEST_S_B_]);
    if(getNB(_WEST_N_B_) != _NO_NEIGHBOUR_ )
        MPI_Irecv(pack_recv_wnb, len_1pt, MPI_DOUBLE, getNB(_WEST_N_B_), 0, RHEA_3DCOMM, &req_r[_WEST_N_B_]);
    if(getNB(_WEST_S_F_) != _NO_NEIGHBOUR_ )
        MPI_Irecv(pack_recv_wsf, len_1pt, MPI_DOUBLE, getNB(_WEST_S_F_), 0, RHEA_3DCOMM, &req_r[_WEST_S_F_]);
    if(getNB(_WEST_N_F_) != _NO_NEIGHBOUR_ )
        MPI_Irecv(pack_recv_wnf, len_1pt, MPI_DOUBLE, getNB(_WEST_N_F_), 0, RHEA_3DCOMM, &req_r[_WEST_N_F_]);

    if(getNB(_EAST_S_B_) != _NO_NEIGHBOUR_ )
        MPI_Irecv(pack_recv_esb, len_1pt, MPI_DOUBLE, getNB(_EAST_S_B_), 0, RHEA_3DCOMM, &req_r[_EAST_S_B_]);
    if(getNB(_EAST_N_B_) != _NO_NEIGHBOUR_ )
        MPI_Irecv(pack_recv_enb, len_1pt, MPI_DOUBLE, getNB(_EAST_N_B_), 0, RHEA_3DCOMM, &req_r[_EAST_N_B_]);
    if(getNB(_EAST_S_F_) != _NO_NEIGHBOUR_ )
        MPI_Irecv(pack_recv_esf, len_1pt, MPI_DOUBLE, getNB(_EAST_S_F_), 0, RHEA_3DCOMM, &req_r[_EAST_S_F_]);
    if(getNB(_EAST_N_F_) != _NO_NEIGHBOUR_ )
        MPI_Irecv(pack_recv_enf, len_1pt, MPI_DOUBLE, getNB(_EAST_N_F_), 0, RHEA_3DCOMM, &req_r[_EAST_N_F_]);

// Wait for the communications to be done

    if(getNB(_WEST_) != _NO_NEIGHBOUR_ ){
        MPI_Wait(&req_r[_WEST_], &stat_r[_WEST_]);
        MPI_Wait(&req_s[_WEST_], &stat_s[_WEST_]);
    }
    if(getNB(_EAST_) != _NO_NEIGHBOUR_ ){
        MPI_Wait(&req_r[_EAST_], &stat_r[_EAST_]);
        MPI_Wait(&req_s[_EAST_], &stat_s[_EAST_]);
    }
    if(getNB(_SOUTH_) != _NO_NEIGHBOUR_ ){
        MPI_Wait(&req_r[_SOUTH_], &stat_r[_SOUTH_]);
        MPI_Wait(&req_s[_SOUTH_], &stat_s[_SOUTH_]);
    }
    if(getNB(_NORTH_) != _NO_NEIGHBOUR_ ){
        MPI_Wait(&req_r[_NORTH_], &stat_r[_NORTH_]);
        MPI_Wait(&req_s[_NORTH_], &stat_s[_NORTH_]);
    }
    if(getNB(_BACK_) != _NO_NEIGHBOUR_ ){
        MPI_Wait(&req_r[_BACK_], &stat_r[_BACK_]);
        MPI_Wait(&req_s[_BACK_], &stat_s[_BACK_]);
    }
    if(getNB(_FRONT_) != _NO_NEIGHBOUR_ ){
        MPI_Wait(&req_r[_FRONT_], &stat_r[_FRONT_]);
        MPI_Wait(&req_s[_FRONT_], &stat_s[_FRONT_]);
    }
//2nd level

    if(getNB(_WEST_S_) != _NO_NEIGHBOUR_ ){
        MPI_Wait(&req_r[_WEST_S_], &stat_r[_WEST_S_]);
        MPI_Wait(&req_s[_WEST_S_], &stat_s[_WEST_S_]);
    }
    if(getNB(_WEST_N_) != _NO_NEIGHBOUR_ ){
        MPI_Wait(&req_r[_WEST_N_], &stat_r[_WEST_N_]);
        MPI_Wait(&req_s[_WEST_N_], &stat_s[_WEST_N_]);
    }
    if(getNB(_WEST_B_) != _NO_NEIGHBOUR_ ){
        MPI_Wait(&req_r[_WEST_B_], &stat_r[_WEST_B_]);
        MPI_Wait(&req_s[_WEST_B_], &stat_s[_WEST_B_]);
    }
    if(getNB(_WEST_F_) != _NO_NEIGHBOUR_ ){
        MPI_Wait(&req_r[_WEST_F_], &stat_r[_WEST_F_]);
        MPI_Wait(&req_s[_WEST_F_], &stat_s[_WEST_F_]);
    }

    if(getNB(_EAST_S_) != _NO_NEIGHBOUR_ ){
        MPI_Wait(&req_r[_EAST_S_], &stat_r[_EAST_S_]);
        MPI_Wait(&req_s[_EAST_S_], &stat_s[_EAST_S_]);
    }
    if(getNB(_EAST_N_) != _NO_NEIGHBOUR_ ){
        MPI_Wait(&req_r[_EAST_N_], &stat_r[_EAST_N_]);
        MPI_Wait(&req_s[_EAST_N_], &stat_s[_EAST_N_]);
    }
    if(getNB(_EAST_B_) != _NO_NEIGHBOUR_ ){
        MPI_Wait(&req_r[_EAST_B_], &stat_r[_EAST_B_]);
        MPI_Wait(&req_s[_EAST_B_], &stat_s[_EAST_B_]);
    }
    if(getNB(_EAST_F_) != _NO_NEIGHBOUR_ ){
        MPI_Wait(&req_r[_EAST_F_], &stat_r[_EAST_F_]);
        MPI_Wait(&req_s[_EAST_F_], &stat_s[_EAST_F_]);
    }

    if(getNB(_SOUTH_B_) != _NO_NEIGHBOUR_ ){
        MPI_Wait(&req_r[_SOUTH_B_], &stat_r[_SOUTH_B_]);
        MPI_Wait(&req_s[_SOUTH_B_], &stat_s[_SOUTH_B_]);
    }
    if(getNB(_SOUTH_F_) != _NO_NEIGHBOUR_ ){
        MPI_Wait(&req_r[_SOUTH_F_], &stat_r[_SOUTH_F_]);
        MPI_Wait(&req_s[_SOUTH_F_], &stat_s[_SOUTH_F_]);
    }
    if(getNB(_NORTH_B_) != _NO_NEIGHBOUR_ ){
        MPI_Wait(&req_r[_NORTH_B_], &stat_r[_NORTH_B_]);
        MPI_Wait(&req_s[_NORTH_B_], &stat_s[_NORTH_B_]);
    }
    if(getNB(_NORTH_F_) != _NO_NEIGHBOUR_ ){
        MPI_Wait(&req_r[_NORTH_F_], &stat_r[_NORTH_F_]);
        MPI_Wait(&req_s[_NORTH_F_], &stat_s[_NORTH_F_]);
    }


//3er Level

    if(getNB(_WEST_S_B_) != _NO_NEIGHBOUR_ ){
        MPI_Wait(&req_r[_WEST_S_B_], &stat_r[_WEST_S_B_]);
        MPI_Wait(&req_s[_WEST_S_B_], &stat_s[_WEST_S_B_]);
    }
    if(getNB(_WEST_N_B_) != _NO_NEIGHBOUR_ ){
        MPI_Wait(&req_r[_WEST_N_B_], &stat_r[_WEST_N_B_]);
        MPI_Wait(&req_s[_WEST_N_B_], &stat_s[_WEST_N_B_]);
    }
    if(getNB(_WEST_S_F_) != _NO_NEIGHBOUR_ ){
        MPI_Wait(&req_r[_WEST_S_F_], &stat_r[_WEST_S_F_]);
        MPI_Wait(&req_s[_WEST_S_F_], &stat_s[_WEST_S_F_]);
    }
    if(getNB(_WEST_N_F_) != _NO_NEIGHBOUR_ ){
        MPI_Wait(&req_r[_WEST_N_F_], &stat_r[_WEST_N_F_]);
        MPI_Wait(&req_s[_WEST_N_F_], &stat_s[_WEST_N_F_]);
    }

    if(getNB(_EAST_S_B_) != _NO_NEIGHBOUR_ ){
        MPI_Wait(&req_r[_EAST_S_B_], &stat_r[_EAST_S_B_]);
        MPI_Wait(&req_s[_EAST_S_B_], &stat_s[_EAST_S_B_]);
    }
    if(getNB(_EAST_N_B_) != _NO_NEIGHBOUR_ ){
        MPI_Wait(&req_r[_EAST_N_B_], &stat_r[_EAST_N_B_]);
        MPI_Wait(&req_s[_EAST_N_B_], &stat_s[_EAST_N_B_]);
    }
    if(getNB(_EAST_S_F_) != _NO_NEIGHBOUR_ ){
        MPI_Wait(&req_r[_EAST_S_F_], &stat_r[_EAST_S_F_]);
        MPI_Wait(&req_s[_EAST_S_F_], &stat_s[_EAST_S_F_]);
    }
    if(getNB(_EAST_N_F_) != _NO_NEIGHBOUR_ ){
        MPI_Wait(&req_r[_EAST_N_F_], &stat_r[_EAST_N_F_]);
        MPI_Wait(&req_s[_EAST_N_F_], &stat_s[_EAST_N_F_]);
    }


}


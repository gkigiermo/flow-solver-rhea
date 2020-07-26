#include "comm_scheme.h"

comm_scheme::comm_scheme(domain* dom, int nprocsx, int nprocsy, int nprocsz)
{

    MPI_Comm_dup(MPI_COMM_WORLD, &RHEA_3DCOMM);
    MPI_Comm_rank(RHEA_3DCOMM, &rank);
    MPI_Comm_size(RHEA_3DCOMM, &np);


    npx= nprocsx;
    npy= nprocsy;
    npz= nprocsz;

    if( np != npx*npy*npz ){

        cout<<"Mismatched in the domain partition and number of procs launched"<<endl;
        MPI_Abort(RHEA_3DCOMM,0);
    }

    cout<<"testing 1"<<endl;
    MPI_Barrier(MPI_COMM_WORLD);
    int cellsx=dom->getGNx();
    int localNx= cellsx/npx;
    int divx = cellsx%npx;

    if(rank%npx < divx) 
        lNx=localNx+1;
    else
        lNx=localNx;


    lNx = lNx + 2;
    cout<<"testing 2"<<endl;
    int cellsy=dom->getGNy();
    int localNy= cellsy/npy;
    int divy = cellsy%npy;

    int plane_rank= rank%(npx*npy);

    if(plane_rank/npx < divy) 
        lNy=localNy+1;
    else
        lNy=localNy;


    lNy = lNy + 2;
    cout<<"testing 3"<<endl;
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


    //Finding neighbours

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
        else{
            neighb[_WEST_] = _NO_NEIGHBOUR_; 
        }
    }
    //EAST
    if(rank%npx != (npx - 1) ){
        neighb[_EAST_] = rank + 1;
    }
    else{
        if( periodic_x == 1 ){
            neighb[_EAST_] = rank - npx + 1;
        }
        else{
            neighb[_EAST_] =  _NO_NEIGHBOUR_;
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
        }
        else{
            neighb[_SOUTH_] = _NO_NEIGHBOUR_;
        }
    }

    //NORTH
    if(plane_rank/npx != npy - 1){
        neighb[_NORTH_] = rank + npx;
    }
    else{
        if(periodic_y == 1) {
            neighb[_NORTH_] = rank - (npy-1)*npx;
        }
        else{
            neighb[_NORTH_] = _NO_NEIGHBOUR_;
        }
    }

    //BACK
    int periodic_z;
    if(dom->getBoco(_BACK_) == _PERIODIC_)
        periodic_z = 1;
    else
        periodic_z = 0;

    plane_rank = rank/(npx*npy);
    if(plane_rank != 0){
        neighb[_BACK_] = rank - npx*npy;
    }
    else{
        if(periodic_z == 1 ){
            neighb[_BACK_] = rank + (npx*npy)*(npz-1);
        }
        else{
            neighb[_BACK_] = _NO_NEIGHBOUR_;
        }
    }

    //FRONT;
    if(plane_rank != npz-1){
        neighb[_FRONT_] = rank + npx*npy;
    }
    else{
        if(periodic_z == 1){
            neighb[_FRONT_] = rank - (npx*npy)*(npz-1);
        }
        else{
            neighb[_FRONT_] = _NO_NEIGHBOUR_;
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

}

void comm_scheme::printCommScheme(int proc)
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

void comm_scheme::printCommSchemeToFile(int proc)
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
            myfile<<" X : WEST  "<<neighb[_WEST_] <<" EAST  "<<neighb[_EAST_]<<endl;
            myfile<<" Y : SOUTH "<<neighb[_SOUTH_]<<" NORTH "<<neighb[_NORTH_]<<endl;
            myfile<<" Z : BACK  "<<neighb[_BACK_] <<" FRONT "<<neighb[_FRONT_]<<endl;
            myfile<<"--------"<<endl;
            myfile<<" Internal iterators: "<<endl;
            myfile<<" INI_X  "<<iter_common[_INNER_][_INIX_]<<" END_X "<<iter_common[_INNER_][_ENDX_]<<endl;
            myfile<<" INI_Y  "<<iter_common[_INNER_][_INIY_]<<" END_Y "<<iter_common[_INNER_][_ENDY_]<<endl;
            myfile<<" INI_Z  "<<iter_common[_INNER_][_INIZ_]<<" END_Z "<<iter_common[_INNER_][_ENDZ_]<<endl;
            myfile<<" Total domain iterators: "<<endl;
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
void comm_scheme::create_common_iters()
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
void comm_scheme::create_basic_bound_iters()
{
    // By default non halos are created
    for(int bd=_WEST_; bd<= _FRONT_; bd++)
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


/* This part assumes halo size 1*/
/* This rarely will be used */
void comm_scheme::create_halo_iters()
{
    // By default non halos are created
    for(int bd=_WEST_; bd<= _FRONT_; bd++)
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
    if(neighb[_SOUTH_]!= _NO_NEIGHBOUR_ ) {
        iter_halo[_SOUTH_][_INIX_] =   1;
        iter_halo[_SOUTH_][_ENDX_] = lNx-2;
        iter_halo[_SOUTH_][_INIY_] =   0;
        iter_halo[_SOUTH_][_ENDY_] =   0;
        iter_halo[_SOUTH_][_INIZ_] =   1;
        iter_halo[_SOUTH_][_ENDZ_] = lNz-2;
    }

    //NORTH
    if(neighb[_NORTH_]!= _NO_NEIGHBOUR_ ) {
        iter_halo[_NORTH_][_INIX_] =   1;
        iter_halo[_NORTH_][_ENDX_] =   lNx-2;
        iter_halo[_NORTH_][_INIY_] =   lNy-1;
        iter_halo[_NORTH_][_ENDY_] =   lNy-1;
        iter_halo[_NORTH_][_INIZ_] =   1;
        iter_halo[_NORTH_][_ENDZ_] = lNz-2;
    }

    //BACK
    if(neighb[_BACK_]!= _NO_NEIGHBOUR_ ) {
        iter_halo[_BACK_][_INIX_] =   1;
        iter_halo[_BACK_][_ENDX_] =   lNx-2;
        iter_halo[_BACK_][_INIY_] =   1;
        iter_halo[_BACK_][_ENDY_] =   lNy-2;
        iter_halo[_BACK_][_INIZ_] =   0;
        iter_halo[_BACK_][_ENDZ_] =   0;
    }

    //FRONT
    if(neighb[_FRONT_]!= _NO_NEIGHBOUR_ ) {
        iter_halo[_FRONT_][_INIX_] =   1;
        iter_halo[_FRONT_][_ENDX_] =   lNx-2;
        iter_halo[_FRONT_][_INIY_] =   1;
        iter_halo[_FRONT_][_ENDY_] =   lNy-2;
        iter_halo[_FRONT_][_INIZ_] =   lNz-1;
        iter_halo[_FRONT_][_ENDZ_] =   lNz-1;
    }

}


/* This part assumes halo size 1*/
void comm_scheme::create_toRecv_iters()
{

    for(int bd=_WEST_; bd<= _FRONT_; bd++)
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
    if(neighb[_SOUTH_]!= _NO_NEIGHBOUR_ ) {
        iter_toRecv[_SOUTH_][_INIX_] =   1;
        iter_toRecv[_SOUTH_][_ENDX_] = lNx-2;
        iter_toRecv[_SOUTH_][_INIY_] =   0;
        iter_toRecv[_SOUTH_][_ENDY_] =   0;
        iter_toRecv[_SOUTH_][_INIZ_] =   1;
        iter_toRecv[_SOUTH_][_ENDZ_] = lNz-2;
    }

    //NORTH
    if(neighb[_NORTH_]!= _NO_NEIGHBOUR_ ) {
        iter_toRecv[_NORTH_][_INIX_] =   1;
        iter_toRecv[_NORTH_][_ENDX_] =   lNx-2;
        iter_toRecv[_NORTH_][_INIY_] =   lNy-1;
        iter_toRecv[_NORTH_][_ENDY_] =   lNy-1;
        iter_toRecv[_NORTH_][_INIZ_] =   1;
        iter_toRecv[_NORTH_][_ENDZ_] = lNz-2;
    }

    //BACK
    if(neighb[_BACK_]!= _NO_NEIGHBOUR_ ) {
        iter_toRecv[_BACK_][_INIX_] =   1;
        iter_toRecv[_BACK_][_ENDX_] =   lNx-2;
        iter_toRecv[_BACK_][_INIY_] =   1;
        iter_toRecv[_BACK_][_ENDY_] =   lNy-2;
        iter_toRecv[_BACK_][_INIZ_] =   0;
        iter_toRecv[_BACK_][_ENDZ_] =   0;
    }

    //FRONT
    if(neighb[_FRONT_]!= _NO_NEIGHBOUR_ ) {
        iter_toRecv[_FRONT_][_INIX_] =   1;
        iter_toRecv[_FRONT_][_ENDX_] =   lNx-2;
        iter_toRecv[_FRONT_][_INIY_] =   1;
        iter_toRecv[_FRONT_][_ENDY_] =   lNy-2;
        iter_toRecv[_FRONT_][_INIZ_] =   lNz-1;
        iter_toRecv[_FRONT_][_ENDZ_] =   lNz-1;
    }

}

/* This part assumes halo size 1*/
void comm_scheme::create_toSend_iters()
{

    for(int bd=_WEST_; bd<= _FRONT_; bd++)
    {
        iter_toSend[bd][_INIX_] =   0;
        iter_toSend[bd][_ENDX_] =  -1;
        iter_toSend[bd][_INIY_] =   0;
        iter_toSend[bd][_ENDY_] =  -1;
        iter_toSend[bd][_INIZ_] =   0;
        iter_toSend[bd][_ENDZ_] =  -1;
    }

    //WEST
    if(neighb[_WEST_]!= _NO_NEIGHBOUR_ ) {
        iter_toSend[_WEST_][_INIX_] =   1;// <= this
        iter_toSend[_WEST_][_ENDX_] =   1;
        iter_toSend[_WEST_][_INIY_] =   1;
        iter_toSend[_WEST_][_ENDY_] = lNy-2;
        iter_toSend[_WEST_][_INIZ_] =   1;
        iter_toSend[_WEST_][_ENDZ_] = lNz-2;
    }

    //EAST
    if(neighb[_EAST_]!= _NO_NEIGHBOUR_ ) {
        iter_toSend[_EAST_][_INIX_] = lNx-2; // <=this
        iter_toSend[_EAST_][_ENDX_] = lNx-2;
        iter_toSend[_EAST_][_INIY_] =   1;
        iter_toSend[_EAST_][_ENDY_] = lNy-2;
        iter_toSend[_EAST_][_INIZ_] =   1;
        iter_toSend[_EAST_][_ENDZ_] = lNz-2;
    }

    //SOUTH
    if(neighb[_SOUTH_]!= _NO_NEIGHBOUR_ ) {
        iter_toSend[_SOUTH_][_INIX_] =   1;
        iter_toSend[_SOUTH_][_ENDX_] = lNx-2;
        iter_toSend[_SOUTH_][_INIY_] =   1; // <=this
        iter_toSend[_SOUTH_][_ENDY_] =   1;
        iter_toSend[_SOUTH_][_INIZ_] =   1;
        iter_toSend[_SOUTH_][_ENDZ_] = lNz-2;
    }

    //NORTH
    if(neighb[_NORTH_]!= _NO_NEIGHBOUR_ ) {
        iter_toSend[_NORTH_][_INIX_] =   1;
        iter_toSend[_NORTH_][_ENDX_] =   lNx-2;
        iter_toSend[_NORTH_][_INIY_] =   lNy-2; //<=this
        iter_toSend[_NORTH_][_ENDY_] =   lNy-2;
        iter_toSend[_NORTH_][_INIZ_] =   1;
        iter_toSend[_NORTH_][_ENDZ_] = lNz-2;
    }

    //BACK
    if(neighb[_BACK_]!= _NO_NEIGHBOUR_ ) {
        iter_toSend[_BACK_][_INIX_] =   1;
        iter_toSend[_BACK_][_ENDX_] =   lNx-2;
        iter_toSend[_BACK_][_INIY_] =   1;
        iter_toSend[_BACK_][_ENDY_] =   lNy-2;
        iter_toSend[_BACK_][_INIZ_] =   1; //<= this
        iter_toSend[_BACK_][_ENDZ_] =   1;
    }

    //FRONT
    if(neighb[_FRONT_]!= _NO_NEIGHBOUR_) {
        iter_toSend[_FRONT_][_INIX_] =   1;
        iter_toSend[_FRONT_][_ENDX_] =   lNx-2;
        iter_toSend[_FRONT_][_INIY_] =   1;
        iter_toSend[_FRONT_][_ENDY_] =   lNy-2;
        iter_toSend[_FRONT_][_INIZ_] =   lNz-2; //<= this
        iter_toSend[_FRONT_][_ENDZ_] =   lNz-2;
    }

}


void comm_scheme::create_global_iters()
{
    
    iter_glob_ind[_INIX_] = iter_common[_INNER_][_INIX_] + offx - 1;
    iter_glob_ind[_ENDX_] = iter_common[_INNER_][_ENDX_] + offx + 1;
    iter_glob_ind[_INIY_] = iter_common[_INNER_][_INIY_] + offy - 1;
    iter_glob_ind[_ENDY_] = iter_common[_INNER_][_ENDY_] + offy + 1;
    iter_glob_ind[_INIZ_] = iter_common[_INNER_][_INIZ_] + offz - 1;
    iter_glob_ind[_ENDZ_] = iter_common[_INNER_][_ENDZ_] + offz + 1;


    cout<<" Globals "<<rank<<" offz "<<offz<<endl;
    cout<<" "<<iter_glob_ind[_INIX_] <<endl;
    cout<<" "<<iter_glob_ind[_ENDX_] <<endl;
    cout<<" "<<iter_glob_ind[_INIY_] <<endl;
    cout<<" "<<iter_glob_ind[_ENDY_] <<endl;
    cout<<" "<<iter_glob_ind[_INIZ_] <<endl;
    cout<<" "<<iter_glob_ind[_ENDZ_] <<endl;





}

void comm_scheme::create_comm_arrays()
{
    //Size of the faces

    len_xy = (lNx-2)*(lNy-2);
    len_xz = (lNx-2)*(lNz-2);
    len_yz = (lNy-2)*(lNz-2);

    cout<<" lenxy "<<len_xy<<" lenxz "<<len_xz<<" lenyz "<<len_yz<<endl;
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

void comm_scheme::update(double *vec)
{
    pack(vec);
    halo_exchange();
    unpack(vec);
}


void comm_scheme::pack(double *vec)
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


void comm_scheme::unpack(double *vec)
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

void comm_scheme :: halo_exchange()
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


    if(getNB(_WEST_) != _NO_NEIGHBOUR_ )
        MPI_Wait(&req_r[_WEST_], &stat_r[_WEST_]);
    if(getNB(_EAST_) != _NO_NEIGHBOUR_ )
        MPI_Wait(&req_r[_EAST_], &stat_r[_EAST_]);
    if(getNB(_SOUTH_) != _NO_NEIGHBOUR_ )
        MPI_Wait(&req_r[_SOUTH_], &stat_r[_SOUTH_]);
    if(getNB(_NORTH_) != _NO_NEIGHBOUR_ )
        MPI_Wait(&req_r[_NORTH_], &stat_r[_NORTH_]);
    if(getNB(_BACK_) != _NO_NEIGHBOUR_ )
        MPI_Wait(&req_r[_BACK_], &stat_r[_BACK_]);
    if(getNB(_FRONT_) != _NO_NEIGHBOUR_ )
        MPI_Wait(&req_r[_FRONT_], &stat_r[_FRONT_]);



}

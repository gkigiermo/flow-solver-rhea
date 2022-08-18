#include "ManagerHDF5.hpp"

using namespace std;

ManagerHDF5::ManagerHDF5(ParallelTopology* topo,  char const* outputName,bool _gen_xdmf)
{
    myTopo = topo;
    sprintf(outname,"%s",outputName);

    gen_xdmf = _gen_xdmf;

    //To Generate the hdf5 file
    num_dims = 3;

    dim_gsize[0] = (topo->getMesh())->getGNz() + 2; //  _lNz_;
    dim_gsize[1] = (topo->getMesh())->getGNy() + 2;
    dim_gsize[2] = (topo->getMesh())->getGNx() + 2;

	dim_lsize[0] = topo->lenslabz; //  _lNz_;
    dim_lsize[1] = topo->lenslaby;
    dim_lsize[2] = topo->lenslabx;

	dim_offset[0] = topo->offslab_z; //  _lNz_;
    dim_offset[1] = topo->offslab_y;
    dim_offset[2] = topo->offslab_x;

    stride[0] = 1;
    stride[1] = 1;
    stride[2] = 1;

    block[0] = 1;
    block[1] = 1;
    block[2] = 1;


    slab_vec = new double[topo->lenslab]; 

    _lNx_ = topo->getlNx();
    _lNy_ = topo->getlNy();
    _lNz_ = topo->getlNz();

    num_dim1D = 1;
    array_1D[0] = 1;

}

void ManagerHDF5::addField(DistributedArray* newfield)
{
    fieldList.push_back(newfield);

}

void ManagerHDF5::printOnScreen()
{
    if(myTopo->getRank() == 0){

        cout<<"Fields on the output "<<outname<<endl;

        for(std::list<DistributedArray*>::iterator it=fieldList.begin();it!=fieldList.end();++it)
        {

            cout<<(*it)->printName()<<endl;
        }
    }
}

void ManagerHDF5::write(int it ) //, double time, bool xdmf_file)
{
    char filename[100];

    sprintf(filename,"%s_%d.h5",outname,it);

    hid_t fa_plist_id = H5Pcreate(H5P_FILE_ACCESS);

    MPI_Comm comm  = MPI_COMM_WORLD;
    MPI_Info info  = MPI_INFO_NULL;
 
    H5Pset_fapl_mpio(fa_plist_id, comm , info);

    hid_t       file_id;   /* file identifier */
    herr_t      status;


    file_id = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, fa_plist_id);

    if(file_id < 0 )
    {
        cout<<" Error creating the file "<<filename<<endl;
        MPI_Abort(MPI_COMM_WORLD,0);
    }

    hid_t fa_plist_id2 = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(fa_plist_id2, H5FD_MPIO_COLLECTIVE);


    hid_t filespace_id = H5Screate_simple(num_dims, dim_gsize, NULL); /*global data space*/
    hid_t memspace_id  = H5Screate_simple(num_dims, dim_lsize, NULL); /*local data space*/

     
    for(std::list<DistributedArray*>::iterator it=fieldList.begin();it!=fieldList.end();++it)
    {

        hid_t dataset_id1 = H5Dcreate2(file_id,(*it)->printName(), H5T_NATIVE_DOUBLE, filespace_id, 
                H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

        status = H5Sselect_hyperslab(filespace_id,H5S_SELECT_SET,dim_offset,stride,dim_lsize,block);

        int l=0;
        for(int k = myTopo->iter_slab[_INIZ_]; k <= myTopo->iter_slab[_ENDZ_]; k++)
            for(int j =  myTopo->iter_slab[_INIY_]; j <=  myTopo->iter_slab[_ENDY_]; j++)
                for(int i =  myTopo->iter_slab[_INIX_]; i <=  myTopo->iter_slab[_ENDX_]; i++){
                    slab_vec[l] =  (*it)->vector[I1D(i,j,k)];
                    l=l+1;
                }

        status = H5Dwrite(dataset_id1, H5T_NATIVE_DOUBLE, memspace_id, filespace_id, fa_plist_id2,
                slab_vec);


        status = H5Dclose(dataset_id1);

    }

   status = H5Sclose(filespace_id);
   status = H5Sclose(memspace_id);

   if(status < 0)
   {
       cout<<" ERROR Writing HDF5 file"<<endl;
       MPI_Abort(MPI_COMM_WORLD,666);
   }
   hid_t attspace_id = H5Screate_simple(num_dim1D, array_1D, NULL);

   char attribute_name[100];


   map<string, int>::iterator it_i;
   for (it_i = iattrib.begin(); it_i != iattrib.end(); it_i++)
   {
 
       sprintf(attribute_name,"%s",(it_i->first).c_str());
       hid_t attri0 = H5Acreate(file_id,attribute_name,H5T_NATIVE_INT, attspace_id,H5P_DEFAULT,H5P_DEFAULT);
       status = H5Awrite(attri0, H5T_NATIVE_INT, &(it_i->second));
       status = H5Aclose(attri0);
   }

   map<string, double>::iterator it_d;
   for (it_d = dattrib.begin(); it_d != dattrib.end(); it_d++)
   {
    
       sprintf(attribute_name,"%s",(it_d->first).c_str());
       hid_t attri1 = H5Acreate(file_id,attribute_name,H5T_NATIVE_DOUBLE, attspace_id,H5P_DEFAULT,H5P_DEFAULT);
       status = H5Awrite(attri1, H5T_NATIVE_DOUBLE, &(it_d->second));
       status = H5Aclose(attri1);
   }

   status = H5Sclose(attspace_id);


   status = H5Fclose(file_id); 

   if( gen_xdmf ) {
	   if(myTopo->getRank() == 0 ) {
		   char filename2[100];
		   char filename3[100];


		   sprintf(filename2,"%s_%d.xdmf",outname,it);
		   sprintf(filename3,"%s_%d",outname,it);

		   ofstream outfile(filename2);

		   outfile << "<?xml version='1.0' ?>" << endl;
		   outfile << "<!DOCTYPE Xdmf SYSTEM 'Xdmf.dtd' []>" << endl;
		   outfile << "<Xdmf Version='2.0'>" << endl;
		   outfile << "  <Domain> " << endl;
		   outfile << "    <Grid Name='"<<filename3<<"' GridType='Uniform'>" << endl;
		   outfile << "      <Topology TopologyType='3DSMesh' Dimensions='"<<dim_gsize[0]<<" "<<dim_gsize[1]<<" "<<dim_gsize[2]<<"'/>" << endl;
		   outfile << "      <Geometry GeometryType='X_Y_Z'>" << endl;
		   outfile << "        <DataItem Name='x' Dimensions='"<<dim_gsize[0]<<" "<<dim_gsize[1]<<" "<<dim_gsize[2]<<"' NumberType='Float' Precision='8' Format='HDF'>" << endl;
		   outfile << "            "<<filename<<":/x" << endl;
		   outfile << "        </DataItem>" << endl;
		   outfile << "        <DataItem Name='y' Dimensions='"<<dim_gsize[0]<<" "<<dim_gsize[1]<<" "<<dim_gsize[2]<<"' NumberType='Float' Precision='8' Format='HDF'>" << endl;
		   outfile << "            "<<filename<<":/y" << endl;
		   outfile << "        </DataItem>" << endl;
		   outfile << "        <DataItem Name='z' Dimensions='"<<dim_gsize[0]<<" "<<dim_gsize[1]<<" "<<dim_gsize[2]<<"' NumberType='Float' Precision='8' Format='HDF'>" << endl;
		   outfile << "            "<<filename<<":/z" << endl;
		   outfile << "        </DataItem>" << endl;
		   outfile << "      </Geometry>" << endl;
		   for(std::list<DistributedArray*>::iterator it=fieldList.begin();it!=fieldList.end();++it) {
			   if( strcmp((*it)->printName(),"x") != 0 && strcmp((*it)->printName(),"y") != 0 && strcmp((*it)->printName(),"z") != 0 ) {
				   outfile << "      <Attribute Name='"<<(*it)->printName()<<"' AttributeType='Scalar' Center='Node'>" << endl;
				   outfile << "        <DataItem Dimensions='"<<dim_gsize[0]<<" "<<dim_gsize[1]<<" "<<dim_gsize[2]<<"' NumberType='Float' Precision='8' Format='HDF'>" << endl;
				   outfile << "          "<<filename<<":/"<<(*it)->printName() << endl;
				   outfile << "        </DataItem>" << endl;
				   outfile << "      </Attribute>" << endl;
			   }
		   }
		   outfile << "    </Grid>" << endl;
		   outfile << "  </Domain>" << endl;
		   outfile << "</Xdmf>" << endl;

		   outfile.close();



	   }
   }


}

void ManagerHDF5::read( char const* inputName)
{
    char filename[100];
    
    sprintf(filename,"%s",inputName);

    hid_t fa_plist_id = H5Pcreate(H5P_FILE_ACCESS);

    MPI_Comm comm  = MPI_COMM_WORLD;
    MPI_Info info  = MPI_INFO_NULL;
 
    H5Pset_fapl_mpio(fa_plist_id, comm , info);


    hid_t       file_id;   /* file identifier */
    herr_t      status;

    file_id = H5Fopen(filename, H5F_ACC_RDWR, fa_plist_id);

    if(file_id < 0 )
    {
        cout<<" Error opening the file "<<filename<<endl;
        MPI_Abort(MPI_COMM_WORLD,0);
    }


    hid_t fa_plist_id2 = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(fa_plist_id2, H5FD_MPIO_COLLECTIVE);

    hid_t filespace_id = H5Screate_simple(num_dims, dim_gsize, NULL);
    hid_t memspace_id  = H5Screate_simple(num_dims, dim_lsize, NULL);

     


    for(std::list<DistributedArray*>::iterator it=fieldList.begin();it!=fieldList.end();++it)
    {
    
        hid_t dataset_id1 = H5Dopen(file_id,(*it)->printName(), H5P_DEFAULT); 

        status = H5Sselect_hyperslab(filespace_id,H5S_SELECT_SET,dim_offset,stride,dim_lsize,block);

        status = H5Dread(dataset_id1, H5T_NATIVE_DOUBLE, memspace_id, filespace_id, fa_plist_id2,
                slab_vec);

        int l=0;
        for(int k = myTopo->iter_slab[_INIZ_]; k <= myTopo->iter_slab[_ENDZ_]; k++)
            for(int j =  myTopo->iter_slab[_INIY_]; j <=  myTopo->iter_slab[_ENDY_]; j++)
                for(int i =  myTopo->iter_slab[_INIX_]; i <=  myTopo->iter_slab[_ENDX_]; i++){
                    (*it)->vector[I1D(i,j,k)] = slab_vec[l];
                    l=l+1;
                }

        (*it)->update();

        status = H5Dclose(dataset_id1);

    }

   status = H5Sclose(filespace_id);
   status = H5Sclose(memspace_id);

   if(status < 0)
   {
       cout<<" ERROR Reading HDF5 file"<<endl;
       MPI_Abort(MPI_COMM_WORLD,666);
   }


   hid_t attspace_id = H5Screate_simple(num_dim1D, array_1D, NULL);

   char attribute_name[100];


   map<string, int>::iterator it_i;
   for (it_i = iattrib.begin(); it_i != iattrib.end(); it_i++)
   {
       sprintf(attribute_name,"%s",(it_i->first).c_str());
       hid_t attri0 = H5Aopen(file_id,attribute_name,H5P_DEFAULT);
       status = H5Aread(attri0, H5T_NATIVE_INT, &(it_i->second));
       status = H5Aclose(attri0);
   }

   map<string, double>::iterator it_d;
   for (it_d = dattrib.begin(); it_d != dattrib.end(); it_d++)
   {
       sprintf(attribute_name,"%s",(it_d->first).c_str());
       //hid_t attri1 = H5Aopen(file_id,attribute_name,H5T_NATIVE_DOUBLE);
       hid_t attri1 = H5Aopen(file_id,attribute_name,H5P_DEFAULT);
       status = H5Aread(attri1, H5T_NATIVE_DOUBLE, &(it_d->second));
       status = H5Aclose(attri1);
   }

   status = H5Sclose(attspace_id);

   status = H5Fclose(file_id); 

}

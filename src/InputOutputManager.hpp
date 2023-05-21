#ifndef _INPUT_OUTPUT_MANAGER_
#define _INPUT_OUTPUT_MANAGER_

////////// INCLUDES //////////
#include <cmath>
#include <iostream>
#include <hdf5.h>
#include <list>
#include <string.h>
#include <limits>
#include <map>
#include "MacroParameters.hpp"
#include "ParallelTopology.hpp"
#include "DistributedArray.hpp"

class WriteReadHDF5 {

    public:
        WriteReadHDF5(){};
        WriteReadHDF5(ParallelTopology*,char const*,bool);
        void addField(DistributedArray*);
        void printOnScreen();
        void write(int);
        void read(char const*);

        void addAttributeDouble(std::string str){ double dval=0.0; dattrib[str]=dval;};
        void addAttributeInt(std::string str){ int ival =0; iattrib[str] = ival; };

        void setAttribute(std::string str,double dval){ dattrib[str]=dval;};
        void setAttribute(std::string str,int ival){ iattrib[str] = ival; };

        double getAttributeDouble(std::string str){ return dattrib[str];};
        int    getAttributeInt(std::string str){ return iattrib[str];};


    private:
        
        ParallelTopology* myTopo;
        char outname[50];

        std::list<DistributedArray*> fieldList;

        bool gen_xdmf;

        int num_dims;
        hsize_t dim_gsize[3];               // dataset dimensions
        hsize_t dim_lsize[3];               // dataset dimensions
        hsize_t dim_offset[3];
        hsize_t stride[3];
        hsize_t block[3];

        double *slab_vec;

        int _lNx_;
        int _lNy_;
        int _lNz_;

        // Variables to store the attribute

        hsize_t array_1D[1];
        int num_dim1D;


        //Attribute list

	std::map<std::string,double> dattrib;
	std::map<std::string,int> iattrib;


};

//////////// TemporalPointProbe CLASS //////////
class TemporalPointProbe {

    public:

        ////////// CONSTRUCTORS & DESTRUCTOR //////////
        TemporalPointProbe();					/// Default constructor
        TemporalPointProbe(const int &x_position_, const int &y_position_, const int &z_position_, const std::string output_file_name_, ComputationalDomain* mesh_, ParallelTopology* topo_);			     /// Parametrized constructor
        virtual ~TemporalPointProbe();				/// Destructor

	////////// GET FUNCTIONS //////////
        inline double getPositionX() { return( x_position ); };
        inline double getPositionY() { return( y_position ); };
        inline double getPositionZ() { return( z_position ); };
        inline double getTime() { return( time ); };
        inline std::string getOutputFileName() { return( output_file_name ); };
        inline int getGlobalOwnerRank() { return( global_owner_rank ); };
        inline int getLocalIndexI() { return( i_local_index ); };
        inline int getLocalIndexJ() { return( j_local_index ); };
        inline int getLocalIndexK() { return( k_local_index ); };

	////////// SET FUNCTIONS //////////
        inline void setPositionX(double x_position_) { x_position = x_position_; };
        inline void setPositionY(double y_position_) { y_position = y_position_; };
        inline void setPositionZ(double z_position_) { z_position = z_position_; };
        inline void setTime(double time_) { time = time_; };
        inline void setOutputFileName(std::string output_file_name_) { output_file_name = output_file_name_; };
        inline void setGlobalOwnerRank(int global_owner_rank_) { global_owner_rank = global_owner_rank_; };
        inline void setLocalIndexI(int i_local_index_) { i_local_index = i_local_index_; };
        inline void setLocalIndexJ(int j_local_index_) { j_local_index = j_local_index_; };
        inline void setLocalIndexK(int k_local_index_) { k_local_index = k_local_index_; };

	////////// PROBE METHODS //////////
	
	/// Locate closest grid point to probe
        virtual void locateClosestGridPointToProbe();

    protected:

        ////////// PROBE PARAMETERS //////////
	
        /// Probe properties 
        double x_position;         	 	 	/// Position in x [m]
        double y_position;         	 	 	/// Position in y [m]
        double z_position;         	 	 	/// Position in z [m]
        double time = -1.0;         	 	 	/// Time [s]
	std::string output_file_name;			/// Output file name

        /// Computational parameters
	int global_owner_rank = -1;			/// Rank of processor owning the probe
	int i_local_index = -1;				/// Local index i
	int j_local_index = -1;				/// Local index j
	int k_local_index = -1;				/// Local index k

	////////// COMPUTATIONAL DOMtimeLEL TOPOLOGY //////////
        ComputationalDomain *mesh;			/// Computational domain
        ParallelTopology *topo;				/// Parallel topology

    private:

};

#endif /*_INPUT_OUTPUT_MANAGER_*/

EXECUTABLE   = RHEA.exe
MAIN         = myRHEA.cpp
PROJECT_PATH = $(RHEA_PATH)
SRC_DIR      = $(PROJECT_PATH)/src
CXX          = mpicxx
# CPU FLAGS
CXXFLAGS     = -O3 -Wall -std=c++0x -Wno-unknown-pragmas -I$(PROJECT_PATH)
# CPU-GPU FLAGS
#CXXFLAGS     = -fast -acc -ta=tesla:managed -Minfo=accel -O3 -Wall -std=c++0x -I$(PROJECT_PATH)
#CXXFLAGS     = -fast -acc -ta=tesla -Minfo=accel -O3 -Wall -std=c++0x -I$(PROJECT_PATH)
# UBUNTU - LINUX
INC_LIB_YAML =
INC_DIR_YAML =
INC_LIB_HDF5 = -L/usr/lib/x86_64-linux-gnu/hdf5/openmpi 
INC_DIR_HDF5 = -I/usr/include/hdf5/openmpi
# MAC - OS X
#INC_LIB_YAML = -L/usr/local/lib
#INC_DIR_YAML = -I/usr/local/include
#INC_LIB_HDF5 =
#INC_DIR_HDF5 =
LDFLAGS      = -lyaml-cpp -lhdf5



# !! THE LINES BELOW SHOULD NOT BE MODIFIED !! #

OBJS = $(SRC_DIR)/*.cpp
INC_LIB = $(INC_LIB_YAML) $(INC_LIB_HDF5)
INC_DIR = $(INC_DIR_YAML) $(INC_DIR_HDF5)

$(EXECUTABLE): $(OBJS)
	$(CXX) $(MAIN) $(CXXFLAGS) $(OBJS) -o $@ $(INC_LIB) $(INC_DIR) $(LDFLAGS)

.PHONY: clean
clean:
	$(RM) $(EXECUTABLE) *.o


EXECUTABLE   = RHEA.exe
MAIN         = myRHEA.cpp
SRC_DIR      = ./src
CXX          = mpic++
CXXFLAGS     = -O3 -Wall -std=c++0x
INC_LIB_YAML = -L/usr/local/lib
INC_DIR_YAML = -I/usr/local/include
INC_LIB_HDF5 = 
INC_DIR_HDF5 = 
LDFLAGS      = -lyaml-cpp -lhdf5



# !! THE LINES BELOW SHOULD NOT BE MODIFIED !! #

OBJS = $(SRC_DIR)/*.cpp
INC_LIB = $(INC_LIB_YAML) $(INC_LIB_HDF5)
INC_DIR = $(INC_DIR_YAML) $(INC_DIR_HDF5)

$(EXECUTABLE): $(OBJS)
	$(CXX) $(MAIN) $(CXXFLAGS) $(OBJS) -o $@ $(INC_LIB) $(INC_DIR) $(LDFLAGS)

.PHONY: clean
clean:
	$(RM) $(EXECUTABLE)


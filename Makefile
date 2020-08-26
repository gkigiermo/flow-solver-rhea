TARGET       = RHEA.exe
SRC_DIRS     = ./src
CXX          = mpic++
CXXFLAGS     = -O3 -Wall -std=c++0x
INC_LIB_YAML = -L/usr/local/lib
INC_DIR_YAML = -I/usr/local/include
INC_LIB_HDF5 = -L/usr/local/lib
INC_DIR_HDF5 = -I/usr/local/include
LDFLAGS      = -lyaml-cpp -lhdf5



SRCS = $(shell find $(SRC_DIRS) -name *.cpp)
OBJS = $(addsuffix .o,$(basename $(SRCS)))

INC_LIB = $(INC_LIB_YAML) $(INC_LIB_HDF5)
INC_DIR = $(INC_DIR_YAML) $(INC_DIR_HDF5)

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) $(OBJS) -o $@ $(INC_LIB) $(INC_DIR) $(LDFLAGS)

.PHONY: clean
clean:
	$(RM) $(TARGET) $(OBJS)

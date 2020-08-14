CC=mpic++
#CFLAGS=-O3 -L/usr/local/Cellar/hdf5/1.12.0_1/lib/ -lhdf5_cpp -I/usr/local/Cellar/hdf5/1.12.0_1/include 
CFLAGS=-O3 -lhdf5  
OBJ=src/*.cpp
OUT= test1

all:
#	$(CC) heat_h5_hyper.cpp $(OBJ) -o $(OUT) $(CFLAGS)
	$(CC) FlowSolverRHEA.cpp -L/usr/local/lib -I/usr/local/include -std=c++0x $(OBJ) -o $(OUT) $(CFLAGS) -lyaml-cpp

clean:
	rm  src/*.o *.info


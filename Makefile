CC=mpic++
#CFLAGS=-O3 -L/usr/local/Cellar/hdf5/1.12.0_1/lib/ -lhdf5_cpp -I/usr/local/Cellar/hdf5/1.12.0_1/include 
CFLAGS=-O3 -L/usr/local/Cellar/hdf5-mpi/1.12.0_1/lib/ -lhdf5 -I/usr/local/Cellar/hdf5-mpi/1.12.0_1/include   
OBJ=src/*.cpp
OUT= test1

all:
	$(CC) heat_h5_hyper.cpp $(OBJ) -o $(OUT) $(CFLAGS)

clean:
	rm  src/*.o *.info


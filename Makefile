CC=mpic++
#CFLAGS=-O3 -L/usr/local/Cellar/hdf5/1.12.0_1/lib/ -lhdf5_cpp -I/usr/local/Cellar/hdf5/1.12.0_1/include 
CFLAGS=-O3 -lhdf5  
OBJ=src/*.cpp
OUT=RHEA.exe
GUIFLAGS=-Wall -I/Users/goyarzun/DEV/apps/yaml-cpp/include -L/Users/goyarzun/DEV/apps/yaml-cpp/build5/ -lyaml-cpp


all:
#	$(CC) heat_h5_hyper.cpp $(OBJ) -o $(OUT) $(CFLAGS)
#	$(CC) myRHEA.cpp -L/usr/local/lib -I/usr/local/include -std=c++0x $(OBJ) -o $(OUT) $(CFLAGS) -lyaml-cpp
#	$(CC) myRHEA.cpp $(GUIFLAGS) -std=c++11 $(OBJ) -o $(OUT) $(CFLAGS)
	$(CC) myRHEA.cpp $(OBJ) -L/usr/lib/x86_64-linux-gnu/hdf5/openmpi -I/usr/include/hdf5/openmpi -o $(OUT) $(CFLAGS) -lyaml-cpp

clean:
	rm $(OUT)


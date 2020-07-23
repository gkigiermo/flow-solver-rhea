CC=mpic++
CFLAGS=-O3 
OBJ=src/*.cpp
OUT= test1

all:
	$(CC) heat.cpp $(OBJ) -o $(OUT) $(CFLAGS)

clean:
	rm  src/*.o *.info


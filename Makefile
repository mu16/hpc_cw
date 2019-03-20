default: all

main.o: main.cpp
	mpicxx -std=c++11 -Wall -O2 -o main.o -c main.cpp

model.o: Model.cpp Model.h
	mpicxx -std=c++11 -Wall -o model.o -c Model.cpp

burgers.o: burgers.cpp burgers.h Model.h
	mpicxx -std=c++11 -Wall -o burgers.o -c burgers.cpp

compile: main.o  model.o burgers.o
	mpicxx -o my_prog main.o  model.o burgers.o -O3 -ffast-math -funroll-loops -march=native -ftree-vectorize

.PHONY: clean # Specify that ’clean’ is not a real file
	target

diff: compile
	mpiexec -np  5 my_prog 0 0 0 1 5 1
	
advx: compile
	mpiexec -np  5 my_prog 1 0 0 0 5 1
	
advy: compile
	mpiexec -np 16 my_prog 0 1 0 0 4 4

burg: compile
	mpiexec -np  36 my_prog 1.0 0.5 1.0 0.02 6 6

clean:
	rm -f *.o my_prog  

all: diff clean

CXXFILES=r_dice.c
CXXFLAGS=-O3 -o DiCE -fopenmp -Wall -Wconversion -fPIC
LIBS=-lm

all:
	g++ -c svm.cpp -fopenmp
	gcc -c dice.c xdrfile.c xdrfile_trr.c xdrfile_xtc.c -fopenmp -lm
	g++ $(CXXFLAGS) *.o

clean:
	rm -f prog *.o

test_flag:
	./DiCE -f1 examples/input1.pdb -f2 examples/input2.pdb -t $(t) 

test_pdb:
	./DiCE examples/input1.pdb examples/input2.pdb

test_trr:
	./DiCE -f1 examples/input1.trr -f2 examples/input2.trr

test_xtc:
	./DiCE -f1 examples/input1.xtc -f2 examples/input2.xtc

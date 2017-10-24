CC=g++-5
CFLAGS= -std=c++11 -g -O3 -larmadillo -llapack -lblas 
OBJ = hartree.o
DEPS = mass.h 
all: craw3.exe

craw3.exe: $(OBJ)
	$(CC) -o craw3.exe $^ $(CFLAGS)

%.o: %.cpp $(DEP)
	$(CC) -c $< $(CFLAGS)

clean:
	rm *.o *.exe

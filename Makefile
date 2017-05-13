# Note that this code only works for odd total number of photon dressing

CC = g++
CFLAGS = -Ofast -std=c++11 
LIBDIR =
EXTRALIBS = -lm
INCLUDEDIR = -I/usr/include/eigen3
OBJECTS = floquet.o dynamics.o parameters.o floquet_bloch_kernel.o
BINARYS =

all: $(OBJECTS) $(BINARYS)

%.o : %.cpp %.h
	$(CC) -c $(CFLAGS) $(INCLUDEDIR) $< -o $@

floquet_bloch_kernel.test : floquet_bloch_kernel.cpp floquet_bloch_kernel.h e_grid.o dynamics.o floquet.o parameters.o
	$(CC) $(CFLAGS) $(INCLUDEDIR) $< e_grid.o dynamics.o floquet.o parameters.o -DTEST -o $@ $(EXTRALIBS)

floquet.test : floquet.cpp floquet.h
	$(CC) $(CFLAGS) $(INCLUDEDIR) $< -DTEST -o $@ $(EXTRALIBS)

dynamics.test : dynamics.cpp dynamics.h floquet.o
	$(CC) $(CFLAGS) $(INCLUDEDIR) $< floquet.o -DTEST -o $@ $(EXTRALIBS)

e_grid.test : e_grid.cpp e_grid.h dynamics.o floquet.o
	$(CC) $(CFLAGS) $(INCLUDEDIR) $< dynamics.o floquet.o -DTEST -o $@ $(EXTRALIBS)

#SFC: Serial Fortran Compiler
#PFC: Parallel Fortran Compiler

# --- macOS --- #
#SFC=gfortran
#PFC=mpif90
#CFLAGS=-O3
#MPIEXEC=mpirun
#INCF=-J

# -- CIRRUS -- #
SFC=ifort
PFC=mpif90
CFLAGS=-O3
MPIEXEC=mpirun
INCF=-module

#Preprocessing
CC=gcc
DEFINESERIAL= -D SERIAL
PREP=-cpp -E -P -o

# PROJECT VARIABLES
EXECUTABLE=imagempi
SERIALSRC=src/pgmio.f90 src/serial.f90
PARALLELSRC= src/pgmio.f90 src/parallel.f90
SERIALOBJ=$(patsubst src/%.f90,binserial/%.o,$(SERIALSRC))
PARALLELOBJ=$(patsubst src/%.f90,bin/%.o,$(PARALLELSRC))
INC=mod

# COMPILATION RULES
# by default create mpi version
all: parallel 

# Files which needs preprocessing (C-style preprocessor used)
# precisiondef must be always compiled to create the correct .mod
# when changing between serial and parallel
.PHONY: bin/precision.o binserial/precision.o

bin/precision.o: src/precision.f90
		$(CC) $< $(PREP) bin/precision.f90
		$(PFC) bin/precision.f90 $(CFLAGS) $(INCF) $(INC) -c -o $@

bin/imagempi.f90: src/imagempi.f90
		$(CC) $< $(PREP) $@

binserial/precision.o: src/precision.f90
		$(CC) $(DEFINESERIAL) $< $(PREP) binserial/precision.f90
		$(SFC)  binserial/precision.f90 $(CFLAGS) $(INCF) $(INC) -c -o $@

binserial/imagempi.f90: src/imagempi.f90
		$(CC) $(DEFINESERIAL) $< $(PREP) $@

# compile other objects
bin/%.o: src/%.f90 Makefile
		$(PFC) $(CFLAGS) $(INCF) $(INC) $< -c -o $@

binserial/%.o: src/%.f90 Makefile
		$(SFC) $(CFLAGS) $(INCF) $(INC) $< -c -o $@

# create executables
serial: binserial/imagempi.f90 binserial/precision.o $(SERIALOBJ) Makefile
		$(SFC) binserial/precision.o $(SERIALOBJ) $< $(CFLAGS) $(INCF) $(INC) -o binserial/$(EXECUTABLE)

parallel: bin/imagempi.f90 bin/precision.o $(PARALLELOBJ) Makefile
		$(PFC) bin/precision.o $(PARALLELOBJ) $< $(CFLAGS) $(INCF) $(INC) -o bin/$(EXECUTABLE)


# OTHER ACTIONS (WITHOUT DEPENDENCIES)
.PHONY: test mpirun run clean

test:
		./test/correctness/ctest.sh bin/imagempi test/inputs test/morarserial

runmpi:
		$(MPIEXEC) -n 8 ./bin/$(EXECUTABLE) test/inputs/edge768x768.pgm 100 10

runserial:
		./binserial/$(EXECUTABLE) test/inputs/edge768x768.pgm 100 10

clean:
		rm -rf bin/*.o bin/*.f90 binserial/*.f90 binserial/*.o mod/*.mod bin/$(EXECUTABLE) binserial/$(EXECUTABLE)

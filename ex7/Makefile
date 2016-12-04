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
SERIALOBJ=$(patsubst src/%.f90,build_serial/%.o,$(SERIALSRC))
PARALLELOBJ=$(patsubst src/%.f90,build/%.o,$(PARALLELSRC))
INC=mod

# Parallel as default version 
all: parallel 

.PHONY: build/precision.o build_serial/precision.o

build/precision.o: src/precision.f90
		$(CC) $< $(PREP) build/precision.f90
		$(PFC) build/precision.f90 $(CFLAGS) $(INCF) $(INC) -c -o $@

build/imagempi.f90: src/imagempi.f90
		$(CC) $< $(PREP) $@

build_serial/precision.o: src/precision.f90
		$(CC) $(DEFINESERIAL) $< $(PREP) build_serial/precision.f90
		$(SFC)  build_serial/precision.f90 $(CFLAGS) $(INCF) $(INC) -c -o $@

build_serial/imagempi.f90: src/imagempi.f90
		$(CC) $(DEFINESERIAL) $< $(PREP) $@

# compile other objects
build/%.o: src/%.f90 Makefile
		$(PFC) $(CFLAGS) $(INCF) $(INC) $< -c -o $@

build_serial/%.o: src/%.f90 Makefile
		$(SFC) $(CFLAGS) $(INCF) $(INC) $< -c -o $@

# create executables
serial: build_serial/imagempi.f90 build_serial/precision.o $(SERIALOBJ) Makefile
		$(SFC) build_serial/precision.o $(SERIALOBJ) $< $(CFLAGS) $(INCF) $(INC) -o build_serial/$(EXECUTABLE)

parallel: build/imagempi.f90 build/precision.o $(PARALLELOBJ) Makefile
		$(PFC) build/precision.o $(PARALLELOBJ) $< $(CFLAGS) $(INCF) $(INC) -o build/$(EXECUTABLE)


# OTHER ACTIONS (WITHOUT DEPENDENCIES)
.PHONY: clean

clean:
		rm -rf build/*.o build/*.f90 build_serial/*.f90 build_serial/*.o mod/*.mod build/$(EXECUTABLE) build_serial/$(EXECUTABLE)

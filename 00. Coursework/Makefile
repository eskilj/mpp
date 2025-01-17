# MODIFICABLE VARIABLES


#CIRRUS
SFC=pgf90
PFC=mpif90
#PFC=vtfort -vt:f90 mpif90
CFLAGS=-O3
MPIEXEC=mpirun
INCF=-module

#Preprocessing
CC=gcc
DEFINESERIAL= -D SERIAL
PREP=-cpp -E -P -o

# PROJECT VARIABLES
EXECUTABLE=image
SERIALSRC=src/pgmio.f90 src/serial
PARALLELSRC= src/pgmio.f90 src/parallel
SERIALOBJ=$(patsubst src/%.f90,binserial/%.o,$(SERIALSRC))
PARALLELOBJ=$(patsubst src/%.f90,bin/%.o,$(PARALLELSRC))
INC=include

# COMPILATION RULES
# by default create mpi version
all: mpiversion

# Files which needs preprocessing (C-style preprocessor used)
# precisiondef must be always compiled to create the correct .mod
# when changing between serialversion and mpiversion
.PHONY: bin/precision.o binserial/precision.o

bin/precisiondef.o: src/precision.f90
		$(CC) $< $(PREP) bin/precision.f90
		$(PFC) bin/precision.f90 $(CFLAGS) $(INCF) $(INC) -c -o $@

bin/image.f90: src/image.f90
		$(CC) $< $(PREP) $@

binserial/precision.o: src/precision.f90
		$(CC) $(DEFINESERIAL) $< $(PREP) binserial/precision.f90
		$(SFC)  binserial/precision.f90 $(CFLAGS) $(INCF) $(INC) -c -o $@

binserial/image.f90: src/image.f90
		$(CC) $(DEFINESERIAL) $< $(PREP) $@

# compile other objects
bin/%.o: src/%.f90 Makefile
		$(PFC) $(CFLAGS) $(INCF) $(INC) $< -c -o $@

binserial/%.o: src/%.f90 Makefile
		$(SFC) $(CFLAGS) $(INCF) $(INC) $< -c -o $@

# create executables
serialversion: binserial/image.f90 binserial/precision.o $(SERIALOBJ) Makefile
		$(SFC) binserial/precision.o $(SERIALOBJ) $< $(CFLAGS) $(INCF) $(INC) -o binserial/$(EXECUTABLE)

mpiversion: bin/invertedges.f90 bin/precision.o $(PARALLELOBJ) Makefile
		$(PFC) bin/precision.o $(PARALLELOBJ) $< $(CFLAGS) $(INCF) $(INC) -o bin/$(EXECUTABLE)


# OTHER ACTIONS (WITHOUT DEPENDENCIES)
.PHONY: test mpirun run clean

test:
		./test/correctness/ctest.sh bin/image test/inputs test/morarserial

runmpi:
		$(MPIEXEC) -n 8 ./bin/$(EXECUTABLE) test/inputs/edge768x768.pgm 100 10

runserial:
		./binserial/$(EXECUTABLE) test/inputs/edge768x768.pgm 100 10

clean:
		rm -rf bin/*.o bin/*.f90 binserial/*.f90 binserial/*.o include/*.mod bin/$(EXECUTABLE) binserial/$(EXECUTABLE)

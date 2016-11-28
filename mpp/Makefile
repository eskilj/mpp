

# MODIFICABLE VARIABLES

## mylaptop openMPI
#SFC=gfortran
#PFC=mpif90
#CFLAGS=-Wall -g -fbounds-check -O3
#MPIEXEC=mpiexec
#INCF=-J

## mylaptop mpichv2
#SFC=gfortran
#PFC=/opt/mpich/bin/mpif90
#CFLAGS=-Wall -g -fbounds-check -O3
#MPIEXEC=/opt/mpich/bin/mpiexec
#INCF=-J


#MORAR gnu
#FCS=gfortran
#FCP=/opt/mpich2-gnu/bin/mpif90
#CFLAGS=-O3
#MPIEXEC=/opt/mpich2-gnu/bin/mpiexec
#INCF=-J

#MORAR pgi
SFC=pgf90
PFC=mpif90
#PFC=vtfort -vt:f90 mpif90
CFLAGS=-O3
MPIEXEC=mpiexec
INCF=-module

#Preprocessing
CC=gcc
DEFINESERIAL= -D SERIALVERSION
PREP=-cpp -E -P -o

# PROJECT VARIABLES
EXECUTABLE=invertedges
SERIALSRC=src/pgmio.f90 src/ieserial.f90
PARALLELSRC= src/pgmio.f90 src/iempi.f90
SERIALOBJ=$(patsubst src/%.f90,binserial/%.o,$(SERIALSRC))
PARALLELOBJ=$(patsubst src/%.f90,bin/%.o,$(PARALLELSRC))
INC=include

# COMPILATION RULES
# by default create mpi version
all: mpiversion

# Files which needs preprocessing (C-style preprocessor used)
# precisiondef must be always compiled to create the correct .mod
# when changing between serialversion and mpiversion
.PHONY: bin/precisiondef.o binserial/precisiondef.o

bin/precisiondef.o: src/precisiondef.F90
	$(CC) $< $(PREP) bin/precisiondef.f90
	$(PFC) bin/precisiondef.f90 $(CFLAGS) $(INCF) $(INC) -c -o $@

bin/invertedges.f90: src/invertedges.F90
	$(CC) $< $(PREP) $@

binserial/precisiondef.o: src/precisiondef.F90
	$(CC) $(DEFINESERIAL) $< $(PREP) binserial/precisiondef.f90
	$(SFC)  binserial/precisiondef.f90 $(CFLAGS) $(INCF) $(INC) -c -o $@

binserial/invertedges.f90: src/invertedges.F90
	$(CC) $(DEFINESERIAL) $< $(PREP) $@



# compile other objects
bin/%.o: src/%.f90 Makefile
	$(PFC) $(CFLAGS) $(INCF) $(INC) $< -c -o $@

binserial/%.o: src/%.f90 Makefile
	$(SFC) $(CFLAGS) $(INCF) $(INC) $< -c -o $@


# create executables
serialversion: binserial/invertedges.f90 binserial/precisiondef.o $(SERIALOBJ) Makefile
	$(SFC) binserial/precisiondef.o $(SERIALOBJ) $< $(CFLAGS) $(INCF) $(INC) -o binserial/$(EXECUTABLE)

mpiversion: bin/invertedges.f90 bin/precisiondef.o $(PARALLELOBJ) Makefile
	$(PFC) bin/precisiondef.o $(PARALLELOBJ) $< $(CFLAGS) $(INCF) $(INC) -o bin/$(EXECUTABLE)


# OTHER ACTIONS (WITHOUT DEPENDENCIES)
.PHONY: test mpiexec run clean

test:
	./test/correctness/ctest.sh bin/invertedges test/inputs test/morarserial

runmpi:
	$(MPIEXEC) -n 8 ./bin/$(EXECUTABLE) test/inputs/edge768x768.pgm 100 10

runserial:
	./binserial/$(EXECUTABLE) test/inputs/edge768x768.pgm 100 10

clean:
	rm -rf bin/*.o bin/*.f90 binserial/*.f90 binserial/*.o include/*.mod bin/$(EXECUTABLE) binserial/$(EXECUTABLE)

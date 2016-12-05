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

VPATH=src
BUILD_DIR = build

# PROJECT VARIABLES
EXECUTABLE=imagempi
SERIALSRC=src/pgmio.f90 src/serial.f90
PARALLELSRC= $(VPATH)/pgmio.f90 src/parallel.f90
SERIALOBJ=$(patsubst src/%.f90,build_serial/%.o,$(SERIALSRC))
PARALLELOBJ=$(patsubst src/%.f90,build/%.o,$(PARALLELSRC))
INC=mod

info :
	@echo "The following are possible:"
	@echo "  make parallel"
	@echo "  make serial"
	@echo "  make clean"

.PHONY: $(BUILD_DIR)/precision.o build_serial/precision.o

$(BUILD_DIR)/precision.o: src/precision.f90
		$(CC) $< $(PREP) $(BUILD_DIR)/precision.f90
		$(PFC) $(BUILD_DIR)/precision.f90 $(CFLAGS) $(INCF) $(INC) -c -o $@

$(BUILD_DIR)/imagempi.f90: src/imagempi.f90
		$(CC) $< $(PREP) $@ | $(BUILD_DIR)

build_serial/precision.o: src/precision.f90
		$(CC) $(DEFINESERIAL) $< $(PREP) build_serial/precision.f90
		$(SFC)  build_serial/precision.f90 $(CFLAGS) $(INCF) $(INC) -c -o $@

build_serial/imagempi.f90: src/imagempi.f90
		$(CC) $(DEFINESERIAL) $< $(PREP) $@

# Compile
$(BUILD_DIR)/%.o: src/%.f90 Makefile
		$(PFC) $(CFLAGS) $(INCF) $(INC) $< -c -o $@

build_serial/%.o: src/%.f90 Makefile
		$(SFC) $(CFLAGS) $(INCF) $(INC) $< -c -o $@

# Make the programs
serial: build_serial/imagempi.f90 build_serial/precision.o $(SERIALOBJ) Makefile
		$(SFC) build_serial/precision.o $(SERIALOBJ) $< $(CFLAGS) $(INCF) $(INC) -o build_serial/$(EXECUTABLE)

parallel: $(BUILD_DIR)/imagempi.f90 $(BUILD_DIR)/precision.o $(PARALLELOBJ) Makefile
		$(PFC) $(BUILD_DIR)/precision.o $(PARALLELOBJ) $< $(CFLAGS) $(INCF) $(INC) -o $(BUILD_DIR)/$(EXECUTABLE)

$(BUILD_DIR):
	mkdir $(BUILD_DIR)

# OTHER ACTIONS (WITHOUT DEPENDENCIES)
.PHONY: clean

clean:
		rm -rf $(BUILD_DIR)/*.o $(BUILD_DIR)/*.f90 build_serial/*.f90 build_serial/*.o mod/*.mod build/$(EXECUTABLE) build_serial/$(EXECUTABLE)

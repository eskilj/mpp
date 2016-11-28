MF=		Makefile

FC=		mpif90
FFLAGS=	-O3

LFLAGS=

EXE=	imagempi

SRC= \
		precisiondef.f90 \
		iempi.f90 \
		ieserial.f90 \
		invertedges.f90 \
		pgmio.f90

#
# No need to edit below this line
#

.SUFFIXES:
.SUFFIXES: .f90 .o

OBJ=	$(SRC:.f90=.o)

.f90.o:
		$(FC) $(FFLAGS) -c $<

all:	$(EXE)

$(EXE):	$(OBJ)
		$(FC) $(FFLAGS) -o $@ $(OBJ) $(LFLAGS)

$(OBJ):	$(MF)

imagempi.o:		pgmio.o

clean:
		rm -f $(OBJ) $(EXE) core

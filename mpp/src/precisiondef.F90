MODULE precisiondef
#ifndef SERIALVERSION
  use mpi
#endif
  implicit none
  integer, parameter :: REALNUMBER = kind(1.0e0)
#ifndef SERIALVERSION
  integer, parameter :: MPI_REALNUMBER = MPI_REAL
#endif

END MODULE precisiondef


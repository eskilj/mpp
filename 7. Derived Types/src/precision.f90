MODULE precision
#ifndef SERIAL
  use mpi
#endif
  implicit none
  integer, parameter :: REALNUMBER = kind(1.0e0)
#ifndef SERIAL
  integer, parameter :: MPI_REALNUMBER = MPI_REAL
#endif

END MODULE precision


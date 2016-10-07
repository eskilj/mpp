program hello

  use mpi

  implicit none

  integer :: ierr, rank, size

  call MPI_Init(ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
  call MPI_COMM_SIZE( MPI_COMM_WORLD, size, ierr)

  if(rank.eq.0) then
  print *,  "Rank :", rank, "size : ", size 
  endif
  !call MPI_Finalize(ierr)

end program hello

program ping
use mpi
implicit none

  integer :: ierr, rank, size, i, N
  double precision :: t_start, t_stop, t_diff
  integer, dimension(5,5) :: message
  integer, dimension(MPI_STATUS_SIZE) :: status
  N = 1
  message(:,:) = 10
  call MPI_Init(ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
  call MPI_COMM_SIZE( MPI_COMM_WORLD, size, ierr)
  call MPI_Barrier(MPI_COMM_WORLD, ierr)
  t_start = MPI_Wtime()
  do i = 1, N
  if(rank.eq.0) then
   call MPI_SSEND(message, 1, MPI_INTEGER, 1, 0, MPI_COMM_WORLD, ierr)
   !print *, "Loop :", i, "Send :", message
   call MPI_RECV(message, 1, MPI_INTEGER, 1, 0, MPI_COMM_WORLD, status, ierr)
   print *, "Loop :", i, "Recv :", message, "MASTER"
  endif

  if(rank.eq.1) then
    call MPI_RECV(message, 1, MPI_REAL, 0, 0, MPI_COMM_WORLD, status, ierr)
    !print *, "Loop :", i, "Recv :", message
    !message(1:1) = i * 100
    call MPI_SSEND(message, 1, MPI_INTEGER, 0, 0, MPI_COMM_WORLD, ierr)
  endif
  end do

  call MPI_Barrier(MPI_COMM_WORLD, ierr)
  t_stop = MPI_WTime()
  t_diff = t_stop - t_start
  if(rank.eq.0) then
     print *, t_diff
     print *, size
     print *, message
     print *, "arraysize :", sizeof(message)
  endif
  call MPI_Finalize(ierr)

end program ping

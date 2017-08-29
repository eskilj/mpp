module sums
implicit none
contains
  
  function partialSum(rank, size, N)result(part_sum)
    implicit none
    integer :: rank, size, i, N
    real :: p_sum, part_sum
    p_sum = 0.0
    do 10 i = lowerLimit(rank, size, N), upperLimit(rank, size, N)
       p_sum = p_sum + (1/(1 + ((i-0.5)/N)**2))
    10 continue
    part_sum = 4 * p_sum / N
  end function partialSum

  function upperLimit(rank, size, N)result(upper_limit)
    implicit none
    integer :: upper_limit, rank, size, N
    upper_limit = N * (rank + 1) / size
  end function upperLimit

  function lowerLimit(rank, size, N)result(lower_limit)
    implicit none
    integer :: lower_limit, rank, size, N
    lower_limit = 1+(N*(rank)/size)
  end function lowerLimit

end module sums

program pi
use mpi
use sums
implicit none

  integer :: ierr, rank, size, size_max, i, N, message, calc_pi_n, pi_i
  real :: p_sum, t_sum, calc_sum
  double precision :: t_start, t_stop, t_diff
  real, dimension(:), allocatable :: parallel_sums
  integer, dimension(MPI_STATUS_SIZE) :: status
  N = 840
  size_max = 8
  calc_pi_n = 100000

  call MPI_Init(ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
  call MPI_COMM_SIZE( MPI_COMM_WORLD, size, ierr)
  allocate(parallel_sums(size))
  call MPI_Barrier(MPI_COMM_WORLD, ierr)
  t_start = MPI_Wtime()
  do pi_i = 1, calc_pi_n
  if(size.gt.size_max) then
     size = size_max
     !print *, "Please use maximum ", size_max, " proccesses "
  endif

  if(rank.eq.0) then
    t_sum = 0.0
    !do i = 1, size-1
      !call MPI_SSEND(i*100, 1, MPI_INTEGER, i, 0, MPI_COMM_WORLD, ierr)
    !end do

    parallel_sums(1) = partialSum(rank, size, N)

    do i = 1, size-1
      !print *, "RECV ", i
      call MPI_RECV(calc_sum, 1, MPI_REAL, i, 0, MPI_COMM_WORLD, status, ierr)
      !print *, "part sum :", calc_sum, "Rank :", rank
      parallel_sums(i+1) = calc_sum
    end do
    t_sum = SUM(parallel_sums)
    !print *, "RESULTS: ", parallel_sums
    !print *, t_sum
  endif

  if(rank.gt.0) then
      !call MPI_RECV(message, 1, MPI_INTEGER, 0, 0, MPI_COMM_WORLD, status, ierr)
      !print *, "message :", message, "Rank :", rank
      p_sum = partialSum(rank, size, N)
      !print *, "P Sum : ",  p_sum
      call MPI_SSEND(p_sum, 1, MPI_REAL, 0, 0, MPI_COMM_WORLD, ierr)
  endif
  end do
  call MPI_Barrier(MPI_COMM_WORLD, ierr)
  t_stop = MPI_WTime()
  t_diff = t_stop - t_start
  if(rank.eq.0) then
     print *, t_diff
     print *, size
     print *, t_sum
  endif
  deallocate(parallel_sums)
  call MPI_Finalize(ierr)

end program pi

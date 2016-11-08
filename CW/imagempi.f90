!
! A simple solution to the Case Study exercise from the EPCC MPI course.
!
! Communications is done using the sendrecv routine; a proper solution
! would use non-blocking communications (ie some combination of issend/recv
! and ssend/irecv).
!
! Note that the special rank of MPI_PROC_NULL is a "black hole" for
! communications (similar to /dev/null in Unix). The use of this value
! for processes off the edges of the image means we do not need any
! additional logic to ensure that processes at the edges do not attempt
! to send to or receive from invalid ranks (ie rank = -1 and rank = P).
!

program casestudy

  use mpi
  use pgmio
 
  implicit none

  integer :: M, N, MP, NP

  integer, parameter :: P = 4

  integer, parameter :: MAXITER   = 1500
  integer, parameter :: PRINTFREQ =  200

  real, dimension(:,:), allocatable :: new, old, edge, buf, masterbuf

  real :: val, boundaryval

  integer, parameter :: maxlen = 32

  character*(maxlen) :: filename

  integer :: i, j, iter

  integer :: comm, rank, size, ierr, uprank, dnrank

  integer, dimension(MPI_STATUS_SIZE) :: status

  filename = 'edgenew768x768.pgm'

  call pgmsize(filename, M, N)

  MP = M
  NP = N/P

  allocate(new(0:MP+1, 0:NP+1))
  allocate(edge(0:MP+1, 0:NP+1))
  allocate(old(0:MP+1, 0:NP+1))

  allocate(buf(MP, NP))
  allocate(masterbuf(M, N))

  call mpi_init(ierr)

  comm = MPI_COMM_WORLD

  call mpi_comm_size(comm, size, ierr)
  call mpi_comm_rank(comm, rank, ierr)

  if (size .ne. P) then
    if (rank .eq. 0) write(*,*) 'ERROR: size = ', size, ', P = ', P
    call mpi_finalize(ierr)
    stop
  end if

  if (rank .eq. 0) then

    write(*,*) 'Processing ', M, ' x ' , N, ' image on ', P, ' processes'
    write(*,*) 'Number of iterations = ', MAXITER

    write(*,*)
    write(*,*) 'Reading ', filename
    write(*,*)

    call pgmread(filename, masterbuf)

  end if

  call mpi_scatter(masterbuf, MP*NP, MPI_REAL, &
                   buf,       MP*NP, MPI_REAL, &
                   0, comm, ierr)

  do j = 1, NP
    do i = 1, MP

      edge(i,j) = buf(i,j)

    end do
  end do


  do j = 0, NP+1
    do i = 0, MP+1

      old(i,j) = 255.0

    end do
  end do

! Set fixed boundary conditions on the left and right sides

  do j = 1, NP

! compute sawtooth value

    val = boundaryval(j, NP)

    old(0,   j) = 255.0*val
    old(MP+1,j) = 255.0*(1.0-val)

  end do

  uprank = rank + 1
  if (uprank .ge. P) uprank = MPI_PROC_NULL

  dnrank = rank - 1
  if (dnrank .lt. 0) dnrank = MPI_PROC_NULL

  do iter = 1, MAXITER

    if (mod(iter,PRINTFREQ) .eq. 0)  then

      if (rank .eq. 0) write(*,*) 'Iteration ', iter

    end if

! Implement periodic boundary conditions on top and bottom sides

    do i = 1, MP

      old(i,  0) = old(i,NP)
      old(i,NP+1) = old(i,1)

    end do

    call mpi_sendrecv(old(1,NP), MP, MPI_REAL, uprank, 1, &
                      old(1,0),  MP, MPI_REAL, dnrank, 1, &
                      comm, status, ierr)
    call mpi_sendrecv(old(1,1),    MP, MPI_REAL, dnrank, 1, &
                      old(1,NP+1), MP, MPI_REAL, uprank, 1, &
                      comm, status, ierr)

    do j = 1, NP
      do i = 1, MP

        new(i,j) = 0.25*(old(i+1,j)+old(i-1,j)+old(i,j+1)+old(i,j-1) &
                         - edge(i,j))

      end do
    end do

    do j = 1, NP
      do i = 1, MP

        old(i,j) = new(i,j)

      end do
    end do

  end do

  if (rank .eq. 0) then
    write(*,*)
    write(*,*) 'Finished ', iter-1, ' iterations'
  end if

!  Gather data

  do j = 1, NP
    do i = 1, MP

      buf(i,j) = old(i,j)

    end do
  end do

  call mpi_gather(buf,       MP*NP, MPI_REAL, &
                  masterbuf, MP*NP, MPI_REAL, &
                  0, comm, ierr)

  if (rank .eq. 0) then

    filename='imagenew768x768.pgm'
    write(*,*)
    write(*,*) 'Writing ', filename
    call pgmwrite(filename, masterbuf)

  end if

  call mpi_finalize(ierr)

  deallocate(new)
  deallocate(old)
  deallocate(edge)
  deallocate(buf)
  deallocate(masterbuf)

end program casestudy

real function boundaryval(i, m)

  implicit none

  integer :: i, m
  real :: val

  val = 2.0*float(i-1)/float(m-1)
  if (i .ge. m/2+1) val = 2.0-val

  boundaryval = val

end function boundaryval

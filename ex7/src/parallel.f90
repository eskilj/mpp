MODULE parallel

  use mpi
  use precision

  implicit none

  ! Define time datatye
  type timetype
      real(kind=8) :: value
  end type timetype
  
  ! Constants
  integer, parameter :: ROOT = 0
  integer, parameter :: N_DIMS = 2 !Num of dimentions

  ! MPI VARIABLES
  integer :: comm, size, cartcomm, rank, ierr, errorcode
  integer :: n_left, n_right, n_up, n_down !process neighbours
  integer, dimension(8) :: request
  
  ! Problem parameters
  integer, dimension(N_DIMS) :: dims
  integer :: MP, NP ! Local array sizes

  ! MPI NEW DATATYPES
  integer :: MASTER_BLOCK_T, BLOCK_T, V_HALO_T, H_HALO_T
  integer, dimension(:), allocatable :: send_counts, displacements

contains

  !  --------------- GLOBAL MPI -------------------------!

  subroutine par_init()
    ! Initialise MPI
    call MPI_INIT(ierr)
    comm = MPI_COMM_WORLD
    call MPI_COMM_SIZE(comm, size, ierr)
  end subroutine par_init

  subroutine par_finalize()
    ! Free the used resources, finalize MPI
    call MPI_TYPE_FREE(MASTER_BLOCK_T,ierr)
    call MPI_TYPE_FREE(BLOCK_T,ierr)
    call MPI_TYPE_FREE(H_HALO_T,ierr)
    call MPI_TYPE_FREE(V_HALO_T,ierr)
    call MPI_Finalize(ierr)
  end subroutine par_finalize

  logical function par_isroot()
    ! Return true when executen in the ROOT process
    par_isroot = rank == ROOT
    return
  end function par_isroot

  !  --------------- GLOBAL MPI METHODS -------------------------!

  subroutine par_domain_decomposition_2D(nx,ny,npx,npy)

    integer, intent(in) :: nx, ny
    integer, intent(out) :: npx, npy
    integer :: y_dir, x_dir, disp
    character(len=100) :: message
    logical, dimension(N_DIMS) :: periods
    logical :: reorder
    
    ! Create cartesian topology 2D
    dims(:) = 0
    periods(:) = .false.
    reorder = .true.

    call MPI_DIMS_CREATE(size, N_DIMS, dims, ierr)
    call MPI_CART_CREATE(comm, N_DIMS, (/dims(2), dims(1) /), periods, reorder, cartcomm, ierr)
    call MPI_COMM_RANK(cartcomm, rank, ierr)

    y_dir = 0
    x_dir = 1
    disp = 1

    ! Get neighbours
    call MPI_CART_SHIFT(cartcomm, y_dir, disp, n_down, n_up, ierr)
    call MPI_CART_SHIFT(cartcomm, x_dir, disp, n_left, n_right, ierr)

    ! Compute the new array dimensions
    if ((mod(nx, dims(1)) .ne. 0) .or. (mod(ny, dims(2)) .ne. 0)) then
        write(message,*) nx," is not divisible in ",dims(1)," parts"
        call exit_all("Could not decompose domain!"//message)
    end if

    npx = nx/dims(1)
    npy = ny/dims(2)
    MP = npx
    NP = npy

    ! Print out the information about the domain decomposition
    write(message,'(I2,A3,I2,A19,I4,A3,I4)') &
        dims(1), " x ", dims(2), "with each block of", MP, " x ", NP
    call print_once("Domain decomposed in:")
    call print_once(" ->  Grid of size: "//message)

    ! Create the derived datatypes
    call create_types()
  end subroutine

  subroutine create_types()
    ! Create all the derived datatypes used in the program, they are:
    ! Block type, master block type, vertical halo and horitzontal halo
    integer, dimension(N_DIMS) :: sizes, subsizes, starts
    integer(kind=mpi_address_kind) :: start, extent, lb, realextent
    integer :: AllocateStatus, i, base, LONG_BLOCK_T
    
    ! Block type: space in local arrays, inside halos
    sizes    = (/ MP+2, NP+2 /)
    subsizes = (/ MP , NP /)
    starts(:)   = 1
    call MPI_TYPE_CREATE_SUBARRAY(N_DIMS, sizes, subsizes, starts, &
             MPI_ORDER_FORTRAN, MPI_REALNUMBER, BLOCK_T, ierr)
    
    ! Master block type: distribution unit from master to working
    ! processes, needs a extent resize and the proper send_counts and
    ! displacements in order to be accessed iteratively.
    sizes    = (/ MP*dims(1), NP*dims(2) /)
    subsizes = (/ MP, NP /)
    starts(:)   = 0
    call MPI_TYPE_CREATE_SUBARRAY(N_DIMS, sizes, subsizes, starts, &
             MPI_ORDER_FORTRAN, MPI_REALNUMBER, LONG_BLOCK_T, ierr)
    call MPI_TYPE_GET_EXTENT(MPI_REALNUMBER, lb, realextent, ierr)

    start = 0
    extent = MP*realextent
    call MPI_TYPE_CREATE_RESIZED(LONG_BLOCK_T, start, extent, &
             MASTER_BLOCK_T,ierr)
  
    allocate(send_counts(size), displacements(size), STAT=AllocateStatus)
    if(allocateStatus .ne. 0) call exit_all("*** NOT enough memory ***")
    
    base = 1
    do i= 1, size
        send_counts(i) = 1
        displacements(i) = (base-1) + mod(i-1,dims(1))
        if(mod(i,dims(1))==0) base = base + NP * dims(1)
    end do

    ! HALO VECTORS (HALO-HALO INTERSECTIONS ARE NOT NEDDED)
    ! Horitzontal halo data is contiguous, defined to maintain
    ! code cohesion in all halo swaps.
    call MPI_TYPE_VECTOR(NP, 1 , MP+2, MPI_REALNUMBER, V_HALO_T, ierr)
    call MPI_TYPE_VECTOR(1 , MP, MP  , MPI_REALNUMBER, H_HALO_T, ierr)
    
    ! COMMIT NEW MPI DATATYPES
    call MPI_TYPE_COMMIT(MASTER_BLOCK_T, ierr)
    call MPI_TYPE_COMMIT(BLOCK_T,  ierr)
    call MPI_TYPE_COMMIT(V_HALO_T, ierr)
    call MPI_TYPE_COMMIT(H_HALO_T, ierr)
  end subroutine

  subroutine par_scatter(source, dest)
    real(kind=REALNUMBER), dimension(:,:), intent(in) :: source
    real(kind=REALNUMBER), dimension(:,:), intent(out) :: dest
    call MPI_Scatterv(source, send_counts, displacements, MASTER_BLOCK_T, &
                      dest, MP*NP, MPI_REALNUMBER, 0, cartcomm,ierr)
  end subroutine par_scatter

  subroutine par_Gather(source, dest)
    real(kind=REALNUMBER), dimension(0:,0:), intent(in) :: source
    real(kind=REALNUMBER), dimension(:,:), intent(out) :: dest
    call MPI_GATHERV(source, 1, BLOCK_T, dest, send_counts, displacements, &
                     MASTER_BLOCK_T, 0, cartcomm, ierr)
  end subroutine


  subroutine par_HalosSwap(old)
    ! Non-blocking send and reveive of all the halos, afterwards the
    ! par_WaitHalos() routine should be called to ensure these communications
    ! are completed.
    real(kind=REALNUMBER), dimension(0:,0:), intent(in) :: old
   
    call MPI_Issend(old(MP,1),1, V_HALO_T, n_right ,0,cartcomm,request(1),ierr)
    call MPI_Issend(old(1,1) ,1, V_HALO_T, n_left   ,0,cartcomm,request(3),ierr)
    call MPI_Issend(old(1,1),1, H_HALO_T, n_down,0,cartcomm,request(7),ierr)
    call MPI_Issend(old(1,NP) ,1, H_HALO_T, n_up ,0,cartcomm,request(5),ierr)
    
    call MPI_Irecv(old(MP+1,1),1, V_HALO_T, n_right ,0,cartcomm,request(4),ierr)
    call MPI_Irecv(old(0,1)   ,1, V_HALO_T, n_left ,0,cartcomm,request(2),ierr)
    call MPI_Irecv(old(1,0),1, H_HALO_T, n_down, 0,cartcomm,request(6),ierr)
    call MPI_Irecv(old(1,NP+1)   ,1, H_HALO_T, n_up ,0,cartcomm,request(8),ierr)
    
  end subroutine par_HalosSwap

  subroutine par_WaitHalos()
    integer, dimension(MPI_STATUS_SIZE,8) :: status
    call MPI_Waitall(8,request,status,ierr)
  end subroutine par_WaitHalos

  !Calculate maximum change across the image
  real(kind=REALNUMBER) function par_calc_max_diff(new, old)

    real(kind=REALNUMBER), dimension(0:,0:), intent(in) :: new, old
    real(kind=REALNUMBER) :: max_diff

    max_diff = maxval(abs( new(1:MP, 1:NP) - old(1:MP, 1:NP) ))
    call MPI_ALLREDUCE(max_diff, par_calc_max_diff, 1, MPI_REALNUMBER, MPI_MAX, cartcomm, ierr)

  end function par_calc_max_diff

  real(kind=REALNUMBER) function par_calc_ave(new, num_pixels)
    ! Calculate average the average pixel value by finding the local sum and ALLREDUCE
    real(kind=REALNUMBER), dimension(0:,0:), intent(in) :: new
    real(kind=REALNUMBER) :: localsum, totalsum
    integer, intent(in) :: num_pixels

    localsum = real(sum(real(new(1:MP,1:NP),kind=8)), kind=REALNUMBER)
    call MPI_ALLREDUCE(localsum, totalsum, 1, MPI_REALNUMBER, MPI_SUM, cartcomm, ierr)

    par_calc_ave = totalsum / num_pixels
  end function par_calc_ave


  ! -----------------------------------------------------!
  ! Set of helper routines which are not related to the  !
  ! message passing model but they depent on the num of  !
  ! processors and/or the MPI library.                   !
  ! -----------------------------------------------------!

  type(timetype) function get_time()
    get_time%value = MPI_WTIME()
    return
  end function get_time

  real function time_diff(time_start,time_finish)
    type(timetype), intent(in) :: time_start, time_finish
    time_diff = real(time_finish%value - time_start%value)
    return
  end function time_diff

  subroutine exit_all(message)
    character(*), intent(in) :: message
    write(*,*) "Error in process", rank, ":", message
    call MPI_ABORT(comm, 2, ierr)
  end subroutine exit_all
  
  subroutine print_once(message)
    character(*), intent(in) :: message
    if (par_isroot()) then
      write(*,*) message
    end if
  end subroutine print_once

  subroutine print_all(message)
    character(*), intent(in) :: message
    character(len=12) :: pn
    write(pn,'(A7,I4)') "Process ",rank
    write(*,*) pn,": ", message
  end subroutine print_all

END MODULE parallel


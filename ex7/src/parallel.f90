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
  integer :: FULL_WINDOW, INNER_WINDOW, HALO_V, HALO_H
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
    call MPI_TYPE_FREE(FULL_WINDOW,ierr)
    call MPI_TYPE_FREE(INNER_WINDOW,ierr)
    call MPI_TYPE_FREE(HALO_H,ierr)
    call MPI_TYPE_FREE(HALO_V,ierr)
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
    
    write *, "Grid size :", dims(1), " by ", dims(2)

    ! Create the derived datatypes
    call create_types()
  end subroutine

  subroutine create_types()
    ! Create all the derived datatypes used in the program, they are:
    ! Block type, master block type, vertical halo and horitzontal halo
    integer, dimension(N_DIMS) :: sizes, subsizes, starts
    integer(kind=mpi_address_kind) :: start, extent, lb, realextent
    integer :: AllocateStatus, i, base, LONG_INNER_WINDOW
    
    ! Block type: space in local arrays, inside halos
    sizes    = (/ MP+2, NP+2 /)
    subsizes = (/ MP , NP /)
    starts(:)   = 1
    call MPI_TYPE_CREATE_SUBARRAY(N_DIMS, sizes, subsizes, starts, &
             MPI_ORDER_FORTRAN, MPI_REALNUMBER, INNER_WINDOW, ierr)
    
    ! Master block type: distribution unit from master to working
    ! processes, needs a extent resize and the proper send_counts and
    ! displacements in order to be accessed iteratively.
    sizes    = (/ MP*dims(1), NP*dims(2) /)
    subsizes = (/ MP, NP /)
    starts(:)   = 0
    call MPI_TYPE_CREATE_SUBARRAY(N_DIMS, sizes, subsizes, starts, &
             MPI_ORDER_FORTRAN, MPI_REALNUMBER, LONG_INNER_WINDOW, ierr)
    call MPI_TYPE_GET_EXTENT(MPI_REALNUMBER, lb, realextent, ierr)

    start = 0
    extent = MP*realextent
    call MPI_TYPE_CREATE_RESIZED(LONG_INNER_WINDOW, start, extent, &
             FULL_WINDOW,ierr)
  
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
    call MPI_TYPE_VECTOR(NP, 1 , MP+2, MPI_REALNUMBER, HALO_V, ierr)
    call MPI_TYPE_VECTOR(1 , MP, MP  , MPI_REALNUMBER, HALO_H, ierr)
    
    ! COMMIT NEW MPI DATATYPES
    call MPI_TYPE_COMMIT(FULL_WINDOW, ierr)
    call MPI_TYPE_COMMIT(INNER_WINDOW,  ierr)
    call MPI_TYPE_COMMIT(HALO_V, ierr)
    call MPI_TYPE_COMMIT(HALO_H, ierr)
  end subroutine

  subroutine par_scatter(source, dest)
    real(kind=REALNUMBER), dimension(:,:), intent(in) :: source
    real(kind=REALNUMBER), dimension(:,:), intent(out) :: dest
    call MPI_Scatterv(source, send_counts, displacements, FULL_WINDOW, &
                      dest, MP*NP, MPI_REALNUMBER, 0, cartcomm,ierr)
  end subroutine par_scatter

  subroutine par_Gather(source, dest)
    real(kind=REALNUMBER), dimension(0:,0:), intent(in) :: source
    real(kind=REALNUMBER), dimension(:,:), intent(out) :: dest
    call MPI_GATHERV(source, 1, INNER_WINDOW, dest, send_counts, displacements, &
                     FULL_WINDOW, 0, cartcomm, ierr)
  end subroutine


  subroutine par_HalosSwap(old)
    ! Non-blocking send and reveive of all the halos, afterwards the
    ! par_WaitHalos() routine should be called to ensure these communications
    ! are completed.
    real(kind=REALNUMBER), dimension(0:,0:), intent(in) :: old
   
    call MPI_Issend(old(MP,1),1, HALO_V, n_right ,0,cartcomm,request(1),ierr)
    call MPI_Issend(old(1,1) ,1, HALO_V, n_left   ,0,cartcomm,request(3),ierr)
    call MPI_Issend(old(1,1),1, HALO_H, n_down,0,cartcomm,request(7),ierr)
    call MPI_Issend(old(1,NP) ,1, HALO_H, n_up ,0,cartcomm,request(5),ierr)
    
    call MPI_Irecv(old(MP+1,1),1, HALO_V, n_right ,0,cartcomm,request(4),ierr)
    call MPI_Irecv(old(0,1)   ,1, HALO_V, n_left ,0,cartcomm,request(2),ierr)
    call MPI_Irecv(old(1,0),1, HALO_H, n_down, 0,cartcomm,request(6),ierr)
    call MPI_Irecv(old(1,NP+1)   ,1, HALO_H, n_up ,0,cartcomm,request(8),ierr)
    
  end subroutine par_HalosSwap

  subroutine par_WaitHalos()
    integer, dimension(MPI_STATUS_SIZE,8) :: status
    call MPI_Waitall(8,request,status,ierr)
  end subroutine par_WaitHalos

  real(kind=REALNUMBER) function par_calc_max_diff(new, old)

    ! Calculate maximum change across the image
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

  function par_print(message)
    character(*), intent(in) :: message
    write(*,*) message
  end function par_print

END MODULE parallel


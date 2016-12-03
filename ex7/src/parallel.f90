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
  
  ! Problem parameters
  integer, dimension(N_DIMS) :: dims
  integer :: MP, NP ! Local array sizes

  ! MPI NEW DATATYPES
  integer :: FULL_WINDOW, INNER_WINDOW, HALO_V, HALO_H
  integer, dimension(:), allocatable :: send_counts, displacements

contains

  !--------------- GLOBAL MPI -------------------------!

  subroutine par_init()
    ! Initialise MPI
    call MPI_INIT(ierr)
    comm = MPI_COMM_WORLD
    call MPI_COMM_SIZE(comm, size, ierr)
  end subroutine par_init

  subroutine par_finalize()
    ! Free the used resources, finalize MPI
    call MPI_TYPE_FREE(FULL_WINDOW, ierr)
    call MPI_TYPE_FREE(INNER_WINDOW, ierr)
    call MPI_TYPE_FREE(HALO_H, ierr)
    call MPI_TYPE_FREE(HALO_V, ierr)
    call MPI_Finalize(ierr)
  end subroutine par_finalize

  subroutine par_abort(message)
    character(*), intent(in) :: message
    write(*,*) "Error in process", rank, ":", message
    call MPI_ABORT(comm, 2, ierr)
  end subroutine par_abort

  logical function par_isroot()
    ! Return true when executen in the ROOT process
    par_isroot = rank == ROOT
  end function par_isroot

  !--------------- MPI V.TOPOLOGY METHODS -------------------------!

  subroutine par_decompose(nx,ny,npx,npy)

    integer, intent(in) :: nx, ny
    integer, intent(out) :: npx, npy
    character(len=100) :: message
    logical, dimension(N_DIMS) :: periods
    logical :: reorder
    
    ! Create virtual topology
    dims(:) = 0
    periods(:) = .false.
    reorder = .true.

    ! Find a well balanced procesor distribution, create the topology and get the rank in cartcomm
    call MPI_DIMS_CREATE(size, N_DIMS, dims, ierr)
    call MPI_CART_CREATE(comm, N_DIMS, (/dims(2), dims(1) /), periods, reorder, cartcomm, ierr)
    call MPI_COMM_RANK(cartcomm, rank, ierr)

    ! Check if decomposition is valid, and calculate the pixel distribution
    if ((mod(nx, dims(1)) .ne. 0) .or. (mod(ny, dims(2)) .ne. 0)) then
        write(message,*) nx," is not divisible in ",dims(1)," parts"
        call par_abort("Could not decompose domain!"//message)
    end if

    npx = nx/dims(1)
    npy = ny/dims(2)
    MP = npx
    NP = npy

    if (par_isroot()) then
      ! Print out the information about the domain decomposition
      write(message,'(I2,A3,I2,A19,I4,A3,I4)') &
        dims(1), " x ", dims(2), "with each block of", MP, " x ", NP
      call par_print("Domain decomposed in:")
      call par_print(" ->  Grid of size: "//message)

    end if

    ! Create the derived datatypes
    call par_get_neighbours()
    call par_create_derived()
  end subroutine par_decompose

  subroutine par_get_neighbours()
    ! Get neighbours using cart_shift in x and y
    integer :: y_dir, x_dir, disp
    y_dir = 0
    x_dir = 1
    disp = 1

    call MPI_CART_SHIFT(cartcomm, y_dir, disp, n_down, n_up, ierr)
    call MPI_CART_SHIFT(cartcomm, x_dir, disp, n_left, n_right, ierr)
  end subroutine par_get_neighbours

  subroutine par_create_derived()

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
    ! processes, needs an extent resize and the proper send_counts and
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
  
    allocate(send_counts(size), displacements(size), STAT=allocateStatus)
    if(allocateStatus .ne. 0) call par_abort("*** NOT enough memory ***")
    
    call MPI_TYPE_COMMIT(FULL_WINDOW, ierr)
    call MPI_TYPE_COMMIT(INNER_WINDOW,  ierr)
    
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
    
    call MPI_TYPE_COMMIT(HALO_V, ierr)
    call MPI_TYPE_COMMIT(HALO_H, ierr)

  end subroutine par_create_derived

  !--------------- MPI COMM METHODS -------------------------!

  subroutine par_scatter(source, dest)
    real(kind=REALNUMBER), dimension(:,:), intent(in) :: source
    real(kind=REALNUMBER), dimension(:,:), intent(out) :: dest
    call MPI_Scatterv(source, send_counts, displacements, FULL_WINDOW, dest, MP*NP, MPI_REALNUMBER, 0, cartcomm, ierr)
  end subroutine par_scatter

  subroutine par_gather(source, dest)
    real(kind=REALNUMBER), dimension(0:,0:), intent(in) :: source
    real(kind=REALNUMBER), dimension(:,:), intent(out) :: dest
    call MPI_GATHERV(source, 1, INNER_WINDOW, dest, send_counts, displacements, FULL_WINDOW, 0, cartcomm, ierr)
  end subroutine par_gather


  subroutine par_swap_halos(old)
    ! Use non-blocking point-to-point communication 
    real(kind=REALNUMBER), dimension(0:,0:), intent(in) :: old
    integer, dimension(MPI_STATUS_SIZE) :: recv_status, send_status
    integer :: request
   
    call MPI_ISSEND(old(MP,1), 1, HALO_V, n_right, 0, cartcomm, request, ierr)
    call MPI_RECV(old(0,1), 1, HALO_V, n_left, 0, cartcomm, recv_status, ierr)
    call MPI_WAIT(request, send_status, ierr)

    call MPI_ISSEND(old(1,1), 1, HALO_V, n_left, 0, cartcomm, request, ierr)
    call MPI_RECV(old(MP+1,1), 1, HALO_V, n_right, 0, cartcomm, recv_status, ierr)
    call MPI_WAIT(request, send_status, ierr)

    call MPI_ISSEND(old(1,1), 1, HALO_H, n_down, 0, cartcomm, request, ierr)
    call MPI_RECV(old(1,NP+1), 1, HALO_H, n_up , 0, cartcomm, recv_status, ierr)  
    call MPI_WAIT(request, send_status, ierr)

    call MPI_ISSEND(old(1,NP), 1, HALO_H, n_up, 0, cartcomm, request, ierr)
    call MPI_RECV(old(1,0),1, HALO_H, n_down, 0, cartcomm, recv_status, ierr)
    call MPI_WAIT(request, send_status, ierr)
    
  end subroutine par_swap_halos

  !--------------- PROGRESS METHODS -------------------------!

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

  !--------------- ADDITIONAL METHODS -------------------------!

  type(timetype) function get_time()
    get_time%value = MPI_WTIME()
    return
  end function get_time

  real function time_diff(time_start,time_finish)
    type(timetype), intent(in) :: time_start, time_finish
    time_diff = real(time_finish%value - time_start%value)
    return
  end function time_diff

  subroutine par_print(message)
    character(*), intent(in) :: message
    write(*,*) message
  end subroutine par_print

END MODULE parallel


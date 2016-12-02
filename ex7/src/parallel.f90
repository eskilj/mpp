MODULE parallel

  use mpi
  use precisiondef

  implicit none

  ! Define time datatye
  type timetype
      real(kind=8) :: value
  end type timetype
  
  ! Constants
  ! IT DOESN'T WORK FOR DIFERENT CONSTANTS, defined to improve readability
  integer, parameter :: root = 0
  integer, parameter :: D = 2 !Num of dimentions

  ! MPI VARIABLES
  integer :: comm, cartcomm, myrank, commsize, ierr, errorcode
  integer :: n_left, n_right, n_up, n_down !process neighbours
  integer, dimension(8) :: request
  
  ! Problem parameters
  integer, dimension(D) :: dims
  integer :: Mp, Np, GM, GN ! Local and global array sizes

  ! MPI NEW DATATYPES
  integer :: MASTER_BLOCK_T, BLOCK_T, V_HALO_T, H_HALO_T
  integer, dimension(:), allocatable :: counts, displs

contains

  logical function MP_ISROOT()
    ! Return true when executen in the root process
    MP_ISROOT = myrank == root
    return
  end function MP_ISROOT

  subroutine MP_init()
    ! Initialize the MPI structures
    call MPI_INIT(ierr)
    comm = MPI_COMM_WORLD
    call MPI_COMM_SIZE(comm, commsize, ierr)
  end subroutine MP_init

  subroutine MP_domain_decomposition_2D(nx,ny,npx,npy)
    ! Creates a 2D cartesian communicator with the appropriate dimensions,
    ! compute the neighbours for each process and computes the MP and NP
    ! At the end calls the routine to create the derived datatypes
    integer, intent(in) :: nx, ny
    integer, intent(out) :: npx, npy
    character(len=100) :: message
    logical, dimension(D) :: periods
    
    ! Create cartesian topology 2D
    dims = (/ 0,0 /)
    periods = (/ .False., .False. /)
    call MPI_DIMS_CREATE(commsize, D, dims, ierr)
    ! dimensions shifted to be consisten with fortran array order
    call MPI_CART_CREATE(comm, D, (/dims(2), dims(1) /), periods, .True., cartcomm, ierr)
    call MPI_COMM_RANK(cartcomm, myrank, ierr)
 
    
    ! Get neighbours
    call MPI_CART_SHIFT(cartcomm,0,1,n_down,n_up,ierr)
    call MPI_CART_SHIFT(cartcomm,1,1,n_left,n_right,ierr)

    ! Compute the new array dimensions
    if(mod(nx,dims(1))/=0) then
        write(message,*) nx," is not divisible in ",dims(1)," parts"
        call exit_all("Could not decompose domain!"//message)
    end if
    if(mod(ny,dims(2))/=0) then
        write(message,*) ny," is not divisible in ",dims(2)," parts"
        call exit_all("Could not decompose domain!"//message)
    end if
    npx = nx/dims(1)
    npy = ny/dims(2)
    Mp = npx
    Np = npy
    GM = nx
    GN = ny

    ! Print out the information about the domain decomposition
    write(message,'(I2,A3,I2,A19,I4,A3,I4)') &
        dims(1), " x ", dims(2), "with each block of", Mp, " x ", Np
    call print_once("Domain decomposed in:")
    call print_once(" ->  Grid of size: "//message)

    ! Create the derived datatypes
    call create_types()
  end subroutine

  subroutine create_types()
    ! Create all the derived datatypes used in the program, they are:
    ! Block type, master block type, vertical halo and horitzontal halo
    integer, dimension(D) :: sizes, subsizes, starts
    integer(kind=mpi_address_kind) :: start, extent, lb, realextent
    integer :: AllocateStatus, i, base, LONG_BLOCK_T
    
    ! Block type: space in local arrays, inside halos
    sizes    = (/ MP+2, NP+2 /)
    subsizes = (/ MP , NP /)
    starts   = (/  1 ,  1 /)
    call MPI_TYPE_CREATE_SUBARRAY(D, sizes, subsizes, starts, &
             MPI_ORDER_FORTRAN, MPI_REALNUMBER, BLOCK_T, ierr)
    
    ! Master block type: distribution unit from master to working
    ! processes, needs a extent resize and the proper counts and
    ! displacements in order to be accessed iteratively.
    sizes    = (/ MP*dims(1), NP*dims(2) /)
    subsizes = (/ MP, NP /)
    starts   = (/  0,  0 /)
    call MPI_TYPE_CREATE_SUBARRAY(D, sizes, subsizes, starts, &
             MPI_ORDER_FORTRAN, MPI_REALNUMBER, LONG_BLOCK_T, ierr)
    call MPI_TYPE_GET_EXTENT(MPI_REALNUMBER, lb, realextent, ierr)
    start = 0
    extent = Mp*realextent
    call MPI_TYPE_CREATE_RESIZED(LONG_BLOCK_T, start, extent, &
             MASTER_BLOCK_T,ierr)
  
    allocate(counts(commsize), STAT=AllocateStatus)
    if(allocateStatus /= 0) call exit_all("*** NOT enough memory ***")
    allocate(displs(commsize), STAT=AllocateStatus)
    if(allocateStatus /= 0) call exit_all("*** NOT enough memory ***")
    
    base = 1
    do i= 1, commsize
        counts(i) = 1
        displs(i) = (base-1) + mod(i-1,dims(1))
        if(mod(i,dims(1))==0) base = base + Np * dims(1)
    end do

    ! HALO VECTORS (HALO-HALO INTERSECTIONS ARE NOT NEDDED)
    ! Horitzontal halo data is contiguous, defined to maintain
    ! code cohesion in all halo swaps.
    call MPI_TYPE_VECTOR(Np, 1 , Mp+2, MPI_REALNUMBER, V_HALO_T, ierr)
    call MPI_TYPE_VECTOR(1 , Mp, Mp  , MPI_REALNUMBER, H_HALO_T, ierr)
    
    ! COMMIT NEW MPI DATATYPES
    call MPI_TYPE_COMMIT(MASTER_BLOCK_T, ierr)
    call MPI_TYPE_COMMIT(BLOCK_T,  ierr)
    call MPI_TYPE_COMMIT(V_HALO_T, ierr)
    call MPI_TYPE_COMMIT(H_HALO_T, ierr)
  end subroutine

  subroutine MP_Scatter(source, dest)
    real(kind=REALNUMBER), dimension(:,:), intent(in) :: source
    real(kind=REALNUMBER), dimension(:,:), intent(out) :: dest
    call MPI_Scatterv(source, counts, displs, MASTER_BLOCK_T, &
                      dest, Mp*Np, MPI_REALNUMBER, 0, cartcomm,ierr)
  end subroutine

  subroutine MP_Gather(source, dest)
    real(kind=REALNUMBER), dimension(0:,0:), intent(in) :: source
    real(kind=REALNUMBER), dimension(:,:), intent(out) :: dest
    call MPI_GATHERV(source, 1, BLOCK_T, dest, counts, displs, &
                     MASTER_BLOCK_T, 0, cartcomm, ierr)
  end subroutine


  subroutine MP_HalosSwap(old)
    ! Non-blocking send and reveive of all the halos, afterwards the
    ! MP_WaitHalos() routine should be called to ensure these communications
    ! are completed.
    real(kind=REALNUMBER), dimension(0:,0:), intent(in) :: old
   
    call MPI_Issend(old(Mp,1),1, V_HALO_T, n_right ,0,cartcomm,request(1),ierr)
    call MPI_Issend(old(1,1) ,1, V_HALO_T, n_left   ,0,cartcomm,request(3),ierr)
    call MPI_Issend(old(1,1),1, H_HALO_T, n_down,0,cartcomm,request(7),ierr)
    call MPI_Issend(old(1,Np) ,1, H_HALO_T, n_up ,0,cartcomm,request(5),ierr)
    
    call MPI_Irecv(old(Mp+1,1),1, V_HALO_T, n_right ,0,cartcomm,request(4),ierr)
    call MPI_Irecv(old(0,1)   ,1, V_HALO_T, n_left ,0,cartcomm,request(2),ierr)
    call MPI_Irecv(old(1,0),1, H_HALO_T, n_down, 0,cartcomm,request(6),ierr)
    call MPI_Irecv(old(1,Np+1)   ,1, H_HALO_T, n_up ,0,cartcomm,request(8),ierr)
    
  end subroutine MP_HalosSwap

  subroutine MP_WaitHalos()
    integer, dimension(MPI_STATUS_SIZE,8) :: stats
    call MPI_Waitall(8,request,stats,ierr)
  end subroutine MP_WaitHalos

  subroutine MP_GetMaxChange(new,old,maxchange)
    ! Compute the local max change and execute a reduce operation to get the
    ! global one.
    real(kind=REALNUMBER), dimension(0:,0:), intent(in) :: new
    real(kind=REALNUMBER), dimension(0:,0:), intent(in) :: old
    real(kind=REALNUMBER), intent(out) :: maxchange
    real(kind=REALNUMBER) :: localmaxchange
    localmaxchange = maxval(abs(new(1:MP,1:NP)-old(1:MP,1:NP)))
    call MPI_ALLREDUCE(localmaxchange,maxchange,1,MPI_REALNUMBER, &
                       MPI_MAX, cartcomm, ierr)
  end subroutine MP_GetMaxChange

  subroutine MP_GetAverage(new, average)
    ! Compute the local sumation of the pixels,
    ! reduce the summation of all processes and
    ! divide by the total number of pixels to get the average
    real(kind=REALNUMBER), dimension(0:,0:), intent(in) :: new
    real(kind=REALNUMBER), intent(out) :: average
    real(kind=REALNUMBER) :: localsum, totalsum
    
    ! The the sum is not make in 8 bytes reals, it differs from serial version
    localsum = real(sum(real(new(1:MP,1:NP),kind=8)),kind=REALNUMBER)
    call MPI_ALLREDUCE(localsum,totalsum,1,MPI_REALNUMBER, &
                       MPI_SUM, cartcomm, ierr)
    
    average = totalsum / (GM*GN)
  end subroutine MP_GetAverage

  subroutine MP_Finalize()
    ! Free the used resources
    call MPI_TYPE_FREE(MASTER_BLOCK_T,ierr)
    call MPI_TYPE_FREE(BLOCK_T,ierr)
    call MPI_TYPE_FREE(H_HALO_T,ierr)
    call MPI_TYPE_FREE(V_HALO_T,ierr)
    call MPI_Finalize(ierr)
  end subroutine MP_Finalize


  ! -----------------------------------------------------!
  ! Set of helper routines which are not related to the  !
  ! message passing model but they depent on the num of  !
  ! processors and/or the MPI library.                   !
  ! -----------------------------------------------------!

  type(timetype) function get_time()
    get_time%value = MPI_WTIME()
    return
  end function get_time

  real function time_diff(tstart,tend)
    type(timetype), intent(in) :: tstart, tend
    time_diff = real(tend%value - tstart%value)
    return
  end function time_diff

  subroutine exit_all(message)
    character(*), intent(in) :: message
    write(*,*) "Error in process", myrank, ":", message
    call MPI_ABORT(comm, 2, ierr)
  end subroutine exit_all
  
  subroutine print_once(message)
    character(*), intent(in) :: message
    if (MP_ISROOT()) then
      write(*,*) message
    end if
  end subroutine print_once

  subroutine print_all(message)
    character(*), intent(in) :: message
    character(len=12) :: pn
    write(pn,'(A7,I4)') "Process ",myrank
    write(*,*) pn,": ", message
  end subroutine print_all

END MODULE parallel


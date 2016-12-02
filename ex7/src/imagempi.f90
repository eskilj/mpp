program imagempi
  use precision
  use pgmio

#ifdef SERIAL
    use serial
#else
    use parallel
#endif

  implicit none

  integer :: AllocateStatus
  integer :: i, j, nx, ny, npx, npy, it = 0
  integer, parameter :: maxlen = 100, MAX_ITER = 2500, PROGRESS_INTERVAL = 100
  character(maxlen) :: filename, outfile,  message
  real(kind=REALNUMBER), dimension(:,:), allocatable :: master, edge, old, new
  real(kind=REALNUMBER), parameter :: DIFF_THRESHOLD = 0.1
  real(kind=REALNUMBER) :: average, max_diff = 1
  type(timetype) :: time_start, time_finish
  
  !  --------------- INITIALIZATION  -------------------------! 
  ! Read parameters
  call getParameters(filename)
  call pgmsize(filename, nx, ny)

  outfile = 'out.pgm'

  ! Initialize MessagePassing Library
  call par_init()
  call par_domain_decomposition_2D(nx, ny, npx, npy)

  allocate(master(nx, ny), edge(npx, npy), old(0:npx+1, 0:npy+1), new(0:npx+1, 0:npy+1))

  if (par_isroot()) call pgmread(filename, master)
  call par_scatter(master, edge)
  old(:,:) = 255
 
  ! EXECUTE INVERT EDGES ALGORITHM
  time_start = get_time()

  do while ((it .lt. MAX_ITER) .and. (max_diff .gt. DIFF_THRESHOLD))
    ! Cange the halos between surrounding processors if such exist
    call par_HalosSwap(old)
    call par_WaitHalos()
    
    ! Compute the local new values not dependent to halos
    do j = 1,npy
      do i = 1,npx
        new(i,j) = 0.25 * (old(i-1,j) + old(i+1,j) + old(i,j-1) + old(i,j+1) - edge(i,j))
      end do
    end do

    ! Perform the necessary reductions every PROGRESS_INTERVAL
    it = it + 1
    if (mod(it,PROGRESS_INTERVAL) == 0) then

      max_diff = par_calc_max_diff(new, old)
      average = par_calc_ave(new, nx*ny)
      
      write(message,'(A10,I5,A17,F6.1)') "Iteration ",it,", pixel average: ", average
      call print_once(message)
    end if

    !Re-assign new to old
    old(1:npx,1:npy) = new(1:npx,1:npy)

  end do
  
  time_finish = get_time()
  write(message,'(A9,I5,A15,F8.3,A8)')"Executed ", it," iterations in ", &
                                       time_diff(time_start,time_finish)," seconds."
  call print_once(message)

  call par_Gather(old, master)
  if (par_isroot()) call pgmwrite(outfile,master)

  !FINALIZE MessagePassing
  call par_finalize()
  
  deallocate(master, edge, new, old)

contains

  subroutine getParameters(filename)
    character(len=100), intent(out) :: filename
    integer :: num_args

    num_args = command_argument_count()
    call get_command_argument(1, filename)
  end subroutine getParameters

end program imagempi

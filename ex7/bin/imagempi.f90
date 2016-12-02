

program imagempi
  use precision
  use pgmio




    use parallel


  implicit none

  integer :: AllocateStatus
  integer :: i, j, nx, ny, npx, npy, it
  integer, parameter :: maxlen = 100, MAX_ITER = 2500, PROGRESS_INTERVAL = 100
  character(maxlen) :: filename, outfile,  message
  real(kind=REALNUMBER), dimension(:,:), allocatable :: image, edge, old, new
  real(kind=REALNUMBER), parameter :: MAX_CHANGE = 0.1
  real(kind=REALNUMBER) :: average, maxchange = 1
  type(timetype) :: time_start, time_finish
  
  !  --------------- INITIALIZATION  -------------------------! 
  ! Read parameters
  call getParameters(filename)
  call pgmsize(filename, nx, ny)

  outfile = 'out.pgm'

  ! Initialize MessagePassing Library
  call par_Init()
  call par_domain_decomposition_2D(nx,ny,npx,npy)

  ! Allocate arrays memory
  allocate(image(nx,ny) , STAT = AllocateStatus )
  if (allocateStatus /= 0) call exit_all(" *** NOT enough memory ***")
  allocate(edge(npx,npy) , STAT = AllocateStatus )
  if (allocateStatus /= 0) call exit_all(" *** NOT enough memory ***")
  allocate(old(0:npx+1,0:npy+1) , STAT = AllocateStatus )
  if (allocateStatus /= 0) call exit_all(" *** NOT enough memory ***")
  allocate(new(0:npx+1,0:npy+1) , STAT = AllocateStatus )
  if (allocateStatus /= 0) call exit_all(" *** NOT enough memory ***")
  call print_all("All arrays allocated successfully!")


  ! Read image and distribute between processes
  time_start = get_time()

  if (par_ISROOT()) call pgmread(filename,image)
  call par_Scatter(image,edge)
  old(:,:) = 255

  time_finish = get_time()
  write(message,*) "Data readed and distributed in ",  &
                   time_diff(time_start,time_finish), " sec."
  call print_once(message)
 
  ! EXECUTE INVERT EDGES ALGORITHM
  time_start = get_time()
  it = 0
  do while (.not. condition(it, maxchange, MAX_ITER))
    ! Cange the halos between surrounding processors if such exist
    call par_HalosSwap(old)
    call par_WaitHalos()
    
    ! Compute the local new values not dependent to halos
    do j = 1,npy
      do i = 1,npx
        new(i,j) = ( old(i-1,j) + old(i+1,j) + old(i,j-1) &
                    + old(i,j+1) - edge(i,j) ) * 0.25
      end do
    enddo

    ! Perform the necessary reductions every PROGRESS_INTERVAL
    it = it + 1
    if (mod(it,PROGRESS_INTERVAL) == 0) then
      if (MAX_ITER == 0 ) then ! Only when a fixed number of iterations is not set
        call par_GetMaxChange(new, old, maxchange)
      end if
      call par_GetAverage(new, average)
      
      write(message,'(A10,I5,A17,F6.1)') "Iteration ",it,", pixel average: ", average
      call print_once(message)
    end if

    ! Loop updates
    old(1:npx,1:npy) = new(1:npx,1:npy)
  end do
  
  time_finish = get_time()
  write(message,'(A9,I5,A15,F8.3,A8)')"Executed ", it," iterations in ", &
                                       time_diff(time_start,time_finish)," seconds."
  call print_once(message)
  
  ! Gather the data again to the root proces and write it to the output file
  time_start = get_time()

  call par_Gather(old,image)
  if (par_ISROOT()) call pgmwrite(outfile,image)
  
  time_finish = get_time()
  write(message,*) "Data gathered and writed in ", time_diff(time_start,time_finish), " sec"
  call print_once(message)
  
  !FINALIZE MessagePassing
  call par_Finalize()
  
  deallocate(image)
  deallocate(edge)
  deallocate(new)
  deallocate(old)

contains

  subroutine getParameters(filename)
    character(len=100), intent(out) :: filename
    integer :: num_args

    num_args = command_argument_count()
    call get_command_argument(1, filename)
  end subroutine getParameters

  logical function condition(it,maxchange,MAX_ITER)
    integer, intent(in) :: it, MAX_ITER
    real(kind=REALNUMBER), intent(in) :: maxchange

    if (MAX_ITER == 0) then !if num iteration not fixed
      condition = maxchange < MAX_CHANGE !Stopping criterion
    else
      condition = MAX_ITER <= it
    end if
    return
  end function condition

end program imagempi
program image
  use precision
  use pgmio

#ifdef SERIAL
    use serial
#else
    use parallel
#endif

  implicit none
  
  real(kind=REALNUMBER), parameter :: STOP_CHANGE = 0.1
  integer :: IT_BTW_REDUCTIONS
  integer :: numiter, AllocateStatus
  integer :: i, j, nx, ny, npx, npy, it
  character(len=100) :: infile, outfile, message
  real(kind=REALNUMBER), dimension(:,:), allocatable :: image, edge, old, new
  real(kind=REALNUMBER) :: average, maxchange = 1
  type(timetype) :: tstart, tend
  
  !  --------------- INITIALIZATION  -------------------------! 
  ! Read parameters
  call getParameters(infile, outfile, numiter, IT_BTW_REDUCTIONS)
  call pgmsize(infile, nx, ny)

  ! Initialize MessagePassing Library
  call MP_Init()
  call MP_domain_decomposition_2D(nx,ny,npx,npy) 

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
  tstart = get_time()

  if (MP_ISROOT()) call pgmread(infile,image)
  call MP_Scatter(image,edge)
  old(:,:) = 255

  tend = get_time()
  write(message,*) "Data readed and distributed in ",  &
                   time_diff(tstart,tend), " sec."
  call print_once(message)
 
  ! EXECUTE INVERT EDGES ALGORITHM
  tstart = get_time()
  it = 0
  do while (.not. condition(it, maxchange, numiter))
    ! Cange the halos between surrounding processors if such exist
    call MP_HalosSwap(old)
    call MP_WaitHalos()
    
    ! Compute the local new values not dependent to halos
    do j = 1,npy
      do i = 1,npx
        new(i,j) = ( old(i-1,j) + old(i+1,j) + old(i,j-1) &
                    + old(i,j+1) - edge(i,j) ) * 0.25
      end do
    enddo

    ! Perform the necessary reductions every IT_BTW_REDUCTIONS
    it = it + 1
    if (mod(it,IT_BTW_REDUCTIONS) == 0) then
      if (numiter == 0 ) then ! Only when a fixed number of iterations is not set
        call MP_GetMaxChange(new, old, maxchange)
      end if
      call MP_GetAverage(new, average)
      
      write(message,'(A10,I5,A17,F6.1)') "Iteration ",it,", pixel average: ", average
      call print_once(message)
    end if

    ! Loop updates
    old(1:npx,1:npy) = new(1:npx,1:npy)
  end do
  
  tend = get_time()
  write(message,'(A9,I5,A15,F8.3,A8)')"Executed ", it," iterations in ", &
                                       time_diff(tstart,tend)," seconds."
  call print_once(message)
  
  ! Gather the data again to the root proces and write it to the output file
  tstart = get_time()

  call MP_Gather(old,image)
  if (MP_ISROOT()) call pgmwrite(outfile,image)
  
  tend = get_time()
  write(message,*) "Data gathered and writed in ", time_diff(tstart,tend), " sec"
  call print_once(message)
  
  !FINALIZE MessagePassing
  call MP_Finalize()
  
  deallocate(image)
  deallocate(edge)
  deallocate(new)
  deallocate(old)

contains

  subroutine getParameters(infile, outfile, numiter, it_btw_red)
    character(len=100), intent(out) :: infile, outfile
    integer, intent(out) :: numiter, it_btw_red
    integer :: num_args
    character(len=100) :: numit_s, it_btw_red_s

    num_args = command_argument_count()
    if (num_args == 3)  then
      outfile = 'output.pgm'
    elseif (num_args == 4) then
      call get_command_argument(4, outfile)
    else
      call exit_all("WRONG ARGUMENTS!")
    end if
    call get_command_argument(1, infile)
    call get_command_argument(2, numit_s)
    call get_command_argument(3, it_btw_red_s)
    read(numit_s,*) numiter
    read(it_btw_red_s,*) it_btw_red
  end subroutine getParameters

  logical function condition(it,maxchange,numiter)
    integer, intent(in) :: it, numiter
    real(kind=REALNUMBER), intent(in) :: maxchange

    if (numiter == 0) then !if num iteration not fixed
      condition = maxchange < STOP_CHANGE !Stopping criterion
    else
      condition = numiter <= it
    end if
    return
  end function condition

end program image

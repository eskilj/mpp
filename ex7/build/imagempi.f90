

program imagempi
  use precision
  use pgmio




    use parallel


  implicit none

  integer :: i, j, nx, ny, npx, npy, array_status, iter = 0
  integer, parameter :: MAXLEN = 100, MAX_ITER = 2500, PROGRESS_INTERVAL = 50
  character(MAXLEN) :: filename, outfile,  message
  real(kind=REALNUMBER), dimension(:,:), allocatable :: master, edge, old, new
  real(kind=REALNUMBER), parameter :: DIFF_THRESHOLD = 0.1
  real(kind=REALNUMBER) :: average, max_diff = 1
  type(timetype) :: time_start, time_finish
  
  !  --------------- INITIALIZATION  -------------------------! 
  
  ! Get program parameter and load image
  call get_params(filename, outfile)
  call pgmsize(filename, nx, ny)

  ! Initialize MPI and virtual topologies
  call par_init()
  call par_decompose(nx, ny, npx, npy)

  ! allocate arrays used for the image processing
  allocate(master(nx, ny), edge(npx, npy), old(0:npx+1, 0:npy+1), new(0:npx+1, 0:npy+1), STAT=array_status)
  call check_allocation(array_status)

  ! read image data into master and distribute
  if (par_isroot()) call pgmread(filename, master)
  call par_scatter(master, edge)
  old(:,:) = 255
 
  !  --------------- IMAGE PROCESSING  -------------------------!

  time_start = get_time()

  do while ((iter .lt. MAX_ITER) .and. (max_diff .gt. DIFF_THRESHOLD))
    
    ! Swap halos
    call par_swap_halos(old)
    
    ! Compute the local new values not dependent to halos
    do j = 1,npy
      do i = 1,npx
        new(i,j) = 0.25 * (old(i-1,j) + old(i+1,j) + old(i,j-1) + old(i,j+1) - edge(i,j))
      end do
    end do

    ! Calc average pixel values and max_diff
    iter = iter + 1
    if (mod(iter, PROGRESS_INTERVAL) == 0) then

      max_diff = par_calc_max_diff(new, old)
      average = par_calc_ave(new, nx*ny)
      
      if (par_isroot()) then 
        write(message,'(A10,I5,A17,F6.1)') "Iter: ", iter, " ==> Average px value: ", average
        call par_print(message)
      end if
    end if

    ! Set old = new
    old(1:npx,1:npy) = new(1:npx,1:npy)

  end do
  
  time_finish = get_time()

  if (par_isroot()) then 
    write(message,'(A9,I5,A15,F8.3,A8)') "Finished image processing in", iter, " iterations. Execution time: ", time_diff(time_start,time_finish)," seconds."
    call par_print(message)
  end if

  !  --------------- GATHER DATA and COMPLETE PROGRAM -------------------------!

  call par_gather(old, master)
  if (par_isroot()) call pgmwrite(outfile,master)

  call par_finalize()
  deallocate(master, edge, new, old)

contains

  subroutine check_allocation(allocation_status)
    integer, intent(in) :: allocation_status
    if (allocation_status .ne. 0) call par_abort("Memory Allocation unsuccessful.")
  end subroutine

  subroutine get_params(filename, outfile)
    character(MAXLEN), intent(out) :: filename, outfile
    call get_command_argument(1, filename)
    call get_command_argument(2, outfile)
  end subroutine get_params

end program imagempi

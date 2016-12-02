MODULE serial
  use precision

  implicit none

  type timetype
     integer :: value
  end type timetype
  integer :: M, N

  contains

  logical function par_isroot()
    par_isroot = .True.
    return
  end function par_isroot

  subroutine par_init()
  end subroutine par_init

  subroutine par_decompose(nx,ny,npx,npy)
    integer, intent(in) :: nx, ny
    integer, intent(out) :: npx, npy
    npx = nx
    npy = ny
    M = nx
    N = ny
  end subroutine

  subroutine par_scatter(source, destination)
    real(kind=REALNUMBER), dimension(:,:), intent(in) :: source
    real(kind=REALNUMBER), dimension(:,:), intent(out) :: destination
    destination = source
  end subroutine

  subroutine par_Gather(source, destination)
    real(kind=REALNUMBER), dimension(0:,0:), intent(in) :: source
    real(kind=REALNUMBER), dimension(:,:), intent(out) :: destination
    destination = source(1:M,1:N)
  end subroutine
 
  subroutine par_HalosSwap(old)
    real(kind=REALNUMBER), dimension(:,:), intent(inout) :: old
    old = old   !to avoid unused warnings
  end subroutine
  
  subroutine par_WaitHalos()
  end subroutine

  subroutine par_calc_max_diff(new, old, maxchange)
    real(kind=REALNUMBER), dimension(0:,0:), intent(in) :: new, old
    real(kind=REALNUMBER), intent(inout) :: maxchange
    maxchange = maxval(abs(new(1:M,1:N)-old(1:M,1:N)))
  end subroutine par_calc_max_diff
  
  subroutine par_calc_ave(new, average)
    real(kind=REALNUMBER), dimension(0:,0:), intent(in) :: new
    real(kind=REALNUMBER), intent(inout) :: average
    real(kind=8) :: accumulate
    
    accumulate = sum(real(new(1:M,1:N),kind=8))
    average = real(accumulate,kind=REALNUMBER) / (M*N)
  end subroutine par_calc_ave

  subroutine par_finalize()
  end subroutine


   ! -----------------------------------------------------!
   ! Set of helper routines which are not related to the  !
   ! message passing model but they depent on the num of  !
   ! processors and/or the existence of the MPI library.  !
   ! -----------------------------------------------------!
 
  type(timetype) function get_time()
    call system_clock(get_time%value)
    return
  end function get_time
  
  real function time_diff(time_start,time_finish)
    type(timetype), intent(in) :: time_start, time_finish
    integer :: clockrate
    call system_clock(count_rate=clockrate)
    time_diff = real(time_finish%value - time_start%value) / clockrate
    return
  end function time_diff

  subroutine par_abort(message)
    character(*), intent(in) :: message
    write(*,*) "Error in process 0 (serial code) :", message
    stop -1
  end subroutine par_abort

  subroutine print_once(message)
    character(*), intent(in) :: message
    write(*,*) message
  end subroutine print_once

  subroutine print_all(message)
    character(*), intent(in) :: message
    write(*,*) "Process 0 (serial code): ", message
  end subroutine print_all

END MODULE serial

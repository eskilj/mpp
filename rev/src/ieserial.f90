MODULE ieserial
  use precisiondef

  implicit none

  type timetype
     integer :: value
  end type timetype
  integer :: M, N

  contains

  logical function MP_ISROOT()
    MP_ISROOT = .True.
    return
  end function MP_ISROOT

  subroutine MP_init()
  end subroutine MP_init

  subroutine MP_domain_decomposition_2D(nx,ny,npx,npy)
    integer, intent(in) :: nx, ny
    integer, intent(out) :: npx, npy
    npx = nx
    npy = ny
    M = nx
    N = ny
  end subroutine

  subroutine MP_Scatter(source, destination)
    real(kind=REALNUMBER), dimension(:,:), intent(in) :: source
    real(kind=REALNUMBER), dimension(:,:), intent(out) :: destination
    destination = source
  end subroutine

  subroutine MP_Gather(source, destination)
    real(kind=REALNUMBER), dimension(0:,0:), intent(in) :: source
    real(kind=REALNUMBER), dimension(:,:), intent(out) :: destination
    destination = source(1:M,1:N)
  end subroutine
 
  subroutine MP_HalosSwap(old)
    real(kind=REALNUMBER), dimension(:,:), intent(inout) :: old
    old = old   !to avoid unused warnings
  end subroutine
  
  subroutine MP_WaitHalos()
  end subroutine

  subroutine MP_GetMaxChange(new, old, maxchange)
    real(kind=REALNUMBER), dimension(0:,0:), intent(in) :: new
    real(kind=REALNUMBER), dimension(0:,0:), intent(in) :: old
    real(kind=REALNUMBER), intent(inout) :: maxchange
    maxchange = maxval(abs(new(1:M,1:N)-old(1:M,1:N)))
  end subroutine MP_GetMaxChange
  
  subroutine MP_GetAverage(new, average)
    real(kind=REALNUMBER), dimension(0:,0:), intent(in) :: new
    real(kind=REALNUMBER), intent(inout) :: average
    real(kind=8) :: accumulate
    
    accumulate = sum(real(new(1:M,1:N),kind=8))
    average = real(accumulate,kind=REALNUMBER) / (M*N)
  end subroutine MP_GetAverage

  subroutine MP_Finalize()
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
  
  real function time_diff(tstart,tend)
    type(timetype), intent(in) :: tstart, tend
    integer :: clockrate
    call system_clock(count_rate=clockrate)
    time_diff = real(tend%value - tstart%value) / clockrate
    return
  end function time_diff

  subroutine exit_all(message)
    character(*), intent(in) :: message
    write(*,*) "Error in process 0 (serial code) :", message
    stop -1
  end subroutine exit_all

  subroutine print_once(message)
    character(*), intent(in) :: message
    write(*,*) message
  end subroutine print_once

  subroutine print_all(message)
    character(*), intent(in) :: message
    write(*,*) "Process 0 (serial code): ", message
  end subroutine print_all

END MODULE ieserial

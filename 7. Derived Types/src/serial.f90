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

  subroutine par_abort(message)
    character(*), intent(in) :: message
    write(*,*) "Error in process 0 (serial code) :", message
    stop -1
  end subroutine par_abort

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

  subroutine par_gather(source, destination)
    real(kind=REALNUMBER), dimension(0:,0:), intent(in) :: source
    real(kind=REALNUMBER), dimension(:,:), intent(out) :: destination
    destination = source(1:M,1:N)
  end subroutine
 
  subroutine par_swap_halos(old)
    real(kind=REALNUMBER), dimension(:,:), intent(inout) :: old
    old = old   !to avoid unused warnings
  end subroutine par_swap_halos
  
  real(kind=REALNUMBER) function par_calc_max_diff(new, old)
    real(kind=REALNUMBER), dimension(0:,0:), intent(in) :: new, old
    par_calc_max_diff = maxval(abs(new(1:M,1:N)-old(1:M,1:N)))
  end function par_calc_max_diff
  
  real(kind=REALNUMBER) function par_calc_ave(new, num_pixels)
    real(kind=REALNUMBER), dimension(0:,0:), intent(in) :: new
    integer, intent(in) :: num_pixels
    real(kind=8) :: pixel_sum
    
    pixel_sum = sum(real(new(1:M,1:N),kind=8))
    par_calc_ave = real(pixel_sum, kind=REALNUMBER) / (num_pixels)
  end function par_calc_ave

  subroutine par_finalize()
  end subroutine


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

  subroutine par_print(message)
    character(*), intent(in) :: message
    write(*,*) message
  end subroutine par_print

END MODULE serial

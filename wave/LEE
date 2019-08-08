
  module constants

    implicit none

    real(8) :: pi
    real(8) :: us

  end module constants
!***************************************************************************************************

  module subprograms

    use constants
    implicit none

  contains

    subroutine set_dt_uw(t,dt_uw)

      !\frac{\partial u_{w}}{\partial t}(t) = dt_uwを求める

      !u_w = a[-sin(t)+t*exp(-t)]
      !dt_uw = a[-cos(t)+(1-t)*exp(-t)]

      real(8), intent(in)   :: t
      real(8), intent(out)  :: dt_uw

      dt_uw = -cos(t) + (1.0d0-t)*exp(-t)

    end subroutine set_dt_uw
  !-------------------------------------------------------------------------------------------------

    subroutine set_b(N,p_pre,p,dt_uw,dx,delta,b)

      !p_pre, p, dt_uwのデータからbを求める

      integer, intent(in)  :: N
      real(8), intent(in)  :: p_pre(1:N-1), p(1:N-1)
      real(8), intent(in)  :: dt_uw
      real(8), intent(in)  :: dx, delta
      real(8), intent(out) :: b(1:N-1)

      !i = 1
      b(1) = 2.0d0*p(1) - p_pre(1) + 4.0d0/3.0d0*delta*dx*dt_uw
      !2 <= i <= N-1
      b(2:N-1) = 2.0d0*p(2:N-1) - p_pre(2:N-1)

    end subroutine set_b
  !-------------------------------------------------------------------------------------------------

    subroutine sol_x(N,A,b,x)

      !A*x = bの解xを求める

      integer, intent(in)   :: N
      real(8), intent(in)   :: A(1:3,1:N-1), b(1:N-1)
      real(8), intent(out)  :: x(1:N-1)

      real(8), allocatable :: tl_A(:), tl_b(:)
      integer :: i

      allocate( tl_A(1:N-1), tl_b(1:N-1) )

      !\tilde{a}_{ii} = tl_A(i)と\tilde{b}_{i} = tl_b(i)を求める
      !- i = 1
      tl_A(1) = A(2,1)
      tl_b(1) = b(1)
      !- i >= 2
      do i = 2, N-1
        tl_A(i) = A(2,i) - A(1,i)/tl_A(i-1)*A(3,i-1)
        tl_b(i) = b(i) - A(1,i)/tl_A(i-1)*tl_b(i-1)
      enddo

      !解xを求める
      x(N-1) = tl_b(N-1)/tl_A(N-1)
      do i = N-2, 1, -1
        x(i) = (tl_b(i)-A(3,i)*x(i+1))/tl_A(i)
      enddo

    end subroutine sol_x
  !-------------------------------------------------------------------------------------------------

    function int2char(x,fmt) result(res)

    !整数xを与えられた書式fmtの文字列に変換して返す

      integer, intent(in) :: x
      character(*), intent(in) :: fmt
      character*30 :: res

      write(res,trim(fmt)) x

    end function int2char
  !-------------------------------------------------------------------------------------------------

  end module subprograms
!***************************************************************************************************

  program main

    !18/04/21-18/04/22,18/05/19のノートを参照

    use constants
    use subprograms
    implicit none

    character*300 :: idname

    integer :: imax
    integer :: div_period
    real(8) :: hat_L, hat_L_div_wavelength
    integer :: num_period
    integer :: output_div_period

    integer :: nmax
    real(8) :: dt
    real(8) :: dx
    real(8) :: delta

    real(8), allocatable :: A(:,:), b(:)

    real(8), allocatable :: x(:)
    real(8), allocatable :: p(:), p_pre(:) !p_pre(:): p(t^{n})

    real(8) :: dt_uw

    real(8) :: t

    integer :: i, n

    pi = 4.0d0*datan(1.0d0)
    us = sqrt(5.0d0/6.0d0)

    !パラメータの入力
    read(*,*) 
    read(*,'(A)') idname
    read(*,*) 
    read(*,*) hat_L_div_wavelength
    read(*,*) 
    read(*,*) imax
    read(*,*) 
    read(*,*) div_period
    read(*,*) 
    read(*,*) num_period
    read(*,*) 
    read(*,*) output_div_period

    !組み立て単位の定義
    nmax = div_period*num_period
    dt = (2.0d0*pi)/dble(div_period)
    hat_L = hat_L_div_wavelength*(2.0d0*pi*us)
    dx = hat_L/dble(imax)
    delta = (5.0d0/6.0d0)*(dt/dx)**2

    write(*,'(A)') 'idname'
    write(*,'(A)') trim(idname)
    write(*,'(A)') 'hat_L_div_wavelength'
    write(*,'(f16.8)') hat_L_div_wavelength
    write(*,'(A)') 'nmax'
    write(*,'(i0)') nmax
    write(*,'(A)') 'dt'
    write(*,'(f16.8)') dt
    write(*,'(A)') 'dx'
    write(*,'(f16.8)') dx
    write(*,'(A)') 'delta'
    write(*,'(f16.8)') delta
    write(*,'(A)') ''

    open(unit=30,status='unknown',file=trim(idname)//'parameter.dat')
    write(30,'(A)') 'idname'
    write(30,'(A)') trim(idname)
    write(30,'(A)') 'hat_L_div_wavelength'
    write(30,'(f16.8)') hat_L_div_wavelength
    write(30,'(A)') 'imax'
    write(30,'(i0)') imax
    write(30,'(A)') 'div_period'
    write(30,'(i0)') div_period
    write(30,'(A)') 'num_period'
    write(30,'(i0)') num_period
    write(30,'(A)') 'output_div_period'
    write(30,'(i0)') output_div_period
    write(30,'(A)') ''
    write(30,'(A)') 'nmax'
    write(30,'(i0)') nmax
    write(30,'(A)') 'dt'
    write(30,'(f16.8)') dt
    write(30,'(A)') 'dx'
    write(30,'(f16.8)') dx
    write(30,'(A)') 'delta'
    write(30,'(f16.8)') delta
    write(30,'(A)') ''
    close(30)

    allocate( A(1:3,1:imax-1), b(1:imax-1) )
    allocate( x(0:imax) )
    allocate( p_pre(0:imax), p(0:imax) )

    !xの定義
    do i = 0, imax
      x(i) = dx*dble(i)
    enddo

    !Aの定義
    A = 0.0d0
    !-対角成分A(2,:)
    A(2,1) = 1.0d0 + 2.0d0/3.0d0*delta
    do i = 2, imax-2
      A(2,i) = 1.0d0 + 2.0d0*delta
    enddo
    A(2,imax-1) = 1.0d0 + 2.0d0/3.0d0*delta
    !-左側帯成分A(1,:)
    A(1,1) = 0.0d0
    do i = 2, imax-2
      A(1,i) = -delta
    enddo
    A(1,imax-1) = -2.0d0/3.0d0*delta
    !-右側帯成分A(3,:)
    A(3,1) = -2.0d0/3.0d0*delta
    do i = 2, imax-2
      A(3,i) = -delta
    enddo
    A(3,imax-1) = 0.0d0

    !初期条件( n = 0, 1での各データを定義できる )
    p_pre(0:imax) = 0.0d0
    p(0:imax) = 0.0d0

    open(unit=20,status='unknown',file=trim(idname)//'p.dat')
    write(20,'(A)') 'variables = "x", "p"'

    do n = 0, nmax-1
      t = dble(n+2)*dt
      !新時刻t^{n+2}でのdt_uwを求める
      call set_dt_uw(t,dt_uw)
      !方程式の非斉次項bを計算
      call set_b(imax,p_pre(1:imax-1),p(1:imax-1),dt_uw,dx,delta,b(1:imax-1))
      !p_preのデータをt^{n}のものからt^{n+1}のものへここで更新しておく
      p_pre(0:imax) = p(0:imax)
      !新時刻t^{n+2}でのp(1:imax-1)の値を求める
      call sol_x(imax,A,b(1:imax-1),p(1:imax-1))
      !新時刻t^{n+2}でのp(0), p(imax)の値を境界条件から求める
      p(0) = ( 4.0d0*p(1) - p(2) + 4.0d0*dx*dt_uw )/3.0d0
      p(imax) = ( 4.0d0*p(imax-1) - p(imax-2) )/3.0d0
      !データを出力
      if( mod(n+2,output_div_period) .eq. 0 ) then
        write(*,'(A,f0.8)') 'now t/period = ', t/(2.0d0*pi)
        write(20,'(A,f0.8,A)') 'zone t = "t/period = ', t/(2.0d0*pi), '"'
        do i = 0, imax
          write(20,'(2E16.8)') x(i), p(i)
        enddo
      endif
    enddo
    close(20)


  end program main
!***************************************************************************************************

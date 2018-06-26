      program main
      implicit none
      integer :: ix, j, ios, pgopen
      real(8), parameter :: hbar = 197.3d0
      real(8), parameter :: mass = 938.6d0 
      real(8), parameter :: a = 6.0d0
      real(8), parameter :: V0 = 35.0d0
      real(8), parameter :: xmin = -20.0d0, xmax = 20.0d0, dx = 0.1d0
      real(8), allocatable :: xa(:)
      real(8) :: k, kp, E, psimax
      complex(8), parameter :: i = (0.0d0,1.0d0)
      complex(8), allocatable, dimension(:) :: psi
      complex(8) :: T, R, C, D

      ix = nint((xmax - xmin)/dx)
      allocate(xa(ix+1),psi(ix+1))
      forall(j=1:ix) xa(j) = xmin + dble(j-1) * dx

      E = 32.0d0
      k = sqrt(2.0d0 * mass * E) / hbar
      kp = sqrt(2.0d0 * mass * (V0 - E)) / hbar
      write(6,*) kp*a

      T = 2.0d0 * i * k * kp * exp(-i*k*a)   &
          / ((k*k-kp*kp)*sinh(kp*a) + 2.0d0*i*k*kp*cosh(kp*a))
      R = (k*k + kp*kp) / (2.0d0*i*k*kp) * exp(i*k*a) &
          * T * sinh(kp*a)
      C = (kp - i * k) / (2.0d0*kp) * exp(i*k*a+kp*a) * T
      D = (kp + i * k) / (2.0d0*kp) * exp(i*k*a-kp*a) * T
      write(6,*) abs(T) ** 2 + abs(R) ** 2

      j = 1
      do
        if (xa(j) > - 0.5d0*a) exit
        psi(j) = T * exp(-i*k*(xa(j)-0.5d0*a))
        j = j + 1
      end do
      do
        if (xa(j) > 0.5d0*a) exit
        psi(j) = C * exp( kp*(xa(j)-0.5d0*a)) &
               + D * exp(-kp*(xa(j)-0.5d0*a))
        j = j + 1
      end do
      do while(j < ix+1)
        psi(j) = exp(-i*k*(xa(j)-0.5d0*a)) &
                + R * exp(i*k*(xa(j)-0.5d0*a))
        j = j + 1
      end do

      open(7,file='wf.dat')
      open(8,file='prob.dat')
      do j=1,ix
        write(7,*) xa(j), dble(psi(j))
        write(8,*) xa(j), abs(psi(j)) ** 2
      end do
      close(7)
      close(8)

      ios = pgopen('/xserv')
      call pgenv(real(xmin),real(xmax),-2.0,2.0,0,1)
      call pglab('x','psi','tunneling')
      call pgsci(2)
      call pgline(ix,real(xa),real(dble(psi)),0,1)
      call pgsci(4)
      call pgline(ix,real(xa),real(aimag(psi)),0,1)
      call pgend

      stop
      end program


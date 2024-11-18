program heat

   implicit none

   integer :: jmax, nmax, jmap, ajm, n, ns
   real :: alph, s, tmax, delx, delt, t
   real, dimension(101) :: x, tn, te
   real :: lamn
   real, parameter :: pi = 3.1415927
   integer :: j, aj, m
   real :: dte, an, dm, dn

   ! I/O variables
   character(len=100) :: data_dir = 'data/'
   logical :: dir_exists
   integer :: io_status
   character(len=100) :: error_msg

   ! Create data directory
   inquire(file=trim(data_dir), exist=dir_exists)
   if (.not. dir_exists) then
      call execute_command_line('mkdir -p '//trim(data_dir), exitstat=io_status)
      if (io_status /= 0) error stop 'Failed to create data directory'
   end if

   ! Open files with error checking
   call open_files()

   ! Initialize grid and parameters
   call initialize_parameters()

   ! Time integration
   time_loop: do while (t < tmax)
      ! Set boundary conditions
      tn(1) = 0.0
      tn(jmax) = 0.0

      if (mod(n, 100) == 1) then
         ! Output numerical solution
         call output_solution(tn, t, n, 'Numerical')

         ! Calculate and output analytical solution
         call calculate_analytical_solution()
         call output_solution(te, t, n, 'Analytical')

         write(20, *) '--'
         ns = ns + 1
      end if

      ! Update temperature (FTCS scheme)
      call update_temperature()

      ! Update time
      n = n + 1
      t = t + delt

      if (n >= nmax) exit time_loop
   end do time_loop

   ! Final output
   call output_solution(tn, t, n, 'Numerical')
   call calculate_analytical_solution()
   call output_solution(te, t, n, 'Analytical')

   write(20, '(A,I5)') 'Number of output time steps =', ns

   ! Close files
   close(10)
   close(20)
   close(30)

contains
   subroutine open_files()
      open(unit=10, file=trim(data_dir)//'diff2.in', status='old', iostat=io_status, iomsg=error_msg)
      if (io_status /= 0) error stop 'Failed to open input file: '//trim(error_msg)

      open(unit=20, file=trim(data_dir)//'diff2.out', iostat=io_status, iomsg=error_msg)
      if (io_status /= 0) error stop 'Failed to open output file: '//trim(error_msg)

      open(unit=30, file=trim(data_dir)//'diff2.da', form='unformatted', iostat=io_status, iomsg=error_msg)
      if (io_status /= 0) error stop 'Failed to open data file: '//trim(error_msg)

      read(10, *, iostat=io_status) jmax, nmax, alph, s, tmax
      if (io_status /= 0) error stop 'Failed to read input parameters'
   end subroutine open_files

   subroutine initialize_parameters()

      jmap = jmax - 1
      ajm = jmap
      delx = 1./real(jmap)
      delt = delx*delx*s/alph

      write(30) jmax
      write(20, '(A,I5,A,I5,A,F8.2)') ' jmax=', jmax, ' nmax=', nmax, ' tmax=', tmax
      write(20, '(A,F5.3,A,E10.3,A,E10.3,A,E10.3,//)') &
         ' s=', s, ' alpha=', alph, ' dt=', delt, ' dx=', delx
      write(20, '(A,F5.3,//)') ' FTCS (Explicit) Scheme, s=', s

      ! in initialize_parameters
      if (s > 0.5) then
         error stop 'Stability condition violated: s must be <= 0.5'
      end if

      ! Initialize grid and initial conditions
      do j=1,jmax
         aj = j - 1
         x(j) = delx*aj
      end do

      do j=1,jmax
         if (x(j) < 0.5) then
            tn(j)=2.*x(j)
         else
            tn(j)=2.-2.*x(j)
         endif
      end do

      write(30) x(1:jmax)

      n = 1
      ns = 1
      t = 0.0

   end subroutine initialize_parameters

   subroutine calculate_analytical_solution()

      integer, parameter :: M_MAX = 100

      do j=1, jmax
         dte = 0.

         do m=1, M_MAX
            dm = m
            dn = 2.*dm-1.
            lamn = dn*pi/1.
            an = 8.*(-1.)**((dn-1.)/2.) / (dn**2*pi**2)
            dte = an*sin(lamn*x(j))*exp(-alph*lamn*lamn*t) + dte
         end do

         te(j) = dte
      end do

   end subroutine calculate_analytical_solution

   subroutine update_temperature()

      do j=2,jmap
         tn(j) = (1.-2.*s)*tn(j) + s*(tn(j-1) + tn(j+1))
      end do

   end subroutine update_temperature

   subroutine output_solution(temp, time, step, solution_type)
      real, intent(in) :: temp(:), time
      integer, intent(in) :: step
      character(len=*), intent(in) :: solution_type

      write(30) step, time, temp(1:jmax)
      ! write(20, '(A,F5.0,A,11F6.2)') ' Time=', time, ' Temp=', (temp(j), j=1,jmax,10)
      write(20, '(A,A,A,F5.0,A,11F6.2)') ' ', solution_type, ' at Time=', &
         time, ' Temp=', (temp(j), j=1,jmax,10)

   end subroutine output_solution

end program heat



program heat

   implicit none

   integer :: jmax, nmax, jmap, ajm, n, ns, j, aj
   real :: alph, s, tmax, delx, delt, t
   real, dimension(101) :: x, tn
   real, parameter :: pi = 3.1415927

   character(len=100) :: data_dir = 'data/'
   logical :: dir_exists

   ! Input/Output variables
   integer :: io_status
   character(len=100) :: error_msg

   ! Create data directory if it doesn't exist
   inquire(file=trim(data_dir), exist=dir_exists)
   if (.not. dir_exists) then
      call execute_command_line('mkdir -p '//trim(data_dir), exitstat=io_status)
      if (io_status /= 0) error stop 'Failed to create data directory'
   end if


   ! Open files with error checking
   open(unit=10, file=trim(data_dir)//'diff.in', status='old', iostat=io_status, iomsg=error_msg)
   if (io_status /= 0) error stop 'Failed to open input file: '//trim(error_msg)

   open(unit=20, file=trim(data_dir)//'diff.out', iostat=io_status, iomsg=error_msg)
   if (io_status /= 0) error stop 'Failed to open output file: '//trim(error_msg)

   open(unit=30, file=trim(data_dir)//'diff.da', form='unformatted', iostat=io_status, iomsg=error_msg)
   if (io_status /= 0) error stop 'Failed to open'

   ! Read input parameters
   read(10, *, iostat=io_status) jmax, nmax, alph, s, tmax
   if (io_status /= 0) error stop 'Failed to read input parameters'


   jmap = jmax - 1
   ajm = jmap
   delx = 1./jmap
   delt = delx*delx*s/alph

   ! Write parameters to output file
   write(30) jmax
   write(20, '(A,I5,A,I5,A,F8.2)') ' jmax=', jmax, ' nmax=', nmax, ' tmax=', tmax
   write(20, '(A,F5.3,A,E10.3,A,E10.3,A,E10.3,//)') &
      ' s=', s, ' alpha=', alph, ' dt=', delt, ' dx=', delx
   write(20, '(A,F5.3,//)') ' FTCS (Explicit) Scheme, s=', s

   ! Initialize spatial grid and initial conditions
   do j=1, jmax
      aj = j - 1
      x(j) = delx*aj
   end do

   write(30) x(1:jmax)

   do j=1,jmax
      if (x(j) < 0.5) then
         tn(j)=2.*x(j)
      else
         tn(j)=2.-2.*x(j)
      endif
   end do

   ! Time integration
   n = 1
   ns = 1
   t = 0.

   time_loop: do while (t < tmax)
      ! Set boundary conditions
      tn(1) = 0.
      tn(jmax) = 0.


      if (mod(n, 100) == 1) then
         call output_results(tn, t, n, jmax, ns)
      end if

      ! Update temperature
      do j=2, jmap
         tn(j) = (1.0-2.0*s)*tn(j) + s*(tn(j-1) + tn(j+1))
      end do

      n = n + 1
      t = t + delt

      if (n >= nmax) exit time_loop

   end do time_loop

   ! Final output
   call output_results(tn, t, n, jmax, ns)
   write(20, '(A,I5)') 'Number of output time steps =', ns

   ! Close files
   close(10)
   close(20)
   close(30)

contains
   subroutine output_results(temp, time, step, grid_size, out_steps)
      real, intent(in) :: temp(:), time
      integer, intent(in) :: step, grid_size
      integer, intent(inout) :: out_steps

      write(30) step, time, temp(1:grid_size)
      write(20, '(A,F5.0,A,11F6.2)') ' Time=', time, ' Temp=', (temp(j), j=1,grid_size,10)
      out_steps = out_steps + 1

   end subroutine output_results

end program heat


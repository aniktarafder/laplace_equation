!*************************************************
! Laplace Serial Fortran Version
!
! Temperature is initially 0.0
! Boundaries are as follows:
!
!            T=linear     100
!    0 +-------------------+  100
!      |                   |
!      |                   |
!      |                   |
! T=0  |                   |  T=linear
!      |                   |
!      |                   |
!      |                   |
!      +-------------------+ 
!               T=0        0
!
!  John Urbanic, PSC 2014
!  MOdified by anik 2021
!*************************************************
program serial
      implicit none

      !Size of plate
      integer, parameter             :: nx = 128
      integer, parameter             :: ny = 128
      real,    parameter             :: max_temp_error = 0.01
      integer                        :: i, j, iteration, max_iterations
      double precision, parameter    :: dt=0.0001, last_time=10.0
      real                           :: start_time, stop_time
      real                           :: Lx, Ly, hxi, hyi, time
	  real                           :: temp_err

      double precision               :: temperature(nx+2,ny+2), temperature_last(nx+2,ny+2)

      Lx = 10.0
      Ly = 10.0
	  iteration = 1
      hxi= (nx+1)/Lx      
      hyi= (ny+1)/Ly
      time=0.0      
      temp_err=100.0

      print*, 'Temperature array', size(Temperature)
      print*, 'Maximum iterations [100-4000]?'
      read*,   max_iterations	  
	  
	  
	  
	  
	  
	  
      call cpu_time(start_time)      !Fortran timer

      call initialize(temperature_last, nx, ny)	  
      call initialize(temperature, nx, ny)	  
!      call print_field(temperature_last, iteration, nx+2, ny+2, hxi, hyi, time)
!	  stop
      !do until error is minimal or until maximum steps

!      do while ( temp_err > max_temp_error .and. time <= last_time)
      do while ( temp_err > max_temp_error .and. iteration <= max_iterations)

         do j=2, ny+1
            do i=2,nx+1
               temperature(i,j)=0.25*(temperature_last(i+1,j)+temperature_last(i-1,j)+ &
                                      temperature_last(i,j+1)+temperature_last(i,j-1) )
            enddo
         enddo

         temp_err=0.0

         !copy grid to old grid for next iteration and find max change
         do j=2,ny+1
            do i=2,nx+1
               temp_err = max( abs(temperature(i,j) - temperature_last(i,j)), temp_err )
               temperature_last(i,j) = temperature(i,j)
            enddo
         enddo
         time = dt*iteration
         !periodically print test values
         if( mod(iteration,100).eq.0 ) then
            call track_progress(temperature, iteration, nx, ny)
			call print_field(temperature, iteration, nx+2, ny+2, hxi, hyi, time)
         endif
         
!         time = dt*iteration
         iteration = iteration+1
      enddo

      call cpu_time(stop_time)

      print*, 'Max error at iteration ', iteration-1, ' was ',temp_err
      print*, 'Total time was ',stop_time-start_time, ' seconds.'

end program serial


! initialize plate and boundery conditions
! temp_last is used to start first iteration
subroutine initialize( temperature_last, nx, ny)
      implicit none

      integer            :: nx
      integer            :: ny
      integer            :: i,j

      double precision   :: temperature_last(nx+2, ny+2)

      temperature_last = 0.0

      !these boundary conditions never change throughout run

      !set left side to 0 and right to linear increase
      do j=1,ny+2
         temperature_last(1,j)    = 0.0
         temperature_last(nx+2,j) = (100.0/(ny+1)) * (j-1) ! Though at j=1, T=100 but at j=ny+2 T is not 0. But the Top BCs will take care of it as we will put zero at all top nodes.
      enddo

      !set bottom to 0 and top to linear increase
      do i=1,nx+2
         temperature_last(i,ny+2) = ((100.0)/(nx+1)) * (i-1)
         temperature_last(i,1)    = 0.0
      enddo

end subroutine initialize


!print diagonal in bottom corner where most action is
subroutine track_progress(temperature, iteration, nx, ny)
      implicit none

      integer             :: nx
      integer             :: ny
      integer             :: i,iteration

      double precision    :: temperature(nx+2, ny+2)

      print *, '---------- Iteration number: ', iteration, ' ---------------'
      do i=5,0,-1
         write (*,'("("i4,",",i4,"):",f6.2,"  ")',advance='no') &
                   nx+1-i,ny+1-i,temperature(nx+1-i,ny+1-i)
      enddo
      print *
end subroutine track_progress



subroutine print_field(temperature, iteration, nxp2, nyp2, hxi, hyi, time)
      implicit none

      integer             :: nxp2
      integer             :: nyp2
      integer             :: i,j,iteration
	  logical             :: fext
	  real                :: hxi, hyi, time

      double precision    :: temperature(nxp2, nyp2)

      INQUIRE(FILE='Temperature_field.dat', EXIST=fext)
      if(.not.fext) then

         open(30,file='Temperature_field.dat',form='formatted')

         write(30,*)'VARIABLES = "X","Y","Temperature"'      

         write(30,*)'ZONE ','T= "',time,'" I=',nxp2,' J=',nyp2,&
        ' F=BLOCK',', SOLUTIONTIME=',time
     

         write(30,'(8e16.8)')(((DBLE(i-1))/hxi,i=1,nxp2),j=1,nyp2)
         write(30,'(8e16.8)')(((DBLE(j-1))/hyi,i=1,nxp2),j=1,nyp2)   
         write(30,'(8e16.8)')temperature


       close(30)
      else
       open(30,file='Temperature_field.dat',form='formatted',status='old', &
               position='append')
      write(30,*) 'ZONE ','T="',time,'" I=',nxp2,' J=',nyp2, &
       ' F=BLOCK, ','VARSHARELIST=([1-2]=1)',', SOLUTIONTIME=',time
     
         write(30,'(8e16.8)')temperature 
	  
       close(30)
      endif    	  
	  
	  

end subroutine print_field

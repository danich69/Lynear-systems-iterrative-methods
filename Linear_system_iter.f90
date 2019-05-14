module Linear_system_iter
	implicit none
    
!    -----------------------------------------------------
!    |   Module contains 3 iterrative methods            |
!    |   for solving linear systems.                     |
!    |   All subroutines have the same input parameters: |
!    |   a -- system square matrix (real)                |
!    |   b -- column vector (real)                       |
!    |   x -- solution vector (real)                     |
!    |   Sys_size -- size of system (integer)            |
!	 -----------------------------------------------------

	real(8), parameter :: eps=1e-8                                      ! precision
	contains
	
	subroutine Jacobi(a,b,x_new)                                        ! Jacobi method (system matrix, vector, solution vector)
	
		integer :: i, j, Sys_size 		                                ! indexes and system size
		real(8), dimension(:,:) :: a							        ! system matrix
		real(8), dimension(:) :: b	                                    ! column vector
		real(8), dimension( size(b) , size(b) ) :: z, d, d_inv			! matrixes used for method
		real(8), dimension( size(b) ) ::  x_new, x_old, g		        ! vectors for iterations
		
        Sys_size = size(b) 
        
		do i=1, Sys_size                                                ! check if system matrix is diagonally mdomaint matrix
			if( SUM(abs(a(i,:))) > 2*abs(a(i,i))) then
				write(*,*) '--------------------WARNING!--------------------'
				write(*,*) 'Inpit matrix is not a diagonally dominant matrix'
				write(*,*) '--------------------WARNING!--------------------'
				write(*,*)
                exit
			endif
		enddo
		
		x_new=0.0                                                       ! first initialization
		x_old=1.0
		d=0.0
		d_inv=0.0
		
		do i=1, Sys_size                                                
			d(i,i)=a(i,i)                                               ! d is diagonal matrix of a
			d_inv(i,i)=1.0/d(i,i)                                       ! d_inv is inverted matrix of d
		enddo
		
		z=matmul(d_inv,(d-a))                                           ! z=d_inv*(d-a)
		g=matmul(d_inv,b)                                               ! g=d_inv*b
        
		do while(norm2(x_new-x_old)>eps)                                ! iterative method
			
            x_old = x_new			
            x_new = matmul(z,x_old) + g                                 ! new solution calculated usinf previous
            
		enddo
		
	end subroutine Jacobi	
	
	subroutine Zeidel(a,b,x_new,Sys_size)                               ! Zeidel method
	
		integer :: i, j, Sys_size 	                                    ! indexes and system size
		real(8) :: interim                                              ! intermediate sum
		real(8), dimension(:,:) :: a							        ! system matrix
		real(8), dimension(:) :: b			                            ! column vector
		real(8), dimension(Sys_size) ::  x_new, x_old		            ! new end previous vectors
		
		do i=1, Sys_size
			if( SUM(abs(a(i,:))) > 2*abs(a(i,i))) then                       ! check if system matrix is diagonally mdomaint matrix
				write(*,*) '--------------------WARNING!--------------------'
				write(*,*) 'Inpit matrix is not a diagonally dominant matrix'
				write(*,*) '--------------------WARNING!--------------------'
				write(*,*)
                exit
			endif
		enddo
		
		x_old=0.0                                                       ! first initialization
		x_new=0.0
		
		do while (.TRUE.)                                               ! cycle with postpredict
		
			x_old = x_new
			
			do i=1,Sys_size                                             ! iteration cycle
			
				interim=0.0
			
				do j=1,i-1
					interim = interim + a(i,j)*x_new(j)                 ! lower matrix with new x 
				enddo
				
				do j=i+1,Sys_size
					interim = interim + a(i,j)*x_old(j)                 ! upper matrix with precious x
				enddo
				
				x_new(i) = ( b(i) - interim )/a(i,i)                    ! calculation of new x
				
			enddo
				
			if (norm2(x_new-x_old) < eps) then                          ! postpredict
				exit
			endif
			
		enddo
		
	end subroutine Zeidel
	
	
	subroutine Relax(a,b,x,Sys_size)                                    ! relaxation method
	
		integer :: i, j, Sys_size                                       ! indexes and system size
		real(8) :: maximum		                                        ! maximum element of residual vector
		real(8), dimension(:,:) :: a							        ! system matrix
		real(8), dimension(:) :: b	                                    ! column matrix
		real(8), dimension(Sys_size,Sys_size) :: P			            ! intermediate matrix
		real(8), dimension(Sys_size) ::  x, Q		                    ! residual vector and solution vector
		
		do i=1, Sys_size
			if( SUM(abs(a(i,:))) > 2*abs(a(i,i))) then
				write(*,*) '--------------------WARNING!--------------------'
				write(*,*) 'Inpit matrix is not a diagonally dominant matrix'
				write(*,*) '--------------------WARNING!--------------------'
				write(*,*)
                exit
			endif
		enddo
		
		X=0                                                             ! first initialization
		maximum=1
		
		do i=1, Sys_size                                                ! calculation of P matrix and Q vector
        
			do j=1, Sys_size
				P(i,j) = -a(i,j)/a(i,i)
			enddo
            
			Q(i) = b(i)/a(i,i)
            
		enddo
		
		do while (.TRUE.)                                               ! algorithm realisation
			
			if ( abs(maximum) < eps) then
				exit
			endif
			
			j = maxloc(abs(Q), dim=1)
			
			maximum=Q(j)
			
			x(j) = x(j) + maximum
			
			do i=1,Sys_size
				Q(i) = Q(i) + P(i,j)*maximum
			enddo
			
		enddo
		
	end subroutine Relax	
		
end module Linear_system_iter

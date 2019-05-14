program in_out
	use Linear_system_iter
	implicit none
	
!	-----------------------------------------------------
!	|	Program solves linear equation system			|
!	|			with iterration methods					|
!	|													|
!	|	Method types and keys:							|
!	|													|
!	|		Jacobi method -- Jacobi						|
!	|		Zeidel method -- Zeidel 					|
!	|		Successive over-relaxation -- Relax			|
!	|													|
!	----------------------------------------------------

	real(8), dimension(:,:), allocatable :: A							!system amtrix
	real(8), dimension(:), allocatable :: B, X							!column vector and solution vector
	integer Sys_size 													!square matrix size
	integer i,j,k														!string
	character*14 method_type											!Method key
	
	
	call getarg(1, method_type)											!Console method key input
	
	do while(method_type/='Jacobi' .and. method_type/='Zeidel' .and. method_type/='Relax')	!Corectness check
		write(*,*) 'Wrong input! Try Again!'
			write(*,*)
			write(*,*) 'Method types and keys:				'	
			write(*,*)
			write(*,*) 'Jacobi method -- Jacobi				'		
			write(*,*) 'Zeidel method -- Zeidel 			'		
			write(*,*) 'Successive over-relaxation -- Relax'
			write(*,*)
		read(*,*) method_type
	enddo
	
	
	open ( unit=1, file='data.dat', status='old', action='read' )
	open ( unit=3, file='result.dat',status='REPLACE')

	read(1,'(2x,I6)') Sys_size											!read sharp symbol and matrix size
	
	allocate(A(Sys_size,Sys_size))
	allocate(B(Sys_size),X(Sys_size))
		
	read(1,*) A															!read system matrix
	A=TRANSPOSE(A)
	
	do i=1, Sys_size
		read(1,*) B(i)													!read column vector
	enddo
	
	select case (method_type)											!call selected method
		case ('Jacobi')
			call Jacobi(A,B,X,Sys_size)
		case ('Zeidel')
			call Zeidel(A,B,X,Sys_size)
		case ('Relax')
			call Relax(A,B,X,Sys_size)
	end select
	
	write(3,55) Sys_size												!print system size
	55 format('#',I5)													
	
	do j=1,Sys_size
				write(3,*) X(j)											!write solution vector
	enddo
	
	write(*,*) 'Norm of residual vector is:', norm2( matmul(A,X)-B )	!console output norm of residal vector
	
	close(1)
	close(3)
	
	deallocate(A,B,X)
	
end program in_out

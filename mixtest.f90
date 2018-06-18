program mixers

implicit none
	integer					:: i, j, t
	integer					:: ceil, grid, tim
	real					:: pressure
	integer, dimension(:,:), allocatable :: arrin, delta
	
	
write(*,*) "Grid size?"
read(*,*) grid
write(*,*) "Timesteps?"
read(*,*) tim

allocate(arrin(grid,grid))
allocate(delta(grid,grid))

arrin = 0
delta = 0
pressure = 0.5

do i = 1, grid, 1
	
	do j =  1, grid, 1
		
		write(*,*) i, j
		read(*,*) arrin(i,j)
	
	end do
	
end do

ceil = maxval(arrin)

write(*,*) "max value:", ceil
call printer(arrin,grid)
write(*,*) 

do t = 1, tim, 1

do i = 1, grid, 1
	
	do j = 1, grid, 1
		
		delta = 0
		
		if (i .lt. grid) then
	
			delta(i,j) = delta(i,j) + abs(floor(pressure*(arrin(i,j) - arrin(i+1,j))))
			
		end if
		
		if (i .ne. 1) then
		
			delta(i,j) = delta(i,j) + abs(floor(pressure*(arrin(i,j) - arrin(i-1,j))))
		
		end if
		
		if (j .lt. grid) then
			
			delta(i,j) = delta(i,j) + abs(floor(pressure*(arrin(i,j) - arrin(i,j+1))))
			
		end if
		
		if (j .ne. 1) then
			
			delta(i,j) = delta(i,j) + abs(floor(pressure*(arrin(i,j) - arrin(i,j-1))))
			
		end if
		
			arrin(i,j) = arrin(i,j) + delta(i,j)
		
		if (arrin(i,j) .gt. ceil) then
			arrin(i,j) = ceil
		end if
		
	end do
	
end do

call printer(arrin,grid)
write(*,*) 

end do


end program

subroutine printer(arr,big)
	
! Printing subroutine for arrays of two dimensions

implicit none
	integer,dimension(big,big),intent(in)			:: arr	! Input matix
	integer,intent(in)							:: big	! Size of matrix
	integer										:: i, j	! Looping integers

! Format statement
50 format(5g11.5)		

! Writing loop
do i =1,big
	write(*,50)(arr(i,j),j=1,big)
end do

end subroutine

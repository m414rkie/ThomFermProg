program diffuse

implicit none
	integer							:: dime
	real, allocatable				:: arrin(:,:)
	real, allocatable				:: delta(:,:)
	integer							:: i, j, t
	real							:: diffco, tot
	
	
write(*,*) "input dime:"
read(*,*) dime

allocate(arrin(dime,dime))
allocate(delta(dime,dime)) 

do i = 1, dime, 1
	do j = 1, dime, 1
		write(*,*) "input coordinate", i, j
		read(*,*) arrin(i,j)
	end do
end do

tot = sum(arrin)

write(*,*) "total:", tot
write(*,*) "timesteps?"
read(*,*) t


	
diffco = 0.1
delta = 0.0

do t = 1, t, 1

delta = 0.0

do i = 1, dime, 1
	
	do j = 1, dime, 1
		
		if ((arrin(i,j) .gt. arrin(i+1,j)) .and. (i .lt. dime)) then
	
			delta(i+1,j) = delta(i+1,j) +  diffco*(arrin(i,j) - arrin(i+1,j))
			delta(i,j) = delta(i,j) - diffco*(arrin(i,j) - arrin(i+1,j))
		
		end if
		
		if ((arrin(i,j) .gt. arrin(i-1,j)) .and. (i .ne. 1)) then
		
			delta(i-1,j) = delta(i-1,j) + diffco*(arrin(i,j) - arrin(i-1,j))
			delta(i,j) = delta(i,j) - diffco*(arrin(i,j) - arrin(i-1,j))
		
		end if
		
		if ((arrin(i,j) .gt. arrin(i,j+1)) .and. (j .lt. dime)) then
			
			delta(i,j+1) = delta(i,j+1) + diffco*(arrin(i,j) - arrin(i,j+1))
			delta(i,j) = delta(i,j) - diffco*(arrin(i,j) - arrin(i,j+1))
			
		end if
		
		if ((arrin(i,j) .gt. arrin(i,j-1)) .and. (j .ne. 1)) then
			
			delta(i,j-1) = delta(i,j-1) + diffco*(arrin(i,j) - arrin(i,j-1))
			delta(i,j) = delta(i,j) - diffco*(arrin(i,j) - arrin(i,j-1))
			
		end if
		
	end do
	
end do





arrin = arrin + delta

tot = sum(arrin)
write(*,*) tot
call printer(arrin,dime)
write(*,*) 

end do


end program

subroutine printer(arr,big)
	
! Printing subroutine for arrays of two dimensions

implicit none
	real,dimension(big,big),intent(in)			:: arr	! Input matix
	integer,intent(in)							:: big	! Size of matrix
	integer										:: i, j	! Looping integers

! Format statement
50 format(5g11.5)		

! Writing loop
do i =1,big
	write(*,50)(arr(i,j),j=1,big)
end do

end subroutine


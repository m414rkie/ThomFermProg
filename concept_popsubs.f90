subroutine fishdist(arrout)

! Subroutine for populating fish layer. Currently an even layer is created using coral 'biomass' 
! as an initializer with a multiplier.

use globalvars

implicit none
	real,dimension(grid,grid)			 	:: arrin				! Input array
	real,dimension(grid,grid),intent(out)	:: arrout				! Output array
	real									:: coraltot, fishtot	! Summed values of the arrays
	real									:: coralfishmult = 1.5	! Multiplies coral 'biomass' 

! Initializing 
coraltot = sum(coral)
fishtot = coralfishmult*coraltot

! Distribution across grid
arrout = fishtot/(grid**2)

write(*,*) "Populating the initial fish layer."

end subroutine

subroutine hppop(arrin)

! Generates a population of coral and algae. Algae is represented by 0 in a gridpoint. 

use globalvars

implicit none
	real,dimension(grid,grid)			    :: arrin			! Input arrays
	integer									:: i, j				! Looping integers
	
write(*,*) "Populating the initial coral layer."
! Generates a seed for use in populating layers

		call system_clock(count=clock)
		seed = clock + 8*(/(i-1,i=1,randall)/)
		call random_seed(put=seed)
		call random_number(arrin)

! Removes coral from the grid based on the percent cover of the total grid. This also determines the 
! amount of algae.
!	do i = 1, grid, 1
!		
!		do j = 1, grid, 1
!		
!			if (arrin(i,j) .lt. (1.0-percentcover)) then
!				arrin(i,j) = 0.0
!			end if
!		
!		end do
!		
!	end do

where (arrin .lt. (1.0-percentcover)) arrin = 0.0

end subroutine

subroutine tightcluster(arrin)
	
! Generates additional coral as circular clusters with a user-input grid radius at random points on the grid.
! Will likely be used as a base to generate more linear distributions.
! Does not impose coral growth on algae areas
	
use globalvars

	real,dimension(grid,grid)				:: arrin					! Input array
	integer									:: i,j,x,y, k				! Looping integers
	integer									:: counter					! Lowers the value of cluster as it gets further from center
	real									:: temp(1:2)				! Holds randomly generated numbers 
	integer									:: coordinate(1:2)			! Holds x,y coordinates of center of cluster
	real									:: tightclustermult = 2.0	! Determines the increase in coral in cluster
	real									:: disttrail, far			! Spreads the increase across the cluster.
																		!  interacts with counter to linearly decrease the 
																		!  increase in coral with distance from center


! Initializations
disttrail = tightclustermult/real(distance)
counter = 0

!call random_seed(size=randall)
call system_clock(count_rate=clock)
seed = clock + 34*(/(i-1,i=1,randall)/)	

do k=1, clusnum, 1

		call random_seed(put=seed)
		call random_number(temp)
		
		seed = seed*2 - 6
		
coordinate = floor(grid*temp)

x = coordinate(1)
y = coordinate(2)	
	
write(*,*) "Cluster at:", coordinate(1), coordinate(2)

	arrin(x,y) = arrin(x,y)*tightclustermult

	do j = -distance, distance, 1
	
		do i = -distance, distance, 1
		
			far = (distance - counter)
			if (far .lt. 0) then
				far = 0
			end if
		
			if ((y+j) .le. grid) then
				arrin(x,y+j) = arrin(x,y+j) + arrin(x,y+j)*disttrail*real(far)
			end if
	
			if ((y-j) .gt. 0) then
				arrin(x,y-j) = arrin(x,y-j) + arrin(x,y-j)*disttrail*real(far)
			end if	

			if ((x+i) .le. grid) then
				arrin(x+i,y) = arrin(x+i,y) + arrin(x+i,y)*disttrail*real(far)	
			end if
		
			if ((x-i) .gt. 0) then
				arrin(x-i,y) = arrin(x-i,y) + arrin(x-i,y)*disttrail*real(far)
			end if

			if (((x+j) .le. grid) .and. ((y+i) .le. grid) .and. ((x+j) .ge. 0) .and. ((y+i) .ge. 0)) then
				arrin(x+j,y+i) = arrin(x+j,y+i) + arrin(x+j,y+i)*disttrail*real(far)*0.707
			end if

			counter = counter + 1

		end do

		counter = 0

	end do

end do
	
end subroutine
		
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		
subroutine newcoral

use globalvars

implicit none
	real			:: coraltot, area
	real			:: avgcoral
	real			:: temp
	real			:: coord(1:2)
	integer			:: x, y, i
	
	
coraltot  = sum(coral)
area 	  = float(grid)**2
avgcoral  = coraltot/area

if (avgcoral .ge. threshold) then
	
	call system_clock(count=clock)
	seed = clock + 3*(/(i-1,i=1,randall)/)
	call random_seed(put=seed)
	call random_number(temp)

	if (temp .ge. 0.7) then

		seed = seed*2
		call random_seed(put=seed)
		call random_number(coord)
				
		x = floor(grid*coord(1))
		y = floor(grid*coord(2))
		
		numnew = numnew + 1
		
		write(*,*) "New coral growth at:", x, y

		if (coral(x,y) .eq. 0.0) then
			coral(x,y) = 1.2
		else
			x = x+1 ; y = y+1
				
				if (x .gt. grid) then
					x = x - floor(0.5*coord(1))
				end if	
			
				if (y .gt. grid) then
					y = y - floor(0.5*coord(2))
				end if
			
			coral(x,y) = 1.2
					
		end if
		
	end if

end if
		
end subroutine	
		
subroutine bacteriapop

use globalvars
use functions

implicit none
	real 				:: avgspec, ran
	real				:: coord(2)
	real				:: average, area
	integer				:: i, j
	real				:: deviation
		
write(*,*) "Populating initial Bacteria layer."

write(*,*) "Average number of bacteria species?"
read(*,*) avgspec

write(*,*) "Maximum number of bacteria species?"
read(*,*) maxspec

write(*,*) "Deviation of species?"
read(*,*) deviation

avgpop = 1000.0
area = (2*float(grid))**2

bacteria%totalpop = int(avgpop)
bacteria%numspecies = 1

	call random_seed(size=randall)
	call system_clock(count=clock)
	seed = clock + 4*(/(i-1,i=1,randall)/)
	call random_seed(put=seed)


do i = 1, 2*grid, 1
	
	do j = 1, 2*grid, 1
	
			call random_number(ran)
			
			ran = floor(maxspec*ran + 3.0*deviation*(1.0-ran))
			
			bacteria(i,j)%numspecies = (int(ran))

	end do

end do

average = sum(bacteria%numspecies)/area
write(*,*) average


open(unit=14,file="bactlayer.dat",status="replace",position="append")
	
	do i = 1, 2*grid, 1
		do j = 1, 2*grid, 1
			write(14,*) i, j, bacteria(i,j)%numspecies, bacteria(i,j)%totalpop
		end do
	end do

close(14)

end subroutine



	
	
	
	
	
	
	
	
	
	
	




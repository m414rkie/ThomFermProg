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
		seed = clock + 34*(/(i-1,i=1,randall)/)
		call random_seed(put=seed)
		call random_number(arrin)

! Removes coral from the grid based on the percent cover of the total grid. This also determines the 
! amount of algae.
	do i = 1, grid, 1
		
		do j = 1, grid, 1
		
			if (arrin(i,j) .lt. (1.0-percentcover)) then
				arrin(i,j) = 0.0
			end if
		
		end do
		
	end do

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
	real									:: disttrail				! Spreads the increase across the cluster.
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
		
		seed = seed*2 - 600
		
coordinate = floor(grid*temp)

x = coordinate(1)
y = coordinate(2)	
	
write(*,*) "Cluster at:", x, y

if (arrin(x,y) .ne. 0.0) then
	arrin(x,y) = arrin(x,y)*tightclustermult*0.7
end if


	do j = x, x + distance, 1
		
		if ((arrin(j,y) .ne. 0.0) .and. (j .le. grid)) then
		
			arrin(x,y+j) = arrin(x,y+j) + arrin(x,y+j)*disttrail*real(distance-counter)
			arrin(x,y-j) = arrin(x,y-j) + arrin(x,y-j)*disttrail*real(distance-counter)
		
		end if

		counter = counter + 1

	end do
	
	counter = 0
	
	do i = y, y + distance, 1
	
		if ((arrin(x,i) .ne. 0.0) .and. (i .le.grid)) then
		
			arrin(x+i,y) = arrin(x+i,y) + arrin(x+i,y)*disttrail*real(distance-counter)
			arrin(x-i,y) = arrin(x-i,y) + arrin(x-i,y)*disttrail*real(distance-counter)
	
		end if
		
		counter = counter + 1		
		
	end do
	
end do
	
end subroutine
		
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		
subroutine newcoral

use globalvars

implicit none
	real			:: coraltot, area
	real			:: avgcoral, threshold
	real			:: temp
	real			:: coord(1:2)
	integer			:: x, y, i
	
	
coraltot  = sum(coral)
area 	  = float(grid)**2
avgcoral  = coraltot/area
threshold = 3.0



if (avgcoral .ge. threshold) then
	
	call system_clock(count=clock)
	seed = clock + 34*(/(i-1,i=1,randall)/)
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
	real 				:: avgspec
	real				:: coord(2)
	real				:: average, area
	integer				:: i, j
	
write(*,*) "Populating initial Bacteria layer."
avgspec = 100.0
avgpop = 1000.0
area = (2*float(grid))**2
average = 0.0

bacteria%totalpop = int(avgpop)
bacteria%numspecies = 1

	call random_seed(size=randall)
	call system_clock(count=clock)
	seed = clock + 34*(/(i-1,i=1,randall)/)

do while (average .lt. avgspec)
	
			call random_seed(put=seed)
			call random_number(coord)
			
			coord = 2*grid*coord
			
			bacteria(floor(coord(1)),floor(coord(2)))%numspecies = bacteria(floor(coord(1)),floor(coord(2)))%numspecies + 1
			
			average = float(sum(bacteria%numspecies))/area

			seed = seed*3 - 50
end do

open(unit=14,file="bactlayer.dat",status="replace",position="append")
	
	do i = 1, 2*grid, 1
		do j = 1, 2*grid, 1
			write(14,*) i, j, bacteria(i,j)%numspecies, bacteria(i,j)%totalpop
		end do
	end do

close(14)

end subroutine



	
	
	
	
	
	
	
	
	
	
	




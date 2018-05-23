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

end subroutine

subroutine hppop(arrin)

! Generates a population of coral and algae. Algae is represented by 0 in a gridpoint. 

use globalvars

implicit none
	real,dimension(grid,grid)			    :: arrin			! Input arrays
	integer									:: i, j				! Looping integers
	
! Populates layer with random numbers between zero and one.	
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
	integer									:: i,j,x,y					! Looping integers
	integer									:: counter					! Lowers the value of cluster as it gets further from center
	real									:: temp(1:2)				! Holds randomly generated numbers 
	integer									:: coordinate(1:2)			! Holds x,y coordinates of center of cluster
	real									:: tightclustermult = 2.5	! Determines the increase in coral in cluster
	real									:: disttrail				! Spreads the increase across the cluster.
																		!  interacts with counter to linearly decrease the 
																		!  increase in coral with distance from center


! Initializations
disttrail = tightclustermult/real(distance)
counter = 0

! Random  grid point generation initialization
call random_seed(put=seed)
call random_number(temp)
coordinate = grid*floor(temp)

x = coordinate(1)
y = coordinate(2)	
	
	do j = x, x + distance, 1
		
		arrin(x,y+j) = arrin(x,y+j) + arrin(x,y+j)*disttrail*(distance-counter)
		arrin(x,y-j) = arrin(x,y-j) + arrin(x,y-j)*disttrail*(distance-counter)
		
		counter = counter + 1
	
	end do
	
	counter = 0
	
	do i = y, y + distance, 1
	
		arrin(x+i,y) = arrin(x+i,y) + arrin(x+i,y)*disttrail*(distance-counter)
		arrin(x-i,y) = arrin(x-i,y) + arrin(x-i,y)*disttrail*(distance-counter)
		
		counter = counter + 1
		
	end do
	
end subroutine
		
		
		
		
		




	
	
	
	
	
	
	
	
	
	
	




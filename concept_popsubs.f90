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
coordinate = floor(grid*temp)

x = coordinate(1)
y = coordinate(2)	
	
if (arrin(x,y) .ne. 0.0) then
	arrin(x,y) = arrin(x,y)*tightclustermult*0.7
end if


	do j = x, x + distance, 1
		
		if (arrin(j,y) .ne. 0.0) then
		
			arrin(x,y+j) = arrin(x,y+j) + arrin(x,y+j)*disttrail*real(distance-counter)
			arrin(x,y-j) = arrin(x,y-j) + arrin(x,y-j)*disttrail*real(distance-counter)
		
		end if

		counter = counter + 1

	end do
	
	counter = 0
	
	do i = y, y + distance, 1
	
		if (arrin(x,i) .ne. 0.0) then
		
			arrin(x+i,y) = arrin(x+i,y) + arrin(x+i,y)*disttrail*real(distance-counter)
			arrin(x-i,y) = arrin(x-i,y) + arrin(x-i,y)*disttrail*real(distance-counter)
	
		end if
		
		counter = counter + 1		
		
	end do
	
end subroutine
		
subroutine newcoral

use globalvars

implicit none
	real			:: coraltot, area
	real			:: avgcoral, threshold
	real			:: temp
	real			:: coord(1:2)
	integer			:: x, y
	
	
coraltot  = sum(coral)
area 	  = grid**2
avgcoral  = coraltot/area
threshold = 4.0

if (avgcoral .ge. threshold) then
	
	call random_seed(put=seed)
	call random_number(temp)

	if (temp .ge. 0.7) then
		call random_seed(put=seed)
		call random_number(coord)
		x = floor(grid*coord(1))
		y = floor(grid*coord(2))
		numnew = numnew + 1
		
102		if (coral(x,y) .eq. 0.0) then
			coral(x,y) = 1.2
		else
			x = x+1 ; y = y+1
				
				if (x .gt. grid) then
					x = x - grid
				end if	
			
				if (y .gt. grid) then
					y = y - grid
				end if
			
			goto 102
		
		end if
		
	end if

end if
		
end subroutine	
		
subroutine bacteriapop

use globalvars

implicit none
	real		:: avgpop
	integer		:: i, j
	real		:: dist(10), randini(10)
	integer		:: randfin(10)
	real		:: r
	integer		:: in
	
write(*,*) "Please input the average population per (volume) as  (num)e(exp), typicall values would be 1.0e12"
read(*,*) avgpop

allocate(bacteria(4*grid,4*grid))

dist = (/0.96,0.97,0.98,0.99,1.0,1.1,1.2,1.3,1.4,1.5) 

call random_number(randini)
randfin = floor(10.0*randini)

bacteria%totalpop = avgpop

do i = 1, 4*grid, 1

	do j = 1, 4*grid, 1
		
		call random_number(r)
		in = floor(10*r)
		
		bacteria(i,j)%totalpop = bacteria%totalpop*in
		bacteria(i,j)%numspecies = numspec(bacteria(i,j)%totalpop)
		
	end do

end do
		




end subroutine



	
	
	
	
	
	
	
	
	
	
	




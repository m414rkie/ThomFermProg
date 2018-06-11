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
	
write(*,*) "Cluster at:", x, y

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
threshold = 1.0

if (avgcoral .ge. threshold) then
	
	call random_seed(put=seed)
	call random_number(temp)

	if (temp .ge. 0.7) then
		call random_seed(put=seed)
		call random_number(coord)
		x = floor(grid*coord(1))
		y = floor(grid*coord(2))
		numnew = numnew + 1
		write(*,*) "New coral growth at:", x, y
		
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
use functions

implicit none
	real (kind=16)		:: avgpop, avgspec
	integer				:: i, j, k, l
	real				:: randini
	integer				:: randfin(10)
	real				:: r, factor, deviation, average
	integer				:: in, percent, maxspec
	real (kind=16)		:: spectopop(500000), probability(500000)

factor = 2.0
maxspec = 500000
spectopop(1) = factor
avgspec = 350000.0
deviation = 25000.0
average = 350000.0

do k = 2, maxspec, 1
	spectopop(k) = 1.00004**k + spectopop(k-1)
end do

do l = 1, 500000, 1

probability(l) = (1.0/((2.0*pi*(deviation**2))**(1.0/2.0)))*exp(-((float(l)-average)**2)/(2.0*(deviation**2)))

end do

do i = 1, 2*grid, 1
	
	do j = 1, 2*grid
		
		call random_number(randini)

		randini = 500000.0*randini + 200000.0*(1.0-randini)

		bacterialayer(i,j)%numpop = probability(floor(randini))
		
	end do

end do


end subroutine



	
	
	
	
	
	
	
	
	
	
	




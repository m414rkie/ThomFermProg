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

end subroutine

subroutine tightcluster(arrin)
	
! Generates additional coral as circular clusters with a user-input grid radius at random points on the grid.
! Will likely be used as a base to generate more linear distributions.
! Does not impose coral growth on algae areas
	
use globalvars

	real,dimension(grid,grid)				:: arrin					! Input array
	integer									:: i,j,x,y, k				! Looping integers
	integer									:: counter					! Lowers the value of cluster as it gets further from center
	real, allocatable						:: coordinate(:,:)			! Holds x,y coordinates of center of cluster
	real									:: tightclustermult = 2.0	! Determines the increase in coral in cluster
	real									:: disttrail, far			! Spreads the increase across the cluster.
																		!  interacts with counter to linearly decrease the 
																		!  increase in coral with distance from center

allocate(coordinate(2,clusnum))

! Initializations
disttrail = tightclustermult/real(distance)
counter = 0

!call random_seed(size=randall)
call system_clock(count_rate=clock)
seed = clock + 34*(/(i-1,i=1,randall)/)	
call random_seed(put=seed)
call random_number(coordinate)

coordinate = grid*coordinate

do k = 1, clusnum, 1

x = floor(coordinate(1,k)) + 1
y = floor(coordinate(2,k)) + 1	
	
write(*,*) "Cluster at:", x, y

	arrin(x,y) = arrin(x,y) + tightclustermult

	do j = -distance, distance, 1
	
		do i = -distance, distance, 1
		
			far = (distance - counter)
			
			if (far .lt. 0) then
				far = 0
			end if
		
			if ((y+j) .le. grid) then
				arrin(x,y+j) = arrin(x,y+j) + tightclustermult*disttrail*real(far)
			end if
	
			if ((y-j) .gt. 0) then
				arrin(x,y-j) = arrin(x,y-j) + tightclustermult*disttrail*real(far)
			end if	

			if ((x+i) .le. grid) then
				arrin(x+i,y) = arrin(x+i,y) + tightclustermult*disttrail*real(far)	
			end if
		
			if ((x-i) .gt. 0) then
				arrin(x-i,y) = arrin(x-i,y) + tightclustermult*disttrail*real(far)
			end if

			if (((x+j) .le. grid) .and. ((y+i) .le. grid) .and. ((x+j) .ge. 0) .and. ((y+i) .ge. 0)) then
				arrin(x+j,y+i) = arrin(x+j,y+i) +tightclustermult*disttrail*real(far)*0.707
			end if

			counter = counter + 1

		end do

		counter = 0

	end do

end do
	
where (arrin .lt. (1.0-percentcover)) arrin = 0.0
	
end subroutine
		
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		
subroutine newcoral

use globalvars

implicit none
	real					:: coraltot, area
	real					:: avgcoral
	real					:: temp
	real					:: coord
	integer					:: x, y, i, j, l, c
	integer,allocatable		:: algaeloc(:,:)
	
	
coraltot  = sum(coral)
area 	  = float(grid)**2
avgcoral  = coraltot/area
c = 0

check = (coral .eq. 0.0)

l = count(check)
write(*,*) "l", l

allocate(algaeloc(2,l))
algaeloc = 0

do i = 1, grid, 1
	
	do j = 1, grid, 1
		
		if (coral(i,j) .eq. 0.0) then
			
			c = c + 1
			algaeloc(1,c) = i
			algaeloc(2,c) = j
		
		end if
		
	end do
	
end do

write(*,*) "c", c

if (avgcoral .ge. threshold) then
	
	call system_clock(count=clock)
	seed = clock + 3*(/(i-1,i=1,randall)/)
	call random_seed(put=seed)
	call random_number(temp)

	if (temp .ge. 0.4) then

		seed = seed*2
		call random_seed(put=seed)
		call random_number(coord)
				
		x = algaeloc(1,floor(l*coord))
		y = algaeloc(2,floor(l*coord))
		
		numnew = numnew + 1
		
		write(*,*) "New coral growth at:", x, y

			coral(x,y) = 1.2
		
	end if

end if
		
deallocate(algaeloc)
		
end subroutine	
		
subroutine bacteriapop

use globalvars
use functions

implicit none
	real 				:: avgspec, ran, minispec
	real				:: coord(2)
	real				:: average, area
	integer				:: i, j
	real				:: deviation
		
write(*,*) "Populating initial Bacteria layer."

write(*,*) "Maximum number of bacteria species?"
read(*,*) maxspec

write(*,*) "Minimum number of species?"
read(*,*) minispec

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
			
			ran = floor(maxspec*ran + minispec*(1.0-ran))
			
			bacteria(i,j)%numspecies = (int(ran))

	end do

end do

average = sum(bacteria%numspecies)/area
write(*,*) "Average umber of species:" ,average


open(unit=14,file="bactlayer.dat",status="replace",position="append")
	
	do i = 1, 2*grid, 1
		do j = 1, 2*grid, 1
			write(14,*) i, j, bacteria(i,j)%numspecies, bacteria(i,j)%totalpop
		end do
	end do

close(14)

end subroutine



	
	
	
	
	
	
	
	
	
	
	




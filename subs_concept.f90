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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine neighborsum(x,y,arr,total)
	
!Sums together the neighboring values of a given array point. Used for testing. 
	
use globalvars

implicit none
	integer, intent(in)						:: x, y									! position of current gridpoint
	real,dimension(grid,grid), intent(in)	:: arr									! array of the moment
	real, intent(out)						:: total								! sum of nearest neighbors
	integer									:: a, b, c, d, e, f, g, h				! These sets of integers are used 
	integer									:: i, j, k, l, m, n, o, p, q, r, s, t	!  to determine the summed value
																					!  of nearest neighbors of input
	
! Initializations
norm = 8.0
total = 0.0
test = 0.0
a = 1 ; b = 1 ; c = 1 ; d = 1
e = 1 ; f = 1 ; g = 1 ; h = 1
i = x+1 ; j = x-1 ; k = y+1 ; l = y-1
m = x+1 ; n = x-1 ; o = x+1 ; p = x-1
q = y+1 ; r = y-1 ; s = y-1 ; t = y+1 

		! Logic statements to catch out-of-bounds
		if ((x+1) .gt. grid) then
			a = 0 ; i = 1
			norm = norm - 1.0
		end if
			
		if ((x-1) .lt. 1) then
			b = 0 ; j = 1
			norm = norm - 1.0
		end if 
		
		if ((y+1) .gt. grid) then
			c = 0 ; k = 1
			norm = norm - 1.0
		end if
	 
		if ((y-1) .lt. 1) then 
			d = 0 ; l = 1
			norm = norm - 1.0
		end if
 
 		if (((x+1) .gt. grid) .or. ((y+1) .gt. grid)) then
			e = 0 ; m = 1 ; q = 1
			norm = norm - 1.0
		end if
 
 		if (((x-1) .lt. 1) .or. ((y-1) .lt. 1)) then
			f = 0 ; n = 1 ; r = 1
			norm = norm - 1.0
		end if
		 
		if (((x+1) .gt. grid) .or. ((y-1) .lt. 1)) then
			g = 0 ; o = 1 ; s = 1
			norm = norm - 1.0
		end if

		if (((x-1) .lt. 1) .or. ((y+1) .gt. grid)) then
			h = 0 ; p = 1 ; t = 1
			norm = norm -1.0
		end if
				
		total = (a*arr(i,y)+b*arr(j,y)+c*arr(x,k)+d*arr(x,l)+e*arr(m,q)+ &
				f*arr(n,r)+g*arr(o,s)+h*arr(p,t))/norm

end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine growth(x,y,arrin,arrout)
	
! Grows the input grid location based on value and neighbors.
! 
!(straight percentage at this point). 
	
use globalvars
use functions
	
implicit none
	integer, intent(in)							:: x, y					! Input coordinates
	real,dimension(grid,grid), intent(in) 		:: arrin				! Input array
	real,dimension(grid,grid), intent(out)		:: arrout				! Output array
	real										:: fishpop				! Fish population of a gridpoint
	real										:: growpercent = 1.1	! Flat percentage growth for coral
	
! Grows coral	
arrout(x,y) = arrin(x,y)*growpercent

! Initializations
fishlocal = 0.0
fishpop = fish(x,y)

	
	! Checks neighboring gridpoints for fish population and grows faster with neighbors.
	! Also checks boundaries to help with edge effects.
	if(holding(x,y) .eq. 0.0) then
		fishlocal = fishlocal - 1.0
	else
		fishlocal = fishlocal + holding(x,y)*0.1	
	end if
	
	if ((holding(x-1,y) .eq. 0.0) .and. (x .gt. 1)) then
		fishlocal = fishlocal + 0.1
	end if
	
	if ((holding(x+1,y) .eq. 0.0) .and. (x .lt. grid)) then
		fishlocal = fishlocal + 0.1
	end if

	if ((holding(x,y+1) .eq. 0.0) .and. (y .lt. grid)) then
		fishlocal = fishlocal + 0.1
	end if
	
	if ((holding(x,y-1) .eq. 0.0) .and. (y .gt. 1)) then
		fishlocal = fishlocal + 0.1
	end if

	if ((holding(x+1,y+1) .eq. 0.0) .and. (x .lt. grid) .and. (y .lt. grid)) then
		fishlocal = fishlocal + 0.1
	end if

	if ((holding(x-1,y-1) .eq. 0.0) .and. (x .gt. 1) .and. (y .gt. 1)) then
		fishlocal = fishlocal + 0.1
	end if
	
	if ((holding(x+1,y-1) .eq. 0.0) .and. (x .lt. grid) .and. (y .gt. 1)) then
		fishlocal = fishlocal + 0.1
	end if
	
	if ((holding(x-1,y+1) .eq. 0.0) .and. (x .gt. 1) .and. (y .lt. grid)) then
		fishlocal = fishlocal + 0.1
	end if


! Finalizes the population growth of fish, faster with more coral.
fish(x,y) = fish(x,y) + fishdelta(sum(coral),fish(x,y))

end subroutine
	
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
subroutine decay(x,y,arrin)
	
! Represents the algae killing the coral

use globalvars
	
implicit none
	integer,intent(in)			:: x, y			! Input coordinates
	real,dimension(grid,grid)	:: arrin		! Input array
	real						:: algcount		! Amount of algae near input coordinates
	real						:: decayconst	! Percent of coral loss per nearby algae 
	
	! Initializations
	algcount = 0.0
	decayconst = 0.08
	
	
	! Checks for algae around the input gridpoint and out-of-bounds
	if ((holding(x-1,y) .eq. 0.0) .and. (x .gt. 1)) then
		algcount = algcount + 1.0
	end if
	
	if ((holding(x+1,y) .eq. 0.0) .and. (x .lt. grid)) then
		algcount = algcount + 1.0
	end if

	if ((holding(x,y+1) .eq. 0.0) .and. (y .lt. grid)) then
		algcount = algcount + 1.0
	end if
	
	if ((holding(x,y-1) .eq. 0.0) .and. (y .gt. 1)) then
		algcount = algcount + 1.0
	end if

	if ((holding(x+1,y+1) .eq. 0.0) .and. (x .lt. grid) .and. (y .lt. grid)) then
		algcount = algcount + 1.0
	end if

	if ((holding(x-1,y-1) .eq. 0.0) .and. (x .gt. 1) .and. (y .gt. 1)) then
		algcount = algcount + 1.0
	end if
	
	if ((holding(x+1,y-1) .eq. 0.0) .and. (x .lt. grid) .and. (y .gt. 1)) then
		algcount = algcount + 1.0
	end if
	
	if ((holding(x-1,y+1) .eq. 0.0) .and. (x .gt. 1) .and. (y .lt. grid)) then
		algcount = algcount + 1.0
	end if
	
	! Calls the fish layer to reduce the effect of algae on coral
	call fishinteraction(decayconst,x,y)
	
	! Coral being eaten.
	if (decayconst .lt. 0.0) then
		decayconst = 0.0
	end if
	
	! Coral less the algae eating it
	arrin(x,y) = arrin(x,y) - decayconst*algcount

	! Resets negative values to zero
	if (arrin(x,y) .le. 0.1) then
		arrin(x,y) = 0.0
	end if

end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine fishinteraction(modify,i,j)

! Interaction of fish with algae layer. Lessens the pressure of algae against the fish.

use globalvars

	real						:: modify				! Input variable to be modified
	integer,intent(in)			:: i, j					! Looping integers
	real						:: fisheat = 0.05		! Lowers input variable based on how much nearby fish eat the algae
	
	
	! Checks for fish around algae and lowers the amount of coral destroyed by the algae
	if(fish(i,j) .ne. 0.0) then
		modify = modify - fish(i,j)*fisheat
	end if
	
	if ((fish(i-1,j) .ne. 0.0) .and. (i .gt. 1)) then
		modify = modify - fish(i-1,j)*fisheat
	end if
	
	if ((fish(i+1,j) .eq. 0.0) .and. (i .lt. grid)) then
		modify = modify - fish(i+1,j)*fisheat
	end if

	if ((fish(i,j+1) .eq. 0.0) .and. (j .lt. grid)) then
		modify = modify - fish(i,j+1)*fisheat
	end if
	
	if ((fish(i,j-1) .eq. 0.0) .and. (j .gt. 1)) then
		modify = modify - fish(i,j-1)*fisheat
	end if

	if ((fish(i+1,j+1) .eq. 0.0) .and. (i .lt. grid) .and. (j .lt. grid)) then
		modify = modify - fish(i+1,j+1)*fisheat
	end if

	if ((fish(i-1,j-1) .eq. 0.0) .and. (i .gt. 1) .and. (j .gt. 1)) then
		modify = modify - fish(i-1,j-1)*fisheat
	end if
	
	if ((fish(i+1,j-1) .eq. 0.0) .and. (i .lt. grid) .and. (j .gt. 1)) then
		fmodify = modify - fish(i+1,j-1)*fisheat
	end if
	
	if ((fish(i-1,j+1) .eq. 0.0) .and. (i .gt. 1) .and. (j .lt. grid)) then
		modify = modify - fish(i-1,j+1)*fisheat
	end if

end subroutine
	
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
subroutine printtofile(arrin)

! Prints to a file in x-y-z format

use globalvars

implicit none
	real,dimension(grid,grid),intent(in)		:: arrin
	integer 									:: i, j
	
open(unit=11,file=filename,status="replace",position="append")

do i = 1, grid, 1
	do j = 1, grid, 1
		write(11,*) i, j, arrin(i,j)
	end do
end do

close(11)

end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
subroutine kgrid

! Finds the carrying capacity of each gridpoint for the bacteria layer

use globalvars

implicit none
	integer							:: i, j								! Looping integers
	real							:: algaemod, coralmod, barriermod	! Variables for varying the carrying capacity
	real,dimension(2*grid,2*grid)	:: kdelta							! Change in carrying capacity


! Initializations 
kbact = avgpop
kdelta = 0.0
algaemod = 1.2
coralmod = 0.8
barriermod = 1.5

! Loops for initial set up, no barrier interaction
do i = 1, grid, 1
	
	do j = 1, grid, 1
	
		if (coral(i,j) .ne. 0) then 
			kbact(2*i-1,2*j-1) = avgpop*coralmod
			kbact(2*i-1,2*j) = avgpop*coralmod
			kbact(2*i,2*j) = avgpop*coralmod
			kbact(2*i,2*j-1) = avgpop*coralmod
		else if (coral(i,j) .eq. 0) then
			kbact(2*i-1,2*j-1) = avgpop*algaemod
			kbact(2*i,2*j) = avgpop*algaemod
			kbact(2*i,2*j-1) = avgpop*algaemod
			kbact(2*i-1,2*j) = avgpop*algaemod
		end if
			
	end do
	
end do

! Loops to determine barrier interaction
do i = 1, 2*grid, 1
	
	do j = 1, 2*grid, 1

		if ((i .lt. 2*grid) .and. (kbact(i,j) .ne. kbact(i+1,j))) then
			kdelta(i,j) = avgpop*barriermod - avgpop
		end if
		
		if ((j .gt. 2*grid) .and. (kbact(i,j) .ne. kbact(i,j+1))) then
			kdelta(i,j) = avgpop*barriermod - avgpop
		end if
	
	end do
	
end do

! Final updating of the layer
kbact = kbact + kdelta

	open(unit=17,file="kbact.dat",position="append",status="replace")
	
	do i = 1, 2*grid, 1
		do j = 1, 2*grid, 1	
			write(17,*) i,j,kbact(i,j)
		end do
	end do
	
	close(17)

end subroutine
	
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine bactgrow

! Grows the bacteria layer, both totalpop and species

use globalvars
use functions
	
implicit none
	integer		:: i, j						! Looping integers
	real		:: groperc, delbactpop		! determines new species and change in bacteria population 
	
! Initializations
delbactpop = 0.0
groperc = 0.0
	
do i = 1, 2*grid, 1
	
	do j = 1, 2*grid, 1
		
		! Finds change in population
		delbactpop = floor(bacgrowth(real(bacteria(i,j)%totalpop),real(bacteria(i,j)%numspecies),kbact(i,j)))
		
		! Determines how many new species show up
		groperc = delbactpop/real(bacteria(i,j)%totalpop)
		
		if (groperc .ge. 0.02) then
			bacteria(i,j)%numspecies = bacteria(i,j)%numspecies + 4
		else if ((groperc .gt. 0.01) .and. (groperc .lt. 0.2) .and. (groperc .gt. 0.0)) then
			bacteria(i,j)%numspecies = bacteria(i,j)%numspecies + 2
		else if ((groperc .gt. -0.02) .and.  (groperc .lt. 0.0))then
			bacteria(i,j)%numspecies = bacteria(i,j)%numspecies - 2
		else if (groperc .lt. -0.2) then
			bacteria(i,j)%numspecies = bacteria(i,j)%numspecies - 4
		end if

		bacteria(i,j)%totalpop = bacteria(i,j)%totalpop + int(delbactpop)
		
	end do

end do

! Trims the layer 
where (bacteria%numspecies .gt. maxspec) bacteria%numspecies = maxspec
where (bacteria%numspecies .lt. 0) bacteria%numspecies = 0
	
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
subroutine diffuse

! Diffusion subroutine. Primary driver is population pressure

use globalvars

implicit none
	real, dimension(2*grid,2*grid)				:: delta		! Holds overall change from diffusion
	integer										:: i, j			! Looping integers
	real										:: diffco		! Diffusion coefficient
	
! Initializations
diffco = 0.01
delta = 0.0
	
! Working loops
do i = 1, 2*grid, 1
	
	do j = 1, 2*grid, 1
		
		if ((i .lt. 2*grid) .and. (bacteria(i,j)%totalpop .gt. bacteria(i+1,j)%totalpop)) then
	
			delta(i+1,j) = delta(i+1,j) +  diffco*(bacteria(i,j)%totalpop - bacteria(i+1,j)%totalpop)
			delta(i,j) = delta(i,j) - diffco*(bacteria(i,j)%totalpop - bacteria(i+1,j)%totalpop)
		
		end if
		
		if ((i .gt. 1) .and. (bacteria(i,j)%totalpop .gt. bacteria(i-1,j)%totalpop)) then
		
			delta(i-1,j) = delta(i-1,j) + diffco*(bacteria(i,j)%totalpop - bacteria(i-1,j)%totalpop)
			delta(i,j) = delta(i,j) - diffco*(bacteria(i,j)%totalpop - bacteria(i-1,j)%totalpop)
		
		end if
		
		if ((j .lt. 2*grid) .and. (bacteria(i,j)%totalpop .gt. bacteria(i,j+1)%totalpop)) then
			
			delta(i,j+1) = delta(i,j+1) + diffco*(bacteria(i,j)%totalpop - bacteria(i,j+1)%totalpop)
			delta(i,j) = delta(i,j) - diffco*(bacteria(i,j)%totalpop - bacteria(i,j+1)%totalpop)
			
		end if
		
		if ((j .gt. 1) .and. (bacteria(i,j)%totalpop .gt. bacteria(i,j-1)%totalpop)) then
			
			delta(i,j-1) = delta(i,j-1) + diffco*(bacteria(i,j)%totalpop - bacteria(i,j-1)%totalpop)
			delta(i,j) = delta(i,j) - diffco*(bacteria(i,j)%totalpop - bacteria(i,j-1)%totalpop)
			
		end if
		
	end do
	
end do

! Final update and trim
bacteria%totalpop = bacteria%totalpop + int(delta)

where (bacteria%totalpop .lt. 0) bacteria%totalpop = 0

end subroutine
	
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine mixing

! Mixing subroutine for redistributing species

use globalvars

implicit none
	integer									:: i,j 			! Looping integers
	real									:: mixpress		! Pressure for mixing
	integer									:: delta		! Change in number of species
								
mixpress = 0.5

! Working loops
do i = 1, 2*grid, 1
	
	do j = 1, 2*grid, 1
		
		delta = 0
		
		if (i .lt. 2*grid) then
	
			delta = delta + bacteria(i,j)%numspecies - bacteria(i+1,j)%numspecies
			
		end if
		
		if (i .ne. 1) then
		
			delta = delta + bacteria(i,j)%numspecies - bacteria(i-1,j)%numspecies
		
		end if
		
		if (j .lt. 2*grid) then
			
			delta = delta + bacteria(i,j)%numspecies - bacteria(i,j+1)%numspecies
			
		end if
		
		if (j .ne. 1) then
			
			delta = delta + bacteria(i,j)%numspecies - bacteria(i,j-1)%numspecies
			
		end if
		
		if (delta .lt. 0) then
			delta = 0
		end if
		
			bacteria(i,j)%numspecies = bacteria(i,j)%numspecies + delta
		
		
	end do
	
end do

! Trimming 
where (bacteria%numspecies .gt. maxspec) bacteria%numspecies = maxspec

end subroutine




	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
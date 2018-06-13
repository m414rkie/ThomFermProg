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
	real										:: fishpop
	real										:: growpercent = 1.1
	
	arrout(x,y) = arrin(x,y)*growpercent
	
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
fish(x,y) = fish(x,y) + fishdelta(fishpop)
!write(*,*) fishpop, fishdelta(fishpop), fishlocal

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
	decayconst = 0.25
	
	
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
	arrin(x,y) = arrin(x,y) - decayconst*algcount

	! Resets negative values to zero
	if (arrin(x,y) .le. 0.0) then
		arrin(x,y) = 0.0
	end if

end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine fishinteraction(modify,i,j)

! Interaction of fish with algae layer. 

use globalvars

	real						:: modify				! Input variable to be modified
	integer,intent(in)			:: i, j					! Looping integers
	real						:: fisheat = 0.08		! Lowers input variable based on how much nearby fish eat the algae
	
	
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
	

	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
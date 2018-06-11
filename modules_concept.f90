module globalvars

implicit none

type bacteriavar
	sequence
	integer (kind=16)	:: totalpop
	integer 			:: numspecies
end type bacteriavar

	integer									:: grid									! Array size
	real									:: norm, nearsum, test					! Variables for interactions
	real, allocatable						:: holding(:,:), coral(:,:), fish(:,:)  ! Layer names
	type (bacteriavar) , allocatable		:: bacteria(:,:)						! Layer names
	integer, allocatable					:: seed(:)								! Random number holding array
	integer									:: clock, distance						! System time and radial distance of coral clusters
	character*20							:: filename								! Changes for what is being put into the file
	real									:: percentcover							! Percent of grid to have coral on it 'groundcover'
	real									:: fishconst, fishlocal, fgrowfact
	integer									:: numnew = 0
	real									:: popconstant
	real									:: pi = acos(-1.0)
	
end module


module functions

contains

real function fishdelta(input)

use globalvars

implicit none
	real		:: input
	
	fishdelta = -fgrowfact*(input - fishconst) + fishlocal
	
end function fishdelta

real function bacgrowth(totalpop,specpop,carry)

use globalvars

implicit none
	real		:: totalpop, specpop, carry
	real		:: rate	
	
	bacgrowth = rate*(1.0 - (totalpop/carry))*specpop
	
end function bacgrowth

integer (kind=16) function numspec(pop)

use globalvars

implicit none
	integer	(kind=16)	:: pop
	real (kind=16)		:: summed, variable
	integer		:: i
	
summed  = 0.0
variable = 1.00004
popconstant = 1.0
	
do i = 1, 500000, 1

	summed = summed + popconstant*variable**i
	
	if (summed .gt. pop) then
		numspec = i-1
		goto 102
	end if
	
102 end do

end function











end module
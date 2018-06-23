module globalvars

implicit none

type bacteriavar
	sequence
	integer				:: totalpop
	integer 			:: numspecies
end type bacteriavar

	integer									:: grid, clusnum						! Array size
	real									:: norm, nearsum, test					! Variables for interactions
	real, allocatable						:: holding(:,:), coral(:,:), fish(:,:)  ! Layer names
	real, allocatable						:: kbact(:,:)							! Holds carrying capacity for bacteria
	type (bacteriavar) , allocatable		:: bacteria(:,:)						! Layer names
	integer, allocatable					:: seed(:)								! Random number holding array
	integer									:: clock, distance						! System time and radial distance of coral clusters
	character*20							:: filename								! Changes for what is being put into the file
	real									:: percentcover							! Percent of grid to have coral on it 'groundcover'
	real									:: fishlocal, fgrowfact
	integer									:: numnew = 0
	real									:: popconstant
	real									:: pi = acos(-1.0)
	integer									:: randall = 12
	real									:: avgpop, threshold
	integer									:: maxspec
	logical,allocatable						:: check(:,:)

	
end module


module functions

contains

real function fishdelta(input,pop)

use globalvars

implicit none
	real		:: input, pop
	
	fishdelta = -fgrowfact*(input - pop)
	
end function fishdelta

real function bacgrowth(totalpop,specpop,carry)

use globalvars

implicit none
	real		:: totalpop, specpop, carry
	real		:: rate	
	
rate = 0.75
bacgrowth = 0.0

	bacgrowth = rate*(1.0 - (real(totalpop)/real(carry)))*real(specpop)
	
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
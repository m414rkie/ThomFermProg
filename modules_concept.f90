module globalvars

type bacteria
	sequence
	integer	:: totalpop
	integer :: numspecies
end type bacteria

implicit none
	integer					:: grid									! Array size
	real					:: norm, nearsum, test					! Variables for interactions
	real, allocatable		:: holding(:,:), coral(:,:), fish(:,:)  ! Layer names
	real (bacteria) , allocatable		:: bacteria(:,:)						! Layer names
	integer, allocatable	:: seed(:)								! Random number holding array
	integer					:: clock, distance						! System time and radial distance of coral clusters
	character*20			:: filename								! Changes for what is being put into the file
	real					:: percentcover							! Percent of grid to have coral on it 'groundcover'
	real					:: fishconst, fishlocal, fgrowfact
	integer					:: numnew = 0
	
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

use global vars

implicit none
	real		:: totalpop, specpop, carry
	real		:: rate	
	
	bacgrowth = rate*(1.0 - (totalpop/carry)*specpop
	
end fnction bacgrowth

end module
module globalvars

implicit none
	integer					:: grid									! Array size
	real					:: norm, nearsum, test					! Variables for interactions
	real, allocatable		:: holding(:,:), coral(:,:), fish(:,:)  ! Layer names
	integer, allocatable	:: seed(:)								! Random number holding array
	integer					:: clock, distance						! System time and radial distance of coral clusters
	character*20			:: filename								! Changes for what is being put into the file
	real					:: percentcover							! Percent of grid to have coral on it 'groundcover'

	
end module

module globalvars

! Global variables and constants.

	real		(kind=8)						:: alpha, beta, gamma, sigma, zeta, xi
	real		(kind=8)     					:: alphal, betal, gammal, sigmal, bl
	real		(kind=8)					    :: alphau, betau, gammau, sigmau, bu
	integer					      		        :: n, checkarray, m, exitcondition
	real		(kind=8)					    :: dv, mass, matchoice, rhodv
	real		(kind=8)						:: pfin, rho0
	real		(kind=8)						:: momfin, momin, momentumcurrent
	real		(kind=8)						:: temp, rho, rhobar
	real		(kind=8)						:: pot, p0
	real		(kind=8)						:: kinout
	real		(kind=8)						:: pi = acos(-1.0)
	real  		(kind=8), allocatable			:: epsrho(:,:), pressure(:,:), enerray(:,:)
	character(1)								:: densenum, densenumparse, mattype, mattypeparse
	character(50)								:: pressdense, pressrho, rhorange
	real		(kind=8)						:: singlenergy, neutchem

end module
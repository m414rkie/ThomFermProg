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

end module


module functions

contains

real (kind=8) function potunmixed(x,y)

use globalvars

implicit none
	real (kind=8)	:: x, y, a, b
	real (kind=8)	:: blo, glo, sl
	
	
	blo(a,b) = betal*((abs(b-a)/momentumcurrent)**2)
	glo(a,b) = gammal*(momfin/abs(b-a))
	sl		 = sigmal*(((2.0*rhobar)/(rho0))**(2.0/3.0))


potunmixed = (alphal - blo(x,y)  - sl + glo(x,y))*(x**2)*(y**2)

end function potunmixed


real (kind=8) function potmixed(x,y)

use globalvars

implicit none
	real (kind=8)	:: x, y, a, b
	real (kind=8) 	:: bup, gup, su


	bup(a,b) = betau*(abs(b-a)/momentumcurrent)**2
	gup(a,b) = gammau*(momfin/abs(b-a))
	su		 = sigmau*(((2.0*rhobar)/(rho0))**(2.0/3.0))


potmixed = (alphau - bup(x,y)  - su + gup(x,y))*(x**2)*(y**2)

end function potmixed

real (kind=8) function kinet(in)

use globalvars

implicit none
	real (kind=8)	:: in
	
kinet = (in**2)/(2.0*mass)

end function kinet

end module










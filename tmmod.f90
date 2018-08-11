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
	real		(kind=8)						:: exclusion
	real (kind=8)								:: massmu = 105.66/197.329
	
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

	
	if (abs(y-x) .lt. 0.415) then
		potunmixed = 0.0
		goto 21
	end if

potunmixed = (alphal - blo(x,y)  - sl + glo(x,y))*(x**2)*(y**2)

21 end function potunmixed

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

real (kind=8) function potmixed(x,y)

use globalvars

implicit none
	real (kind=8)	:: x, y, a, b
	real (kind=8) 	:: bup, gup, su


	bup(a,b) = betau*(abs(b-a)/momentumcurrent)**2
	gup(a,b) = gammau*(momfin/abs(b-a))
	su		 = sigmau*(((2.0*rhobar)/(rho0))**(2.0/3.0))

	if (abs(y-x) .lt. 0.415) then
		potmixed = 0.0
		goto 22
	end if

potmixed = (alphau - bup(x,y)  - su + gup(x,y))*(x**2)*(y**2)

22 end function potmixed

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

real (kind=8) function kinet(in)

use globalvars

implicit none
	real (kind=8)	:: in
	
kinet = (in**2)/(2.0*mass)

end function kinet

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

real (kind=8) function rhotomom(rhoin)

use globalvars

implicit none
	real (kind=8)		::rhoin
	
rhotomom = (3.0*(pi**2)*rhoin)**(1.0/3.0)

end function

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

real(kind=8) function momtorho(mominf)

use globalvars

implicit none
	real(kind=8)		:: mominf
	
momtorho = (mominf**3)/(3.0*(pi**2))

end function

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

real(kind=8) function electronen(mominf)

implicit none
	real(kind=8)		:: masse = 0.511/197.329
	real(kind=8)		:: mominf
	
electronen = ((mominf**2) + (masse)**2)**(1.0/2.0)

end function

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

real(kind=8) function muonmom(chemin)

use globalvars

implicit none
	real (kind=8)		:: chemin

	
muonmom = sqrt((chemin**2) - (massmu**2))

end function

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

real(kind=8) function muonferm(chemenin)

use globalvars

implicit none
	real (kind=8)		:: chemenin
	
muonferm = sqrt((massmu**2)+(chemenin**2))

end function

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

real (kind=8) function potforchem(x,y)

use globalvars

implicit none
	real (kind=8)	:: x, y, a, b
	real (kind=8)	:: blo, glo, sl
	

	blo(a,b) = betal*((abs(b-a)/momentumcurrent)**2)
	glo(a,b) = gammal*(momfin/abs(b-a))
	sl		 = sigmal*(((2.0*rhobar)/(rho0))**(2.0/3.0))

	
	if (abs(y-x) .lt. 0.415) then
		potforchem = 0.0
		goto 21
	end if

potforchem = (alphal - blo(x,y)  - sl + glo(x,y))

21 end function potforchem

end module









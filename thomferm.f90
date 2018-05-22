PROGRAM msenergy
!
! This program is an attempt to replicate the numerical analysis of the Thomas-Fermi nucleon interarction
! energies as performed by Myers and Swiatecki.
! This program version assumes pre neutron matter.
!
! History
!	1.00a - Initial construction of program began
!	1.00b - Kinetic and Potential implemented
!	1.00c - Iteration of temp, rho implemented
!	1.00d - Implementation of pressure, epsilon
!	1.00e - User interface implemented
!
! Version   Programmer         Date       Description
! -------   ----------         --------   --------------------
! 1.00e     Jon Parsons        10-29-17		Snr. Thesis
!
! IN args/commons              Units      Description
! -----------------            -----      ----------------------
! input variables              units      variable purpose
!
! OUT args/commons             Units      Description
! ----------------             -----      ----------------------------
! rho 						  unitless	  density of the system as a ratio with standard density
! energy					   MeV		  average energy per nucleon
!
! Special requirements:
! 
!
!
! ------------------ Variable declarations ------------

! ------------------ Constant declarations ------------

! ------------------ Input/Output files ---------------

! Energy per Nucleon - Stores: rho, energy

! ------------------ Code -----------------------------

use globalvars

implicit none
	integer			:: i
	real (kind=8)	:: energy, momdv, rholow
	real (kind=8)	:: epsil, rhofinal

! Constant declarations
momfin = 1.36
rho0 = 0.161
mass = 939.57/197.329
temp = 37.679/197.329

	
alpha = 3.60928
beta  = 0.37597
gamma = 0.21329
sigma = 1.33677
xi    = 0.44003
zeta  = 0.59778

! l suffix for unmixed matter, u suffix for mixed matter
! For current version only unmixed matter is used (neutron matter)
alphau = 0.5*(1.0 + xi)*alpha
betau  = 0.5*(1.0+zeta)*beta
gammau = 0.5*(1.0+zeta)*gamma
sigmau = 0.5*(1.0+zeta)*sigma

alphal = 0.5*(1.0 - xi)*alpha
betal  = 0.5*(1.0-zeta)*beta
gammal = 0.5*(1.0-zeta)*gamma
sigmal = 0.5*(1.0-zeta)*sigma

! Step variables, n for integration, rhodv chosen for best agreement with incompressibility (derivatives highly
! reliant on step size)
n = 800
rhodv = 0.00833333333333333
! Step variables, n for integration, m for graph output
n = 800
m = 93

call input

if (densenumparse .eq. "S") then
	
	momentumcurrent = (3.0*(pi**2)*rho)**(1.0/3.0)
	rhobar = rho
	dv = momentumcurrent/float(n)
	
	call potential
	call kinetic
	
	energy = (1.0/rho)*(2.0/((2.0*pi)**3))*(kinout + 0.5*pot)
	epsil = 197.329*(energy+mass)*rho
	
	! Screen output
	write(*,*) "Density: ", rho
	write(*,*) "Energy per Nucleon: ", energy*197.329
	write(*,*) "Epsilon: ", epsil

else 



! Output loop

allocate(enerray(4,m))
do i = 1, m, 1

! Iteration of variables

rho = rholow + rhodv*float(i)
rhobar = rho

momentumcurrent = (3.0*(pi**2)*rho)**(1.0/3.0)
dv = momentumcurrent/float(n)

rho = rhodv*float(i)
rhobar = rho
momentumcurrent = (3.0*(pi**2)*rho)**(1.0/3.0)
dv = momentumcurrent/float(n)

call potential
call kinetic

! Output finalization
energy = (1.0/rho)*(2.0/((2.0*pi)**3))*(kinout + 0.5*pot)
epsil = (energy+mass)*rho*197.329
enerray(1,i) = rho
enerray(2,i) = epsil
enerray(3,i) = energy
enerray(4,i) = energy*rho

! Write statements

write(13,*) rho/rho0, energy*197.329


end do

! File Closing statements
close(unit=13)

call press(enerray)

call incompressibility(enerray)

end if

end program
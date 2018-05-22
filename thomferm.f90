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


! User interface 
! Inputs
50 write(*,*) "Symmetric matter or Neutron matter? (S/N)"
read(*,*) mattype

call chartoup(mattype,mattypeparse)

! Parse first input
if ((mattypeparse .ne. "S") .and. (mattypeparse .ne. "N")) then
	write(*,*) "Please enter 'S' for symmetric of 'N' for neutron matter"
	goto 50
end if

if (mattypeparse .eq. "S") then
	
	matchoice = 1.0
	
	else
	
	matchoice = 0.0

end if

51 write(*,*) "Would you like a single density or a range? (S/R)"
read(*,*) densenum

call chartoup(densenum,densenumparse)

! Parse second input
if ((densenumparse .ne. "S") .and. (densenumparse .ne. "R")) then
	write(*,*) "Please enter 'S' for a single density or 'R' for a range."
	goto 51
end if

! Further inputs 

if (densenumparse .eq. "S") then

! For a single density
	write(*,*) "Please enter the density desired."
	read(*,*) rho
	
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

! For a range of densities
	write(*,*) "Please enter the lower bound of density desired:"
	read(*,*) rholow
	write(*,*) "Please enter the upper bound:"
	read(*,*) rhofinal
	
	m = ceiling((rhofinal - rholow)/rhodv)

	
	write(rhorange,'(2f5.2)') rholow, rhofinal

	rhodv = (rhofinal - rholow)/m

	! Filename logic, standard (-std) or neutron (-neut) matter names
	if (mattypeparse .eq. "S") then
	
		open(unit=13,file=trim(rhorange)//"energypernucleonstd.dat",position="append",status="replace")
		pressrho   = trim(rhorange)//"Pressurebyrhostd.dat"
		pressdense = trim(rhorange)//"Pressurebyepsilonstd.dat"
		
	else
	
		open(unit=13,file=trim(rhorange)//"energypernucleonneut.dat",position="append",status="replace")
		pressrho	 = trim(rhorange)//"Pressurebyrhoneut.dat"
		pressdense   = trim(rhorange)//"Pressurebyepsilonneut.dat"

		open(unit=13,file="energypernucleonstd.dat",position="append",status="replace")
		pressrho   = "Pressurebyrhostd.dat"
		pressdense = "Pressurebyepsilonstd.dat"
		comprho	   = "Incompressbyrhostd.dat" 
		
	else
	
		open(unit=13,file="energypernucleonneut.dat",position="append",status="replace")
		pressrho	 = "Pressurebyrhoneut.dat"
		pressdense   = "Pressurebyepsilonneut.dat"
		comprho		 = "Incompressbyrhoneut.dat"
		
	end if

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
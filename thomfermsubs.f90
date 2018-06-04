subroutine boolequad(arr,booleint)

! Integration subroutine, uses boolean quadrature method (5-pt)

use globalvars

implicit none
	real (kind=8), dimension(n)		  :: arr		! Dummy name for array
	real (kind=8), intent(out)		  :: booleint	! Final value for this subroutine
	integer							  :: i			! Looping integer

! Initialization
booleint = 0.0

! Loop for composite Boole's rule. Weights included.
do i = 4, n, 4
	
	booleint = booleint + dv*(2.0/45.0)*(7.0*arr(i)+32.0*arr(i-1)+12.0*arr(i-2)+32.0*arr(i-3)+7.0*arr(i-4))

end do

end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine kinetic

! Calcalates the kinetic energy, calls subroutine boolequad

use globalvars

implicit none
	integer								:: i			! loop integer
	real (kind=8)						:: kin, x, p1	! equation and variable variables
	real (kind=8), allocatable			:: kinarr(:)	! array for holding values for use in boolequad

! Kinetic energy formula, original int(x^2)d3.
! as implemented the d3 is converted to 4pix^4 resulting in formula kin(x)
kin(x) = (x**4)

allocate(kinarr(n), stat=checkarray)
	if(checkarray .ne. 0) stop "kinarr failed"
	
do i = 1, n, 1
	p1 = float(i)*momentumcurrent/n
	kinarr(i) = kin(p1)
end do

call boolequad(kinarr,kinout)

kinout = ((4.0*pi)/(2.0*mass))*kinout
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine potential

! Calculates the potential
! Formulas derived from the thomas-fermi model as interpreted by Myers and Swiatecki
! Calls subroutine boolequad
use globalvars
use functions

implicit none
	integer 					:: i, j							! Looping variables
	real (kind=8), allocatable 	:: upperarr1(:), upperarr2(:)	! Arrays for first and second integration
	real (kind=8), allocatable 	:: lowerarr1(:), lowerarr2(:)	! Arrays for first and second integration
	real (kind=8), allocatable 	:: botharr(:)					! Summed array for output
	real (kind=8)				:: mom1, mom2					! Iterated variables
	real (kind=8)				:: uquad, lquad					! Temp. holding variables for integration
	real (kind=8)				:: dv1, dv2						! Step sizes for the iterated variables
	real (kind=8)				:: pref, prefu					! Prefixes for the explicit portion of the pot. full equation


! Allocation statements, includes mixed and unmixed.
allocate(upperarr1(n))
allocate(upperarr2(n))
allocate(lowerarr1(n))
allocate(lowerarr2(n))
allocate(botharr(n))

! Initialization statements
upperarr1 = 0.0
upperarr2 = 0.0
lowerarr1 = 0.0
lowerarr2 = 0.0
botharr   = 0.0
dv1 = momentumcurrent/float(n)
dv2 = momentumcurrent/float(n)

! Loops that iterate p1 and p2
! I loop iterates p1, J loop iterates p2
do i = 1, n, 1
	
	mom1 = dv1*float(i)

	do j = 1, n, 1
	
		mom2 = dv2*float(j)
		
		if (abs(mom1-mom2) .gt. 0.415) then
		
			
			upperarr1(j) = (matchoice)*potmixed(mom1,mom2)
			lowerarr1(j) = potunmixed(mom1,mom2)
		
		end if
		
	end do

	call boolequad(upperarr1,uquad)
	upperarr2(i) = uquad
	call boolequad(lowerarr1,lquad)
	lowerarr2(i) = lquad
	
	if(exitcondition .eq. 1) then
		lquad = singlenergy
		exit
	else 
		goto 102
	end if

102 end do

! Combines the mixed and pure matter results for use in integration

botharr = lowerarr2 + (matchoice)*upperarr2

call boolequad(botharr,pot)

! Finalization of results
pot = -temp*(2.0/rho0)*((4.0*pi)**2)*(2.0/((2.0*pi)**3))*pot


end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine chartoup(stringin,stringout)

! converts text input to upper case

implicit none
	character(*)					:: stringin
	character(len(stringin))		:: stringout
	integer							:: i, j

do i = 1, len(stringin), 1
	j = iachar(stringin(i:i))
		if(j .ge. iachar("a") .and. j .le. iachar("z")) then
			stringout(i:i) = achar(iachar(stringin(i:i))-32)
		else
			stringout(i:i) = stringin(i:i)
		end if
end do

end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine press(arrin)

! Calculates the pressure as a derivative of energy per nucleon. 

use globalvars

implicit none
	real (kind=8),dimension(4,m), intent(in)	:: arrin
	real (kind=8)								:: derivative
	integer	(kind=8)							:: i, step

open(unit=11,file=trim(pressrho),position="append",status="replace")
open(unit=14,file=trim(pressdense),position="append",status="replace")

step = 4
		
do i = 1, m-step, 1

	derivative = (arrin(3,i+step) - arrin(3,i))/(float(step)*rhodv)
	
	write(11,*) arrin(1,i)/rho0, (arrin(1,i)**2)*derivative*197.329
	write(14,*) arrin(2,i), (arrin(1,i)**2)*derivative*197.329

end do
	
close(11)
close(14)

end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine incompressibility(arrin)

! Calculates the compressibility of the system as a second derivative of the energy.

use globalvars

implicit none
	real (kind=8),dimension(4,m), intent(in)				:: arrin
	real (kind=8)											:: compress
	integer (kind=8), dimension(4)							:: valarray
	integer (kind=8)										:: minim

valarray = minloc(arrin,dim=2)
minim = valarray(3)

compress = (arrin(4,minim+1)-2.0*arrin(4,minim)+arrin(4,minim-1))/((arrin(1,minim+1)-arrin(1,minim))**2)

write(*,*)"Density: ", arrin(1,minim)/rho0, "Compressibility :", 9.0*(arrin(1,minim)**2)*compress*197.329

end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine input

! Parses input

use globalvars

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
	end if
	
end if

end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine chempot

use globalvars

implicit none
	real (kind=8)		:: pfp, pfe
	real (kind=8)		:: protchem, elecchem

!write(*,*) "Input Density to use in chemical potential test:"
!read(*,*) rho

pfp = (3.0*(pi**2)*rho/2.0)**(1.0/3.0)

pfe = pfp
exitcondition = 1
matchoice = 0.0
rhobar = rho

momentumcurrent = pfp
call potential
write(*,*) "A"
protchem = -temp*(2.0/rho0)*((4.0*pi)**2)*(2.0/((2.0*pi)**3))*singlenergy*197.329 + ((pfp**2)/(2.0*mass*197.329))

elecchem = ((pfe**2) + (0.511)**2)**(1.0/2.0)

neutchem = protchem + elecchem

write(*,*) "Proton potential:", protchem
write(*,*) "Electron potential:", elecchem
write(*,*) "Proton density:", rho/2.0
write(*,*) "Proton Fermi momentum:", pfp

call root

end subroutine

subroutine root 

use globalvars

implicit none
	real (kind=8)		:: step, pfn
	real (kind=8)		:: neutsingle, neutrho
	real (kind=8)		:: func, derive, derivepls1, deriveorig
	
matchoice = 0.0
step = 0.01
pfn = 0.0
101 momentumcurrent = pfn
rho = (pfn**3)/(3.0*(pi**2))
rhobar = rho
	
call potential	
neutsingle = -temp*(2.0/rho0)*((4.0*pi)**2)*(2.0/((2.0*pi)**3))*singlenergy*197.329 + ((pfn**2)/(2.0*mass*197.329))
	
if(abs(neutchem - neutsingle) .le. 0.001) then
	write(*,*) "Neutron potential:", neutchem
	write(*,*) "Neutron single energy:", neutsingle
	write(*,*) "Neutron density:", rho
	write(*,*) "Neutron Fermi momentum:", pfn
else
	deriveorig = neutsingle
	momentumcurrent = pfn + step 
	call potential
	derivepls1 = -temp*(2.0/rho0)*((4.0*pi)**2)*(2.0/((2.0*pi)**3))*singlenergy*197.329 + ((momentumcurrent**2)/(2.0*mass*197.329))
	derive = (derivepls1 - deriveorig)/step
	pfn = pfn + step !- deriveorig/derive
	goto 101
end if

end subroutine

	














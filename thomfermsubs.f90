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

implicit none
	real (kind=8)				:: au, gup, su, bup				! Fitted coefficients, mixed matter
	real (kind=8)				:: al, glo, sl, blo				! Fitted coefficients, unmixed matter
	real (kind=8)				:: x, y							! Eqn variables
	integer 					:: i, j							! Looping variables
	real (kind=8), allocatable 	:: upperarr1(:), upperarr2(:)	! Arrays for first and second integration
	real (kind=8), allocatable 	:: lowerarr1(:), lowerarr2(:)	! Arrays for first and second integration
	real (kind=8), allocatable 	:: botharr(:)					! Summed array for output
	real (kind=8)				:: mom1, mom2					! Iterated variables
	real (kind=8)				:: uquad, lquad					! Temp. holding variables for integration
	real (kind=8)				:: dv1, dv2						! Step sizes for the iterated variables
	real (kind=8)				:: pref, prefu					! Prefixes for the explicit portion of the pot. full equation

! Equations used, l for unmixed matter, u for mixed
bup(x,y) = betau*(abs(y-x)/momentumcurrent)**2
blo(x,y) = betal*((abs(y-x)/momentumcurrent)**2)
gup(x,y) = gammau*(momfin/abs(y-x))
glo(x,y) = gammal*(momfin/abs(y-x))
sl		 = sigmal*(((2.0*rhobar)/(rho0))**(2.0/3.0))
su		 = sigmau*(((2.0*rhobar)/(rho0))**(2.0/3.0))


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

			upperarr1(j) = (matchoice)*(alphau - bup(mom1,mom2)  - su + gup(mom1,mom2))*(mom1**2)*(mom2**2)
			lowerarr1(j) = (alphal - blo(mom1,mom2)  - sl + glo(mom1,mom2))*(mom1**2)*(mom2**2)
		
		end if
		
	end do
	
	call boolequad(upperarr1,uquad)
	upperarr2(i) = uquad
	call boolequad(lowerarr1,lquad)
	lowerarr2(i) = lquad
	
end do

! Combines the mixed and pure matter results for use in integration

botharr = lowerarr2 + (matchoice)*upperarr2

botharr = lowerarr2 + upperarr2

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

write(*,*) arrin(1,minim)/rho0, 9.0*(arrin(1,minim)**2)*compress*197.329

end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine input

! Parses input

use globalvars

























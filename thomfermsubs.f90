subroutine boolequad(arr,booleint,upto)

! Integration subroutine, uses boolean quadrature method (5-pt)

use globalvars

implicit none
	real (kind=8), dimension(n)		  :: arr		! Dummy name for array
	integer							  :: upto		! Array input size
	real (kind=8), intent(out)		  :: booleint	! Final value for this subroutine
	integer							  :: i			! Looping integer

! Initialization
booleint = 0.0

! Loop for composite Boole's rule. Weights included.
do i = 2, upto, 1
	
	booleint = booleint + 0.5*dv*(arr(i) + arr(i-1))
	
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

call boolequad(kinarr,kinout,n)

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
		
		if (abs(mom1-mom2) .gt. exclusion) then
		
			upperarr1(j) = (matchoice)*potmixed(mom1,mom2)
			lowerarr1(j) = potunmixed(mom1,mom2)
		
		end if
		
	end do

	call boolequad(upperarr1,uquad,n)
	upperarr2(i) = uquad
	call boolequad(lowerarr1,lquad,n)
	lowerarr2(i) = lquad

end do

! Combines the mixed and pure matter results for use in integration

botharr = lowerarr2 + (matchoice)*upperarr2

call boolequad(botharr,pot,n)

! Finalization of results
pot = -temp*(2.0/rho0)*((4.0*pi)**2)*(2.0/((2.0*pi)**3))*pot


end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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


subroutine chempot

use globalvars
use functions

implicit none
	real (kind=8)					:: pfp, pfe, pfn, pfx, pfy
	real (kind=8)					:: pfm, rhomu, muchem, charge
	real (kind=8)					:: protchem, elecchem, neutchem, elsechem
	integer							:: i, j, arrval, intcount, k
	real (kind=8)					:: protelectot, rhoprot
	real (kind=8) 					:: neutdense, totdense
	real (kind=8), allocatable		:: protchemarr(:)
	real (kind=8)					:: mucheck, elecrho

!write(*,*) "Input density due to protons:"
!read(*,*) rhoprot
!dv = 0.0001
protchem = 0.0
elecchem = 0.0
neutchem = 0.0

open(unit=13,file="reldensities.dat",status="replace",position="append")
write(13,*) "density	proton	electron	neutron		muon"

do j = 1, 200, 1

protchem = 0.0
elecchem = 0.0
neutchem = 0.0
rhomu = 0.0


pfe = float(j)*0.01

elecrho = momtorho(pfe)
elecchem = electronen(pfe)
muchem = elecchem
pfm = muonmom(muchem)
rhomu = momtorho(pfm)

if (pfm .eq. 0.0) then
	write(*,*) "MU NOT!!!"
	muchem = 0.0
	pfm = 0.0
	rhomu = 0.0
end if


rhoprot = elecrho + rhomu
rho = rhoprot
rhobar = rho

pfp = rhotomom(rhoprot)

momentumcurrent = pfp

arrval = 1000!ceiling(pfp/dv)
allocate(protchemarr(arrval))

protchemarr = 0.0

do i = 1, 1000, 1!arrval, 1
 pfx = float(i)*(pfp/1000.0)!dv*float(i)
 
	protchemarr(i) = potunmixed(pfp,pfx)

end do

call boolequad(protchemarr,protchem,arrval)

protchem = -temp*(2.0/rho0)*((4.0*pi))*(2.0/((2.0*pi)**3))*protchem + kinet(pfp)

protelectot = protchem + elecchem

neutchem = protelectot

pfn = 0.5*pfp 
call root(pfn,neutchem)

neutdense =  momtorho(pfn)

charge = rhoprot - elecrho - rhomu

write(*,*) "Charge: ", charge

totdense = (neutdense + rhoprot + rhomu + elecrho)

write(*,*) "Density:", totdense

write(13,*) totdense/rho0, rhoprot/totdense, elecrho/totdense, neutdense/totdense, rhomu/totdense

deallocate(protchemarr)

end do

end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine root(finchem,pcheck) 

use globalvars
use functions

implicit none
	real (kind=8)					:: finchem, pcheck
	integer							:: arrsize, i, nit
	real (kind=8)					:: pfnx, delta, dvp
	real (kind=8)					:: neutsingle, neutrho
	real (kind=8)					:: diff
	real (kind=8), allocatable		:: neutchemarr(:)

nit = 0
!dvp = 0.0001
delta = 0.0001
neutsingle = 0.0 
diff = 1.0

do while ((abs(diff) .gt. 0.001).or.(pcheck*neutsingle .lt. 0.0))

finchem = finchem + delta

nit = nit + 1

if (finchem .gt. 5.0) then
	write(*,*) "No root"
	exit
end if

!arrsize = ceiling(finchem/dvp)
arrsize = 1000
momentumcurrent = finchem

allocate(neutchemarr(arrsize))

neutchemarr = 0.0
neutsingle = 0.0

rho = momtorho(finchem)
rhobar = rho
momentumcurrent = finchem

do i = 1, 1000, 1!arrsize, 1
		
	pfnx = float(i)*(finchem/1000.0)!dvp*float(i)
	
	neutchemarr(i) = potunmixed(finchem,pfnx)	

end do

call boolequad(neutchemarr,neutsingle,arrsize)

neutsingle =  -temp*(2.0/rho0)*((4.0*pi))*(2.0/((2.0*pi)**3))*neutsingle + kinet(finchem)

diff = abs(pcheck)-abs(neutsingle)

deallocate(neutchemarr)
	
end do

write(*,*) "Input: ", pcheck, "Found: ", neutsingle, "delta: ", diff
write(*,*) nit

end subroutine

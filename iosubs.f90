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
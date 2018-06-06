PROGRAM concept
!
! This is a proof of concept for use in a program that will model
! interaction of coral reefs
!
! History
! 1.00a - Initial program construction
! 1.00b - Coral and algae interaction
! 1.00c - Fish layer added
! 1.00d - Coral/Algae layer population routines
!
! Version   Programmer         Date       Description
! -------   ----------         --------   --------------------
! 1.00d     Jon Parsons        3-19-18	  Proof of concept
!
! IN args/commons              Units      Description
! -----------------            -----      ----------------------
! input variables              units      variable purpose
!
! OUT args/commons             Units      Description
! ----------------             -----      ----------------------------
! output variables             units      variable purpose
! 
!
! Special requirements:
! * Module file -  modules_concept.f90
! * Subroutines -  subs_concept.f90
!				   concept_popsubs.f90
!
! ------------------ Code -----------------------------

use globalvars

implicit none
	integer				:: i, j, t, k, n ! Looping integers; n is random seed holder
	integer				:: numtime, clusnum	! Number of timesteps and clusters of coral

50 format ("coraltime",1i2,".dat")		


! User input 
write(*,*) "Enter the dimension of the grid (square):"
read(*,*) grid
write(*,*) "Dimension: ", grid
write(*,*) "Enter the number of time steps :"
read(*,*) numtime
write(*,*) "Enter percentage of bed with coral:"
read(*,*) percentcover
write(*,*) "Number of coral clusters?"		! Currently 
read(*,*) clusnum
write(*,*) "Please input distance for the tightly clustered coral clusters:"
read(*,*) distance

! Allocation statements, error checking coming soon.
allocate(coral(grid,grid))
allocate(holding(grid,grid))
allocate(fish(grid,grid))

! Initializing grids
coral = 0.0
holding = 0.0
fgrowfact = 0.25
fishconst = 5.0

! Generates a seed for use in populating layers
call random_seed(size=n)
allocate(seed(n))

call system_clock(count=clock)

seed = clock

! Populates the coral/algae layer
call hppop(coral)

do k=1, clusnum, 1
	call tightcluster(coral)
end do

! Increases overal population of coral as each gridpoint will be between zero and one beforehand
coral = 2.0*coral

holding = coral

write(*,*) "0.0 represents pure algae; greater than zero represents coral, higher number represents more coral"
write(*,*) "Files are written as (x,y,z) where z is the population/biomass"

! Initial disposition of coral/algae layer. 
filename = "coralini.dat"
call printtofile(coral)
call fishdist(fish)

filename = "fishini.dat"
call printtofile(fish)

! Outer loops iterates time, i and j iterate x and y respectively
do t = 1, numtime, 1

	write(*,*) "timestep", t
		
		do i = 1, grid, 1
	
			do j = 1, grid, 1
		
				call neighborsum(i,j,holding,nearsum)
				call growth(i,j,coral,coral)
				call decay(i,j,coral)
				call newcoral
	
			end do
	
		end do

		if (mod(t,5) .eq. 0) then
			write(filename,50) t
			call printtofile(coral)
		end if

 
 		holding = coral

end do

write(*,*) numnew

! Prints the final coral/algae layer after the number of timesteps is reached.
filename = "coralfin.dat"
call printtofile(coral)

filename = "fishfin.dat"
call printtofile(fish)

 
end program
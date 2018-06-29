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
! 1.00d - Bacteria layer added
!
! Version   Programmer         Date       Description
! -------   ----------         --------   --------------------
! 1.00e     Jon Parsons        6-23-18	  Proof of concept
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
	integer					:: i, j, t ! Looping integers; n is random seed holder
	integer					:: numtime ! Number of timesteps and clusters of coral

! Format statements
50 format ("coraltime",1i2,".dat")		


! User input 
write(*,*) "Enter the dimension of the grid (square):"
read(*,*) grid
write(*,*) "Dimension: ", grid
write(*,*) "Enter the number of time steps :"
read(*,*) numtime
write(*,*) "Enter percentage of bed with coral:"
read(*,*) percentcover
write(*,*) "Number of coral clusters?"
read(*,*) clusnum
write(*,*) "Please input distance for the tightly clustered coral clusters:"
read(*,*) distance
write(*,*) "New coral threshold?"
read(*,*) threshold

! Allocation statements, error checking coming soon.
allocate(coral(grid,grid))
allocate(holding(grid,grid))
allocate(fish(grid,grid))
allocate(check(grid,grid))
allocate(bacteria(2*grid,2*grid))
allocate(kbact(2*grid,2*grid))
allocate(seed(randall))

! Initializing grids
coral = 0.0
holding = 0.0
fgrowfact = 0.25
bacteria%totalpop = 0
bacteria%numspecies = 0
sharkmod = 0.0
hunger = 0.0	

! Populates the coral/algae layer
call hppop(coral)
call tightcluster(coral)

! Increases overal population of coral as each gridpoint will be between zero and one beforehand
coral = 2.0*coral

holding = coral

write(*,*) "0.0 represents pure algae; greater than zero represents coral, higher number represents more coral"
write(*,*) "Files are written as (x,y,z) where z is the population/biomass"

! Initial disposition of coral/algae layer. 
filename = "coralini.dat"
call printtofile(coral)

call fishdist(fish)

! Initial disposition of fish layer
filename = "fishini.dat"
call printtofile(fish)

! Populating initital bacteria layer
call bacteriapop
! Outer loops iterates time, i and j iterate x and y respectively
do t = 1, numtime, 1

	write(*,*) "timestep", t
	
	call shark
		
		do i = 1, grid, 1
	
			do j = 1, grid, 1
		
				!call neighborsum(i,j,holding,nearsum)
				call growth(i,j,coral,coral)
				call decay(i,j,coral)				
	
			end do
	
		end do
		
		call newcoral
		call kgrid
		call bactgrow
		call diffuse
		call mixing	
		if (mod(t,5) .eq. 0) then
			write(filename,50) t
			call printtofile(coral)
		end if

 
 		holding = coral

end do

write(*,*) "Total number of new coral growths:", numnew

! Print statements for final layer after the number of timesteps is reached.
filename = "coralfin.dat"
call printtofile(coral)

filename = "fishfin.dat"
call printtofile(fish)

open(unit=14,file="bactlayerfin.dat",status="replace",position="append")
	
	do i = 1, 2*grid, 1
		do j = 1, 2*grid, 1
			write(14,*) i, j, bacteria(i,j)%numspecies, bacteria(i,j)%totalpop
		end do
	end do

close(14)

deallocate(coral)
deallocate(holding)
deallocate(fish)
deallocate(check)
deallocate(bacteria)
deallocate(kbact)
deallocate(seed)

 
end program
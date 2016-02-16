!!REORDER(4): solnData

subroutine Simulation_initBlock(blockID)
  use Simulation_data, ONLY : rhoOut,rhoIn,P,rcloud, &
                              sim_xctr,sim_yctr,sim_zctr,&
                              sim_gamma,sim_gascon 

  use Grid_interface, ONLY : Grid_getBlkIndexLimits, Grid_getBlkPtr, &
    Grid_getDeltas, Grid_releaseBlkPtr, Grid_getBlkBoundBox


  implicit none

#include "constants.h"
#include "Flash.h"

  ! compute the maximum length of a vector in each coordinate direction 
  ! (including guardcells)
  
  integer, intent(in) :: blockID
  
  real, pointer, dimension(:,:,:,:) :: solnData
  
  real,dimension(MDIM) :: size, coord
  
  integer         i, j, k, n, imax, jmax, kmax, jlo
  integer         Nint, ii, jj, kk
  real            delx, xx, dely, yy, delz, zz, velocity, distinv
  real            xdist, ydist, zdist, dist
  real            xxmin, xxmax, yymin, yymax, zzmin, zzmax,ek
  real            Nintinv,Nintinv1
 ! real 		  sim_gascon,sim_gamma
  integer, dimension(MDIM) :: guard
  integer, dimension(LOW:HIGH,MDIM) :: blkLimits, blkLimitsGC
  real, dimension(LOW:HIGH,MDIM) :: bndBox
  real, dimension(MDIM) :: delta

!==========================================================================
!               Initialize scalar quantities we will need.

  call Grid_getBlkPtr(blockID,solnData,CENTER)
  call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)
  call Grid_getBlkBoundBox(blockID,bndBox)
  call Grid_getDeltas(blockID,delta)
!  sim_gascon=8.2544E7
!  sim_gamma=1.0001
!  print *, "begin Simulation_initBlock"
  imax = blkLimitsGC(HIGH,IAXIS)-blkLimitsGC(LOW,IAXIS)+1
  jmax = blkLimitsGC(HIGH,JAXIS)-blkLimitsGC(LOW,JAXIS)+1
  kmax = blkLimitsGC(HIGH,KAXIS)-blkLimitsGC(LOW,KAXIS)+1
  guard = blkLimits(LOW,:)-blkLimitsGC(LOW,:)

  ! Coordinates of the edges of the block

  xxmax = bndBox(HIGH,IAXIS)
  xxmin = bndBox(LOW,IAXIS)
  yymax = bndBox(HIGH,JAXIS)
  yymin = bndBox(LOW,JAXIS)
  zzmax = bndBox(HIGH,KAXIS)
  zzmin = bndBox(LOW,KAXIS)

  ! Cell size
  
  delx = delta(IAXIS)
  dely = delta(JAXIS)
  delz = delta(KAXIS)
!-------------------------------------------------------------------------------

!               Loop over cells in the block.  For each, compute the physical
!               position of its left and right edge and its center as well as
!               its physical width.  Then decide whether it is inside the
!               initial radius or outside and initialize the hydro variables
!               appropriately.

  Nint    = 7
  Nintinv = 1./float(Nint)
  Nintinv1= 1./(float(Nint)-1.)
  
  do k = 1, kmax
    !zdist = k - sim_zctr 
    do j = 1, jmax
       ! ydist = j - sim_yctr 
        do i = 1, imax
	  !  xdist = i - sim_xctr
           do kk = 0, (Nint-1)*K3D
              zz    = zzmin + delz*(real(k-guard(KAXIS)-1)+kk*Nintinv1)
              zdist = (zz - sim_zctr) * K3D
              do jj = 0, (Nint-1)*K3D
                 yy    = yymin + dely*(real(j-guard(JAXIS)-1)+jj*Nintinv1)
                 ydist = (yy - sim_yctr) * K3D
                 do ii = 0, Nint-1
                    xx    = xxmin + delx*(real(i-guard(IAXIS)-1)+ii*Nintinv1)
                    xdist = (xx - sim_xctr)*K3D
                    dist    = sqrt( xdist**2 + ydist**2 + zdist**2 )
		    !print *,"(i,j,k),dist:",i,j,k,dist
		    if (dist<rcloud) then 
			solnData(DENS_VAR,i,j,k) = rhoIn 
		    else
			solnData(DENS_VAR,i,j,k) = rhoOut              
	            endif
                    
                ! enddo
              !enddo
           !enddo
           
           solnData(PRES_VAR,i,j,k) = P 
           solnData(TEMP_VAR,i,j,k) = solnData(PRES_VAR,i,j,k) /(solnData(DENS_VAR,i,j,k)*sim_gascon)
           solnData(VELX_VAR,i,j,k) = 0.0
           solnData(VELY_VAR,i,j,k) = 0.0
           solnData(VELZ_VAR,i,j,k) = 0.0 
           
        enddo
     enddo
  enddo
  
  !-------------------------------------------------------------------------------
  
  !               Initialize the nuclear abundances.  These are not of interest
  !               in this problem, so we set them to 1. everywhere.
  
  do n = 1, NSPECIES
     do k = 1, kmax
        do j = 1, jmax
           do i = 1, imax
              solnData(SPECIES_BEGIN+n-1,i,j,k) = 1.
           enddo
        enddo
     enddo
  enddo
  
  !               Compute the gas energy and set the gamma-values needed for
  !               the equation of state.
  
  do k = 1, kmax
     do j = 1, jmax
        do i = 1, imax
           
           solnData(GAME_VAR,i,j,k) = sim_gamma
           solnData(GAMC_VAR,i,j,k) = sim_gamma
           
           ek = 0.5 * (solnData(VELX_VAR,i,j,k)**2 + & 
                &                    solnData(VELY_VAR,i,j,k)**2 + & 
                &                    solnData(VELZ_VAR,i,j,k)**2)
           solnData(EINT_VAR,i,j,k) = solnData(PRES_VAR,i,j,k) / & 
                &                                    (solnData(GAME_VAR,i,j,k)-1.)
           solnData(EINT_VAR,i,j,k) = solnData(EINT_VAR,i,j,k) / & 
                &                                    solnData(DENS_VAR,i,j,k)
           solnData(EINT_VAR,i,j,k) = solnData(EINT_VAR,i,j,k)
           solnData(ENER_VAR,i,j,k) = solnData(EINT_VAR,i,j,k) + ek
           
        enddo
     enddo
  enddo
  
  !===============================================================================
  
  call Grid_releaseBlkPtr(blockID,solnData)

end subroutine Simulation_initBlock



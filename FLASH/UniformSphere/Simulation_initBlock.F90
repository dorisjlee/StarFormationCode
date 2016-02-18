!!REORDER(4): solnData
subroutine Simulation_initBlock(blockID,solnData)
#include "Flash.h"
#include "constants.h"
#include "Flash.h"
#include "Eos.h"

    
  use Simulation_data, ONLY : rhoOut,rhoIn,P,rcloud, &
                              sim_xctr,sim_yctr,sim_zctr,&
                              sim_gamma,sim_gascon
  !use Grid_interface, ONLY : Grid_getBlkIndexLimits, &
  !  Grid_getCellCoords, Grid_putPointData,
  use Grid_interface, ONLY : Grid_getBlkIndexLimits, Grid_getBlkPtr, &
    Grid_getDeltas, Grid_releaseBlkPtr, Grid_getBlkBoundBox,Grid_getCellCoords, Grid_putPointData
  use Eos_interface, ONLY : Eos, Eos_wrapped

  implicit none
  integer, intent(in) :: blockID
  real, pointer, dimension(:,:,:,:) :: solnData

  real,dimension(MDIM) :: size, coord

  integer         i, j, k, n, imax, jmax, kmax, jlo
  integer         Nint, ii, jj, kk
  real ::         xx, yy, zz, r,ek
  real ::         lPosn0, lPosn 
  real,allocatable, dimension(:) ::xCoord,yCoord,zCoord
 ! real                   sim_gascon,sim_gamma
  integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC
  real, dimension(LOW:HIGH,MDIM) :: bndBox
  real, dimension(MDIM) :: delta
  integer :: sizeX,sizeY,sizeZ
  integer, dimension(MDIM) :: axis
  real :: rhoZone, velxZone, velyZone, velzZone, presZone, & 
       eintZone, enerZone, ekinZone, gameZone, gamcZone

#ifdef SIMULATION_TWO_MATERIALS
  real, dimension(EOS_NUM) :: eosData
  real, dimension(NSPECIES) :: mfrac
#endif

#ifdef FLASH_3T
  real :: peleZone, eeleZone
  real :: pionZone, eionZone
  real :: pradZone, eradZone
#endif
  
  logical :: gcell = .true.

  
  ! get the integer index information for the current block
  call Grid_getBlkPtr(blockID,solnData,CENTER)
  call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)
  call Grid_getBlkBoundBox(blockID,bndBox)
  call Grid_getDeltas(blockID,delta)  
  
  sizeX = blkLimitsGC(HIGH,IAXIS)
  sizeY = blkLimitsGC(HIGH,JAXIS)
  sizeZ = blkLimitsGC(HIGH,KAXIS)
  allocate(xCoord(sizeX))
  allocate(yCoord(sizeY))
  allocate(zCoord(sizeZ))
  xCoord = 0.0
  yCoord = 0.0
  zCoord = 0.0

  if (NDIM == 3) call Grid_getCellCoords&
                      (KAXIS, blockId, CENTER,gcell, zCoord, sizeZ)
  if (NDIM >= 2) call Grid_getCellCoords&
                      (JAXIS, blockId, CENTER,gcell, yCoord, sizeY)
  call Grid_getCellCoords(IAXIS, blockId, CENTER, gcell, xCoord, sizeX)

!------------------------------------------------------------------------------

! Loop over cells in the block.  For each, compute the physical position of 
! its left and right edge and its center as well as its physical width.  
! Then decide which side of the initial discontinuity it is on and initialize 
! the hydro variables appropriately.


  do k = blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS)
     zz = zCoord(k) -sim_zctr
     do j = blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)
        yy = yCoord(j)-sim_yctr
        do i = blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)
           xx  = xCoord(i)-sim_xctr
!           print *,"(i,j,k):",i,j,k
!           print *,"(xx,yy,zz):",xx,yy,zz
!           print *,"zCoord,yCoord,xCenter:",zCoord(k),yCoord(j),xCenter(i)
           r = sqrt(xx**2+yy**2+zz**2)
           if (r  <= rcloud) then
!              solnData(DENS_VAR,i,j,k)  = rhoIn
               rhoZone = rhoIn
           else 
!              solnData(DENS_VAR,i,j,k)  = rhoOut 
               rhoZone = rhoOut 
           endif
           presZone = P
           velxZone = 0.0
           velyZone = 0.0
           velzZone = 0.0
!	   solnData(PRES_VAR,i,j,k) = 59.6525 !P
          ! solnData(TEMP_VAR,i,j,k) = solnData(PRES_VAR,i,j,k) /(solnData(DENS_VAR,i,j,k)*sim_gascon)
 !          solnData(VELX_VAR,i,j,k) = 0.0
 !          solnData(VELY_VAR,i,j,k) = 0.0
 !          solnData(VELZ_VAR,i,j,k) = 0.0
       
           axis(IAXIS) = i
           axis(JAXIS) = j
           axis(KAXIS) = k

           ekinZone = 0.5 * (velxZone**2 + &
                velyZone**2 + &
                velzZone**2)

#ifdef SIMULATION_TWO_MATERIALS
           eosData(EOS_DENS) = rhoZone
           eosData(EOS_PRES) = presZone
           eosData(EOS_TEMP) = 1.0e8
           call Eos(MODE_DENS_PRES, 1, eosData, mfrac)
           eintZone = eosData(EOS_EINT)
           gameZone = 1.0+eosData(EOS_PRES)/eosData(EOS_DENS)/eosData(EOS_EINT)
           gamcZone = eosData(EOS_GAMC)
#else
           eintZone = presZone / (sim_gamma-1.)
           eintZone = eintZone / rhoZone
           gameZone = sim_gamma
           gamcZone = sim_gamma
#endif
           enerZone = eintZone + ekinZone

        call Grid_putPointData(blockId, CENTER, DENS_VAR, EXTERIOR, axis, rhoZone)
           call Grid_putPointData(blockId, CENTER, PRES_VAR, EXTERIOR, axis, presZone)
           call Grid_putPointData(blockId, CENTER, VELX_VAR, EXTERIOR, axis, velxZone)
           call Grid_putPointData(blockId, CENTER, VELY_VAR, EXTERIOR, axis, velyZone)
           call Grid_putPointData(blockId, CENTER, VELZ_VAR, EXTERIOR, axis, velzZone)

#ifdef ENER_VAR
           call Grid_putPointData(blockId, CENTER, ENER_VAR, EXTERIOR, axis, enerZone)
#endif
#ifdef EINT_VAR
           call Grid_putPointData(blockId, CENTER, EINT_VAR, EXTERIOR, axis, eintZone)
#endif
#ifdef GAME_VAR
           call Grid_putPointData(blockId, CENTER, GAME_VAR, EXTERIOR, axis, gameZone)
#endif
#ifdef GAMC_VAR
           call Grid_putPointData(blockId, CENTER, GAMC_VAR, EXTERIOR, axis, gamcZone)
#endif


        enddo
     enddo
  enddo
           ! store the variables in the current zone via Grid put methods
           ! data is put stored one cell at a time with these calls to Grid_putData           

           call Grid_putPointData(blockId, CENTER, DENS_VAR, EXTERIOR, axis, rhoZone)
           call Grid_putPointData(blockId, CENTER, PRES_VAR, EXTERIOR, axis, presZone)
           call Grid_putPointData(blockId, CENTER, VELX_VAR, EXTERIOR, axis, velxZone)
           call Grid_putPointData(blockId, CENTER, VELY_VAR, EXTERIOR, axis, velyZone)
           call Grid_putPointData(blockId, CENTER, VELZ_VAR, EXTERIOR, axis, velzZone)

!#ifdef ENER_VAR
          ! call Grid_putPointData(blockId, CENTER, ENER_VAR, EXTERIOR, axis, enerZone)   
!#endif
!#ifdef EINT_VAR
!          ! call Grid_putPointData(blockId, CENTER, EINT_VAR, EXTERIOR, axis, eintZone)   
!#endif
!#ifdef GAME_VAR          
          ! call Grid_putPointData(blockId, CENTER, GAME_VAR, EXTERIOR, axis, gameZone)
!#endif
!#ifdef GAMC_VAR
!           call Grid_putPointData(blockId, CENTER, GAMC_VAR, EXTERIOR, axis, gamcZone)
!#endif
!#ifdef TEMP_VAR
!           call Grid_putPointData(blockId, CENTER, TEMP_VAR, EXTERIOR, axis, 1.e-10)
!# endif

!! Cleanup!  Must deallocate arrays

  deallocate(xCoord)
  deallocate(yCoord)
  deallocate(zCoord)
  call Grid_releaseBlkPtr(blockID,solnData)
 
  return
end subroutine Simulation_initBlock

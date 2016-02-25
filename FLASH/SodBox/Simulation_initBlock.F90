subroutine Simulation_initBlock(blockID)

#include "constants.h"
#include "Flash.h"
#include "Eos.h"

  use Simulation_data, ONLY: sim_posn, sim_xCos, sim_yCos, sim_zCos, &    
     &  sim_rhoLeft,  sim_pLeft, sim_uLeft, sim_rhoRight, sim_pRight, sim_uRight, &
     &  sim_smallX, sim_gamma, sim_smallP

#ifdef FLASH_3T
  use Simulation_data, ONLY : &
       sim_pionLeft, sim_peleLeft, sim_pradLeft, &
       sim_pionRight, sim_peleRight, sim_pradRight
#endif
     
  use Grid_interface, ONLY : Grid_getBlkIndexLimits, &
    Grid_getCellCoords, Grid_putPointData
  use Eos_interface, ONLY : Eos, Eos_wrapped


  implicit none

  ! compute the maximum length of a vector in each coordinate direction 
  ! (including guardcells)
  
  integer, intent(in) :: blockID
  

  integer :: i, j, k, n
  integer :: iMax, jMax, kMax
  


  real :: xx, yy,  zz, rr!xxL, xxR
  
  real :: lPosn0, lPosn
  

  real,allocatable, dimension(:) ::xCenter,xLeft,xRight,yCoord,zCoord

  integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC
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
  call Grid_getBlkIndexLimits(blockId,blkLimits,blkLimitsGC)
  
  sizeX = blkLimitsGC(HIGH,IAXIS)
  sizeY = blkLimitsGC(HIGH,JAXIS)
  sizeZ = blkLimitsGC(HIGH,KAXIS)
  allocate(xLeft(sizeX))
  allocate(xRight(sizeX))
  allocate(xCenter(sizeX))
  allocate(yCoord(sizeY))
  allocate(zCoord(sizeZ))
  xCenter = 0.0
  xLeft = 0.0
  xRight = 0.0
  yCoord = 0.0
  zCoord = 0.0

  if (NDIM == 3) call Grid_getCellCoords&
                      (KAXIS, blockId, CENTER,gcell, zCoord, sizeZ)
  if (NDIM >= 2) call Grid_getCellCoords&
                      (JAXIS, blockId, CENTER,gcell, yCoord, sizeY)

  call Grid_getCellCoords(IAXIS, blockId, LEFT_EDGE, gcell, xLeft, sizeX)
  call Grid_getCellCoords(IAXIS, blockId, CENTER, gcell, xCenter, sizeX)
  call Grid_getCellCoords(IAXIS, blockId, RIGHT_EDGE, gcell, xRight, sizeX)

!------------------------------------------------------------------------------

! Loop over cells in the block.  For each, compute the physical position of 
! its left and right edge and its center as well as its physical width.  
! Then decide which side of the initial discontinuity it is on and initialize 
! the hydro variables appropriately.


  do k = blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS)
     ! get the coordinates of the cell center in the z-direction
     zz = zCoord(k)
     do j = blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)
        ! get the coordinates of the cell center in the y-direction
        yy = yCoord(j)
        ! The position of the shock in the current yz-row.
       ! lPosn = lPosn0 - yy*sim_yCos/sim_xCos
        do i = blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)
           ! get the cell center, left, and right positions in x
           xx  = xCenter(i)
     !      xxL = xLeft(i)
      !     xxR = xRight(i)
           rr = xx**2 + yy**2 + zz**2
!           print *,"xx,yy,zz: ", xx,yy,zz
!           print *,"rr: ", rr
           if (rr <= 0.25) then
!               print *,"Inside"
               rhoZone =sim_rhoRight 
           else
!               print *,"Outisde"
               rhoZone =sim_rhoLeft
           endif
    
           presZone = sim_pRight
!           rhoZone = sim_rhoRight
           velxZone = 0.0
           velyZone = 0.0
           velzZone = 0.0
           axis(IAXIS) = i
           axis(JAXIS) = j
           axis(KAXIS) = k

           !put in default mass fraction values of all species
           if (NSPECIES > 0) then
              call Grid_putPointData(blockID, CENTER, SPECIES_BEGIN, EXTERIOR, &
                   axis, 1.0e0-(NSPECIES-1)*sim_smallX)


              !if there is only 1 species, this loop will not execute
              do n = SPECIES_BEGIN+1,SPECIES_END
                 call Grid_putPointData(blockID, CENTER, n, EXTERIOR, &
                      axis, sim_smallX)
              enddo
           end if

           ! Compute the gas energy and set the gamma-values needed for the equation of 
           ! state.
           ekinZone = 0.5 * (velxZone**2 + & 
                velyZone**2 + & 
                velzZone**2)
           
           eintZone = presZone / (sim_gamma-1.)
           eintZone = eintZone / rhoZone
           gameZone = sim_gamma
           gamcZone = sim_gamma
           enerZone = eintZone + ekinZone
           enerZone = max(enerZone, sim_smallP)

           ! store the variables in the current zone via Grid put methods
           ! data is put stored one cell at a time with these calls to Grid_putData           


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
#ifdef TEMP_VAR
# ifdef SIMULATION_TWO_MATERIALS
           call Grid_putPointData(blockId, CENTER, TEMP_VAR, EXTERIOR, axis, eosData(EOS_TEMP))
# else
           call Grid_putPointData(blockId, CENTER, TEMP_VAR, EXTERIOR, axis, 1.e-10)
# endif
#endif

#ifdef SIMULATION_TWO_MATERIALS
           call Grid_putPointData(blockID, CENTER, LEFT_SPEC, EXTERIOR, &
                   axis, mfrac(LEFT_SPEC-SPECIES_BEGIN+1) )
           call Grid_putPointData(blockID, CENTER, RGHT_SPEC, EXTERIOR, &
                   axis, mfrac(RGHT_SPEC-SPECIES_BEGIN+1) )
#endif

#ifdef FLASH_3T
           ! We must now compute the internal energy from the pressure
           ! for the ions, electrons, and radiation field:
           
           ! Electrons...
           eeleZone = peleZone / (sim_gamma - 1.0) / rhoZone
           eionZone = pionZone / (sim_gamma - 1.0) / rhoZone
           eradZone = 3.0 * pradZone / rhoZone
           
           call Grid_putPointData(blockId, CENTER, EELE_VAR, EXTERIOR, axis, eeleZone)
           call Grid_putPointData(blockId, CENTER, EION_VAR, EXTERIOR, axis, eionZone)
           call Grid_putPointData(blockId, CENTER, ERAD_VAR, EXTERIOR, axis, eradZone)
#endif
        enddo
     enddo
  enddo

! #ifdef EELE_VAR
!   call Eos_wrapped(MODE_DENS_EI_SCATTER,blkLimits,blockId)
! #endif

!   do k = blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS)
!      do j = blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)
!         do i = blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)
!            axis(IAXIS) = i
!            axis(JAXIS) = j
!            axis(KAXIS) = k
! #ifdef ERAD_VAR
!            call Grid_putPointData(blockId, CENTER, ERAD_VAR, EXTERIOR, axis, 0.0  )   
! #endif
! #ifdef E3_VAR
!            call Grid_putPointData(blockId, CENTER, E3_VAR,   EXTERIOR, axis, 0.0  )   
! #endif

! #ifdef PRAD_VAR
!            call Grid_putPointData(blockId, CENTER, PRAD_VAR, EXTERIOR, axis, 0.0  )   
! #endif
! #ifdef TRAD_VAR
!            call Grid_putPointData(blockId, CENTER, TRAD_VAR, EXTERIOR, axis, 0.0  )   
! #endif
!         enddo
!      enddo
!   enddo

!! Cleanup!  Must deallocate arrays

  deallocate(xLeft)
  deallocate(xRight)
  deallocate(xCenter)
  deallocate(yCoord)
  deallocate(zCoord)

 
  return
end subroutine Simulation_initBlock

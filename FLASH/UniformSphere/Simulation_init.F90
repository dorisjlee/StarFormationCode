!! SYNOPSIS
!!
!!  Simulation_init()
!!
!! DESCRIPTION
!!
!!  Initializes all the data specified in Simulation_data.
!!  It calls RuntimeParameters_get rotuine for initialization.
!!  Initializes initial conditions for Sedov Spherical Explosion 
!!  problem.
!!
!!***

subroutine Simulation_init()
!  use Simulation_data 
!  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
!  use Logfile_interface,           ONLY : Logfile_stamp
  use Simulation_data
  use Driver_interface, ONLY : Driver_getMype
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use Logfile_interface, ONLY : Logfile_stampMessage
  use PhysicalConstants_interface, ONLY:  PhysicalConstants_get
  implicit none

#include "constants.h"
#include "Flash.h"
  !real, save :: sim_gamma, sim_gascon
  logical :: threadBlockListBuild, threadWithinBlockBuild  

  call RuntimeParameters_get("rhoOut",rhoOut)
  call RuntimeParameters_get("rhoIn",rhoIn)
  call RuntimeParameters_get("P", P)
  call RuntimeParameters_get("rcloud",rcloud)
  call RuntimeParameters_get("sim_xctr",sim_xctr)
  call RuntimeParameters_get("sim_yctr",sim_yctr)
  call RuntimeParameters_get("sim_zctr",sim_zctr)
!  call RuntimeParameters_get("sim_gamma",sim_gamma)
!  call RuntimeParameters_get("sim_gascon",sim_gascon)
  print *, "Finished calling RuntimeParameters_get inside Simulation_init"
  !Any calculation that need to be only done once at the beginning of the run  
!  sim_gamma = 1.0001 
!  sim_gascon =8.2544E7
end subroutine Simulation_init

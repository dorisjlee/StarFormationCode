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
  use Simulation_data 
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use Logfile_interface,           ONLY : Logfile_stamp

  implicit none

#include "constants.h"
#include "Flash.h"

  logical :: threadBlockListBuild, threadWithinBlockBuild  

  call RuntimeParameters_get('rhoOut',rhoOut)
  call RuntimeParameters_get('rhoIn',rhoIn)
  call RuntimeParameters_get('P', P)

  !Any calculation that need to be only done once at the beginning of the run  
end subroutine Simulation_init

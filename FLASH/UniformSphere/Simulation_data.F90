!!****if* source/Simulation/SimulationMain/Sedov/Simulation_data
!!
!! NAME
!!
!!  Simulation_data
!!
!! SYNOPSIS
!!
!!  use Simulation_data 
!!
!!  DESCRIPTION
!!
!!  Stores the local data for Simulation setup: Sedov
!!  
!! PARAMETERS
!!
!!
!!
!!***

module Simulation_data

  implicit none
#include "constants.h"

  !! *** Runtime Parameters *** !!

  real, save     :: rhoOut, rhoIn,P,rcloud
  !! *** Variables pertaining to this Simulation *** !!

  !integer, parameter                  :: sim_nProfile = 1000
  !real   , save                       :: sim_inSubZones, sim_inSubzm1
  !real, dimension(sim_nProfile), save :: sim_rProf, sim_rhoProf, sim_pProf
  !logical, save :: sim_threadBlockList = .false.
  
end module Simulation_data

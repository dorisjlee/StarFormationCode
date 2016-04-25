
module Simulation_data
#include "Flash.h"
  implicit none

  integer,save  :: sim_globalMe
  real, save :: sim_rhoLeft, sim_rhoRight, sim_pLeft, sim_pRight, fattening_factor 
  real, save :: sim_uLeft, sim_uRight, sim_xAngle, sim_yAngle, sim_posn
  real, save :: sim_gamma, sim_smallP, sim_smallX
  logical, save :: sim_testInitialized
  real, save    :: sim_sink_x, sim_sink_y, sim_sink_z, sim_sink_vx, sim_sink_vy, sim_sink_vz, sim_sink_mass
  !! *** Variables pertaining to Simulation Setup 'Sod' *** !!
  real, save :: sim_xCos, sim_yCos, sim_zCos
  logical, save :: sim_gCell

  integer, save :: sim_meshMe
end module Simulation_data



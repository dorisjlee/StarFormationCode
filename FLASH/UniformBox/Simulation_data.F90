module Simulation_data
  implicit none
#include "constants.h"
  !! *** Runtime Parameters *** !!
  real, save     :: rhoOut, rhoIn,P,rcloud
  real, save     :: sim_gascon,sim_gamma,sim_xctr,sim_yctr,sim_zctr
  !integer, parameter                  :: sim_nProfile = 1000
  !real   , save                       :: sim_inSubZones, sim_inSubzm1
  !real, dimension(sim_nProfile), save :: sim_rProf, sim_rhoProf, sim_pProf
  !logical, save :: sim_threadBlockList = .false.
end module Simulation_data

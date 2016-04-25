!!****if* source/Simulation/SimulationMain/Sod/Simulation_init

subroutine Simulation_init()
  
  use Simulation_data
  use Driver_interface, ONLY : Driver_getMype, Driver_abortFlash
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use Logfile_interface, ONLY : Logfile_stamp
  use Particles_sinkData
  use pt_sinkInterface, only: pt_sinkCreateParticle, pt_sinkGatherGlobal
  implicit none
#include "constants.h"
#include "Flash.h"
#include "Particles.h"
  
  integer :: myPE, pno, blockID
  real    :: pt
  logical :: restart

  ! do nothing on restart
  call RuntimeParameters_get("restart", restart)
  if (restart) return

  call Driver_getMype(GLOBAL_COMM, myPE)
  sim_globalMe = myPE
  call RuntimeParameters_get('smallp', sim_smallP)
  call RuntimeParameters_get('smallx', sim_smallX) 
  call RuntimeParameters_get('fattening_factor',fattening_factor ) 
  call RuntimeParameters_get('gamma', sim_gamma)
  
  call RuntimeParameters_get('sim_rhoLeft', sim_rhoLeft)
  call RuntimeParameters_get('sim_rhoRight', sim_rhoRight)
  
  call RuntimeParameters_get('sim_pLeft', sim_pLeft)
  call RuntimeParameters_get('sim_pRight', sim_pRight)
  
  call RuntimeParameters_get('sim_uLeft', sim_uLeft)
  call RuntimeParameters_get('sim_uRight', sim_uRight)
  
  call RuntimeParameters_get('sim_xangle', sim_xAngle)
  call RuntimeParameters_get('sim_yangle', sim_yAngle)
  
  call RuntimeParameters_get('sim_posn', sim_posn)


  call Logfile_stamp( "initializing Sod problem",  &
       "[Simulation_init]")
     
! place initial sink particle

  if (sim_globalMe == MASTER_PE) then

    call RuntimeParameters_get("sim_sink_x", sim_sink_x)
    call RuntimeParameters_get("sim_sink_y", sim_sink_y)
    call RuntimeParameters_get("sim_sink_z", sim_sink_z)
    call RuntimeParameters_get("sim_sink_vx", sim_sink_vx)
    call RuntimeParameters_get("sim_sink_vy", sim_sink_vy)
    call RuntimeParameters_get("sim_sink_vz", sim_sink_vz)
    call RuntimeParameters_get("sim_sink_mass", sim_sink_mass)

    blockID = 1
    pt = 0.0

    pno = pt_sinkCreateParticle(sim_sink_x, sim_sink_y, sim_sink_z, pt, blockID, sim_globalMe)

    particles_local(VELX_PART_PROP, 1) = sim_sink_vx
    particles_local(VELY_PART_PROP, 1) = sim_sink_vy
    particles_local(VELZ_PART_PROP, 1) = sim_sink_vz
    particles_local(MASS_PART_PROP, 1) = sim_sink_mass

    write(*,'(A,4(1X,ES16.9),3I8)') "initial sink particle created (x, y, z, pt, blockID, MyPE, tag): ", &
      & sim_sink_x, sim_sink_y, sim_sink_z, pt, blockID, sim_globalMe, int(particles_local(TAG_PART_PROP,pno))

  endif

 !Looks like tolerance level set for mass, px,py,pz, seems to be specific to the mom test
 !call RuntimeParameters_get("sim_massTol", sim_massTol)
 ! call RuntimeParameters_get("sim_momXTol", sim_momXTol)
 ! call RuntimeParameters_get("sim_momYTol", sim_momYTol)
 ! call RuntimeParameters_get("sim_momZTol", sim_momZTol)


  call pt_sinkGatherGlobal()

  sim_testInitialized = .FALSE.


end subroutine Simulation_init

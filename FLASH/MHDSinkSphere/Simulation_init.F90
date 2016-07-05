!!****if* source/Simulation/SimulationMain/Sod/Simulation_init

subroutine Simulation_init()
  
  use Simulation_data
  use Driver_interface, ONLY : Driver_getMype, Driver_abortFlash
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use Logfile_interface, ONLY : Logfile_stamp
  implicit none
#include "constants.h"
#include "Flash.h"
  
  integer :: myPE, pno, blockID
  real    :: pt
  logical :: restart

  ! do nothing on restart
  call RuntimeParameters_get("restart", restart)

  call Driver_getMype(GLOBAL_COMM, myPE)
  sim_globalMe = myPE
  call RuntimeParameters_get('smallp', sim_smallP)
  call RuntimeParameters_get('smallx', sim_smallX) 
  call RuntimeParameters_get('fattening_factor',fattening_factor ) 
  call RuntimeParameters_get('beta_param',beta_param) !rotation beta
  call RuntimeParameters_get('gamma', sim_gamma)
  call RuntimeParameters_get('alpha',alpha) ! thermal-to-magnetic pressure ratio  

  call Logfile_stamp( "initializing Sod problem",  &
       "[Simulation_init]")
     
! place initial sink particle
!if (restart) then
!	  print *,"placing sink particle on restart"
!	  if (sim_globalMe == MASTER_PE) then
!
!	    call RuntimeParameters_get("sim_sink_x", sim_sink_x)
!	    call RuntimeParameters_get("sim_sink_y", sim_sink_y)
!	    call RuntimeParameters_get("sim_sink_z", sim_sink_z)
!	    call RuntimeParameters_get("sim_sink_vx", sim_sink_vx)
!	    call RuntimeParameters_get("sim_sink_vy", sim_sink_vy)
!	    call RuntimeParameters_get("sim_sink_vz", sim_sink_vz)
!	    call RuntimeParameters_get("sim_sink_mass", sim_sink_mass)
!
!	    blockID = 1
!	    pt = 0.0
!
!	    pno = pt_sinkCreateParticle(sim_sink_x, sim_sink_y, sim_sink_z, pt, blockID, sim_globalMe)
!
!	    particles_local(VELX_PART_PROP, 1) = sim_sink_vx
!	    particles_local(VELY_PART_PROP, 1) = sim_sink_vy
!	    particles_local(VELZ_PART_PROP, 1) = sim_sink_vz
!	    particles_local(MASS_PART_PROP, 1) = sim_sink_mass
!
!	    write(*,'(A,4(1X,ES16.9),3I8)') "initial sink particle created (x, y, z, pt, blockID, MyPE, tag): ", &
!	      & sim_sink_x, sim_sink_y, sim_sink_z, pt, blockID, sim_globalMe, int(particles_local(TAG_PART_PROP,pno))
!
!	  endif
!	  call pt_sinkGatherGlobal()

!	  sim_testInitialized = .FALSE.
!endif

end subroutine Simulation_init

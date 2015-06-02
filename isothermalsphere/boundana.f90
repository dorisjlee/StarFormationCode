!############################################################
!############################################################
!############################################################
!############################################################
subroutine boundana(x,u,dx,ibound,ncell)
  use amr_parameters, ONLY: dp,ndim,nvector,boxlen
  use hydro_parameters, ONLY: nvar,boundary_var,gamma
  implicit none
  integer ::ibound                        ! Index of boundary region
  integer ::ncell                         ! Number of active cells
  real(dp)::dx                            ! Cell size
  real(dp),dimension(1:nvector,1:nvar)::u ! Conservative variables
  real(dp),dimension(1:nvector,1:ndim)::x ! Cell center position.
  !================================================================
  ! This routine generates boundary conditions for RAMSES.
  ! Positions are in user units:
  ! x(i,1:3) are in [0,boxlen]**ndim.
  ! U is the conservative variable vector. Conventions are here:
  ! U(i,1): d, U(i,2:ndim+1): d.u,d.v,d.w and U(i,ndim+2): E.
  ! U is in user units.
  ! ibound is the index of the boundary region defined in the namelist.
  !================================================================

  integer::ivar,i
  integer ::nn                            ! Number of cells
  real(dp),dimension(1:nvector,1:nvar),save::q   ! Primitive variables

  real  xl,xr,xc,yl,yr,yc,zr,zl,zc,rr
  do ivar=1,nvar
     do i=1,ncell
        u(i,ivar)=boundary_var(ibound,ivar)
     end do
  end do

  ! Add here, if you wish, some user-defined boudary conditions
  do i=1,nn
     xl=x(i,1)-0.5*dx-boxlen/2.0
     xr=x(i,1)+0.5*dx-boxlen/2.0
     xc=x(i,1)-boxlen/2.0
     yl=x(i,2)-0.5*dx-boxlen/2.0
     yr=x(i,2)+0.5*dx-boxlen/2.0
     yc=x(i,2)-boxlen/2.0
     zl=x(i,3)-0.5*dx-boxlen/2.0
     zr=x(i,3)+0.5*dx-boxlen/2.0
     zc=x(i,3)-boxlen/2.0
     rr=sqrt(xc**2+yc**2+zc**2)
     !Enforce V=0 of the cloud at the radius at all r
     if (rr .EQ. 1.0E10) then
	 print *,"At rmax!"
         !q(i,2)=0.0
         !q(i,3)=0.0
         !q(i,4)=0.0	
	 U(i,2)=0.0d0
         U(i,3)=0.0d0
         U(i,4)=0.0d0
     endif
   enddo
  ! Convert primitive to conservative variables
  ! density -> density
!  u(1:nn,1)=q(1:nn,1)
  ! velocity -> momentum
  ! kinetic energy
  ! pressure -> total fluid energy
!  u(1:nn,ndim+2)=u(1:nn,ndim+2)+q(1:nn,ndim+2)/(gamma-1.0d0)



end subroutine boundana

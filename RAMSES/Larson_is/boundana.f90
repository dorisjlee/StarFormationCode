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
  real(dp),dimension(1:nvector,1:nvar),save::q   ! Primitive variables
  real  xl,xr,xc,yl,yr,yc,zr,zl,zc,rr
  print *, "Dealing with Boundary Conditions"
  do ivar=1,nvar
     do i=1,ncell
        u(i,ivar)=boundary_var(ibound,ivar)
     end do
  end do
  !!print *,"U:",u
  print *,"ncell:",ncell
  ! Add here, if you wish, some user-defined boudary conditions
  do i=1,32
     !print *,"in here"
     xc=x(i,1)-boxlen/2.0
     yc=x(i,2)-boxlen/2.0
     zc=x(i,3)-boxlen/2.0
     rr=sqrt(xc**2+yc**2+zc**2)
     !Enforce V=0 of the cloud at the radius at all r
     !print *,i,"xc: ",xc
     !print *,i,"yc: ",yc
     !print *,i,"zc: ",zc
     !print *,"before if"
     print *,"rr:",rr

     if ((rr .LE. 1.1E10) .AND. (rr .GE. 0.9E10)) then
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
!  u(1:ncell,1)=q(1:ncell,1)
  ! velocity -> momentum
  ! kinetic energy
  ! pressure -> total fluid energy
!  u(1:ncell,ndim+2)=u(1:ncell,ndim+2)+q(1:ncell,ndim+2)/(gamma-1.0d0)



end subroutine boundana

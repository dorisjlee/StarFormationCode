subroutine condinit(x,u,dx,nn)
  use amr_parameters
  use hydro_parameters
  implicit none
  integer ::nn                            ! Number of cells
  real(dp)::dx                            ! Cell size
  real(dp),dimension(1:nvector,1:nvar)::u ! Conservative variables
  real(dp),dimension(1:nvector,1:ndim)::x ! Cell center position.
  !================================================================
  ! This routine generates initial conditions for RAMSES.
  ! Positions are in user units:
  ! x(i,1:3) are in [0,boxlen]**ndim.
  ! U is the conservative variable vector. Conventions are here:
  ! U(i,1): d, U(i,2:ndim+1): d.u,d.v,d.w and U(i,ndim+2): E.
  ! Q is the primitive variable vector. Conventions are here:
  ! Q(i,1): d, Q(i,2:ndim+1):u,v,w and Q(i,ndim+2): P.
  ! If nvar >= ndim+3, remaining variables are treated as passive
  ! scalars in the hydro solver.
  ! U(:,:) and Q(:,:) are in user units.
  !================================================================
  integer::ivar
  real(dp),dimension(1:nvector,1:nvar),save::q   ! Primitive variables
  real rmax, rho0,P0,xl,xr,xc,yl,yr,yc,zr,zl,zc,rr
  integer i
  !print *, "Inside condinit.f90"
  ! Call built-in initial condition generator
  !call region_condinit(x,q,dx,nn)
  q(1:nn,:)=0.
  !print *,"nn:",nn
  ! Add here, if you wish, some user-defined initial conditions
  do i=1,nn !looping through all grid numbers
    !Defining the left, right, and center positions of each cells.
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
     !print *,"rr_init: ",rr
     !G=1 for self gravity
     rmax=3.27
     rho0=0.02806
     P0= 0.0359
     !Density
     IF (rr .LE. rmax) THEN
        !PRINT *,"Inside Box"
	q(i,1)=rho0
     ELSE !the rest of the box
        !PRINT *,"Outside box"
        !PRINT *, "Radius: ",rr
        q(i,1)=1.0E-6 !few orders of magnitude less dense (~approx 0?)
     END IF
     !Initially static cloud
     q(i,2)=0.0      ! Velocity x
     q(i,3)=0.0      ! Velocity y
     q(i,4)=0.0      ! Velocity z
     !Pressure 
     IF (rr .LE. rmax) THEN     
	q(i,5)=P0
     ELSE
	q(i,5)=P0
     END IF
  end do 
  !Convert primitive to conservative variables
  ! density -> density
  u(1:nn,1)=q(1:nn,1)
  ! velocity -> momentum
  u(1:nn,2)=q(1:nn,1)*q(1:nn,2)
#if NDIM>1
  u(1:nn,3)=q(1:nn,1)*q(1:nn,3)
#endif
#if NDIM>2
  u(1:nn,4)=q(1:nn,1)*q(1:nn,4)
#endif
  ! kinetic energy
  u(1:nn,ndim+2)=0.0d0
  u(1:nn,ndim+2)=u(1:nn,ndim+2)+0.5*q(1:nn,1)*q(1:nn,2)**2
#if NDIM>1
  u(1:nn,ndim+2)=u(1:nn,ndim+2)+0.5*q(1:nn,1)*q(1:nn,3)**2
#endif
#if NDIM>2
  u(1:nn,ndim+2)=u(1:nn,ndim+2)+0.5*q(1:nn,1)*q(1:nn,4)**2
#endif
  ! thermal pressure -> total fluid energy
  u(1:nn,ndim+2)=u(1:nn,ndim+2)+q(1:nn,ndim+2)/(gamma-1.0d0)
#if NENER>0
  ! radiative pressure -> radiative energy
  ! radiative energy -> total fluid energy
  do ivar=1,nener
     u(1:nn,ndim+2+ivar)=q(1:nn,ndim+2+ivar)/(gamma_rad(ivar)-1.0d0)
     u(1:nn,ndim+2)=u(1:nn,ndim+2)+u(1:nn,ndim+2+ivar)
  enddo
#endif
#if NVAR>NDIM+2+NENER
  ! passive scalars
  do ivar=ndim+3+nener,nvar
     u(1:nn,ivar)=q(1:nn,1)*q(1:nn,ivar)
  end do
#endif

end subroutine condinit
!================================================================
!================================================================
!================================================================
!================================================================
subroutine velana(x,v,dx,t,ncell)
  use amr_parameters
  use hydro_parameters  
  implicit none
  integer ::ncell                         ! Size of input arrays
  real(dp)::dx                            ! Cell size
  real(dp)::t                             ! Current time
  real(dp),dimension(1:nvector,1:3)::v    ! Velocity field
  real(dp),dimension(1:nvector,1:ndim)::x ! Cell center position.
  !================================================================
  ! This routine computes the user defined velocity fields.
  ! x(i,1:ndim) are cell center position in [0,boxlen] (user units).
  ! v(i,1:3) is the imposed 3-velocity in user units.
  !================================================================
  integer::i
  real(dp)::xx,yy,zz,vx,vy,vz,rr,tt,rmax,margin

  ! Add here, if you wish, some user-defined initial conditions
  
  rmax=1.07483E10
  margin=1.0E-2
  do i=1,ncell
     xx=x(i,1)
     yy=x(i,2)
     zz=x(i,3)
     rr=sqrt(xx**2+yy**2+zz**2)
     IF (rr<rmax+margin .AND. rr>rmax-margin) THEN
     	v(i,1)=0.0
     	v(i,2)=0.0
     	v(i,3)=0.0
     end IF
  end do
end subroutine velana

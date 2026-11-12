!***********************************************************************
! comm.f90 -- shared data for the ILW-WENO solver
!
!   u_t = ( f(u) )_xx  ( + ( f(u) )_yy in two dimensions ),  f(u) = u^m
!
!   iorder = 6 : WENO-ZQ-522  (6-point stencil, 3 sub-stencils of 4
!                points, mapped weights, ideal weights -2/15, 19/15, -2/15)
!   iorder = 8 : WENO-ZQ      (8-point stencil, 3 sub-stencils of 6
!                points, mapped weights, ideal weights -9/56, 37/28, -9/56)
!
!   ischeme = 3  : SSP-RK3              (baseline, 3 expensive evaluations)
!   ischeme = 4  : classical RK4        (baseline, 4 expensive evaluations)
!   ischeme = 13 : ILW3  (improved Lax-Wendroff on the SSP-RK3 tableau)
!   ischeme = 14 : ILW4  (improved Lax-Wendroff on the classical RK4 tableau)
!
!   nq = 2, 4 or 6 : order of the cheap component-wise operator C
!***********************************************************************
MODULE comm
  IMPLICIT NONE
  SAVE
  INTEGER, PARAMETER :: md   = 4      ! halo width (>= 4 for the 8th-order stencil)
  INTEGER, PARAMETER :: nmax = 4000   ! upper bound for nx

  INTEGER :: idim, iorder, nq, icase, ischeme, ibc, ipos, iops, ieq
  INTEGER :: iblflux      ! 0: Buckley-Leverett flux without gravity, 1: with
  REAL(8) :: ubcl, ubcr   ! Dirichlet data at x = xl, xr (ibc = 3)
  REAL(8) :: epsbl        ! diffusion coefficient of the BL and degenerate cases
  INTEGER :: nx, ny, mt, nlev, nx0, ny0, nsrc
  INTEGER :: jb0, jb1     ! bounds of the transverse direction (0:0 in 1D)
  REAL(8) :: dx, dy, cdx, cdy, cfl, tend, tnum, dt, pmem, tini
  REAL(8) :: xl, xr, yl, yr, aam
  REAL(8) :: kconv, sconv     ! convection term  (b(u) = kconv*u^sconv)_x
  REAL(8) :: alf              ! Lax-Friedrichs speed for the convective flux
  INTEGER :: isgn             ! upwind direction of b'(u): 1, -1, or 0 (LF)

  REAL(8), ALLOCATABLE :: u0(:,:)      ! U^n
  REAL(8), ALLOCATABLE :: s1(:,:), s2(:,:), s3(:,:)       ! stage states
  REAL(8), ALLOCATABLE :: r1(:,:), r2(:,:), r3(:,:), r4(:,:)  ! residuals
  REAL(8), ALLOCATABLE :: Lop(:,:)     ! expensive operator at U^n = r0
  REAL(8), ALLOCATABLE :: Cop(:,:)     ! cheap operator at U^n
  REAL(8), ALLOCATABLE :: dc(:,:)      ! frozen correction d = L(U^n) - C(U^n)
  REAL(8), ALLOCATABLE :: fw(:,:)      ! flux field f(u) = u^m
  REAL(8), ALLOCATABLE :: fwb(:,:)     ! convective flux field b(u)
  REAL(8), ALLOCATABLE :: uw(:,:)      ! copy of the state (for the LF splitting)
  REAL(8), ALLOCATABLE :: fb1d(:)

  REAL(8), ALLOCATABLE :: f1d(:), fh1d(:)

  REAL(8) :: amat(0:4,0:4), bvec(0:4), cvec(0:4)
END MODULE comm

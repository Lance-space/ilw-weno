!***********************************************************************
!  Example 6.5 with scheme ILW3-WENO6 : standalone program
!  Compile:  f90 /optimize:4 /real_size:64 thisfile.f90
!            ifort -O3 -r8 thisfile.f90 -o run
!  Run:      run   (writes results.csv and profile.plt)
!***********************************************************************
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

!***********************************************************************
! para.f90 -- parameters, allocation, Butcher tableaux, grid
!
!   The parameters are read from the namelist file  input.nml
!   (all of them have a default, so an absent file is harmless).
!***********************************************************************
MODULE para
  USE comm
  IMPLICIT NONE
CONTAINS

  SUBROUTINE Canshu()
    INTEGER :: ios

    ! ---------------- defaults ----------------
    idim    = 1        ! 1 or 2 space dimensions
    iorder  = 6        ! spatial order of the WENO operator: 6 or 8
    nq      = 4        ! order of the cheap centred derivative: 2, 4 or 6
    icase   = 1        ! test case
    ischeme = 13       ! 3, 4 (RK baselines) or 13, 14 (ILW)
    nlev    = 1        ! number of grid levels of a convergence study
    nx0     = 200      ! coarsest grid
    ny0     = 100
    cfl     = 0.4d0
    tend    = -1.d0    ! < 0 : use the default of the test case
    tini    = -1.d0   ! < 0 : use the default of the test case
    pmem    = 2.0d0    ! exponent m of the porous medium equation
    ipos    = 1        ! 1: f(u) = max(u,0)^m, 0: f(u) = u^m
    ieq     = -1       ! -1: from icase; 1: g=u^m, 2: g=u, 3: BL, 4: strongly degenerate
    epsbl   = 0.01d0   ! diffusion coefficient of cases 4, 5, 8, 9
    ubcl    = 0.d0
    ubcr    = 0.d0
    iblflux = 0        ! 0: BL flux without gravity, 1: with
    nsrc    = -1       ! -1: decide from icase, 0: no source, 1: source (control runs)
    iops    = 0        ! 0: WENO operator L,  1: centred difference of order nq in place of L
    kconv   = 0.0d0    ! convection  (kconv*u^sconv)_x
    sconv   = 1.0d0
    ibc     = 0        ! 0: default of the test case, 1: periodic, 2: zero, 3: transmissive
    IF (ios == 0) THEN
    END IF

    IF (nx0 > nmax) THEN
       WRITE(*,*) 'Canshu: nx0 too large, reset to nmax'
       nx0 = nmax
    END IF
    IF (ny0 > nmax) ny0 = nmax
  END SUBROUTINE Canshu

  !---------------------------------------------------------------------
  !  Butcher tableaux (entries with j > l only)
  !    SSP-RK3 : stages 0, 1, 1/2 ;  b = 1/6, 1/6, 2/3
  !    RK4     : stages 0, 1/2, 1/2, 1 ;  b = 1/6, 1/3, 1/3, 1/6
  !---------------------------------------------------------------------
  SUBROUTINE SetTableau()
    INTEGER :: j, l
    amat = 0.d0; bvec = 0.d0; cvec = 0.d0
    IF (mt == 3) THEN
       amat(2,1) = 1.d0
       amat(3,1) = 0.25d0
       amat(3,2) = 0.25d0
       bvec(1)   = 1.d0/6.d0
       bvec(2)   = 1.d0/6.d0
       bvec(3)   = 2.d0/3.d0
    ELSE
       amat(2,1) = 0.5d0
       amat(3,1) = 0.d0
       amat(3,2) = 0.5d0
       amat(4,1) = 0.d0
       amat(4,2) = 0.d0
       amat(4,3) = 1.d0
       bvec(1)   = 1.d0/6.d0
       bvec(2)   = 1.d0/3.d0
       bvec(3)   = 1.d0/3.d0
       bvec(4)   = 1.d0/6.d0
    END IF
    DO j = 1, mt
       cvec(j) = 0.d0
       DO l = 1, j-1
          cvec(j) = cvec(j) + amat(j,l)
       END DO
    END DO
  END SUBROUTINE SetTableau

  SUBROUTINE SetScheme()
    IF (ischeme == 3 .OR. ischeme == 13) THEN
       mt = 3
    ELSE
       mt = 4
    END IF
    CALL SetTableau()
  END SUBROUTINE SetScheme

  SUBROUTINE Alloc()
    IF (ALLOCATED(u0)) THEN
       DEALLOCATE(u0, s1, s2, s3, r1, r2, r3, r4, Lop, Cop, dc, fw, fwb, uw)
    END IF
    ALLOCATE( u0(-md:nx+md, jb0:jb1) )
    ALLOCATE( s1(-md:nx+md, jb0:jb1) )
    ALLOCATE( s2(-md:nx+md, jb0:jb1) )
    ALLOCATE( s3(-md:nx+md, jb0:jb1) )
    ALLOCATE( r1(-md:nx+md, jb0:jb1) )
    ALLOCATE( r2(-md:nx+md, jb0:jb1) )
    ALLOCATE( r3(-md:nx+md, jb0:jb1) )
    ALLOCATE( r4(-md:nx+md, jb0:jb1) )
    ALLOCATE( Lop(-md:nx+md, jb0:jb1) )
    ALLOCATE( Cop(-md:nx+md, jb0:jb1) )
    ALLOCATE( dc(-md:nx+md, jb0:jb1) )
    ALLOCATE( fw(-md:nx+md, jb0:jb1) )
    ALLOCATE( fwb(-md:nx+md, jb0:jb1) )
    ALLOCATE( uw(-md:nx+md, jb0:jb1) )
    u0 = 0.d0; s1 = 0.d0; s2 = 0.d0; s3 = 0.d0
    r1 = 0.d0; r2 = 0.d0; r3 = 0.d0; r4 = 0.d0
    Lop = 0.d0; Cop = 0.d0; dc = 0.d0; fw = 0.d0; fwb = 0.d0; uw = 0.d0
    IF (.NOT. ALLOCATED(f1d)) THEN
       ALLOCATE( f1d(-md:nmax+md), fh1d(-md:nmax+md), fb1d(-md:nmax+md) )
    END IF
  END SUBROUTINE Alloc

  SUBROUTINE Dealloc()
    IF (ALLOCATED(u0)) THEN
       DEALLOCATE(u0, s1, s2, s3, r1, r2, r3, r4, Lop, Cop, dc, fw, fwb, uw)
    END IF
  END SUBROUTINE Dealloc

  SUBROUTINE SetGrid()
    IF (idim == 1) THEN
       ny = 0
       dy = 1.d0
       jb0 = 0
       jb1 = 0
    ELSE
       dy = (yr - yl)/REAL(ny, 8)
       jb0 = -md
       jb1 = ny + md
    END IF
    dx  = (xr - xl)/REAL(nx, 8)
    cdx = 1.d0/dx
    cdy = 1.d0/dy
  END SUBROUTINE SetGrid

END MODULE para
!***********************************************************************
! funcs.f90 -- flux, source, exact solutions, initial data, norms
!
!   Physical equation :  u_t = ( f(u) )_xx            f(u) = u^m
!                (2D) :  u_t = ( f(u) )_xx + ( f(u) )_yy
!
!   icase = 1 : linear heat equation, m = 1
!               1D : u = exp(-t) sin(x)          on [-pi,pi], periodic
!               2D : u = exp(-2t) sin(x+y)       on [-pi,pi]^2, periodic
!   icase = 2 : manufactured smooth PME solution, m = 2
!               1D : u = 1 + 0.5 exp(-t) sin(x)  on [-pi,pi], periodic,
!                    with the analytic source term Src
!   icase = 3 : Barenblatt (Zeldovich-Kompaneets) solution, m = pmem
!               1D : u = t^-a (1 - a(m-1)x^2/(2m t^2b))^1/(m-1)
!               2D : u = t^-a (1 - a(m-1)r^2/(4m t^2b))^1/(m-1)
!               1/a = d/(d(m-1)+2) -> decay,  b = 1/(d(m-1)+2) -> spreading
!   icase = 4 : two Barenblatt pulses at x = -2.5 and x = +2.5 (1D)
!***********************************************************************
MODULE funcs
  USE comm
  IMPLICIT NONE
CONTAINS

  ! flux  f(u) = u^m   (ipos = 1 : f = max(u,0)^m, positivity-safe)
  ! integer exponents are expanded: the flux evaluation is then cheaper
  ! than the reconstruction in both operators, as it is in practice.
  FUNCTION Flu(x)
    REAL(8) :: Flu, x, y, uu
    IF (ipos == 1) THEN
       y = MAX(x, 0.d0)
    ELSE
       y = x
    END IF
    SELECT CASE (ieq)
    CASE (2)                      ! linear diffusion, g = u
       Flu = y
       RETURN
    CASE (3)                      ! Buckley-Leverett, g' = 4 eps u (1-u)
       uu = MAX(0.d0, MIN(1.d0, x))
       Flu = epsbl*((uu*uu)*2.d0 - (uu*uu*uu)*(4.d0/3.d0))
       RETURN
    CASE (4)                      ! strongly degenerate, nu = 0 for |u| <= 0.25
       IF (x < -0.25d0) THEN
          Flu = epsbl*(x + 0.25d0)
       ELSE IF (x > 0.25d0) THEN
          Flu = epsbl*(x - 0.25d0)
       ELSE
          Flu = 0.d0
       END IF
       RETURN
    END SELECT
    IF (pmem == 1.d0) THEN
       Flu = y
    ELSE IF (pmem == 2.d0) THEN
       Flu = y*y
    ELSE IF (pmem == 3.d0) THEN
       Flu = y*y*y
    ELSE IF (pmem == 4.d0) THEN
       Flu = (y*y)*(y*y)
    ELSE IF (pmem == 5.d0) THEN
       Flu = (y*y)*(y*y)*y
    ELSE IF (pmem == 6.d0) THEN
       Flu = (y*y*y)*(y*y*y)
    ELSE IF (pmem == 8.d0) THEN
       y = y*y; y = y*y; Flu = y*y
    ELSE
       Flu = y**pmem
    END IF
  END FUNCTION Flu

  ! convective flux  b(u) = kconv * u^sconv   (kconv = 0 : no convection)
  FUNCTION FluB(x)
    REAL(8) :: FluB, x, y
    IF (kconv == 0.d0) THEN
       FluB = 0.d0
       RETURN
    END IF
    IF (ipos == 1) THEN
       y = MAX(x, 0.d0)
    ELSE
       y = x
    END IF
    IF (sconv == 1.d0) THEN
       FluB = kconv*y
    ELSE IF (sconv == 2.d0) THEN
       FluB = kconv*y*y
    ELSE
       FluB = kconv*y**sconv
    END IF
  END FUNCTION FluB

  ! derivative of the diffusive flux g'(u), analytically for each family
  FUNCTION DgD(x)
    REAL(8) :: DgD, x, uu, y
    SELECT CASE (ieq)
    CASE (2)
       DgD = 1.d0
       RETURN
    CASE (3)
       uu = MAX(0.d0, MIN(1.d0, x))
       DgD = epsbl*4.d0*uu*(1.d0 - uu)
       RETURN
    CASE (4)
       IF (ABS(x) > 0.25d0) THEN
          DgD = epsbl
       ELSE
          DgD = 0.d0
       END IF
       RETURN
    END SELECT
    IF (ipos == 1) THEN
       y = MAX(x, 0.d0)
    ELSE
       y = x
    END IF
    IF (pmem == 1.d0) THEN
       DgD = 1.d0
    ELSE
       DgD = pmem*y**(pmem - 1.d0)
    END IF
  END FUNCTION DgD

  ! derivative of the convective flux b'(u)
  FUNCTION DbD(x)
    REAL(8) :: DbD, x, uu, dd, ff, df, gg, dg
    SELECT CASE (ieq)
    CASE (3)
       ! analytic derivative of the Buckley-Leverett flux, with and without gravity
       uu = MAX(0.d0, MIN(1.d0, x))
       dd = uu*uu + (1.d0 - uu)**2
       ff = uu*uu/dd
       df = (2.d0*uu*dd - uu*uu*(4.d0*uu - 2.d0))/dd**2
       IF (iblflux == 1) THEN
          gg = 1.d0 - (1.d0 - uu)**2
          dg = 2.d0*(1.d0 - uu)
          df = df*gg + ff*dg
       END IF
       DbD = df
       RETURN
    CASE (4)
       DbD = 2.d0*x
       RETURN
    END SELECT
    IF (kconv == 0.d0) THEN
       DbD = 0.d0
       RETURN
    END IF
    IF (sconv == 1.d0) THEN
       DbD = kconv
    ELSE
       DbD = kconv*sconv*MAX(x, 1.d-12)**(sconv - 1.d0)
    END IF
  END FUNCTION DbD

  ! exponent and coefficient of the Barenblatt solution
  SUBROUTINE BarenPar(d, m, aa, bb, cc)
    REAL(8), INTENT(IN)  :: d, m
    REAL(8), INTENT(OUT) :: aa, bb, cc
    aa = d/(d*(m-1.d0) + 2.d0)
    bb = 1.d0/(d*(m-1.d0) + 2.d0)
    cc = aa*(m-1.d0)/(2.d0*m*d)
  END SUBROUTINE BarenPar

  FUNCTION Ubaren(x, y, t)
    REAL(8) :: Ubaren, x, y, t
    REAL(8) :: aa, bb, cc, s
    IF (idim == 1) THEN
       CALL BarenPar(1.d0, pmem, aa, bb, cc)
       s = 1.d0 - cc*x*x/t**(2.d0*bb)
    ELSE
       CALL BarenPar(2.d0, pmem, aa, bb, cc)
       s = 1.d0 - cc*(x*x + y*y)/t**(2.d0*bb)
    END IF
    IF (s <= 0.d0) THEN
       Ubaren = 0.d0
    ELSE
       Ubaren = t**(-aa)*s**(1.d0/(pmem-1.d0))
    END IF
  END FUNCTION Ubaren

  ! exact solution of the test case
  FUNCTION Uex(x, y, t)
    REAL(8) :: Uex, x, y, t
    IF (icase == 1) THEN
       IF (idim == 1) THEN
          Uex = EXP(-t)*SIN(x)
       ELSE
          Uex = EXP(-2.d0*t)*SIN(x + y)
       END IF
    ELSE IF (icase == 2) THEN
       Uex = 1.d0 + 0.5d0*EXP(-t)*SIN(x)
    ELSE IF (icase == 5) THEN
       ! convected Barenblatt solution: u_t + (u)_x = (u^m)_xx
       Uex = Ubaren(x - kconv*(t - tini), y, t)
    ELSE IF (icase == 6) THEN
       Uex = 1.d0 + 0.5d0*EXP(-t)*SIN(x)
    ELSE IF (icase == 7) THEN
       Uex = 1.d0 + 0.4d0*EXP(-t)*SIN(x)*SIN(y)
    ELSE IF (icase == 8) THEN
       Uex = 1.d0 + 0.4d0*SIN(x)*SIN(y)
    ELSE
       Uex = Ubaren(x, y, t)
    END IF
  END FUNCTION Uex

  ! analytic source term, icase = 2 only (u = 1 + a e^-t sin x, m = 2)
  FUNCTION Src(x, y, t)
    REAL(8) :: Src, x, y, t
    REAL(8) :: a, e, u, ut, ux, uy, uxx
    a = 0.5d0
    e = EXP(-t)
    IF (icase == 8) THEN
       ! stationary 2D solution: u_t = 0, S = -(u^2)_xx - (u^2)_yy
       a = 0.4d0
       u   = 1.d0 + a*SIN(x)*SIN(y)
       ux  =  a*COS(x)*SIN(y)
       uy  =  a*SIN(x)*COS(y)
       uxx = -a*SIN(x)*SIN(y)
       Src = -2.d0*( ux*ux + uy*uy + u*(2.d0*uxx) )
    ELSE IF (icase == 7) THEN
       ! u_t = (u^2)_xx + (u^2)_yy ,  u = 1 + a e^-t sin x sin y
       a = 0.4d0
       u   = 1.d0 + a*e*SIN(x)*SIN(y)
       ut  = -a*e*SIN(x)*SIN(y)
       ux  =  a*e*COS(x)*SIN(y)
       uy  =  a*e*SIN(x)*COS(y)
       uxx = -a*e*SIN(x)*SIN(y)
       Src = ut - 2.d0*( ux*ux + uy*uy + u*(uxx - a*e*SIN(x)*SIN(y)) )
    ELSE
       ! one-dimensional manufactured solutions (icase = 2 and 6)
       u   = 1.d0 + a*e*SIN(x)
       ut  = -a*e*SIN(x)
       ux  =  a*e*COS(x)
       uxx = -a*e*SIN(x)
       IF (icase == 6) THEN
          Src = ut + u*ux - 2.d0*(ux*ux + u*uxx)
       ELSE
          Src = ut - 2.d0*(ux*ux + u*uxx)
       END IF
    END IF
  END FUNCTION Src

  ! computational domain of the test case
  SUBROUTINE Domain()
    REAL(8), PARAMETER :: pi = 3.14159265358979323846d0
    REAL(8) :: tenduser, tiuser, pmuser
    INTEGER :: ibcuser
    tenduser = tend;  tiuser = tini;  ibcuser = ibc;  pmuser = pmem
    IF (icase == 1) THEN
       pmem = 1.d0
       ipos = 0                       ! the solution changes sign: f(u) = u
       xl = -pi; xr = pi
       yl = -pi; yr = pi
       ibc = 1
       tend = 2.d0
       tini = 0.d0
    ELSE IF (icase == 2) THEN
       pmem = 2.d0
       ipos = 0                       ! the manufactured solution is positive
       xl = -pi; xr = pi
       yl = 0.d0; yr = 1.d0
       ibc = 1
       tend = 1.d0
       tini = 0.d0
    ELSE IF (icase == 3) THEN
       IF (idim == 1) THEN
          xl = -6.d0; xr = 6.d0
       ELSE
          xl = -8.d0; xr = 8.d0
       END IF
       yl = xl; yr = xr
       ibc = 2
       tend = 2.d0
       tini = 1.d0
    ELSE IF (icase == 5) THEN
       ! convected Barenblatt: the pulse moves right by kconv*(tend-tini)
       xl = -6.d0; xr = 9.d0
       yl = 0.d0; yr = 1.d0
       ibc = 2
       tend = 2.d0
       tini = 1.d0
       kconv = 1.d0
       sconv = 1.d0
    ELSE IF (icase == 6) THEN
       ! manufactured convection-diffusion, periodic
       xl = -pi; xr = pi
       yl = 0.d0; yr = 1.d0
       ibc = 1
       tend = 1.d0
       tini = 0.d0
       kconv = 0.5d0
       sconv = 2.d0
    ELSE IF (icase == 20) THEN
       ! example 3 of the reference set: two boxes for the PME with m = 6
       pmem = 6.d0
       ieq = 1
       ipos = 1
       kconv = 0.d0
       xl = -6.d0; xr = 6.d0
       yl = 0.d0; yr = 1.d0
       ibc = 1
       tend = 0.8d0
       tini = 0.d0
    ELSE IF (icase == 21) THEN
       ! example 7: 2D PME, two Gaussian caps, T = 0.5
       pmem = 2.d0
       ieq = 1
       ipos = 1
       kconv = 0.d0
       xl = -10.d0; xr = 10.d0
       yl = -10.d0; yr = 10.d0
       ibc = 1
       tend = 0.5d0
       tini = 0.d0
    ELSE IF (icase == 22) THEN
       ! example 8: 2D Buckley-Leverett, unit disc of height one
       ieq = 3
       epsbl = 0.01d0
       iblflux = 1
       ipos = 0
       kconv = 1.d0
       sconv = 1.d0
       xl = -1.5d0; xr = 1.5d0
       yl = -1.5d0; yr = 1.5d0
       ibc = 1
       tend = 0.5d0
       tini = 0.d0
    ELSE IF (icase == 23) THEN
       ! example 9: 2D strongly degenerate, two discs of opposite sign
       ieq = 4
       epsbl = 0.1d0
       ipos = 0
       kconv = 1.d0
       sconv = 2.d0
       xl = -1.5d0; xr = 1.5d0
       yl = -1.5d0; yr = 1.5d0
       ibc = 1
       tend = 0.5d0
       tini = 0.d0
    ELSE IF (icase == 24) THEN
       ! example 5: 1D strongly degenerate, two square pulses
       ieq = 4
       epsbl = 0.1d0
       ipos = 0
       kconv = 1.d0
       sconv = 2.d0
       xl = -2.d0; xr = 2.d0
       yl = 0.d0; yr = 1.d0
       ibc = 1
       tend = 0.7d0
       tini = 0.d0
    ELSE IF (icase == 25) THEN
       ! example 4: Buckley-Leverett initial-boundary value problem on [0,1]
       ieq = 3
       epsbl = 0.01d0
       iblflux = 0
       ipos = 0
       kconv = 1.d0
       sconv = 1.d0
       xl = 0.d0; xr = 1.d0
       yl = 0.d0; yr = 1.d0
       ibc = 3
       ubcl = 1.d0
       ubcr = 0.d0
       tend = 0.2d0
       tini = 0.d0
    ELSE IF (icase == 26) THEN
       ! example 4, Riemann problem
       ieq = 3
       epsbl = 0.01d0
       iblflux = 1
       ipos = 0
       kconv = 1.d0
       sconv = 1.d0
       xl = 0.d0; xr = 1.d0
       yl = 0.d0; yr = 1.d0
       ibc = 3
       ubcl = 0.d0
       ubcr = 1.d0
       tend = 0.2d0
       tini = 0.d0
    ELSE IF (icase == 8) THEN
       ! stationary two-dimensional test, periodic
       pmem = 2.d0
       ipos = 0
       kconv = 0.d0
       xl = -pi; xr = pi
       yl = -pi; yr = pi
       ibc = 1
       tend = 1.d0
       tini = 0.d0
    ELSE IF (icase == 7) THEN
       ! two-dimensional manufactured porous medium solution, periodic
       pmem = 2.d0
       ipos = 0
       kconv = 0.d0
       xl = -pi; xr = pi
       yl = -pi; yr = pi
       ibc = 1
       tend = 1.d0
       tini = 0.d0
    ELSE
       xl = -8.d0; xr = 8.d0
       yl = 0.d0; yr = 1.d0
       ibc = 2
       tend = 1.5d0
       tini = 0.2d0
    END IF
    ! values given in input.nml take precedence when they are positive
    IF (tenduser > 0.d0) tend = tenduser
    IF (tiuser   > 0.d0) tini = tiuser
    IF (ibcuser  > 0)    ibc  = ibcuser
    IF (icase > 2 .AND. pmuser > 0.d0) pmem = pmuser
  END SUBROUTINE Domain

  ! initial condition
  SUBROUTINE Init()
    INTEGER :: i, j
    REAL(8) :: x, y, s
    DO j = jb0, jb1
       DO i = -md, nx+md
          x = xl + REAL(i, 8)*dx
          y = yl + REAL(j, 8)*dy
          IF (icase == 20) THEN
             ! two boxes: 1 on (-3,-1), 2 on (0,3), zero elsewhere
             IF (x > -3.d0 .AND. x < -1.d0) THEN
                u0(i,j) = 1.d0
             ELSE IF (x > 0.d0 .AND. x < 3.d0) THEN
                u0(i,j) = 2.d0
             ELSE
                u0(i,j) = 0.d0
             END IF
          ELSE IF (icase == 21) THEN
             ! two Gaussian caps, exp(1/(r^2-6)) inside r^2 < 6 about (2,-2) and (-2,2)
             IF ((x - 2.d0)**2 + (y + 2.d0)**2 < 6.d0) THEN
                u0(i,j) = EXP(1.d0/((x - 2.d0)**2 + (y + 2.d0)**2 - 6.d0))
             ELSE IF ((x + 2.d0)**2 + (y - 2.d0)**2 < 6.d0) THEN
                u0(i,j) = EXP(1.d0/((x + 2.d0)**2 + (y - 2.d0)**2 - 6.d0))
             ELSE
                u0(i,j) = 0.d0
             END IF
          ELSE IF (icase == 22) THEN
             ! unit disc of height one
             IF (x*x + y*y < 0.5d0) THEN
                u0(i,j) = 1.d0
             ELSE
                u0(i,j) = 0.d0
             END IF
          ELSE IF (icase == 23) THEN
             ! two discs of radius 0.4, opposite signs
             IF ((x + 0.5d0)**2 + (y + 0.5d0)**2 < 0.16d0) THEN
                u0(i,j) = 1.d0
             ELSE IF ((x - 0.5d0)**2 + (y - 0.5d0)**2 < 0.16d0) THEN
                u0(i,j) = -1.d0
             ELSE
                u0(i,j) = 0.d0
             END IF
          ELSE IF (icase == 24) THEN
             ! two square pulses of height +1 and -1
             IF (ABS(x + 0.70710678d0) < 0.4d0) THEN
                u0(i,j) = 1.d0
             ELSE IF (ABS(x - 0.70710678d0) < 0.4d0) THEN
                u0(i,j) = -1.d0
             ELSE
                u0(i,j) = 0.d0
             END IF
          ELSE IF (icase == 25) THEN
             ! 1 - x on [0,1/3], zero to the right
             IF (x <= 1.d0/3.d0) THEN
                u0(i,j) = 1.d0 - x
             ELSE
                u0(i,j) = 0.d0
             END IF
          ELSE IF (icase == 26) THEN
             ! Riemann data: 0 on the left, 1 on the right
             IF (x < 1.d0 - 0.70710678d0) THEN
                u0(i,j) = 0.d0
             ELSE
                u0(i,j) = 1.d0
             END IF
          ELSE IF (icase == 4) THEN
             s = MAX(Ubaren(x - 2.5d0, y, tini), Ubaren(x + 2.5d0, y, tini))
             u0(i,j) = s
          ELSE IF (icase == 3) THEN
             u0(i,j) = Ubaren(x, y, tini)
          ELSE IF (icase == 5) THEN
             u0(i,j) = Ubaren(x, y, tini)
          ELSE
             u0(i,j) = Uex(x, y, tini)
          END IF
       END DO
    END DO
    tnum = tini
  END SUBROUTINE Init

  ! error norms with respect to the exact solution at t = tend
  SUBROUTINE ErrNorm(u, e1, einf)
    REAL(8), INTENT(IN)  :: u(-md:nx+md, jb0:jb1)
    REAL(8), INTENT(OUT) :: e1, einf
    INTEGER :: i, j
    REAL(8) :: x, y, d
    e1 = 0.d0; einf = 0.d0
    DO j = 0, ny
       DO i = 0, nx
          x = xl + REAL(i, 8)*dx
          y = yl + REAL(j, 8)*dy
          d = ABS(u(i,j) - Uex(x, y, tend))
          e1 = e1 + d
          einf = MAX(einf, d)
       END DO
    END DO
    e1 = e1/REAL((nx+1)*(ny+1), 8)
  END SUBROUTINE ErrNorm

  ! parabolic time-step restriction   dt = cfl/aam
  SUBROUTINE Cflc(aamv)
    REAL(8), INTENT(OUT) :: aamv
    INTEGER :: i, j
    REAL(8) :: w, umax, s, umin, uu, cc
    INTEGER :: k
    umax = 0.d0
    DO j = 0, ny
       DO i = 0, nx
          umax = MAX(umax, u0(i,j))
       END DO
    END DO
    IF (ieq == 3 .OR. ieq == 4) THEN
       ! the diffusive coefficient is not a power of u in these families; the
       ! maxima are taken over the admissible range of the case, not over the
       ! current state, so that the time step of the two schemes is the same
       IF (ieq == 3) THEN
          umin = 0.d0
          umax = 1.d0
       ELSE
          umin = -1.d0
          umax = 1.d0
       END IF
       w = 0.d0
       DO k = 0, 10
          uu = umin + (umax - umin)*REAL(k, 8)/10.d0
          w = MAX(w, ABS(DgD(uu)))
       END DO
    ELSE IF (pmem == 1.d0) THEN
       w = 1.d0
    ELSE
       w = pmem*MAX(umax, 0.d0)**(pmem - 1.d0)
    END IF
    s = 1.d0/dx**2
    IF (idim == 2) s = s + 1.d0/dy**2
    aamv = MAX(w*s, 1.d-30)
    ! the convective restriction is added only when a convective flux is present,
    ! so that the diffusion-only cases keep the very same time step as before
    IF (kconv /= 0.d0 .OR. ieq == 3 .OR. ieq == 4) THEN
       IF (ieq == 3) THEN
          umin = 0.d0
          umax = 1.d0
       ELSE IF (ieq == 4) THEN
          umin = -1.d0
          umax = 1.d0
       ELSE
          umin = 1.d30
          umax = -1.d30
          DO j = 0, ny
             DO i = 0, nx
                umin = MIN(umin, u0(i,j))
                umax = MAX(umax, u0(i,j))
             END DO
          END DO
       END IF
       cc = 0.d0
       DO k = 0, 10
          uu = umin + (umax - umin)*REAL(k, 8)/10.d0
          cc = MAX(cc, ABS(DbD(uu)))
       END DO
       aamv = aamv + cc/dx
       IF (idim == 2) aamv = aamv + cc/dy
    END IF
  END SUBROUTINE Cflc

  ! total mass (a conserved quantity of the PME)
  FUNCTION TotMass(u)
    REAL(8) :: TotMass
    REAL(8), INTENT(IN) :: u(-md:nx+md, jb0:jb1)
    INTEGER :: i, j
    TotMass = 0.d0
    DO j = 0, ny
       DO i = 0, nx
          TotMass = TotMass + u(i,j)
       END DO
    END DO
    TotMass = TotMass*dx
    IF (idim == 2) TotMass = TotMass*dy
  END FUNCTION TotMass

END MODULE funcs

!***********************************************************************
! ops.f90 -- the two approximations of ( f(u) )_xx
!
!   kind = 1 : L , the expensive WENO operator: mapped WENO-ZQ
!              reconstruction of the flux at the half points, nonlinear
!              weights, 6th or 8th order (iorder)
!   kind = 2 : C , the cheap component-wise operator: centred second
!              difference of f(u) at the grid points, fixed coefficients,
!              no reconstruction and no nonlinear weights (order nq)
!
!   Both approximate the same quantity, hence
!      d = L - C
!   is O(dx^nq) in smooth regions and carries the local limiting
!   information (the nonlinear weights) near a free boundary.
!
!   fh(i) approximates f_x(x_{i+1/2});  ( fh(i) - fh(i-1) )/dx^2 is the
!   (2k)-th order approximation of (f(u))_xx.  All the numerical
!   constants below are those of the published WENO-ZQ-522 / WENO-ZQ
!   schemes of the RK+WENO reference code (weno6.f90, ghyz.f90, ...).
!***********************************************************************
MODULE ops
  USE comm
  USE funcs
  IMPLICIT NONE
CONTAINS

  !---------------------------------------------------------------------
  !  one operator evaluation: boundary conditions, flux field, operator
  !---------------------------------------------------------------------
  SUBROUTINE EvalOp(kind, u, out, tt)
    INTEGER, INTENT(IN)  :: kind
    REAL(8), INTENT(INOUT) :: u(-md:nx+md, jb0:jb1)
    REAL(8), INTENT(OUT) :: out(-md:nx+md, jb0:jb1)
    REAL(8), INTENT(IN)  :: tt
    INTEGER :: i, j
    REAL(8) :: x, y, db, bmin, bmax

    CALL Bcond(u)
    CALL FluxField(u, fw)
    IF (kconv /= 0.d0) THEN
       alf = 0.d0
       DO j = 0, ny
          DO i = 0, nx
             alf = MAX(alf, ABS(u(i,j)))
          END DO
       END DO
       IF (ieq == 3 .OR. ieq == 4) THEN
          alf = 0.d0
          DO j = 0, ny
             DO i = 0, nx
                alf = MAX(alf, ABS(DbD(u(i,j))))
             END DO
          END DO
       ELSE IF (sconv == 1.d0) THEN
          alf = ABS(kconv)
       ELSE IF (sconv == 2.d0) THEN
          alf = ABS(kconv)*2.d0*alf
       ELSE
          alf = ABS(kconv)*ABS(sconv)*alf**(sconv - 1.d0)
       END IF
       alf = MAX(alf, 1.d-12)
       ! sign of b'(u):  +1 upwind from the left, -1 from the right,
       ! 0 mixed -> Lax-Friedrichs with the global speed
       bmin = 1.d30
       bmax = -1.d30
       DO j = 0, ny
          DO i = 0, nx
             db = DbD(u(i,j))
             bmin = MIN(bmin, db)
             bmax = MAX(bmax, db)
          END DO
       END DO
       isgn = 0
       IF (bmin > 0.d0) isgn = 1
       IF (bmax < 0.d0) isgn = -1
    END IF
    IF (kind == 1 .AND. iops == 0) THEN
       CALL OpLField(out)
    ELSE
       CALL OpCField(out)
    END IF

    IF (nsrc == 1) THEN
       DO j = 0, ny
          DO i = 0, nx
             x = xl + REAL(i, 8)*dx
             y = yl + REAL(j, 8)*dy
             out(i,j) = out(i,j) + Src(x, y, tt)
          END DO
       END DO
    END IF
  END SUBROUTINE EvalOp

  !---------------------------------------------------------------------
  !  boundary conditions on the halo
  !    ibc = 1 : periodic (the point nx duplicates the point 0)
  !    ibc = 2 : zero outside the domain (compact-support problems)
  !    ibc = 3 : transmissive (copy of the nearest interior point)
  !---------------------------------------------------------------------
  SUBROUTINE Bcond(u)
    REAL(8), INTENT(INOUT) :: u(-md:nx+md, jb0:jb1)
    INTEGER :: i, j
    IF (ibc == 1) THEN
       DO j = jb0, jb1
          DO i = 1, md
             u(-i,j)   = u(nx-i,j)
             u(nx+i,j) = u(i,j)
          END DO
       END DO
       IF (idim == 2) THEN
          DO i = -md, nx+md
             DO j = 1, md
                u(i,-j)   = u(i,ny-j)
                u(i,ny+j) = u(i,j)
             END DO
          END DO
       END IF
    ELSE IF (ibc == 3) THEN
       DO j = jb0, jb1
          DO i = 1, md
             u(-i,j)   = ubcl
             u(nx+i,j) = ubcr
          END DO
       END DO
    ELSE IF (ibc == 2) THEN
       DO j = jb0, jb1
          DO i = 1, md
             u(-i,j)   = 0.d0
             u(nx+i,j) = 0.d0
          END DO
       END DO
       IF (idim == 2) THEN
          DO i = -md, nx+md
             DO j = 1, md
                u(i,-j)   = 0.d0
                u(i,ny+j) = 0.d0
             END DO
          END DO
       END IF
    ELSE
       DO j = jb0, jb1
          DO i = 1, md
             u(-i,j)   = u(0,j)
             u(nx+i,j) = u(nx,j)
          END DO
       END DO
       IF (idim == 2) THEN
          DO i = -md, nx+md
             DO j = 1, md
                u(i,-j)   = u(i,0)
                u(i,ny+j) = u(i,ny)
             END DO
          END DO
       END IF
    END IF
  END SUBROUTINE Bcond

  !---------------------------------------------------------------------
  !  flux field f(u)
  !---------------------------------------------------------------------
  SUBROUTINE FluxField(u, f)
    REAL(8), INTENT(IN)  :: u(-md:nx+md, jb0:jb1)
    REAL(8), INTENT(OUT) :: f(-md:nx+md, jb0:jb1)
    INTEGER :: i, j
    DO j = jb0, jb1
       DO i = -md, nx+md
          f(i,j)  = Flu(u(i,j))
          fwb(i,j) = FluB(u(i,j))
          uw(i,j)  = u(i,j)
       END DO
    END DO
  END SUBROUTINE FluxField

  !---------------------------------------------------------------------
  !  L : WENO operator, direction by direction
  !---------------------------------------------------------------------
  SUBROUTINE OpLField(out)
    REAL(8), INTENT(OUT) :: out(-md:nx+md, jb0:jb1)
    INTEGER :: i, j
    DO j = 0, ny
       DO i = -md, nx+md
          f1d(i) = fw(i,j)
       END DO
       CALL Weno1D(nx, f1d, fh1d)
       DO i = 0, nx
          out(i,j) = cdx*cdx*( fh1d(i) - fh1d(i-1) )
       END DO
       IF (kconv /= 0.d0) THEN
          DO i = -md, nx+md
             fb1d(i) = fwb(i,j)
             f1d(i)  = uw(i,j)
          END DO
          CALL Weno5Flux(nx, fb1d, f1d, fh1d, alf, isgn)
          DO i = 0, nx
             out(i,j) = out(i,j) - cdx*( fh1d(i) - fh1d(i-1) )
          END DO
       END IF
    END DO
    IF (idim == 2) THEN
       DO i = 0, nx
          DO j = jb0, jb1
             f1d(j) = fw(i,j)
          END DO
          CALL Weno1D(ny, f1d, fh1d)
          DO j = 0, ny
             out(i,j) = out(i,j) + cdy*cdy*( fh1d(j) - fh1d(j-1) )
          END DO
          IF (kconv /= 0.d0) THEN
             DO j = jb0, jb1
                fb1d(j) = fwb(i,j)
                f1d(j)  = uw(i,j)
             END DO
             CALL Weno5Flux(ny, fb1d, f1d, fh1d, alf, isgn)
             DO j = 0, ny
                out(i,j) = out(i,j) - cdy*( fh1d(j) - fh1d(j-1) )
             END DO
          END IF
       END DO
    END IF
  END SUBROUTINE OpLField

  !---------------------------------------------------------------------
  !  C : cheap centred second difference of f(u), order nq = 2, 4 or 6
  !---------------------------------------------------------------------
  SUBROUTINE OpCField(out)
    REAL(8), INTENT(OUT) :: out(-md:nx+md, jb0:jb1)
    INTEGER :: i, j
    REAL(8) :: ax, ay
    ax = 1.d0/dx**2
    ay = 1.d0/dy**2
    IF (nq == 2) THEN
       DO j = 0, ny
          DO i = 0, nx
             out(i,j) = ax*( fw(i-1,j) - 2.d0*fw(i,j) + fw(i+1,j) )
          END DO
       END DO
    ELSE IF (nq == 4) THEN
       DO j = 0, ny
          DO i = 0, nx
             out(i,j) = ax/12.d0*( -fw(i-2,j) + 16.d0*fw(i-1,j) - 30.d0*fw(i,j) &
                      + 16.d0*fw(i+1,j) - fw(i+2,j) )
          END DO
       END DO
    ELSE IF (nq == 6) THEN
       DO j = 0, ny
          DO i = 0, nx
             out(i,j) = ax/180.d0*( 2.d0*fw(i-3,j) - 27.d0*fw(i-2,j) + 270.d0*fw(i-1,j) &
                      - 490.d0*fw(i,j) + 270.d0*fw(i+1,j) - 27.d0*fw(i+2,j) + 2.d0*fw(i+3,j) )
          END DO
       END DO
    ELSE
       DO j = 0, ny
          DO i = 0, nx
             out(i,j) = ax*( -fw(i-4,j)/560.d0 + 8.d0*fw(i-3,j)/315.d0 - fw(i-2,j)/5.d0 &
                      + 8.d0*fw(i-1,j)/5.d0 - 205.d0*fw(i,j)/72.d0 + 8.d0*fw(i+1,j)/5.d0 &
                      - fw(i+2,j)/5.d0 + 8.d0*fw(i+3,j)/315.d0 - fw(i+4,j)/560.d0 )
          END DO
       END DO
    END IF
    IF (kconv /= 0.d0) THEN
       IF (nq == 2) THEN
          DO j = 0, ny
             DO i = 0, nx
                out(i,j) = out(i,j) - ax*dx*( fwb(i+1,j) - fwb(i-1,j) )/2.d0
             END DO
          END DO
       ELSE IF (nq == 4) THEN
          DO j = 0, ny
             DO i = 0, nx
                out(i,j) = out(i,j) - cdx*( fwb(i-2,j) - 8.d0*fwb(i-1,j) &
                         + 8.d0*fwb(i+1,j) - fwb(i+2,j) )/12.d0
             END DO
          END DO
       ELSE
          DO j = 0, ny
             DO i = 0, nx
                out(i,j) = out(i,j) - cdx*( -fwb(i-3,j) + 9.d0*fwb(i-2,j) - 45.d0*fwb(i-1,j) &
                         + 45.d0*fwb(i+1,j) - 9.d0*fwb(i+2,j) + fwb(i+3,j) )/60.d0
             END DO
          END DO
       END IF
       IF (idim == 2) THEN
          IF (nq == 2) THEN
             DO j = 0, ny
                DO i = 0, nx
                   out(i,j) = out(i,j) - cdy*( fwb(i,j+1) - fwb(i,j-1) )/2.d0
                END DO
             END DO
          ELSE IF (nq == 4) THEN
             DO j = 0, ny
                DO i = 0, nx
                   out(i,j) = out(i,j) - cdy*( fwb(i,j-2) - 8.d0*fwb(i,j-1) &
                            + 8.d0*fwb(i,j+1) - fwb(i,j+2) )/12.d0
                END DO
             END DO
          ELSE
             DO j = 0, ny
                DO i = 0, nx
                   out(i,j) = out(i,j) - cdy*( -fwb(i,j-3) + 9.d0*fwb(i,j-2) - 45.d0*fwb(i,j-1) &
                            + 45.d0*fwb(i,j+1) - 9.d0*fwb(i,j+2) + fwb(i,j+3) )/60.d0
                END DO
             END DO
          END IF
       END IF
    END IF
    IF (idim == 2) THEN
       IF (nq == 2) THEN
          DO j = 0, ny
             DO i = 0, nx
                out(i,j) = out(i,j) + ay*( fw(i,j-1) - 2.d0*fw(i,j) + fw(i,j+1) )
             END DO
          END DO
       ELSE IF (nq == 4) THEN
          DO j = 0, ny
             DO i = 0, nx
                out(i,j) = out(i,j) + ay/12.d0*( -fw(i,j-2) + 16.d0*fw(i,j-1) - 30.d0*fw(i,j) &
                         + 16.d0*fw(i,j+1) - fw(i,j+2) )
             END DO
          END DO
       ELSE IF (nq == 6) THEN
          DO j = 0, ny
             DO i = 0, nx
                out(i,j) = out(i,j) + ay/180.d0*( 2.d0*fw(i,j-3) - 27.d0*fw(i,j-2) &
                         + 270.d0*fw(i,j-1) - 490.d0*fw(i,j) + 270.d0*fw(i,j+1) &
                         - 27.d0*fw(i,j+2) + 2.d0*fw(i,j+3) )
             END DO
          END DO
       ELSE
          DO j = 0, ny
             DO i = 0, nx
                out(i,j) = out(i,j) + ay*( -fw(i,j-4)/560.d0 + 8.d0*fw(i,j-3)/315.d0 &
                         - fw(i,j-2)/5.d0 + 8.d0*fw(i,j-1)/5.d0 - 205.d0*fw(i,j)/72.d0 &
                         + 8.d0*fw(i,j+1)/5.d0 - fw(i,j+2)/5.d0 + 8.d0*fw(i,j+3)/315.d0 &
                         - fw(i,j+4)/560.d0 )
             END DO
          END DO
       END IF
    END IF
  END SUBROUTINE OpCField

  !---------------------------------------------------------------------
  !  WENO reconstruction of the flux at the half points
  !---------------------------------------------------------------------
  SUBROUTINE Weno1D(ns, f, fh)
    INTEGER, INTENT(IN)  :: ns
    REAL(8), INTENT(IN)  :: f(-md:nmax+md)
    REAL(8), INTENT(OUT) :: fh(-md:nmax+md)
    IF (iorder == 6) THEN
       CALL Weno6_1D(ns, f, fh)
    ELSE
       CALL Weno8_1D(ns, f, fh)
    END IF
  END SUBROUTINE Weno1D

  !---------------------------------------------------------------------
  !  sixth order: 6-point stencil, sub-stencils of 4 points
  !---------------------------------------------------------------------
  SUBROUTINE Weno6_1D(ns, f, fh)
    INTEGER, INTENT(IN)  :: ns
    REAL(8), INTENT(IN)  :: f(-md:nmax+md)
    REAL(8), INTENT(OUT) :: fh(-md:nmax+md)
    INTEGER :: i
    REAL(8) :: q1, q2, q3, b1, b2, b3, s1, s2, s3, d1, d2, d3, tt, ep
    REAL(8) :: g1, g2, g3
    ep = 1.d-6
    g1 = -2.d0/15.d0
    g2 = 19.d0/15.d0
    g3 = -2.d0/15.d0
    DO i = -1, ns
       ! candidate values: derivatives of the reconstruction polynomials
       q1 = ( f(i-2) - 3.d0*f(i-1) - 9.d0*f(i) + 11.d0*f(i+1) )/12.d0
       q2 = ( f(i-1) - 15.d0*f(i) + 15.d0*f(i+1) - f(i+2) )/12.d0
       q3 = ( -11.d0*f(i) + 9.d0*f(i+1) + 3.d0*f(i+2) - f(i+3) )/12.d0
       ! WENO-JS smoothness indicators of the three sub-stencils
       b1 = 13.d0/12.d0*( f(i-2) - 3.d0*f(i-1) + 3.d0*f(i) - f(i+1) )**2 &
          + 1.d0/4.d0*( f(i-2) - 5.d0*f(i-1) + 7.d0*f(i) - 3.d0*f(i+1) )**2
       b2 = 13.d0/12.d0*( f(i-1) - 3.d0*f(i) + 3.d0*f(i+1) - f(i+2) )**2 &
          + 1.d0/4.d0*( f(i-1) - f(i) - f(i+1) + f(i+2) )**2
       b3 = 13.d0/12.d0*( f(i) - 3.d0*f(i+1) + 3.d0*f(i+2) - f(i+3) )**2 &
          + 1.d0/4.d0*( -3.d0*f(i) + 7.d0*f(i+1) - 5.d0*f(i+2) + f(i+3) )**2
       ! WENO-Z type nonlinear weights of the positive split weights
       s1 = (1.d0/21.d0)/(ep + b1)**2
       s2 = (19.d0/21.d0)/(ep + b2)**2
       s3 = (1.d0/21.d0)/(ep + b3)**2
       tt = 1.d0/(s1 + s2 + s3)
       s1 = s1*tt; s2 = s2*tt; s3 = s3*tt
       d1 = (4.d0/27.d0)/(ep + b1)**2
       d2 = (19.d0/27.d0)/(ep + b2)**2
       d3 = (4.d0/27.d0)/(ep + b3)**2
       tt = 1.d0/(d1 + d2 + d3)
       d1 = d1*tt; d2 = d2*tt; d3 = d3*tt
       ! split of the negative ideal weights
       s1 = 42.d0/15.d0*s1 - 27.d0/15.d0*d1
       s2 = 42.d0/15.d0*s2 - 27.d0/15.d0*d2
       s3 = 42.d0/15.d0*s3 - 27.d0/15.d0*d3
       ! mapping onto the (negative) ideal weights
       d1 = s1*(g1 + g1**2 - 3.d0*g1*s1 + s1**2)/(g1**2 + s1*(1.d0 - 2.d0*g1))
       d2 = s2*(g2 + g2**2 - 3.d0*g2*s2 + s2**2)/(g2**2 + s2*(1.d0 - 2.d0*g2))
       d3 = s3*(g3 + g3**2 - 3.d0*g3*s3 + s3**2)/(g3**2 + s3*(1.d0 - 2.d0*g3))
       tt = 1.d0/(d1 + d2 + d3)
       s1 = d1*tt; s2 = d2*tt; s3 = d3*tt
       fh(i) = s1*q1 + s2*q2 + s3*q3
    END DO
  END SUBROUTINE Weno6_1D

  !---------------------------------------------------------------------
  !  eighth order: 8-point stencil, sub-stencils of 6 points
  !---------------------------------------------------------------------
  SUBROUTINE Weno8_1D(ns, f, fh)
    INTEGER, INTENT(IN)  :: ns
    REAL(8), INTENT(IN)  :: f(-md:nmax+md)
    REAL(8), INTENT(OUT) :: fh(-md:nmax+md)
    INTEGER :: i
    REAL(8) :: q1, q2, q3, b1, b2, b3, s1, s2, s3, d1, d2, d3, tt, ep
    REAL(8) :: g1, g2, g3
    REAL(8) :: uf3, uf2, uf1, ui, uz1, uz2, uz3, uz4
    ep = 1.d-6
    g1 = -9.d0/56.d0
    g2 = 37.d0/28.d0
    g3 = -9.d0/56.d0
    DO i = -1, ns
       uf3 = f(i-3); uf2 = f(i-2); uf1 = f(i-1); ui  = f(i)
       uz1 = f(i+1); uz2 = f(i+2); uz3 = f(i+3); uz4 = f(i+4)
       q1 = (-2.d0*uf3 + 10.d0*uf2 - 5.d0*uf1 - 205.d0*ui + 215.d0*uz1 - 13.d0*uz2)/180.d0
       q2 = (-2.d0*uf2 + 25.d0*uf1 - 245.d0*ui + 245.d0*uz1 - 25.d0*uz2 + 2.d0*uz3)/180.d0
       q3 = (13.d0*uf1 - 215.d0*ui + 205.d0*uz1 + 5.d0*uz2 - 10.d0*uz3 + 2.d0*uz4)/180.d0
       b1 = (39797.d0/28350.d0)*uf3**2 &
          + (1.d0/907200.d0*(-13781185.d0*uf2 + 29445860.d0*uf1 - 30853430.d0*ui &
          + 15779140.d0*uz1 - 3137393.d0*uz2))*uf3 &
          + (3811447.d0/90720.d0)*uf2**2 &
          + (1.d0/907200.d0*(-166306870.d0*uf1 + 177529300.d0*ui - 92248625.d0*uz1 &
          + 18578440.d0*uz2))*uf2 &
          + (9303521.d0/45360.d0)*uf1**2 &
          + (1.d0/907200.d0*(-407561420.d0*ui + 216748900.d0*uz1 - 44467310.d0*uz2))*uf1 &
          + (11542331.d0/45360.d0)*ui**2 &
          + (1.d0/907200.d0*(-254485870.d0*uz1 + 53678180.d0*uz2))*ui &
          + (229151.d0/2835.d0)*uz1**2 - (6490037.d0/181440.d0)*uz1*uz2 &
          + (1949567.d0/453600.d0)*uz2**2
       b2 = (39797.d0/28350.d0)*uf2**2 &
          + (1.d0/907200.d0*(-12144655.d0*uf1 + 22425980.d0*ui - 20086730.d0*uz1 &
          + 8759260.d0*uz2 - 1500863.d0*uz3))*uf2 &
          + (773023.d0/22680.d0)*uf1**2 &
          + (1.d0/907200.d0*(-119945170.d0*ui + 111450700.d0*uz1 - 49961975.d0*uz2 &
          + 8759260.d0*uz3))*uf1 &
          + (6158981.d0/45360.d0)*ui**2 &
          + (1.d0/907200.d0*(-240204020.d0*uz1 + 111450700.d0*uz2 - 20086730.d0*uz3))*ui &
          + (6158981.d0/45360.d0)*uz1**2 &
          + (1.d0/907200.d0*(-119945170.d0*uz2 + 22425980.d0*uz3))*uz1 &
          + (773023.d0/22680.d0)*uz2**2 - (2428931.d0/181440.d0)*uz3*uz2 &
          + (39797.d0/28350.d0)*uz3**2
       b3 = (1949567.d0/453600.d0)*uf1**2 &
          + (1.d0/907200.d0*(-32450185.d0*ui + 53678180.d0*uz1 - 44467310.d0*uz2 &
          + 18578440.d0*uz3 - 3137393.d0*uz4))*uf1 &
          + (229151.d0/2835.d0)*ui**2 &
          + (1.d0/907200.d0*(-254485870.d0*uz1 + 216748900.d0*uz2 - 92248625.d0*uz3 &
          + 15779140.d0*uz4))*ui &
          + (11542331.d0/45360.d0)*uz1**2 &
          + (1.d0/907200.d0*(-407561420.d0*uz2 + 177529300.d0*uz3 - 30853430.d0*uz4))*uz1 &
          + (9303521.d0/45360.d0)*uz2**2 &
          + (1.d0/907200.d0*(-166306870.d0*uz3 + 29445860.d0*uz4))*uz2 &
          + (3811447.d0/90720.d0)*uz3**2 - (2756237.d0/181440.d0)*uz4*uz3 &
          + (39797.d0/28350.d0)*uz4**2
       s1 = (9.d0/166.d0)/(ep + b1)**2
       s2 = (148.d0/166.d0)/(ep + b2)**2
       s3 = (9.d0/166.d0)/(ep + b3)**2
       tt = 1.d0/(s1 + s2 + s3)
       s1 = s1*tt; s2 = s2*tt; s3 = s3*tt
       d1 = (18.d0/110.d0)/(ep + b1)**2
       d2 = (74.d0/110.d0)/(ep + b2)**2
       d3 = (18.d0/110.d0)/(ep + b3)**2
       tt = 1.d0/(d1 + d2 + d3)
       d1 = d1*tt; d2 = d2*tt; d3 = d3*tt
       s1 = 166.d0/56.d0*s1 - 110.d0/56.d0*d1
       s2 = 166.d0/56.d0*s2 - 110.d0/56.d0*d2
       s3 = 166.d0/56.d0*s3 - 110.d0/56.d0*d3
       d1 = s1*(g1 + g1**2 - 3.d0*g1*s1 + s1**2)/(g1**2 + s1*(1.d0 - 2.d0*g1))
       d2 = s2*(g2 + g2**2 - 3.d0*g2*s2 + s2**2)/(g2**2 + s2*(1.d0 - 2.d0*g2))
       d3 = s3*(g3 + g3**2 - 3.d0*g3*s3 + s3**2)/(g3**2 + s3*(1.d0 - 2.d0*g3))
       tt = 1.d0/(d1 + d2 + d3)
       s1 = d1*tt; s2 = d2*tt; s3 = d3*tt
       fh(i) = s1*q1 + s2*q2 + s3*q3
    END DO
  END SUBROUTINE Weno8_1D


  !---------------------------------------------------------------------
  !  fifth-order WENO (Jiang-Shu) Lax-Friedrichs flux for the convective
  !  term:  F_{i+1/2} = WENO5_left(f+) + WENO5_right(f-),
  !          f+- = ( b(u) +- alf*u )/2 ,   alf = max |b'(u)|
  !---------------------------------------------------------------------
  SUBROUTINE Weno5Flux(ns, f, uu, fh, alpha, isgn)
    INTEGER, INTENT(IN)  :: ns
    REAL(8), INTENT(IN)  :: f(-md:nmax+md), uu(-md:nmax+md)
    REAL(8), INTENT(OUT) :: fh(-md:nmax+md)
    REAL(8), INTENT(IN)  :: alpha
    INTEGER, INTENT(IN)  :: isgn
    INTEGER :: i
    REAL(8) :: fp(-md:nmax+md), fm(-md:nmax+md)
    REAL(8) :: v0, v1, v2, b0, b1, b2, a0, a1, a2, ss, ep
    ep = 1.d-6
    IF (isgn == 0) THEN
       DO i = -md, ns + md
          fp(i) = 0.5d0*( f(i) + alpha*uu(i) )
          fm(i) = 0.5d0*( f(i) - alpha*uu(i) )
       END DO
    ELSE IF (isgn > 0) THEN
       DO i = -md, ns + md
          fp(i) = f(i)
          fm(i) = 0.d0
       END DO
    ELSE
       DO i = -md, ns + md
          fp(i) = 0.d0
          fm(i) = f(i)
       END DO
    END IF
    DO i = -1, ns
       ! ---- left-biased reconstruction of f+ ----
       v0 = ( 2.d0*fp(i-2) - 7.d0*fp(i-1) + 11.d0*fp(i) )/6.d0
       v1 = ( -fp(i-1) + 5.d0*fp(i) + 2.d0*fp(i+1) )/6.d0
       v2 = ( 2.d0*fp(i) + 5.d0*fp(i+1) - fp(i+2) )/6.d0
       b0 = 13.d0/12.d0*( fp(i-2) - 2.d0*fp(i-1) + fp(i) )**2 &
          + 1.d0/4.d0*( fp(i-2) - 4.d0*fp(i-1) + 3.d0*fp(i) )**2
       b1 = 13.d0/12.d0*( fp(i-1) - 2.d0*fp(i) + fp(i+1) )**2 &
          + 1.d0/4.d0*( fp(i-1) - fp(i+1) )**2
       b2 = 13.d0/12.d0*( fp(i) - 2.d0*fp(i+1) + fp(i+2) )**2 &
          + 1.d0/4.d0*( 3.d0*fp(i) - 4.d0*fp(i+1) + fp(i+2) )**2
       a0 = (1.d0/10.d0)/(ep + b0)**2
       a1 = (6.d0/10.d0)/(ep + b1)**2
       a2 = (3.d0/10.d0)/(ep + b2)**2
       ss = 1.d0/(a0 + a1 + a2)
       fh(i) = ss*( a0*v0 + a1*v1 + a2*v2 )
       ! ---- right-biased reconstruction of f- ----
       v0 = ( 2.d0*fm(i+3) - 7.d0*fm(i+2) + 11.d0*fm(i+1) )/6.d0
       v1 = ( -fm(i+2) + 5.d0*fm(i+1) + 2.d0*fm(i) )/6.d0
       v2 = ( 2.d0*fm(i+1) + 5.d0*fm(i) - fm(i-1) )/6.d0
       b0 = 13.d0/12.d0*( fm(i+3) - 2.d0*fm(i+2) + fm(i+1) )**2 &
          + 1.d0/4.d0*( fm(i+3) - 4.d0*fm(i+2) + 3.d0*fm(i+1) )**2
       b1 = 13.d0/12.d0*( fm(i+2) - 2.d0*fm(i+1) + fm(i) )**2 &
          + 1.d0/4.d0*( fm(i+2) - fm(i) )**2
       b2 = 13.d0/12.d0*( fm(i+1) - 2.d0*fm(i) + fm(i-1) )**2 &
          + 1.d0/4.d0*( 3.d0*fm(i+1) - 4.d0*fm(i) + fm(i-1) )**2
       a0 = (3.d0/10.d0)/(ep + b0)**2
       a1 = (6.d0/10.d0)/(ep + b1)**2
       a2 = (1.d0/10.d0)/(ep + b2)**2
       ss = 1.d0/(a0 + a1 + a2)
       fh(i) = fh(i) + ss*( a0*v0 + a1*v1 + a2*v2 )
    END DO
  END SUBROUTINE Weno5Flux

END MODULE ops


!***********************************************************************
! timeint.f90 -- time integration
!
!   ischeme = 3  : SSP-RK3            (baseline, 3 expensive evaluations)
!   ischeme = 4  : classical RK4      (baseline, 4 expensive evaluations)
!   ischeme = 13 : ILW3 -- one-step improved Lax-Wendroff scheme on the
!                  SSP-RK3 tableau: ONE evaluation of the expensive WENO
!                  operator L per step and three of the cheap one C
!   ischeme = 14 : ILW4 -- the same construction on the classical RK4
!                  tableau: one L and four C per step
!
!   Generic form (any s-stage tableau a, b), with G = C:
!
!     U^{n+1} = U^n + dt [ L(U^n) + sum_{j>=2} b_j C(V_j)
!                          - ( sum_{j>=2} b_j ) C(U^n) ]
!     V_j     = U^n + dt [ sum_{l<j} a_jl C(V_l) + sigma_j d ]
!     sigma_j = sum_{l<j} a_jl ,   d = L(U^n) - C(U^n)   (frozen)
!
!   Because d = L - C,  C(V_l) + d  is the frozen approximation of
!   L(V_l): every stage increment has the leading-order term sigma_j L
!   with sigma_j > 0 (the sign condition), and the weights of the cheap
!   derivative in the final update sum to zero (the zero-sum condition).
!   When C = L (d = 0) the scheme reduces exactly to the RK method.
!***********************************************************************
MODULE timeint
  USE comm
  USE ops
  IMPLICIT NONE
CONTAINS

  SUBROUTINE Step()
    INTEGER :: j
    REAL(8) :: sg, bs

    CALL EvalOp(1, u0, Lop, tnum)          ! the single expensive evaluation

    IF (ischeme > 10) THEN                 ! ---------------- ILW ----------------
       CALL EvalOp(2, u0, Cop, tnum)       ! cheap evaluation at U^n
       dc = Lop - Cop                      ! frozen correction d = L - C

       ! stage 2 :  V_2 = U^n + dt sigma_2 ( C(U^n) + d ) = U^n + dt sigma_2 L(U^n)
       sg = amat(2,1)
       s1 = u0 + dt*sg*( Cop + dc )
       CALL EvalOp(2, s1, r1, tnum + cvec(2)*dt)     ! C(V_2)

       IF (mt > 2) THEN
          sg = amat(3,1) + amat(3,2)
          s2 = u0 + dt*( amat(3,1)*Cop + amat(3,2)*r1 + sg*dc )
          CALL EvalOp(2, s2, r2, tnum + cvec(3)*dt)  ! C(V_3)
       END IF
       IF (mt > 3) THEN
          sg = amat(4,1) + amat(4,2) + amat(4,3)
          s3 = u0 + dt*( amat(4,1)*Cop + amat(4,2)*r1 + amat(4,3)*r2 + sg*dc )
          CALL EvalOp(2, s3, r3, tnum + cvec(4)*dt)  ! C(V_4)
       END IF

       bs = 0.d0
       DO j = 1, mt-1
          bs = bs + bvec(j+1)
       END DO
       u0 = u0 + dt*( Lop - bs*Cop )
       u0 = u0 + dt*( bvec(2)*r1 + bvec(3)*r2 + bvec(4)*r3 )
    ELSE                                   ! ---------------- RK ---------------
       r1 = Lop                            ! k_1
       s1 = u0 + dt*amat(2,1)*r1
       CALL EvalOp(1, s1, r2, tnum + cvec(2)*dt)        ! k_2
       IF (mt > 2) THEN
          s2 = u0 + dt*( amat(3,1)*r1 + amat(3,2)*r2 )
          CALL EvalOp(1, s2, r3, tnum + cvec(3)*dt)     ! k_3
       END IF
       IF (mt > 3) THEN
          s3 = u0 + dt*( amat(4,1)*r1 + amat(4,2)*r2 + amat(4,3)*r3 )
          CALL EvalOp(1, s3, r4, tnum + cvec(4)*dt)     ! k_4
       END IF
       u0 = u0 + dt*( bvec(1)*r1 + bvec(2)*r2 + bvec(3)*r3 + bvec(4)*r4 )
    END IF
  END SUBROUTINE Step

END MODULE timeint

!***********************************************************************
! main.f90 -- driver
!
!   Reads input.nml, runs either a single computation (nlev = 1) or a
!   convergence study on a sequence of doubled grids, and writes
!
!      results.csv   (one line per grid: errors, rates, cpu time)
!      profile.plt   (numerical and exact profiles)
!***********************************************************************
PROGRAM ILWWENO
  USE comm
  USE funcs
  USE ops
  USE para
  USE timeint
  IMPLICIT NONE
  INTEGER :: lev, i, j, nstep
  REAL(8) :: e1, einf, rate1, rateinf, cpu0, cpu1, tstart, umin, mass0, mass1
  REAL(8) :: e1p, einfp, x, y

  CALL Canshu()
  ! ---- Example 6.5, ILW3-WENO6 ----
  idim = 1
  iorder = 6
  nq = 4
  ischeme = 13
  icase = 24
  nlev = 1
  nx0 = 400
  cfl = 0.4d0
  tend = 0.7d0

  CALL Domain()
  CALL SetScheme()
  IF (nsrc < 0) THEN                    ! not set in the namelist: decide from icase
     IF (icase == 2 .OR. icase == 6 .OR. icase == 7 .OR. icase == 8) THEN
        nsrc = 1
     ELSE
        nsrc = 0
     END IF
  END IF
  IF (ibc == 0) ibc = 1                  ! safety: a default is always set by Domain
  IF (tend < 0.d0) tend = 2.d0

  WRITE(*,'(a)') '--------------------------------------------------------------'
  WRITE(*,'(a,i2,a)')   ' space dimension      : ', idim, 'D'
  WRITE(*,'(a,f6.2)')   ' equation             : u_t = (u^m)_xx, m = ', pmem
  WRITE(*,'(a,i2)')     ' WENO order           : ', iorder
  WRITE(*,'(a,i2)')     ' cheap operator order : ', nq
  WRITE(*,'(a,i3)')     ' scheme               : ', ischeme
  WRITE(*,'(a,i2)')     ' test case            : ', icase
  WRITE(*,'(a,f8.3)')   ' CFL                  : ', cfl
  WRITE(*,'(a,f8.3)')   ' final time           : ', tend
  WRITE(*,'(a)') '--------------------------------------------------------------'

  CALL CPU_TIME(tstart)
  e1p = 0.d0; einfp = 0.d0
  OPEN(UNIT=20, FILE='results.csv', STATUS='REPLACE', ACTION='WRITE')
  WRITE(20,'(a)') 'level,nx,ny,e1,einf,r1,rinf,nsteps,cpu,umin,mass'

  DO lev = 1, nlev
     IF (idim == 1) THEN
        nx = nx0*2**(lev-1)
        ny = 0
     ELSE
        nx = nx0*2**(lev-1)
        ny = ny0*2**(lev-1)
     END IF
     CALL SetGrid()
     CALL Alloc()
     CALL Init()
     mass0 = TotMass(u0)

     CALL CPU_TIME(cpu0)
     nstep = 0
     DO WHILE (tnum < tend - 1.d-14)
        CALL Cflc(aam)
        dt = cfl/aam
        IF (tnum + dt >= tend) dt = tend - tnum
        CALL Step()
        tnum = tnum + dt
        nstep = nstep + 1
        umin = MINVAL(u0(0:nx,0:ny))
        IF (umin /= umin .OR. ABS(umin) > 1.d+10) THEN
           WRITE(*,*) ' *** divergence at step ', nstep, ' t = ', tnum
           EXIT
        END IF
     END DO
     CALL CPU_TIME(cpu1)

     CALL ErrNorm(u0, e1, einf)
     umin  = MINVAL(u0(0:nx,0:ny))
     mass1 = TotMass(u0)
     IF (lev > 1 .AND. e1 > 0.d0) THEN
        rate1   = LOG(e1p/e1)/LOG(2.d0)
        rateinf = LOG(einfp/einf)/LOG(2.d0)
     ELSE
        rate1 = 0.d0; rateinf = 0.d0
     END IF
     WRITE(*,'(a,i3,a,i6,a,i6,a,es12.4,a,es12.4)') ' level ', lev, '  nx ', nx, &
          '  ny ', ny, '  L1 ', e1, '  Linf ', einf
     IF (lev > 1) THEN
        WRITE(*,'(a,f7.2,a,f7.2,a,i8,a,f9.3,a)') '        rate L1 ', rate1, &
             '  rate Linf ', rateinf, '   steps ', nstep, '   cpu ', cpu1 - cpu0, ' s'
     ELSE
        WRITE(*,'(a,i8,a,f9.3,a)') '        steps ', nstep, '   cpu ', cpu1 - cpu0, ' s'
     END IF
     WRITE(*,'(a,es12.4,a,es12.4,a,es12.4)') '        min u ', umin, '   mass0 ', mass0, &
          '   mass1 ', mass1
     WRITE(20,'(i3,a,i6,a,i6,a,es14.6,a,es14.6,a,f8.3,a,f8.3,a,i8,a,f12.4,a,es12.4,a,es14.7)') &
          lev, ',', nx, ',', ny, ',', e1, ',', einf, ',', rate1, ',', rateinf, ',', nstep, ',', &
          cpu1 - cpu0, ',', umin, ',', mass1
     e1p = e1; einfp = einf

     ! ---------------- profiles of the last run ----------------
     IF (lev == nlev) THEN
        OPEN(UNIT=21, FILE='profile.plt', STATUS='REPLACE', ACTION='WRITE')
        WRITE(21,*) 'TITLE="ILW-WENO"'
        IF (idim == 1) THEN
           WRITE(21,*) 'VARIABLES="x" "u_num" "u_exact"'
           DO i = 0, nx
              x = xl + REAL(i, 8)*dx
              WRITE(21,'(3es16.8)') x, u0(i,0), Uex(x, 0.d0, tend)
           END DO
        ELSE
           WRITE(21,*) 'VARIABLES="x" "y" "u_num" "u_exact"'
           WRITE(21,*) 'ZONE I=', nx+1, ' J=', ny+1
           DO j = 0, ny
              DO i = 0, nx
                 x = xl + REAL(i, 8)*dx
                 y = yl + REAL(j, 8)*dy
                 WRITE(21,'(es16.8,a,es16.8,a,es16.8,a,es16.8)') x, ',', y, ',', &
                      u0(i,j), ',', Uex(x, y, tend)
              END DO
           END DO
        END IF
        CLOSE(21)
     END IF

     CALL Dealloc()
  END DO

  CLOSE(20)
  CALL CPU_TIME(cpu1)
  WRITE(*,'(a,f10.3,a)') ' total cpu ', cpu1 - tstart, ' s'

END PROGRAM ILWWENO



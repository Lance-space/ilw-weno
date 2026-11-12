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

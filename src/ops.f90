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


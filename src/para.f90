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
    NAMELIST /para_nml/ idim, iorder, nq, icase, ischeme, nlev, nx0, ny0, &
                        cfl, tend, tini, pmem, ipos, ibc, kconv, sconv, iops, nsrc, ieq, epsbl, iblflux, ubcl, ubcr

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

    OPEN(UNIT=15, FILE='input.nml', STATUS='OLD', IOSTAT=ios)
    IF (ios == 0) THEN
       READ(15, NML=para_nml, IOSTAT=ios)
       IF (ios /= 0) WRITE(*,*) 'Canshu: error in input.nml, defaults kept'
       CLOSE(15)
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
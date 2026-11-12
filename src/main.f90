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



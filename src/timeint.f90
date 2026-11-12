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

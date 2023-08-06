*> \brief \b SORHR_COL01
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*  Definition:
*  ===========
*
*       SUBROUTINE SORHR_COL01( M, N, MB1, NB1, NB2, RESULT )
*
*       .. Scalar Arguments ..
*       INTEGER           M, N, MB1, NB1, NB2
*       .. Return values ..
*       REAL              RESULT(6)
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> SORHR_COL01 tests SORGTSQR and SORHR_COL using SLATSQR, SGEMQRT.
*> Therefore, SLATSQR (part of SGEQR), SGEMQRT (part of SGEMQR)
*> have to be tested before this test.
*>
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] M
*> \verbatim
*>          M is INTEGER
*>          Number of rows in test matrix.
*> \endverbatim
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          Number of columns in test matrix.
*> \endverbatim
*> \param[in] MB1
*> \verbatim
*>          MB1 is INTEGER
*>          Number of row in row block in an input test matrix.
*> \endverbatim
*>
*> \param[in] NB1
*> \verbatim
*>          NB1 is INTEGER
*>          Number of columns in column block an input test matrix.
*> \endverbatim
*>
*> \param[in] NB2
*> \verbatim
*>          NB2 is INTEGER
*>          Number of columns in column block in an output test matrix.
*> \endverbatim
*>
*> \param[out] RESULT
*> \verbatim
*>          RESULT is REAL array, dimension (6)
*>          Results of each of the six tests below.
*>
*>            A is a m-by-n test input matrix to be factored.
*>            so that A = Q_gr * ( R )
*>                               ( 0 ),
*>
*>            Q_qr is an implicit m-by-m orthogonal Q matrix, the result
*>            of factorization in blocked WY-representation,
*>            stored in SGEQRT output format.
*>
*>            R is a n-by-n upper-triangular matrix,
*>
*>            0 is a (m-n)-by-n zero matrix,
*>
*>            Q is an explicit m-by-m orthogonal matrix Q = Q_gr * I
*>
*>            C is an m-by-n random matrix,
*>
*>            D is an n-by-m random matrix.
*>
*>          The six tests are:
*>
*>          RESULT(1) = |R - (Q**H) * A| / ( eps * m * |A| )
*>            is equivalent to test for | A - Q * R | / (eps * m * |A|),
*>
*>          RESULT(2) = |I - (Q**H) * Q| / ( eps * m ),
*>
*>          RESULT(3) = | Q_qr * C - Q * C | / (eps * m * |C|),
*>
*>          RESULT(4) = | (Q_gr**H) * C - (Q**H) * C | / (eps * m * |C|)
*>
*>          RESULT(5) = | D * Q_qr - D * Q | / (eps * m * |D|)
*>
*>          RESULT(6) = | D * (Q_qr**H) - D * (Q**H) | / (eps * m * |D|),
*>
*>          where:
*>            Q_qr * C, (Q_gr**H) * C, D * Q_qr, D * (Q_qr**H) are
*>            computed using SGEMQRT,
*>
*>            Q * C, (Q**H) * C, D * Q, D * (Q**H)  are
*>            computed using SGEMM.
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee
*> \author Univ. of California Berkeley
*> \author Univ. of Colorado Denver
*> \author NAG Ltd.
*
*> \ingroup single_lin
*
*  =====================================================================
      SUBROUTINE SORHR_COL01( M, N, MB1, NB1, NB2, RESULT )
      IMPLICIT NONE
*
*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      INTEGER           M, N, MB1, NB1, NB2
*     .. Return values ..
      REAL              RESULT(6)
*
*  =====================================================================
*
*     ..
*     .. Local allocatable arrays
      REAL            , ALLOCATABLE ::  A(:,:), AF(:,:), Q(:,:), R(:,:),
     $                   RWORK(:), WORK( : ), T1(:,:), T2(:,:), DIAG(:),
     $                   C(:,:), CF(:,:), D(:,:), DF(:,:)
*
*     .. Parameters ..
      REAL               ONE, ZERO
      PARAMETER          ( ZERO = 0.0E+0, ONE = 1.0E+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            TESTZEROS
      INTEGER            INFO, I, J, K, L, LWORK, NB1_UB, NB2_UB, NRB
      REAL               ANORM, EPS, RESID, CNORM, DNORM
*     ..
*     .. Local Arrays ..
      INTEGER            ISEED( 4 )
      REAL               WORKQUERY( 1 )
*     ..
*     .. External Functions ..
      REAL               SLAMCH, SLANGE, SLANSY
      EXTERNAL           SLAMCH, SLANGE, SLANSY
*     ..
*     .. External Subroutines ..
      EXTERNAL           SLACPY, SLARNV, SLASET, SLATSQR, SORHR_COL,
     $                   SORGTSQR, SSCAL, SGEMM, SGEMQRT, SSYRK
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          CEILING, REAL, MAX, MIN
*     ..
*     .. Scalars in Common ..
      CHARACTER(LEN=32)  SRNAMT
*     ..
*     .. Common blocks ..
      COMMON             / SRMNAMC / SRNAMT
*     ..
*     .. Data statements ..
      DATA ISEED / 1988, 1989, 1990, 1991 /
*
*     TEST MATRICES WITH HALF OF MATRIX BEING ZEROS
*
      TESTZEROS = .FALSE.
*
      EPS = SLAMCH( 'Epsilon' )
      K = MIN( M, N )
      L = MAX( M, N, 1)
*
*     Dynamically allocate local arrays
*
      ALLOCATE ( A(M,N), AF(M,N), Q(L,L), R(M,L), RWORK(L),
     $           C(M,N), CF(M,N),
     $           D(N,M), DF(N,M) )
*
*     Put random numbers into A and copy to AF
*
      DO J = 1, N
         CALL SLARNV( 2, ISEED, M, A( 1, J ) )
      END DO
      IF( TESTZEROS ) THEN
         IF( M.GE.4 ) THEN
            DO J = 1, N
               CALL SLARNV( 2, ISEED, M/2, A( M/4, J ) )
            END DO
         END IF
      END IF
      CALL SLACPY( 'Full', M, N, A, M, AF, M )
*
*     Number of row blocks in SLATSQR
*
      NRB = MAX( 1, CEILING( REAL( M - N ) / REAL( MB1 - N ) ) )
*
      ALLOCATE ( T1( NB1, N * NRB ) )
      ALLOCATE ( T2( NB2, N ) )
      ALLOCATE ( DIAG( N ) )
*
*     Begin determine LWORK for the array WORK and allocate memory.
*
*     SLATSQR requires NB1 to be bounded by N.
*
      NB1_UB = MIN( NB1, N)
*
*     SGEMQRT requires NB2 to be bounded by N.
*
      NB2_UB = MIN( NB2, N)
*
      CALL SLATSQR( M, N, MB1, NB1_UB, AF, M, T1, NB1,
     $              WORKQUERY, -1, INFO )
      LWORK = INT( WORKQUERY( 1 ) )
      CALL SORGTSQR( M, N, MB1, NB1, AF, M, T1, NB1, WORKQUERY, -1,
     $               INFO )

      LWORK = MAX( LWORK, INT( WORKQUERY( 1 ) ) )
*
*     In SGEMQRT, WORK is N*NB2_UB if SIDE = 'L',
*                or  M*NB2_UB if SIDE = 'R'.
*
      LWORK = MAX( LWORK, NB2_UB * N, NB2_UB * M )
*
      ALLOCATE ( WORK( LWORK ) )
*
*     End allocate memory for WORK.
*
*
*     Begin Householder reconstruction routines
*
*     Factor the matrix A in the array AF.
*
      SRNAMT = 'SLATSQR'
      CALL SLATSQR( M, N, MB1, NB1_UB, AF, M, T1, NB1, WORK, LWORK,
     $              INFO )
*
*     Copy the factor R into the array R.
*
      SRNAMT = 'SLACPY'
      CALL SLACPY( 'U', N, N, AF, M, R, M )
*
*     Reconstruct the orthogonal matrix Q.
*
      SRNAMT = 'SORGTSQR'
      CALL SORGTSQR( M, N, MB1, NB1, AF, M, T1, NB1, WORK, LWORK,
     $               INFO )
*
*     Perform the Householder reconstruction, the result is stored
*     the arrays AF and T2.
*
      SRNAMT = 'SORHR_COL'
      CALL SORHR_COL( M, N, NB2, AF, M, T2, NB2, DIAG, INFO )
*
*     Compute the factor R_hr corresponding to the Householder
*     reconstructed Q_hr and place it in the upper triangle of AF to
*     match the Q storage format in SGEQRT. R_hr = R_tsqr * S,
*     this means changing the sign of I-th row of the matrix R_tsqr
*     according to sign of of I-th diagonal element DIAG(I) of the
*     matrix S.
*
      SRNAMT = 'SLACPY'
      CALL SLACPY( 'U', N, N, R, M, AF, M )
*
      DO I = 1, N
         IF( DIAG( I ).EQ.-ONE ) THEN
            CALL SSCAL( N+1-I, -ONE, AF( I, I ), M )
         END IF
      END DO
*
*     End Householder reconstruction routines.
*
*
*     Generate the m-by-m matrix Q
*
      CALL SLASET( 'Full', M, M, ZERO, ONE, Q, M )
*
      SRNAMT = 'SGEMQRT'
      CALL SGEMQRT( 'L', 'N', M, M, K, NB2_UB, AF, M, T2, NB2, Q, M,
     $              WORK, INFO )
*
*     Copy R
*
      CALL SLASET( 'Full', M, N, ZERO, ZERO, R, M )
*
      CALL SLACPY( 'Upper', M, N, AF, M, R, M )
*
*     TEST 1
*     Compute |R - (Q**T)*A| / ( eps * m * |A| ) and store in RESULT(1)
*
      CALL SGEMM( 'T', 'N', M, N, M, -ONE, Q, M, A, M, ONE, R, M )
*
      ANORM = SLANGE( '1', M, N, A, M, RWORK )
      RESID = SLANGE( '1', M, N, R, M, RWORK )
      IF( ANORM.GT.ZERO ) THEN
         RESULT( 1 ) = RESID / ( EPS * MAX( 1, M ) * ANORM )
      ELSE
         RESULT( 1 ) = ZERO
      END IF
*
*     TEST 2
*     Compute |I - (Q**T)*Q| / ( eps * m ) and store in RESULT(2)
*
      CALL SLASET( 'Full', M, M, ZERO, ONE, R, M )
      CALL SSYRK( 'U', 'T', M, M, -ONE, Q, M, ONE, R, M )
      RESID = SLANSY( '1', 'Upper', M, R, M, RWORK )
      RESULT( 2 ) = RESID / ( EPS * MAX( 1, M ) )
*
*     Generate random m-by-n matrix C
*
      DO J = 1, N
         CALL SLARNV( 2, ISEED, M, C( 1, J ) )
      END DO
      CNORM = SLANGE( '1', M, N, C, M, RWORK )
      CALL SLACPY( 'Full', M, N, C, M, CF, M )
*
*     Apply Q to C as Q*C = CF
*
      SRNAMT = 'SGEMQRT'
      CALL SGEMQRT( 'L', 'N', M, N, K, NB2_UB, AF, M, T2, NB2, CF, M,
     $               WORK, INFO )
*
*     TEST 3
*     Compute |CF - Q*C| / ( eps *  m * |C| )
*
      CALL SGEMM( 'N', 'N', M, N, M, -ONE, Q, M, C, M, ONE, CF, M )
      RESID = SLANGE( '1', M, N, CF, M, RWORK )
      IF( CNORM.GT.ZERO ) THEN
         RESULT( 3 ) = RESID / ( EPS * MAX( 1, M ) * CNORM )
      ELSE
         RESULT( 3 ) = ZERO
      END IF
*
*     Copy C into CF again
*
      CALL SLACPY( 'Full', M, N, C, M, CF, M )
*
*     Apply Q to C as (Q**T)*C = CF
*
      SRNAMT = 'SGEMQRT'
      CALL SGEMQRT( 'L', 'T', M, N, K, NB2_UB, AF, M, T2, NB2, CF, M,
     $               WORK, INFO )
*
*     TEST 4
*     Compute |CF - (Q**T)*C| / ( eps * m * |C|)
*
      CALL SGEMM( 'T', 'N', M, N, M, -ONE, Q, M, C, M, ONE, CF, M )
      RESID = SLANGE( '1', M, N, CF, M, RWORK )
      IF( CNORM.GT.ZERO ) THEN
         RESULT( 4 ) = RESID / ( EPS * MAX( 1, M ) * CNORM )
      ELSE
         RESULT( 4 ) = ZERO
      END IF
*
*     Generate random n-by-m matrix D and a copy DF
*
      DO J = 1, M
         CALL SLARNV( 2, ISEED, N, D( 1, J ) )
      END DO
      DNORM = SLANGE( '1', N, M, D, N, RWORK )
      CALL SLACPY( 'Full', N, M, D, N, DF, N )
*
*     Apply Q to D as D*Q = DF
*
      SRNAMT = 'SGEMQRT'
      CALL SGEMQRT( 'R', 'N', N, M, K, NB2_UB, AF, M, T2, NB2, DF, N,
     $               WORK, INFO )
*
*     TEST 5
*     Compute |DF - D*Q| / ( eps * m * |D| )
*
      CALL SGEMM( 'N', 'N', N, M, M, -ONE, D, N, Q, M, ONE, DF, N )
      RESID = SLANGE( '1', N, M, DF, N, RWORK )
      IF( DNORM.GT.ZERO ) THEN
         RESULT( 5 ) = RESID / ( EPS * MAX( 1, M ) * DNORM )
      ELSE
         RESULT( 5 ) = ZERO
      END IF
*
*     Copy D into DF again
*
      CALL SLACPY( 'Full', N, M, D, N, DF, N )
*
*     Apply Q to D as D*QT = DF
*
      SRNAMT = 'SGEMQRT'
      CALL SGEMQRT( 'R', 'T', N, M, K, NB2_UB, AF, M, T2, NB2, DF, N,
     $               WORK, INFO )
*
*     TEST 6
*     Compute |DF - D*(Q**T)| / ( eps * m * |D| )
*
      CALL SGEMM( 'N', 'T', N, M, M, -ONE, D, N, Q, M, ONE, DF, N )
      RESID = SLANGE( '1', N, M, DF, N, RWORK )
      IF( DNORM.GT.ZERO ) THEN
         RESULT( 6 ) = RESID / ( EPS * MAX( 1, M ) * DNORM )
      ELSE
         RESULT( 6 ) = ZERO
      END IF
*
*     Deallocate all arrays
*
      DEALLOCATE ( A, AF, Q, R, RWORK, WORK, T1, T2, DIAG,
     $             C, D, CF, DF )
*
      RETURN
*
*     End of SORHR_COL01
*
      END

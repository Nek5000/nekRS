*> \brief \b CCKCSD
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*  Definition:
*  ===========
*
*       SUBROUTINE CCKCSD( NM, MVAL, PVAL, QVAL, NMATS, ISEED, THRESH,
*                          MMAX, X, XF, U1, U2, V1T, V2T, THETA, IWORK,
*                          WORK, RWORK, NIN, NOUT, INFO )
*
*       .. Scalar Arguments ..
*       INTEGER            INFO, NIN, NM, NMATS, MMAX, NOUT
*       REAL               THRESH
*       ..
*       .. Array Arguments ..
*       INTEGER            ISEED( 4 ), IWORK( * ), MVAL( * ), PVAL( * ),
*      $                   QVAL( * )
*       REAL               RWORK( * ), THETA( * )
*       COMPLEX            U1( * ), U2( * ), V1T( * ), V2T( * ),
*      $                   WORK( * ), X( * ), XF( * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> CCKCSD tests CUNCSD:
*>        the CSD for an M-by-M unitary matrix X partitioned as
*>        [ X11 X12; X21 X22 ]. X11 is P-by-Q.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] NM
*> \verbatim
*>          NM is INTEGER
*>          The number of values of M contained in the vector MVAL.
*> \endverbatim
*>
*> \param[in] MVAL
*> \verbatim
*>          MVAL is INTEGER array, dimension (NM)
*>          The values of the matrix row dimension M.
*> \endverbatim
*>
*> \param[in] PVAL
*> \verbatim
*>          PVAL is INTEGER array, dimension (NM)
*>          The values of the matrix row dimension P.
*> \endverbatim
*>
*> \param[in] QVAL
*> \verbatim
*>          QVAL is INTEGER array, dimension (NM)
*>          The values of the matrix column dimension Q.
*> \endverbatim
*>
*> \param[in] NMATS
*> \verbatim
*>          NMATS is INTEGER
*>          The number of matrix types to be tested for each combination
*>          of matrix dimensions.  If NMATS >= NTYPES (the maximum
*>          number of matrix types), then all the different types are
*>          generated for testing.  If NMATS < NTYPES, another input line
*>          is read to get the numbers of the matrix types to be used.
*> \endverbatim
*>
*> \param[in,out] ISEED
*> \verbatim
*>          ISEED is INTEGER array, dimension (4)
*>          On entry, the seed of the random number generator.  The array
*>          elements should be between 0 and 4095, otherwise they will be
*>          reduced mod 4096, and ISEED(4) must be odd.
*>          On exit, the next seed in the random number sequence after
*>          all the test matrices have been generated.
*> \endverbatim
*>
*> \param[in] THRESH
*> \verbatim
*>          THRESH is REAL
*>          The threshold value for the test ratios.  A result is
*>          included in the output file if RESULT >= THRESH.  To have
*>          every test ratio printed, use THRESH = 0.
*> \endverbatim
*>
*> \param[in] MMAX
*> \verbatim
*>          MMAX is INTEGER
*>          The maximum value permitted for M, used in dimensioning the
*>          work arrays.
*> \endverbatim
*>
*> \param[out] X
*> \verbatim
*>          X is COMPLEX array, dimension (MMAX*MMAX)
*> \endverbatim
*>
*> \param[out] XF
*> \verbatim
*>          XF is COMPLEX array, dimension (MMAX*MMAX)
*> \endverbatim
*>
*> \param[out] U1
*> \verbatim
*>          U1 is COMPLEX array, dimension (MMAX*MMAX)
*> \endverbatim
*>
*> \param[out] U2
*> \verbatim
*>          U2 is COMPLEX array, dimension (MMAX*MMAX)
*> \endverbatim
*>
*> \param[out] V1T
*> \verbatim
*>          V1T is COMPLEX array, dimension (MMAX*MMAX)
*> \endverbatim
*>
*> \param[out] V2T
*> \verbatim
*>          V2T is COMPLEX array, dimension (MMAX*MMAX)
*> \endverbatim
*>
*> \param[out] THETA
*> \verbatim
*>          THETA is REAL array, dimension (MMAX)
*> \endverbatim
*>
*> \param[out] IWORK
*> \verbatim
*>          IWORK is INTEGER array, dimension (MMAX)
*> \endverbatim
*>
*> \param[out] WORK
*> \verbatim
*>          WORK is COMPLEX array
*> \endverbatim
*>
*> \param[out] RWORK
*> \verbatim
*>          RWORK is REAL array
*> \endverbatim
*>
*> \param[in] NIN
*> \verbatim
*>          NIN is INTEGER
*>          The unit number for input.
*> \endverbatim
*>
*> \param[in] NOUT
*> \verbatim
*>          NOUT is INTEGER
*>          The unit number for output.
*> \endverbatim
*>
*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>          = 0 :  successful exit
*>          > 0 :  If CLAROR returns an error code, the absolute value
*>                 of it is returned.
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
*> \ingroup complex_eig
*
*  =====================================================================
      SUBROUTINE CCKCSD( NM, MVAL, PVAL, QVAL, NMATS, ISEED, THRESH,
     $                   MMAX, X, XF, U1, U2, V1T, V2T, THETA, IWORK,
     $                   WORK, RWORK, NIN, NOUT, INFO )
*
*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      INTEGER            INFO, NIN, NM, NMATS, MMAX, NOUT
      REAL               THRESH
*     ..
*     .. Array Arguments ..
      INTEGER            ISEED( 4 ), IWORK( * ), MVAL( * ), PVAL( * ),
     $                   QVAL( * )
      REAL               RWORK( * ), THETA( * )
      COMPLEX            U1( * ), U2( * ), V1T( * ), V2T( * ),
     $                   WORK( * ), X( * ), XF( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      INTEGER            NTESTS
      PARAMETER          ( NTESTS = 15 )
      INTEGER            NTYPES
      PARAMETER          ( NTYPES = 4 )
      REAL               GAPDIGIT, ORTH, REALONE, REALZERO, TEN
      PARAMETER          ( GAPDIGIT = 10.0E0, ORTH = 1.0E-4,
     $                     REALONE = 1.0E0, REALZERO = 0.0E0,
     $                     TEN = 10.0E0 )
      COMPLEX            ONE, ZERO
      PARAMETER          ( ONE = (1.0E0,0.0E0), ZERO = (0.0E0,0.0E0) )
      REAL               PIOVER2
      PARAMETER ( PIOVER2 = 1.57079632679489661923132169163975144210E0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            FIRSTT
      CHARACTER*3        PATH
      INTEGER            I, IINFO, IM, IMAT, J, LDU1, LDU2, LDV1T,
     $                   LDV2T, LDX, LWORK, M, NFAIL, NRUN, NT, P, Q, R
*     ..
*     .. Local Arrays ..
      LOGICAL            DOTYPE( NTYPES )
      REAL               RESULT( NTESTS )
*     ..
*     .. External Subroutines ..
      EXTERNAL           ALAHDG, ALAREQ, ALASUM, CCSDTS, CLACSG, CLAROR,
     $                   CLASET, CSROT
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MIN
*     ..
*     .. External Functions ..
      REAL               SLARAN, SLARND
      EXTERNAL           SLARAN, SLARND
*     ..
*     .. Executable Statements ..
*
*     Initialize constants and the random number seed.
*
      PATH( 1: 3 ) = 'CSD'
      INFO = 0
      NRUN = 0
      NFAIL = 0
      FIRSTT = .TRUE.
      CALL ALAREQ( PATH, NMATS, DOTYPE, NTYPES, NIN, NOUT )
      LDX = MMAX
      LDU1 = MMAX
      LDU2 = MMAX
      LDV1T = MMAX
      LDV2T = MMAX
      LWORK = MMAX*MMAX
*
*     Do for each value of M in MVAL.
*
      DO 30 IM = 1, NM
         M = MVAL( IM )
         P = PVAL( IM )
         Q = QVAL( IM )
*
         DO 20 IMAT = 1, NTYPES
*
*           Do the tests only if DOTYPE( IMAT ) is true.
*
            IF( .NOT.DOTYPE( IMAT ) )
     $         GO TO 20
*
*           Generate X
*
            IF( IMAT.EQ.1 ) THEN
               CALL CLAROR( 'L', 'I', M, M, X, LDX, ISEED, WORK, IINFO )
               IF( M .NE. 0 .AND. IINFO .NE. 0 ) THEN
                  WRITE( NOUT, FMT = 9999 ) M, IINFO
                  INFO = ABS( IINFO )
                  GO TO 20
               END IF
            ELSE IF( IMAT.EQ.2 ) THEN
               R = MIN( P, M-P, Q, M-Q )
               DO I = 1, R
                  THETA(I) = PIOVER2 * SLARND( 1, ISEED )
               END DO
               CALL CLACSG( M, P, Q, THETA, ISEED, X, LDX, WORK )
               DO I = 1, M
                  DO J = 1, M
                     X(I+(J-1)*LDX) = X(I+(J-1)*LDX) +
     $                                ORTH*SLARND(2,ISEED)
                  END DO
               END DO
            ELSE IF( IMAT.EQ.3 ) THEN
               R = MIN( P, M-P, Q, M-Q )
               DO I = 1, R+1
                  THETA(I) = TEN**(-SLARND(1,ISEED)*GAPDIGIT)
               END DO
               DO I = 2, R+1
                  THETA(I) = THETA(I-1) + THETA(I)
               END DO
               DO I = 1, R
                  THETA(I) = PIOVER2 * THETA(I) / THETA(R+1)
               END DO
               CALL CLACSG( M, P, Q, THETA, ISEED, X, LDX, WORK )
            ELSE
               CALL CLASET( 'F', M, M, ZERO, ONE, X, LDX )
               DO I = 1, M
                  J = INT( SLARAN( ISEED ) * M ) + 1
                  IF( J .NE. I ) THEN
                     CALL CSROT( M, X(1+(I-1)*LDX), 1, X(1+(J-1)*LDX),
     $                 1, REALZERO, REALONE )
                  END IF
               END DO
            END IF
*
            NT = 15
*
            CALL CCSDTS( M, P, Q, X, XF, LDX, U1, LDU1, U2, LDU2, V1T,
     $                   LDV1T, V2T, LDV2T, THETA, IWORK, WORK, LWORK,
     $                   RWORK, RESULT )
*
*           Print information about the tests that did not
*           pass the threshold.
*
            DO 10 I = 1, NT
               IF( RESULT( I ).GE.THRESH ) THEN
                  IF( NFAIL.EQ.0 .AND. FIRSTT ) THEN
                     FIRSTT = .FALSE.
                     CALL ALAHDG( NOUT, PATH )
                  END IF
                  WRITE( NOUT, FMT = 9998 )M, P, Q, IMAT, I,
     $               RESULT( I )
                  NFAIL = NFAIL + 1
               END IF
   10       CONTINUE
            NRUN = NRUN + NT
   20    CONTINUE
   30 CONTINUE
*
*     Print a summary of the results.
*
      CALL ALASUM( PATH, NOUT, NFAIL, NRUN, 0 )
*
 9999 FORMAT( ' CLAROR in CCKCSD: M = ', I5, ', INFO = ', I15 )
 9998 FORMAT( ' M=', I4, ' P=', I4, ', Q=', I4, ', type ', I2,
     $      ', test ', I2, ', ratio=', G13.6 )
      RETURN
*
*     End of CCKCSD
*
      END
*
*
*
      SUBROUTINE CLACSG( M, P, Q, THETA, ISEED, X, LDX, WORK )
      IMPLICIT NONE
*
      INTEGER            LDX, M, P, Q
      INTEGER            ISEED( 4 )
      REAL               THETA( * )
      COMPLEX            WORK( * ), X( LDX, * )
*
      COMPLEX            ONE, ZERO
      PARAMETER          ( ONE = (1.0E0,0.0E0), ZERO = (0.0E0,0.0E0) )
*
      INTEGER            I, INFO, R
*
      R = MIN( P, M-P, Q, M-Q )
*
      CALL CLASET( 'Full', M, M, ZERO, ZERO, X, LDX )
*
      DO I = 1, MIN(P,Q)-R
         X(I,I) = ONE
      END DO
      DO I = 1, R
         X(MIN(P,Q)-R+I,MIN(P,Q)-R+I) = CMPLX ( COS(THETA(I)), 0.0E0 )
      END DO
      DO I = 1, MIN(P,M-Q)-R
         X(P-I+1,M-I+1) = -ONE
      END DO
      DO I = 1, R
         X(P-(MIN(P,M-Q)-R)+1-I,M-(MIN(P,M-Q)-R)+1-I) =
     $      CMPLX( -SIN(THETA(R-I+1)), 0.0E0 )
      END DO
      DO I = 1, MIN(M-P,Q)-R
         X(M-I+1,Q-I+1) = ONE
      END DO
      DO I = 1, R
         X(M-(MIN(M-P,Q)-R)+1-I,Q-(MIN(M-P,Q)-R)+1-I) =
     $      CMPLX( SIN(THETA(R-I+1)), 0.0E0 )
      END DO
      DO I = 1, MIN(M-P,M-Q)-R
         X(P+I,Q+I) = ONE
      END DO
      DO I = 1, R
         X(P+(MIN(M-P,M-Q)-R)+I,Q+(MIN(M-P,M-Q)-R)+I) =
     $      CMPLX( COS(THETA(I)), 0.0E0 )
      END DO
      CALL CLAROR( 'Left', 'No init', P, M, X, LDX, ISEED, WORK, INFO )
      CALL CLAROR( 'Left', 'No init', M-P, M, X(P+1,1), LDX,
     $             ISEED, WORK, INFO )
      CALL CLAROR( 'Right', 'No init', M, Q, X, LDX, ISEED,
     $             WORK, INFO )
      CALL CLAROR( 'Right', 'No init', M, M-Q,
     $             X(1,Q+1), LDX, ISEED, WORK, INFO )
*
      END


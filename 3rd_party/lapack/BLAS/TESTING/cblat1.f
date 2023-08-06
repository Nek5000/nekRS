*> \brief \b CBLAT1
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*  Definition:
*  ===========
*
*       PROGRAM CBLAT1
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*>    Test program for the COMPLEX Level 1 BLAS.
*>    Based upon the original BLAS test routine together with:
*>
*>    F06GAF Example Program Text
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
*> \ingroup complex_blas_testing
*
*  =====================================================================
      PROGRAM CBLAT1
*
*  -- Reference BLAS test routine --
*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*  =====================================================================
*
*     .. Parameters ..
      INTEGER          NOUT
      PARAMETER        (NOUT=6)
*     .. Scalars in Common ..
      INTEGER          ICASE, INCX, INCY, MODE, N
      LOGICAL          PASS
*     .. Local Scalars ..
      REAL             SFAC
      INTEGER          IC
*     .. External Subroutines ..
      EXTERNAL         CHECK1, CHECK2, HEADER
*     .. Common blocks ..
      COMMON           /COMBLA/ICASE, N, INCX, INCY, MODE, PASS
*     .. Data statements ..
      DATA             SFAC/9.765625E-4/
*     .. Executable Statements ..
      WRITE (NOUT,99999)
      DO 20 IC = 1, 10
         ICASE = IC
         CALL HEADER
*
*        Initialize PASS, INCX, INCY, and MODE for a new case.
*        The value 9999 for INCX, INCY or MODE will appear in the
*        detailed  output, if any, for cases that do not involve
*        these parameters.
*
         PASS = .TRUE.
         INCX = 9999
         INCY = 9999
         MODE = 9999
         IF (ICASE.LE.5) THEN
            CALL CHECK2(SFAC)
         ELSE IF (ICASE.GE.6) THEN
            CALL CHECK1(SFAC)
         END IF
*        -- Print
         IF (PASS) WRITE (NOUT,99998)
   20 CONTINUE
      STOP
*
99999 FORMAT (' Complex BLAS Test Program Results',/1X)
99998 FORMAT ('                                    ----- PASS -----')
*
*     End of CBLAT1
*
      END
      SUBROUTINE HEADER
*     .. Parameters ..
      INTEGER          NOUT
      PARAMETER        (NOUT=6)
*     .. Scalars in Common ..
      INTEGER          ICASE, INCX, INCY, MODE, N
      LOGICAL          PASS
*     .. Local Arrays ..
      CHARACTER*6      L(10)
*     .. Common blocks ..
      COMMON           /COMBLA/ICASE, N, INCX, INCY, MODE, PASS
*     .. Data statements ..
      DATA             L(1)/'CDOTC '/
      DATA             L(2)/'CDOTU '/
      DATA             L(3)/'CAXPY '/
      DATA             L(4)/'CCOPY '/
      DATA             L(5)/'CSWAP '/
      DATA             L(6)/'SCNRM2'/
      DATA             L(7)/'SCASUM'/
      DATA             L(8)/'CSCAL '/
      DATA             L(9)/'CSSCAL'/
      DATA             L(10)/'ICAMAX'/
*     .. Executable Statements ..
      WRITE (NOUT,99999) ICASE, L(ICASE)
      RETURN
*
99999 FORMAT (/' Test of subprogram number',I3,12X,A6)
*
*     End of HEADER
*
      END
      SUBROUTINE CHECK1(SFAC)
*     .. Parameters ..
      INTEGER           NOUT
      PARAMETER         (NOUT=6)
*     .. Scalar Arguments ..
      REAL              SFAC
*     .. Scalars in Common ..
      INTEGER           ICASE, INCX, INCY, MODE, N
      LOGICAL           PASS
*     .. Local Scalars ..
      COMPLEX           CA
      REAL              SA
      INTEGER           I, IX, J, LEN, NP1
*     .. Local Arrays ..
      COMPLEX           CTRUE5(8,5,2), CTRUE6(8,5,2), CV(8,5,2), CVR(8),
     +                  CX(8), CXR(15), MWPCS(5), MWPCT(5)
      REAL              STRUE2(5), STRUE4(5)
      INTEGER           ITRUE3(5), ITRUEC(5)
*     .. External Functions ..
      REAL              SCASUM, SCNRM2
      INTEGER           ICAMAX
      EXTERNAL          SCASUM, SCNRM2, ICAMAX
*     .. External Subroutines ..
      EXTERNAL          CSCAL, CSSCAL, CTEST, ITEST1, STEST1
*     .. Intrinsic Functions ..
      INTRINSIC         MAX
*     .. Common blocks ..
      COMMON            /COMBLA/ICASE, N, INCX, INCY, MODE, PASS
*     .. Data statements ..
      DATA              SA, CA/0.3E0, (0.4E0,-0.7E0)/
      DATA              ((CV(I,J,1),I=1,8),J=1,5)/(0.1E0,0.1E0),
     +                  (1.0E0,2.0E0), (1.0E0,2.0E0), (1.0E0,2.0E0),
     +                  (1.0E0,2.0E0), (1.0E0,2.0E0), (1.0E0,2.0E0),
     +                  (1.0E0,2.0E0), (0.3E0,-0.4E0), (3.0E0,4.0E0),
     +                  (3.0E0,4.0E0), (3.0E0,4.0E0), (3.0E0,4.0E0),
     +                  (3.0E0,4.0E0), (3.0E0,4.0E0), (3.0E0,4.0E0),
     +                  (0.1E0,-0.3E0), (0.5E0,-0.1E0), (5.0E0,6.0E0),
     +                  (5.0E0,6.0E0), (5.0E0,6.0E0), (5.0E0,6.0E0),
     +                  (5.0E0,6.0E0), (5.0E0,6.0E0), (0.1E0,0.1E0),
     +                  (-0.6E0,0.1E0), (0.1E0,-0.3E0), (7.0E0,8.0E0),
     +                  (7.0E0,8.0E0), (7.0E0,8.0E0), (7.0E0,8.0E0),
     +                  (7.0E0,8.0E0), (0.3E0,0.1E0), (0.5E0,0.0E0),
     +                  (0.0E0,0.5E0), (0.0E0,0.2E0), (2.0E0,3.0E0),
     +                  (2.0E0,3.0E0), (2.0E0,3.0E0), (2.0E0,3.0E0)/
      DATA              ((CV(I,J,2),I=1,8),J=1,5)/(0.1E0,0.1E0),
     +                  (4.0E0,5.0E0), (4.0E0,5.0E0), (4.0E0,5.0E0),
     +                  (4.0E0,5.0E0), (4.0E0,5.0E0), (4.0E0,5.0E0),
     +                  (4.0E0,5.0E0), (0.3E0,-0.4E0), (6.0E0,7.0E0),
     +                  (6.0E0,7.0E0), (6.0E0,7.0E0), (6.0E0,7.0E0),
     +                  (6.0E0,7.0E0), (6.0E0,7.0E0), (6.0E0,7.0E0),
     +                  (0.1E0,-0.3E0), (8.0E0,9.0E0), (0.5E0,-0.1E0),
     +                  (2.0E0,5.0E0), (2.0E0,5.0E0), (2.0E0,5.0E0),
     +                  (2.0E0,5.0E0), (2.0E0,5.0E0), (0.1E0,0.1E0),
     +                  (3.0E0,6.0E0), (-0.6E0,0.1E0), (4.0E0,7.0E0),
     +                  (0.1E0,-0.3E0), (7.0E0,2.0E0), (7.0E0,2.0E0),
     +                  (7.0E0,2.0E0), (0.3E0,0.1E0), (5.0E0,8.0E0),
     +                  (0.5E0,0.0E0), (6.0E0,9.0E0), (0.0E0,0.5E0),
     +                  (8.0E0,3.0E0), (0.0E0,0.2E0), (9.0E0,4.0E0)/
      DATA              CVR/(8.0E0,8.0E0), (-7.0E0,-7.0E0),
     +                  (9.0E0,9.0E0), (5.0E0,5.0E0), (9.0E0,9.0E0),
     +                  (8.0E0,8.0E0), (7.0E0,7.0E0), (7.0E0,7.0E0)/
      DATA              STRUE2/0.0E0, 0.5E0, 0.6E0, 0.7E0, 0.8E0/
      DATA              STRUE4/0.0E0, 0.7E0, 1.0E0, 1.3E0, 1.6E0/
      DATA              ((CTRUE5(I,J,1),I=1,8),J=1,5)/(0.1E0,0.1E0),
     +                  (1.0E0,2.0E0), (1.0E0,2.0E0), (1.0E0,2.0E0),
     +                  (1.0E0,2.0E0), (1.0E0,2.0E0), (1.0E0,2.0E0),
     +                  (1.0E0,2.0E0), (-0.16E0,-0.37E0), (3.0E0,4.0E0),
     +                  (3.0E0,4.0E0), (3.0E0,4.0E0), (3.0E0,4.0E0),
     +                  (3.0E0,4.0E0), (3.0E0,4.0E0), (3.0E0,4.0E0),
     +                  (-0.17E0,-0.19E0), (0.13E0,-0.39E0),
     +                  (5.0E0,6.0E0), (5.0E0,6.0E0), (5.0E0,6.0E0),
     +                  (5.0E0,6.0E0), (5.0E0,6.0E0), (5.0E0,6.0E0),
     +                  (0.11E0,-0.03E0), (-0.17E0,0.46E0),
     +                  (-0.17E0,-0.19E0), (7.0E0,8.0E0), (7.0E0,8.0E0),
     +                  (7.0E0,8.0E0), (7.0E0,8.0E0), (7.0E0,8.0E0),
     +                  (0.19E0,-0.17E0), (0.20E0,-0.35E0),
     +                  (0.35E0,0.20E0), (0.14E0,0.08E0),
     +                  (2.0E0,3.0E0), (2.0E0,3.0E0), (2.0E0,3.0E0),
     +                  (2.0E0,3.0E0)/
      DATA              ((CTRUE5(I,J,2),I=1,8),J=1,5)/(0.1E0,0.1E0),
     +                  (4.0E0,5.0E0), (4.0E0,5.0E0), (4.0E0,5.0E0),
     +                  (4.0E0,5.0E0), (4.0E0,5.0E0), (4.0E0,5.0E0),
     +                  (4.0E0,5.0E0), (-0.16E0,-0.37E0), (6.0E0,7.0E0),
     +                  (6.0E0,7.0E0), (6.0E0,7.0E0), (6.0E0,7.0E0),
     +                  (6.0E0,7.0E0), (6.0E0,7.0E0), (6.0E0,7.0E0),
     +                  (-0.17E0,-0.19E0), (8.0E0,9.0E0),
     +                  (0.13E0,-0.39E0), (2.0E0,5.0E0), (2.0E0,5.0E0),
     +                  (2.0E0,5.0E0), (2.0E0,5.0E0), (2.0E0,5.0E0),
     +                  (0.11E0,-0.03E0), (3.0E0,6.0E0),
     +                  (-0.17E0,0.46E0), (4.0E0,7.0E0),
     +                  (-0.17E0,-0.19E0), (7.0E0,2.0E0), (7.0E0,2.0E0),
     +                  (7.0E0,2.0E0), (0.19E0,-0.17E0), (5.0E0,8.0E0),
     +                  (0.20E0,-0.35E0), (6.0E0,9.0E0),
     +                  (0.35E0,0.20E0), (8.0E0,3.0E0),
     +                  (0.14E0,0.08E0), (9.0E0,4.0E0)/
      DATA              ((CTRUE6(I,J,1),I=1,8),J=1,5)/(0.1E0,0.1E0),
     +                  (1.0E0,2.0E0), (1.0E0,2.0E0), (1.0E0,2.0E0),
     +                  (1.0E0,2.0E0), (1.0E0,2.0E0), (1.0E0,2.0E0),
     +                  (1.0E0,2.0E0), (0.09E0,-0.12E0), (3.0E0,4.0E0),
     +                  (3.0E0,4.0E0), (3.0E0,4.0E0), (3.0E0,4.0E0),
     +                  (3.0E0,4.0E0), (3.0E0,4.0E0), (3.0E0,4.0E0),
     +                  (0.03E0,-0.09E0), (0.15E0,-0.03E0),
     +                  (5.0E0,6.0E0), (5.0E0,6.0E0), (5.0E0,6.0E0),
     +                  (5.0E0,6.0E0), (5.0E0,6.0E0), (5.0E0,6.0E0),
     +                  (0.03E0,0.03E0), (-0.18E0,0.03E0),
     +                  (0.03E0,-0.09E0), (7.0E0,8.0E0), (7.0E0,8.0E0),
     +                  (7.0E0,8.0E0), (7.0E0,8.0E0), (7.0E0,8.0E0),
     +                  (0.09E0,0.03E0), (0.15E0,0.00E0),
     +                  (0.00E0,0.15E0), (0.00E0,0.06E0), (2.0E0,3.0E0),
     +                  (2.0E0,3.0E0), (2.0E0,3.0E0), (2.0E0,3.0E0)/
      DATA              ((CTRUE6(I,J,2),I=1,8),J=1,5)/(0.1E0,0.1E0),
     +                  (4.0E0,5.0E0), (4.0E0,5.0E0), (4.0E0,5.0E0),
     +                  (4.0E0,5.0E0), (4.0E0,5.0E0), (4.0E0,5.0E0),
     +                  (4.0E0,5.0E0), (0.09E0,-0.12E0), (6.0E0,7.0E0),
     +                  (6.0E0,7.0E0), (6.0E0,7.0E0), (6.0E0,7.0E0),
     +                  (6.0E0,7.0E0), (6.0E0,7.0E0), (6.0E0,7.0E0),
     +                  (0.03E0,-0.09E0), (8.0E0,9.0E0),
     +                  (0.15E0,-0.03E0), (2.0E0,5.0E0), (2.0E0,5.0E0),
     +                  (2.0E0,5.0E0), (2.0E0,5.0E0), (2.0E0,5.0E0),
     +                  (0.03E0,0.03E0), (3.0E0,6.0E0),
     +                  (-0.18E0,0.03E0), (4.0E0,7.0E0),
     +                  (0.03E0,-0.09E0), (7.0E0,2.0E0), (7.0E0,2.0E0),
     +                  (7.0E0,2.0E0), (0.09E0,0.03E0), (5.0E0,8.0E0),
     +                  (0.15E0,0.00E0), (6.0E0,9.0E0), (0.00E0,0.15E0),
     +                  (8.0E0,3.0E0), (0.00E0,0.06E0), (9.0E0,4.0E0)/
      DATA              ITRUE3/0, 1, 2, 2, 2/
      DATA              ITRUEC/0, 1, 1, 1, 1/
*     .. Executable Statements ..
      DO 60 INCX = 1, 2
         DO 40 NP1 = 1, 5
            N = NP1 - 1
            LEN = 2*MAX(N,1)
*           .. Set vector arguments ..
            DO 20 I = 1, LEN
               CX(I) = CV(I,NP1,INCX)
   20       CONTINUE
            IF (ICASE.EQ.6) THEN
*              .. SCNRM2 ..
               CALL STEST1(SCNRM2(N,CX,INCX),STRUE2(NP1),STRUE2(NP1),
     +                     SFAC)
            ELSE IF (ICASE.EQ.7) THEN
*              .. SCASUM ..
               CALL STEST1(SCASUM(N,CX,INCX),STRUE4(NP1),STRUE4(NP1),
     +                     SFAC)
            ELSE IF (ICASE.EQ.8) THEN
*              .. CSCAL ..
               CALL CSCAL(N,CA,CX,INCX)
               CALL CTEST(LEN,CX,CTRUE5(1,NP1,INCX),CTRUE5(1,NP1,INCX),
     +                    SFAC)
            ELSE IF (ICASE.EQ.9) THEN
*              .. CSSCAL ..
               CALL CSSCAL(N,SA,CX,INCX)
               CALL CTEST(LEN,CX,CTRUE6(1,NP1,INCX),CTRUE6(1,NP1,INCX),
     +                    SFAC)
            ELSE IF (ICASE.EQ.10) THEN
*              .. ICAMAX ..
               CALL ITEST1(ICAMAX(N,CX,INCX),ITRUE3(NP1))
               DO 160 I = 1, LEN
                  CX(I) = (42.0E0,43.0E0)
  160          CONTINUE
               CALL ITEST1(ICAMAX(N,CX,INCX),ITRUEC(NP1))
            ELSE
               WRITE (NOUT,*) ' Shouldn''t be here in CHECK1'
               STOP
            END IF
*
   40    CONTINUE
         IF (ICASE.EQ.10) THEN
            N = 8
            IX = 1
            DO 180 I = 1, N
               CXR(IX) = CVR(I)
               IX = IX + INCX
  180       CONTINUE
            CALL ITEST1(ICAMAX(N,CXR,INCX),3)
         END IF
   60 CONTINUE
*
      INCX = 1
      IF (ICASE.EQ.8) THEN
*        CSCAL
*        Add a test for alpha equal to zero.
         CA = (0.0E0,0.0E0)
         DO 80 I = 1, 5
            MWPCT(I) = (0.0E0,0.0E0)
            MWPCS(I) = (1.0E0,1.0E0)
   80    CONTINUE
         CALL CSCAL(5,CA,CX,INCX)
         CALL CTEST(5,CX,MWPCT,MWPCS,SFAC)
      ELSE IF (ICASE.EQ.9) THEN
*        CSSCAL
*        Add a test for alpha equal to zero.
         SA = 0.0E0
         DO 100 I = 1, 5
            MWPCT(I) = (0.0E0,0.0E0)
            MWPCS(I) = (1.0E0,1.0E0)
  100    CONTINUE
         CALL CSSCAL(5,SA,CX,INCX)
         CALL CTEST(5,CX,MWPCT,MWPCS,SFAC)
*        Add a test for alpha equal to one.
         SA = 1.0E0
         DO 120 I = 1, 5
            MWPCT(I) = CX(I)
            MWPCS(I) = CX(I)
  120    CONTINUE
         CALL CSSCAL(5,SA,CX,INCX)
         CALL CTEST(5,CX,MWPCT,MWPCS,SFAC)
*        Add a test for alpha equal to minus one.
         SA = -1.0E0
         DO 140 I = 1, 5
            MWPCT(I) = -CX(I)
            MWPCS(I) = -CX(I)
  140    CONTINUE
         CALL CSSCAL(5,SA,CX,INCX)
         CALL CTEST(5,CX,MWPCT,MWPCS,SFAC)
      END IF
      RETURN
*
*     End of CHECK1
*
      END
      SUBROUTINE CHECK2(SFAC)
*     .. Parameters ..
      INTEGER           NOUT
      PARAMETER         (NOUT=6)
*     .. Scalar Arguments ..
      REAL              SFAC
*     .. Scalars in Common ..
      INTEGER           ICASE, INCX, INCY, MODE, N
      LOGICAL           PASS
*     .. Local Scalars ..
      COMPLEX           CA
      INTEGER           I, J, KI, KN, KSIZE, LENX, LENY, LINCX, LINCY,
     +                  MX, MY
*     .. Local Arrays ..
      COMPLEX           CDOT(1), CSIZE1(4), CSIZE2(7,2), CSIZE3(14),
     +                  CT10X(7,4,4), CT10Y(7,4,4), CT6(4,4), CT7(4,4),
     +                  CT8(7,4,4), CTY0(1), CX(7), CX0(1), CX1(7),
     +                  CY(7), CY0(1), CY1(7)
      INTEGER           INCXS(4), INCYS(4), LENS(4,2), NS(4)
*     .. External Functions ..
      COMPLEX           CDOTC, CDOTU
      EXTERNAL          CDOTC, CDOTU
*     .. External Subroutines ..
      EXTERNAL          CAXPY, CCOPY, CSWAP, CTEST
*     .. Intrinsic Functions ..
      INTRINSIC         ABS, MIN
*     .. Common blocks ..
      COMMON            /COMBLA/ICASE, N, INCX, INCY, MODE, PASS
*     .. Data statements ..
      DATA              CA/(0.4E0,-0.7E0)/
      DATA              INCXS/1, 2, -2, -1/
      DATA              INCYS/1, -2, 1, -2/
      DATA              LENS/1, 1, 2, 4, 1, 1, 3, 7/
      DATA              NS/0, 1, 2, 4/
      DATA              CX1/(0.7E0,-0.8E0), (-0.4E0,-0.7E0),
     +                  (-0.1E0,-0.9E0), (0.2E0,-0.8E0),
     +                  (-0.9E0,-0.4E0), (0.1E0,0.4E0), (-0.6E0,0.6E0)/
      DATA              CY1/(0.6E0,-0.6E0), (-0.9E0,0.5E0),
     +                  (0.7E0,-0.6E0), (0.1E0,-0.5E0), (-0.1E0,-0.2E0),
     +                  (-0.5E0,-0.3E0), (0.8E0,-0.7E0)/
      DATA              ((CT8(I,J,1),I=1,7),J=1,4)/(0.6E0,-0.6E0),
     +                  (0.0E0,0.0E0), (0.0E0,0.0E0), (0.0E0,0.0E0),
     +                  (0.0E0,0.0E0), (0.0E0,0.0E0), (0.0E0,0.0E0),
     +                  (0.32E0,-1.41E0), (0.0E0,0.0E0), (0.0E0,0.0E0),
     +                  (0.0E0,0.0E0), (0.0E0,0.0E0), (0.0E0,0.0E0),
     +                  (0.0E0,0.0E0), (0.32E0,-1.41E0),
     +                  (-1.55E0,0.5E0), (0.0E0,0.0E0), (0.0E0,0.0E0),
     +                  (0.0E0,0.0E0), (0.0E0,0.0E0), (0.0E0,0.0E0),
     +                  (0.32E0,-1.41E0), (-1.55E0,0.5E0),
     +                  (0.03E0,-0.89E0), (-0.38E0,-0.96E0),
     +                  (0.0E0,0.0E0), (0.0E0,0.0E0), (0.0E0,0.0E0)/
      DATA              ((CT8(I,J,2),I=1,7),J=1,4)/(0.6E0,-0.6E0),
     +                  (0.0E0,0.0E0), (0.0E0,0.0E0), (0.0E0,0.0E0),
     +                  (0.0E0,0.0E0), (0.0E0,0.0E0), (0.0E0,0.0E0),
     +                  (0.32E0,-1.41E0), (0.0E0,0.0E0), (0.0E0,0.0E0),
     +                  (0.0E0,0.0E0), (0.0E0,0.0E0), (0.0E0,0.0E0),
     +                  (0.0E0,0.0E0), (-0.07E0,-0.89E0),
     +                  (-0.9E0,0.5E0), (0.42E0,-1.41E0), (0.0E0,0.0E0),
     +                  (0.0E0,0.0E0), (0.0E0,0.0E0), (0.0E0,0.0E0),
     +                  (0.78E0,0.06E0), (-0.9E0,0.5E0),
     +                  (0.06E0,-0.13E0), (0.1E0,-0.5E0),
     +                  (-0.77E0,-0.49E0), (-0.5E0,-0.3E0),
     +                  (0.52E0,-1.51E0)/
      DATA              ((CT8(I,J,3),I=1,7),J=1,4)/(0.6E0,-0.6E0),
     +                  (0.0E0,0.0E0), (0.0E0,0.0E0), (0.0E0,0.0E0),
     +                  (0.0E0,0.0E0), (0.0E0,0.0E0), (0.0E0,0.0E0),
     +                  (0.32E0,-1.41E0), (0.0E0,0.0E0), (0.0E0,0.0E0),
     +                  (0.0E0,0.0E0), (0.0E0,0.0E0), (0.0E0,0.0E0),
     +                  (0.0E0,0.0E0), (-0.07E0,-0.89E0),
     +                  (-1.18E0,-0.31E0), (0.0E0,0.0E0), (0.0E0,0.0E0),
     +                  (0.0E0,0.0E0), (0.0E0,0.0E0), (0.0E0,0.0E0),
     +                  (0.78E0,0.06E0), (-1.54E0,0.97E0),
     +                  (0.03E0,-0.89E0), (-0.18E0,-1.31E0),
     +                  (0.0E0,0.0E0), (0.0E0,0.0E0), (0.0E0,0.0E0)/
      DATA              ((CT8(I,J,4),I=1,7),J=1,4)/(0.6E0,-0.6E0),
     +                  (0.0E0,0.0E0), (0.0E0,0.0E0), (0.0E0,0.0E0),
     +                  (0.0E0,0.0E0), (0.0E0,0.0E0), (0.0E0,0.0E0),
     +                  (0.32E0,-1.41E0), (0.0E0,0.0E0), (0.0E0,0.0E0),
     +                  (0.0E0,0.0E0), (0.0E0,0.0E0), (0.0E0,0.0E0),
     +                  (0.0E0,0.0E0), (0.32E0,-1.41E0), (-0.9E0,0.5E0),
     +                  (0.05E0,-0.6E0), (0.0E0,0.0E0), (0.0E0,0.0E0),
     +                  (0.0E0,0.0E0), (0.0E0,0.0E0), (0.32E0,-1.41E0),
     +                  (-0.9E0,0.5E0), (0.05E0,-0.6E0), (0.1E0,-0.5E0),
     +                  (-0.77E0,-0.49E0), (-0.5E0,-0.3E0),
     +                  (0.32E0,-1.16E0)/
      DATA              CT7/(0.0E0,0.0E0), (-0.06E0,-0.90E0),
     +                  (0.65E0,-0.47E0), (-0.34E0,-1.22E0),
     +                  (0.0E0,0.0E0), (-0.06E0,-0.90E0),
     +                  (-0.59E0,-1.46E0), (-1.04E0,-0.04E0),
     +                  (0.0E0,0.0E0), (-0.06E0,-0.90E0),
     +                  (-0.83E0,0.59E0), (0.07E0,-0.37E0),
     +                  (0.0E0,0.0E0), (-0.06E0,-0.90E0),
     +                  (-0.76E0,-1.15E0), (-1.33E0,-1.82E0)/
      DATA              CT6/(0.0E0,0.0E0), (0.90E0,0.06E0),
     +                  (0.91E0,-0.77E0), (1.80E0,-0.10E0),
     +                  (0.0E0,0.0E0), (0.90E0,0.06E0), (1.45E0,0.74E0),
     +                  (0.20E0,0.90E0), (0.0E0,0.0E0), (0.90E0,0.06E0),
     +                  (-0.55E0,0.23E0), (0.83E0,-0.39E0),
     +                  (0.0E0,0.0E0), (0.90E0,0.06E0), (1.04E0,0.79E0),
     +                  (1.95E0,1.22E0)/
      DATA              ((CT10X(I,J,1),I=1,7),J=1,4)/(0.7E0,-0.8E0),
     +                  (0.0E0,0.0E0), (0.0E0,0.0E0), (0.0E0,0.0E0),
     +                  (0.0E0,0.0E0), (0.0E0,0.0E0), (0.0E0,0.0E0),
     +                  (0.6E0,-0.6E0), (0.0E0,0.0E0), (0.0E0,0.0E0),
     +                  (0.0E0,0.0E0), (0.0E0,0.0E0), (0.0E0,0.0E0),
     +                  (0.0E0,0.0E0), (0.6E0,-0.6E0), (-0.9E0,0.5E0),
     +                  (0.0E0,0.0E0), (0.0E0,0.0E0), (0.0E0,0.0E0),
     +                  (0.0E0,0.0E0), (0.0E0,0.0E0), (0.6E0,-0.6E0),
     +                  (-0.9E0,0.5E0), (0.7E0,-0.6E0), (0.1E0,-0.5E0),
     +                  (0.0E0,0.0E0), (0.0E0,0.0E0), (0.0E0,0.0E0)/
      DATA              ((CT10X(I,J,2),I=1,7),J=1,4)/(0.7E0,-0.8E0),
     +                  (0.0E0,0.0E0), (0.0E0,0.0E0), (0.0E0,0.0E0),
     +                  (0.0E0,0.0E0), (0.0E0,0.0E0), (0.0E0,0.0E0),
     +                  (0.6E0,-0.6E0), (0.0E0,0.0E0), (0.0E0,0.0E0),
     +                  (0.0E0,0.0E0), (0.0E0,0.0E0), (0.0E0,0.0E0),
     +                  (0.0E0,0.0E0), (0.7E0,-0.6E0), (-0.4E0,-0.7E0),
     +                  (0.6E0,-0.6E0), (0.0E0,0.0E0), (0.0E0,0.0E0),
     +                  (0.0E0,0.0E0), (0.0E0,0.0E0), (0.8E0,-0.7E0),
     +                  (-0.4E0,-0.7E0), (-0.1E0,-0.2E0),
     +                  (0.2E0,-0.8E0), (0.7E0,-0.6E0), (0.1E0,0.4E0),
     +                  (0.6E0,-0.6E0)/
      DATA              ((CT10X(I,J,3),I=1,7),J=1,4)/(0.7E0,-0.8E0),
     +                  (0.0E0,0.0E0), (0.0E0,0.0E0), (0.0E0,0.0E0),
     +                  (0.0E0,0.0E0), (0.0E0,0.0E0), (0.0E0,0.0E0),
     +                  (0.6E0,-0.6E0), (0.0E0,0.0E0), (0.0E0,0.0E0),
     +                  (0.0E0,0.0E0), (0.0E0,0.0E0), (0.0E0,0.0E0),
     +                  (0.0E0,0.0E0), (-0.9E0,0.5E0), (-0.4E0,-0.7E0),
     +                  (0.6E0,-0.6E0), (0.0E0,0.0E0), (0.0E0,0.0E0),
     +                  (0.0E0,0.0E0), (0.0E0,0.0E0), (0.1E0,-0.5E0),
     +                  (-0.4E0,-0.7E0), (0.7E0,-0.6E0), (0.2E0,-0.8E0),
     +                  (-0.9E0,0.5E0), (0.1E0,0.4E0), (0.6E0,-0.6E0)/
      DATA              ((CT10X(I,J,4),I=1,7),J=1,4)/(0.7E0,-0.8E0),
     +                  (0.0E0,0.0E0), (0.0E0,0.0E0), (0.0E0,0.0E0),
     +                  (0.0E0,0.0E0), (0.0E0,0.0E0), (0.0E0,0.0E0),
     +                  (0.6E0,-0.6E0), (0.0E0,0.0E0), (0.0E0,0.0E0),
     +                  (0.0E0,0.0E0), (0.0E0,0.0E0), (0.0E0,0.0E0),
     +                  (0.0E0,0.0E0), (0.6E0,-0.6E0), (0.7E0,-0.6E0),
     +                  (0.0E0,0.0E0), (0.0E0,0.0E0), (0.0E0,0.0E0),
     +                  (0.0E0,0.0E0), (0.0E0,0.0E0), (0.6E0,-0.6E0),
     +                  (0.7E0,-0.6E0), (-0.1E0,-0.2E0), (0.8E0,-0.7E0),
     +                  (0.0E0,0.0E0), (0.0E0,0.0E0), (0.0E0,0.0E0)/
      DATA              ((CT10Y(I,J,1),I=1,7),J=1,4)/(0.6E0,-0.6E0),
     +                  (0.0E0,0.0E0), (0.0E0,0.0E0), (0.0E0,0.0E0),
     +                  (0.0E0,0.0E0), (0.0E0,0.0E0), (0.0E0,0.0E0),
     +                  (0.7E0,-0.8E0), (0.0E0,0.0E0), (0.0E0,0.0E0),
     +                  (0.0E0,0.0E0), (0.0E0,0.0E0), (0.0E0,0.0E0),
     +                  (0.0E0,0.0E0), (0.7E0,-0.8E0), (-0.4E0,-0.7E0),
     +                  (0.0E0,0.0E0), (0.0E0,0.0E0), (0.0E0,0.0E0),
     +                  (0.0E0,0.0E0), (0.0E0,0.0E0), (0.7E0,-0.8E0),
     +                  (-0.4E0,-0.7E0), (-0.1E0,-0.9E0),
     +                  (0.2E0,-0.8E0), (0.0E0,0.0E0), (0.0E0,0.0E0),
     +                  (0.0E0,0.0E0)/
      DATA              ((CT10Y(I,J,2),I=1,7),J=1,4)/(0.6E0,-0.6E0),
     +                  (0.0E0,0.0E0), (0.0E0,0.0E0), (0.0E0,0.0E0),
     +                  (0.0E0,0.0E0), (0.0E0,0.0E0), (0.0E0,0.0E0),
     +                  (0.7E0,-0.8E0), (0.0E0,0.0E0), (0.0E0,0.0E0),
     +                  (0.0E0,0.0E0), (0.0E0,0.0E0), (0.0E0,0.0E0),
     +                  (0.0E0,0.0E0), (-0.1E0,-0.9E0), (-0.9E0,0.5E0),
     +                  (0.7E0,-0.8E0), (0.0E0,0.0E0), (0.0E0,0.0E0),
     +                  (0.0E0,0.0E0), (0.0E0,0.0E0), (-0.6E0,0.6E0),
     +                  (-0.9E0,0.5E0), (-0.9E0,-0.4E0), (0.1E0,-0.5E0),
     +                  (-0.1E0,-0.9E0), (-0.5E0,-0.3E0),
     +                  (0.7E0,-0.8E0)/
      DATA              ((CT10Y(I,J,3),I=1,7),J=1,4)/(0.6E0,-0.6E0),
     +                  (0.0E0,0.0E0), (0.0E0,0.0E0), (0.0E0,0.0E0),
     +                  (0.0E0,0.0E0), (0.0E0,0.0E0), (0.0E0,0.0E0),
     +                  (0.7E0,-0.8E0), (0.0E0,0.0E0), (0.0E0,0.0E0),
     +                  (0.0E0,0.0E0), (0.0E0,0.0E0), (0.0E0,0.0E0),
     +                  (0.0E0,0.0E0), (-0.1E0,-0.9E0), (0.7E0,-0.8E0),
     +                  (0.0E0,0.0E0), (0.0E0,0.0E0), (0.0E0,0.0E0),
     +                  (0.0E0,0.0E0), (0.0E0,0.0E0), (-0.6E0,0.6E0),
     +                  (-0.9E0,-0.4E0), (-0.1E0,-0.9E0),
     +                  (0.7E0,-0.8E0), (0.0E0,0.0E0), (0.0E0,0.0E0),
     +                  (0.0E0,0.0E0)/
      DATA              ((CT10Y(I,J,4),I=1,7),J=1,4)/(0.6E0,-0.6E0),
     +                  (0.0E0,0.0E0), (0.0E0,0.0E0), (0.0E0,0.0E0),
     +                  (0.0E0,0.0E0), (0.0E0,0.0E0), (0.0E0,0.0E0),
     +                  (0.7E0,-0.8E0), (0.0E0,0.0E0), (0.0E0,0.0E0),
     +                  (0.0E0,0.0E0), (0.0E0,0.0E0), (0.0E0,0.0E0),
     +                  (0.0E0,0.0E0), (0.7E0,-0.8E0), (-0.9E0,0.5E0),
     +                  (-0.4E0,-0.7E0), (0.0E0,0.0E0), (0.0E0,0.0E0),
     +                  (0.0E0,0.0E0), (0.0E0,0.0E0), (0.7E0,-0.8E0),
     +                  (-0.9E0,0.5E0), (-0.4E0,-0.7E0), (0.1E0,-0.5E0),
     +                  (-0.1E0,-0.9E0), (-0.5E0,-0.3E0),
     +                  (0.2E0,-0.8E0)/
      DATA              CSIZE1/(0.0E0,0.0E0), (0.9E0,0.9E0),
     +                  (1.63E0,1.73E0), (2.90E0,2.78E0)/
      DATA              CSIZE3/(0.0E0,0.0E0), (0.0E0,0.0E0),
     +                  (0.0E0,0.0E0), (0.0E0,0.0E0), (0.0E0,0.0E0),
     +                  (0.0E0,0.0E0), (0.0E0,0.0E0), (1.17E0,1.17E0),
     +                  (1.17E0,1.17E0), (1.17E0,1.17E0),
     +                  (1.17E0,1.17E0), (1.17E0,1.17E0),
     +                  (1.17E0,1.17E0), (1.17E0,1.17E0)/
      DATA              CSIZE2/(0.0E0,0.0E0), (0.0E0,0.0E0),
     +                  (0.0E0,0.0E0), (0.0E0,0.0E0), (0.0E0,0.0E0),
     +                  (0.0E0,0.0E0), (0.0E0,0.0E0), (1.54E0,1.54E0),
     +                  (1.54E0,1.54E0), (1.54E0,1.54E0),
     +                  (1.54E0,1.54E0), (1.54E0,1.54E0),
     +                  (1.54E0,1.54E0), (1.54E0,1.54E0)/
*     .. Executable Statements ..
      DO 60 KI = 1, 4
         INCX = INCXS(KI)
         INCY = INCYS(KI)
         MX = ABS(INCX)
         MY = ABS(INCY)
*
         DO 40 KN = 1, 4
            N = NS(KN)
            KSIZE = MIN(2,KN)
            LENX = LENS(KN,MX)
            LENY = LENS(KN,MY)
*           .. initialize all argument arrays ..
            DO 20 I = 1, 7
               CX(I) = CX1(I)
               CY(I) = CY1(I)
   20       CONTINUE
            IF (ICASE.EQ.1) THEN
*              .. CDOTC ..
               CDOT(1) = CDOTC(N,CX,INCX,CY,INCY)
               CALL CTEST(1,CDOT,CT6(KN,KI),CSIZE1(KN),SFAC)
            ELSE IF (ICASE.EQ.2) THEN
*              .. CDOTU ..
               CDOT(1) = CDOTU(N,CX,INCX,CY,INCY)
               CALL CTEST(1,CDOT,CT7(KN,KI),CSIZE1(KN),SFAC)
            ELSE IF (ICASE.EQ.3) THEN
*              .. CAXPY ..
               CALL CAXPY(N,CA,CX,INCX,CY,INCY)
               CALL CTEST(LENY,CY,CT8(1,KN,KI),CSIZE2(1,KSIZE),SFAC)
            ELSE IF (ICASE.EQ.4) THEN
*              .. CCOPY ..
               CALL CCOPY(N,CX,INCX,CY,INCY)
               CALL CTEST(LENY,CY,CT10Y(1,KN,KI),CSIZE3,1.0E0)
               IF (KI.EQ.1) THEN
                  CX0(1) = (42.0E0,43.0E0)
                  CY0(1) = (44.0E0,45.0E0)
                  IF (N.EQ.0) THEN
                     CTY0(1) = CY0(1)
                  ELSE
                     CTY0(1) = CX0(1)
                  END IF
                  LINCX = INCX
                  INCX = 0
                  LINCY = INCY
                  INCY = 0
                  CALL CCOPY(N,CX0,INCX,CY0,INCY)
                  CALL CTEST(1,CY0,CTY0,CSIZE3,1.0E0)
                  INCX = LINCX
                  INCY = LINCY
               END IF
            ELSE IF (ICASE.EQ.5) THEN
*              .. CSWAP ..
               CALL CSWAP(N,CX,INCX,CY,INCY)
               CALL CTEST(LENX,CX,CT10X(1,KN,KI),CSIZE3,1.0E0)
               CALL CTEST(LENY,CY,CT10Y(1,KN,KI),CSIZE3,1.0E0)
            ELSE
               WRITE (NOUT,*) ' Shouldn''t be here in CHECK2'
               STOP
            END IF
*
   40    CONTINUE
   60 CONTINUE
      RETURN
*
*     End of CHECK2
*
      END
      SUBROUTINE STEST(LEN,SCOMP,STRUE,SSIZE,SFAC)
*     ********************************* STEST **************************
*
*     THIS SUBR COMPARES ARRAYS  SCOMP() AND STRUE() OF LENGTH LEN TO
*     SEE IF THE TERM BY TERM DIFFERENCES, MULTIPLIED BY SFAC, ARE
*     NEGLIGIBLE.
*
*     C. L. LAWSON, JPL, 1974 DEC 10
*
*     .. Parameters ..
      INTEGER          NOUT
      REAL             ZERO
      PARAMETER        (NOUT=6, ZERO=0.0E0)
*     .. Scalar Arguments ..
      REAL             SFAC
      INTEGER          LEN
*     .. Array Arguments ..
      REAL             SCOMP(LEN), SSIZE(LEN), STRUE(LEN)
*     .. Scalars in Common ..
      INTEGER          ICASE, INCX, INCY, MODE, N
      LOGICAL          PASS
*     .. Local Scalars ..
      REAL             SD
      INTEGER          I
*     .. External Functions ..
      REAL             SDIFF
      EXTERNAL         SDIFF
*     .. Intrinsic Functions ..
      INTRINSIC        ABS
*     .. Common blocks ..
      COMMON           /COMBLA/ICASE, N, INCX, INCY, MODE, PASS
*     .. Executable Statements ..
*
      DO 40 I = 1, LEN
         SD = SCOMP(I) - STRUE(I)
         IF (ABS(SFAC*SD) .LE. ABS(SSIZE(I))*EPSILON(ZERO))
     +       GO TO 40
*
*                             HERE    SCOMP(I) IS NOT CLOSE TO STRUE(I).
*
         IF ( .NOT. PASS) GO TO 20
*                             PRINT FAIL MESSAGE AND HEADER.
         PASS = .FALSE.
         WRITE (NOUT,99999)
         WRITE (NOUT,99998)
   20    WRITE (NOUT,99997) ICASE, N, INCX, INCY, MODE, I, SCOMP(I),
     +     STRUE(I), SD, SSIZE(I)
   40 CONTINUE
      RETURN
*
99999 FORMAT ('                                       FAIL')
99998 FORMAT (/' CASE  N INCX INCY MODE  I                            ',
     +       ' COMP(I)                             TRUE(I)  DIFFERENCE',
     +       '     SIZE(I)',/1X)
99997 FORMAT (1X,I4,I3,3I5,I3,2E36.8,2E12.4)
*
*     End of STEST
*
      END
      SUBROUTINE STEST1(SCOMP1,STRUE1,SSIZE,SFAC)
*     ************************* STEST1 *****************************
*
*     THIS IS AN INTERFACE SUBROUTINE TO ACCOMMODATE THE FORTRAN
*     REQUIREMENT THAT WHEN A DUMMY ARGUMENT IS AN ARRAY, THE
*     ACTUAL ARGUMENT MUST ALSO BE AN ARRAY OR AN ARRAY ELEMENT.
*
*     C.L. LAWSON, JPL, 1978 DEC 6
*
*     .. Scalar Arguments ..
      REAL              SCOMP1, SFAC, STRUE1
*     .. Array Arguments ..
      REAL              SSIZE(*)
*     .. Local Arrays ..
      REAL              SCOMP(1), STRUE(1)
*     .. External Subroutines ..
      EXTERNAL          STEST
*     .. Executable Statements ..
*
      SCOMP(1) = SCOMP1
      STRUE(1) = STRUE1
      CALL STEST(1,SCOMP,STRUE,SSIZE,SFAC)
*
      RETURN
*
*     End of STEST1
*
      END
      REAL             FUNCTION SDIFF(SA,SB)
*     ********************************* SDIFF **************************
*     COMPUTES DIFFERENCE OF TWO NUMBERS.  C. L. LAWSON, JPL 1974 FEB 15
*
*     .. Scalar Arguments ..
      REAL                            SA, SB
*     .. Executable Statements ..
      SDIFF = SA - SB
      RETURN
*
*     End of SDIFF
*
      END
      SUBROUTINE CTEST(LEN,CCOMP,CTRUE,CSIZE,SFAC)
*     **************************** CTEST *****************************
*
*     C.L. LAWSON, JPL, 1978 DEC 6
*
*     .. Scalar Arguments ..
      REAL             SFAC
      INTEGER          LEN
*     .. Array Arguments ..
      COMPLEX          CCOMP(LEN), CSIZE(LEN), CTRUE(LEN)
*     .. Local Scalars ..
      INTEGER          I
*     .. Local Arrays ..
      REAL             SCOMP(20), SSIZE(20), STRUE(20)
*     .. External Subroutines ..
      EXTERNAL         STEST
*     .. Intrinsic Functions ..
      INTRINSIC        AIMAG, REAL
*     .. Executable Statements ..
      DO 20 I = 1, LEN
         SCOMP(2*I-1) = REAL(CCOMP(I))
         SCOMP(2*I) = AIMAG(CCOMP(I))
         STRUE(2*I-1) = REAL(CTRUE(I))
         STRUE(2*I) = AIMAG(CTRUE(I))
         SSIZE(2*I-1) = REAL(CSIZE(I))
         SSIZE(2*I) = AIMAG(CSIZE(I))
   20 CONTINUE
*
      CALL STEST(2*LEN,SCOMP,STRUE,SSIZE,SFAC)
      RETURN
*
*     End of CTEST
*
      END
      SUBROUTINE ITEST1(ICOMP,ITRUE)
*     ********************************* ITEST1 *************************
*
*     THIS SUBROUTINE COMPARES THE VARIABLES ICOMP AND ITRUE FOR
*     EQUALITY.
*     C. L. LAWSON, JPL, 1974 DEC 10
*
*     .. Parameters ..
      INTEGER           NOUT
      PARAMETER         (NOUT=6)
*     .. Scalar Arguments ..
      INTEGER           ICOMP, ITRUE
*     .. Scalars in Common ..
      INTEGER           ICASE, INCX, INCY, MODE, N
      LOGICAL           PASS
*     .. Local Scalars ..
      INTEGER           ID
*     .. Common blocks ..
      COMMON            /COMBLA/ICASE, N, INCX, INCY, MODE, PASS
*     .. Executable Statements ..
      IF (ICOMP.EQ.ITRUE) GO TO 40
*
*                            HERE ICOMP IS NOT EQUAL TO ITRUE.
*
      IF ( .NOT. PASS) GO TO 20
*                             PRINT FAIL MESSAGE AND HEADER.
      PASS = .FALSE.
      WRITE (NOUT,99999)
      WRITE (NOUT,99998)
   20 ID = ICOMP - ITRUE
      WRITE (NOUT,99997) ICASE, N, INCX, INCY, MODE, ICOMP, ITRUE, ID
   40 CONTINUE
      RETURN
*
99999 FORMAT ('                                       FAIL')
99998 FORMAT (/' CASE  N INCX INCY MODE                               ',
     +       ' COMP                                TRUE     DIFFERENCE',
     +       /1X)
99997 FORMAT (1X,I4,I3,3I5,2I36,I12)
*
*     End of ITEST1
*
      END

C=======================================================================
      SUBROUTINE SRMACHD(I,X)
C.......................................................................
C
C
C                   DOUBLE PRECISION VERSION
C                   ************************
C
C  MACHINE PARAMETERS : TO  ALTER  THIS  SUBROUTINE  FOR  A PARTICULAR 
C  ++++++++++++++++++   ENVIRONMENT, THE DESIRED SET OF DATA STATEMENT
C  SHOULD BE ACTIVATED BY REMOVING THE "C" FROM COLUMN ONE  AND ADDING
C  THE "C" FOR THE TWO LINES AFTER "... VAX FORTRAN (V5) compiler".
C
C  RADIX IS ALWAYS EQUAL TO 2.D0
C  PREC CAN BE FOUND BY CALLING THE ROBETH SUBROUTINE "PRECD"
C  EPMACH IS APPROXIMATELY EQUAL TO THE EXPONENT PART OF PREC
C  EXMIN, XLGMN, YLGMN AND XBIG CAN BE FOUND BY TRIAL AND ERROR
C
      DOUBLE PRECISION X,RADIX,PREC,EXMIN,XLGMN,YLGMN,XBIG,EPMACH
C  for VAX FORTRAN (V5) compiler (VMS)
C       DATA RADIX,PREC,EXMIN,XLGMN,YLGMN,XBIG,EPMACH
C     *  /2.D0,1.38778D-17,-88.722D0,2.939D-39,-88.7227D0,1.7D38,1.D-17/
C  for IBM-PC F77L compiler
C      DATA RADIX,PREC,EXMIN,XLGMN,YLGMN,XBIG,EPMACH
C    *  /2.D0,5.422D-20,-708.D0,1.D-307,-706.591D0,1.D308,1.D-19/
C  for IBM-PC MICROSOFT FORTRAN compiler
C      DATA RADIX,PREC,EXMIN,XLGMN,YLGMN,XBIG,EPMACH
C    *  /2.D0,5.47522D-18,-745.133D0,0.9D-48,-110.629D0,1.D308,1.D-17/
C  for WATCOM F77 Compiler (32 bits version)
      DATA RADIX,PREC,EXMIN,XLGMN,YLGMN,XBIG,EPMACH
     +     /2.,0.1121D-15,-709.782D0,9.74D-290,-718.433D0,1.797D308,
     +     1.0D-17/
C    */2.,0.1121D-15,-709.782D0,0.974D-312,-718.433D0,1.797D308,1.0D-17/
C  for ULTRIX DEC FORTRAN-77 compiler
C     DATA RADIX,PREC,EXMIN,XLGMN,YLGMN,XBIG,EPMACH
C    */2.,1.12133D-16,-744.44D0,2.226D-308,-708.396D0,1.797D308,1.0D-17/
C  for SUN FORTRAN compiler
C     DATA RADIX,PREC,EXMIN,XLGMN,YLGMN,XBIG,EPMACH
C    */2.,1.12133D-16,-745.13D0,0.494D-323,-744.44D0,1.797D308,1.0D-17/
C  for SILICON GRAPHICS MIPS FORTRAN 77 Compiler
C      DATA RADIX,PREC,EXMIN,XLGMN,YLGMN,XBIG,EPMACH
C    */2.,1.12133D-16,-744.04D0,0.758D-323,-743.75D0,1.797D308,1.0D-17/
C  for HP-UX FORTRAN 77 Compiler
C      DATA RADIX,PREC,EXMIN,XLGMN,YLGMN,XBIG,EPMACH
C    */2.,1.12133D-16,-708.396D0,0.1D-308,-709.09D0,1.797D308,1.0D-17/
C  for DEC-ALPHA FORTRAN Compiler (OpenVMS)
C      DATA RADIX,PREC,EXMIN,XLGMN,YLGMN,XBIG,EPMACH
C    */2.,1.1102D-16,-709.782D0,0.1057D-45,-105.863D0,0.898D307,1.0D-17/
C
      IF (I.EQ.1) X=RADIX
      IF (I.EQ.2) X=PREC
      IF (I.EQ.3) X=EXMIN
      IF (I.EQ.4) X=XLGMN
      IF (I.EQ.5) X=YLGMN
      IF (I.EQ.6) X=XBIG
      IF (I.EQ.7) X=EPMACH
      RETURN
      END
C
C-----------------------------------------------------------------------
      double precision FUNCTION SRXEXPD(X)
C.......................................................................
C
C     Extended exp() function
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DATA NCALL,DMIN,DMAX,XBIG/0,0.D0,0.D0,0.D0/
      IF (NCALL.EQ.0) THEN
         CALL SRMACHD(3,DMIN)
         CALL SRMACHD(6,XBIG)
         XBIG=XBIG/10.D0
         DMAX=DLOG(XBIG)
         NCALL=1
      ENDIF
      IF (X.LE.DMIN) THEN
         SRXEXPD=0.D0
      ELSEIF (X.GE.DMAX) THEN
         SRXEXPD=XBIG
      ELSE
         SRXEXPD=DEXP(X)
      ENDIF
      RETURN
      END
C
C-----------------------------------------------------------------------
C
      SUBROUTINE SRINTGRT(F,FARR,N,FEXT,GEXT,A,B,EPSABS,EPSREL,KEY,
     1           LIMIT,RESULT,ABSERR,NEVAL,IER,WORK,IWORK,NPR,PARAM)

C.......................................................................
C
C   AUTHOR :     QUADPACK
C                ADAPTED FOR ROBETH BY A. RANDRIAMIHARISOA
C.......................................................................
C
      implicit double precision(a-h, o-z)
      DOUBLE PRECISION A,ABSERR,B,EPSABS,EPSREL,F,RESULT,FEXT,WORK
      INTEGER IER,KEY,LAST,LIMIT,NEVAL,ALIST,BLIST,ELIST,RLIST
C
      DIMENSION FARR(N),WORK(4*LIMIT),IWORK(LIMIT),PARAM(NPR)
C
      EXTERNAL F,FEXT,GEXT
C
C         LIMIT IS THE MAXIMUM NUMBER OF SUBINTERVALS ALLOWED IN THE
C         SUBDIVISION PROCESS OF QAGE1D. TAKE CARE THAT LIMIT.GE.1.
C
C**** DATA LIMIT/500/
C
C      IF ((EPSABS.LT.0..AND.EPSREL.LT.0.).OR.LIMIT.LE.1
C     1   .OR.LIMIT.GT.500)  CALL MESSGE(500,'INTGRD',1)
      ALIST=1
      BLIST=ALIST+LIMIT
      RLIST=BLIST+LIMIT
      ELIST=RLIST+LIMIT
      CALL SRQAGE1T(F,FARR,N,FEXT,GEXT,A,B,EPSABS,EPSREL,KEY,LIMIT,
     1     RESULT,ABSERR,NEVAL,IER,WORK,WORK(BLIST),
     2     WORK(RLIST),WORK(ELIST),IWORK,LAST,NPR,PARAM)
C
      RETURN
      END
C
C-----------------------------------------------------------------------
C
      SUBROUTINE SRQAGE1T(F,FARR,N,FEXT,GEXT,A,B,EPSABS,EPSREL,KEY,
     *  LIMIT,RESULT,ABSERR,NEVAL,IER,ALIST,BLIST,RLIST,ELIST,
     *  IORD,LAST,NPR,PARAM)
C.......................................................................
C
C   AUTHOR :     QUADPACK
C                ADAPTED FOR ROBETH BY A. RANDRIAMIHARISOA
C.......................................................................
C
       implicit double precision(a-h, o-z)
C      DOUBLE PRECISION A,ABSERR,ALIST,AREA,AREA1,AREA12,AREA2,
C     *  AA1,AA2,B,
C     *  BLIST,BB1,BB2,C,DABS,DEFABS,DEFAB1,DEFAB2,DMAX1,ELIST,EPMACH,
C     *  EPSABS,EPSREL,ERRBND,ERRMAX,ERROR1,ERROR2,ERRO12,ERRSUM,F,OFLOW,
C     *  RESABS,RESULT,RLIST,UFLOW,FEXT
      INTEGER IER,IORD,IROFF1,IROFF2,K,KEY,KEYF,LAST,LIMIT,MAXERR,NEVAL,
     *  NRMAX
C
      DIMENSION ALIST(LIMIT),BLIST(LIMIT),ELIST(LIMIT),IORD(LIMIT),
     *  RLIST(LIMIT),FARR(N),param(npr)
C
      EXTERNAL F,FEXT,GEXT
C
C            LIST OF MAJOR VARIABLES
C            -----------------------
C
C           ALIST     - LIST OF LEFT END POINTS OF ALL SUBINTERVALS
C                       CONSIDERED UP TO NOW
C           BLIST     - LIST OF RIGHT END POINTS OF ALL SUBINTERVALS
C                       CONSIDERED UP TO NOW
C           RLIST(I)  - APPROXIMATION TO THE INTEGRAL OVER
C                      (ALIST(I),BLIST(I))
C           ELIST(I)  - ERROR ESTIMATE APPLYING TO RLIST(I)
C           MAXERR    - POINTER TO THE INTERVAL WITH LARGEST
C                       ERROR ESTIMATE
C           ERRMAX    - ELIST(MAXERR)
C           AREA      - SUM OF THE INTEGRALS OVER THE SUBINTERVALS
C           ERRSUM    - SUM OF THE ERRORS OVER THE SUBINTERVALS
C           ERRBND    - REQUESTED ACCURACY MAX(EPSABS,EPSREL*
C                       ABS(RESULT))
C           *****1    - VARIABLE FOR THE LEFT SUBINTERVAL
C           *****2    - VARIABLE FOR THE RIGHT SUBINTERVAL
C           LAST      - INDEX FOR SUBDIVISION
C
C
C           MACHINE DEPENDENT CONSTANTS
C           ---------------------------
C
C           EPMACH  IS THE LARGEST RELATIVE SPACING.
C           UFLOW  IS THE SMALLEST POSITIVE MAGNITUDE.
C           OFLOW  IS THE LARGEST MAGNITUDE.
C
C***FIRST EXECUTABLE STATEMENTS
      CALL SRMACHD(7,EPMACH)
      CALL SRMACHD(4,UFLOW)
      CALL SRMACHD(6,OFLOW)
C
C           TEST ON VALIDITY OF PARAMETERS
C           ------------------------------
C
      NEVAL = 0
      LAST = 0
      RESULT = 0.0D+00
      ABSERR = 0.0D+00
      ALIST(1) = A
      BLIST(1) = B
      RLIST(1) = 0.0D+00
      ELIST(1) = 0.0D+00
      IORD(1) = 0
      IER=6
C      IF ((EPSABS.LT.0..AND.EPSREL.LT.0.).OR.LIMIT.LE.1
C     1   .OR.LIMIT.GT.500)  CALL MESSGE(500,'QAGE1D',1)
      IER = 0
C
C           FIRST APPROXIMATION TO THE INTEGRAL
C           -----------------------------------
C
      KEYF = KEY
      IF(KEY.LE.0) KEYF = 1
      IF(KEY.GE.7) KEYF = 6
      C = dble(FLOAT(KEYF))
      NEVAL = 0
      IF (KEYF.EQ.1)
     *  CALL SRQ1K15T(F,FARR,N,FEXT,GEXT,A,B,RESULT,ABSERR,DEFABS,
     *  RESABS,NPR,PARAM)

      LAST = 1
      RLIST(1) = RESULT
      ELIST(1) = ABSERR
      IORD(1) = 1
C
C           TEST ON ACCURACY.
C
      ERRBND = DMAX1(EPSABS,EPSREL*DABS(RESULT))
      IF(ABSERR.LE.5.0D+01*EPMACH*DEFABS.AND.ABSERR.GT.ERRBND) IER = 2
C LIMIT==1 MEANS THAT ONLY 1 SUBINTERVAL IS ALLOWED..
      IF(LIMIT.EQ.1) IER = 1
      IF(IER.NE.0.OR.(ABSERR.LE.ERRBND.AND.ABSERR.NE.RESABS)
     *  .OR.ABSERR.EQ.0.0D+00) GO TO 60
C
C           INITIALIZATION
C           --------------
C
C
      ERRMAX = ABSERR
      MAXERR = 1
      AREA = RESULT
      ERRSUM = ABSERR
      NRMAX = 1
      IROFF1 = 0
      IROFF2 = 0
C
C           MAIN DO-LOOP
C           ------------
C
      DO 30 LAST = 2,LIMIT
C
C           BISECT THE SUBINTERVAL WITH THE LARGEST ERROR ESTIMATE.
C
        AA1 = ALIST(MAXERR)
        BB1 = 5.0D-01*(ALIST(MAXERR)+BLIST(MAXERR))
        AA2 = BB1
        BB2 = BLIST(MAXERR)
        IF (KEYF.EQ.1)
     *  CALL SRQ1K15T(F,FARR,N,FEXT,GEXT,AA1,BB1,AREA1,ERROR1,
     *          RESABS,DEFAB1,NPR,PARAM)
        IF (KEYF.EQ.1)
     *  CALL SRQ1K15T(F,FARR,N,FEXT,GEXT,AA2,BB2,AREA2,ERROR2,
     *          RESABS,DEFAB2,NPR,PARAM)
        NEVAL = NEVAL+1
        AREA12 = AREA1+AREA2
        ERRO12 = ERROR1+ERROR2
        ERRSUM = ERRSUM+ERRO12-ERRMAX
        AREA = AREA+AREA12-RLIST(MAXERR)
        IF(DEFAB1.EQ.ERROR1.OR.DEFAB2.EQ.ERROR2) GO TO 5
        IF(DABS(RLIST(MAXERR)-AREA12).LE.1.0D-05*DABS(AREA12)
     *  .AND.ERRO12.GE.9.9D-01*ERRMAX) IROFF1 = IROFF1+1
        IF(LAST.GT.10.AND.ERRO12.GT.ERRMAX) IROFF2 = IROFF2+1
    5   RLIST(MAXERR) = AREA1
        RLIST(LAST) = AREA2
        ERRBND = DMAX1(EPSABS,EPSREL*DABS(AREA))
        IF(ERRSUM.LE.ERRBND) GO TO 8
C
C           TEST FOR ROUNDOFF ERROR AND EVENTUALLY SET ERROR FLAG.
C
        IF(IROFF1.GE.6.OR.IROFF2.GE.20) IER = 2
C
C           SET ERROR FLAG IN THE CASE THAT THE NUMBER OF
C           SUBINTERVALS EQUALS LIMIT.
C
        IF(LAST.EQ.LIMIT) IER = 1
C
C           SET ERROR FLAG IN THE CASE OF BAD INTEGRAND BEHAVIOUR
C           AT A POINT OF THE INTEGRATION RANGE.
C
        IF(DMAX1(DABS(AA1),DABS(BB2)).LE.(1.0D+00+C*1.0D+03*
     *  EPMACH)*(DABS(AA2)+1.0D+04*UFLOW)) IER = 3
C
C           APPEND THE NEWLY-CREATED INTERVALS TO THE LIST.
C
    8   IF(ERROR2.GT.ERROR1) GO TO 10
        ALIST(LAST) = AA2
        BLIST(MAXERR) = BB1
        BLIST(LAST) = BB2
        ELIST(MAXERR) = ERROR1
        ELIST(LAST) = ERROR2
        GO TO 20
   10   ALIST(MAXERR) = AA2
        ALIST(LAST) = AA1
        BLIST(LAST) = BB1
        RLIST(MAXERR) = AREA2
        RLIST(LAST) = AREA1
        ELIST(MAXERR) = ERROR2
        ELIST(LAST) = ERROR1
C
C           CALL SUBROUTINE QSORTD TO MAINTAIN THE DESCENDING ORDERING
C           IN THE LIST OF ERROR ESTIMATES AND SELECT THE SUBINTERVAL
C           WITH THE LARGEST ERROR ESTIMATE (TO BE BISECTED NEXT).
C
   20   CALL SRQSORTD(LIMIT,LAST,MAXERR,ERRMAX,ELIST,IORD,NRMAX)
C***JUMP OUT OF DO-LOOP
        IF(IER.NE.0.OR.ERRSUM.LE.ERRBND) GO TO 40
   30 CONTINUE
C
C           COMPUTE FINAL RESULT.
C           ---------------------
C
   40 RESULT = 0.0D+00
      DO 50 K=1,LAST
        RESULT = RESULT+RLIST(K)
   50 CONTINUE
      ABSERR = ERRSUM
   60 IF(KEYF.NE.1) NEVAL = (10*KEYF+1)*(2*NEVAL+1)
      IF(KEYF.EQ.1) NEVAL = 30*NEVAL+15
C  999 IF (IER.NE.0) CALL MESSGE(400+IER,'QAGE1D',0)
c 999 CONTINUE
      RETURN
      END
C
C-----------------------------------------------------------------------
C
      SUBROUTINE SRQ1K15T
     *  (F,FARR,N,FEXT,GEXT,A,B,RESULT,ABSERR,RESABS,RESASC,NPR,PARAM)
C.......................................................................
C
C   AUTHOR :     QUADPACK
C                ADAPTED FOR ROBETH BY A. RANDRIAMIHARISOA
C.......................................................................
C
      implicit double precision(a-h, o-z)
      DOUBLE PRECISION A,ABSC,ABSERR,B,CENTR,DABS,DHLGTH,DMAX1,DMIN1,
     *  EPMACH,F,FC,FSUM,FVAL1,FVAL2,FV1,FV2,HLGTH,OFLOW,RESABS,RESASC,
     *  RESG,RESK,RESKH,RESULT,UFLOW,WG,WGK,XGK,FEXT
      INTEGER J,JTW,JTWM1
      EXTERNAL F,FEXT,GEXT
C
      DIMENSION FV1(7),FV2(7),WG(4),WGK(8),XGK(8),FARR(N),param(npr)
C
C           THE ABSCISSAE AND WEIGHTS ARE GIVEN FOR THE INTERVAL (-1,1).
C           BECAUSE OF SYMMETRY ONLY THE POSITIVE ABSCISSAE AND THEIR
C           CORRESPONDING WEIGHTS ARE GIVEN.
C
C           XGK    - ABSCISSAE OF THE 15-POINT KRONROD RULE
C                    XGK(2), XGK(4), ...  ABSCISSAE OF THE 7-POINT
C                    GAUSS RULE
C                    XGK(1), XGK(3), ...  ABSCISSAE WHICH ARE OPTIMALLY
C                    ADDED TO THE 7-POINT GAUSS RULE
C
C           WGK    - WEIGHTS OF THE 15-POINT KRONROD RULE
C
C           WG     - WEIGHTS OF THE 7-POINT GAUSS RULE
C
      DATA XGK(1),XGK(2),XGK(3),XGK(4),XGK(5),XGK(6),XGK(7),XGK(8)/
     *     9.914553711208126D-01,   9.491079123427585D-01,
     *     8.648644233597691D-01,   7.415311855993944D-01,
     *     5.860872354676911D-01,   4.058451513773972D-01,
     *     2.077849550078985D-01,   0.0D+00              /
      DATA WGK(1),WGK(2),WGK(3),WGK(4),WGK(5),WGK(6),WGK(7),WGK(8)/
     *     2.293532201052922D-02,   6.309209262997855D-02,
     *     1.047900103222502D-01,   1.406532597155259D-01,
     *     1.690047266392679D-01,   1.903505780647854D-01,
     *     2.044329400752989D-01,   2.094821410847278D-01/
      DATA WG(1),WG(2),WG(3),WG(4)/
     *     1.294849661688697D-01,   2.797053914892767D-01,
     *     3.818300505051189D-01,   4.179591836734694D-01/
C
C
C           LIST OF MAJOR VARIABLES
C           -----------------------
C
C           CENTR  - MID POINT OF THE INTERVAL
C           HLGTH  - HALF-LENGTH OF THE INTERVAL
C           ABSC   - ABSCISSA
C           FVAL*  - FUNCTION VALUE
C           RESG   - RESULT OF THE 7-POINT GAUSS FORMULA
C           RESK   - RESULT OF THE 15-POINT KRONROD FORMULA
C           RESKH  - APPROXIMATION TO THE MEAN VALUE OF F OVER (A,B),
C                    I.E. TO I/(B-A)
C
C           MACHINE DEPENDENT CONSTANTS
C           ---------------------------
C
C           EPMACH IS THE LARGEST RELATIVE SPACING.
C           UFLOW IS THE SMALLEST POSITIVE MAGNITUDE.
C           OFLOW IS THE LARGEST MAGNITUDE.
C
C***FIRST EXECUTABLE STATEMENTS
      CALL SRMACHD(7,EPMACH)
      CALL SRMACHD(4,UFLOW)
      CALL SRMACHD(6,OFLOW)
C
      CENTR = 5.0D-01*(A+B)
      HLGTH = 5.0D-01*(B-A)
      DHLGTH = DABS(HLGTH)
C
C           COMPUTE THE 15-POINT KRONROD APPROXIMATION TO
C           THE INTEGRAL, AND ESTIMATE THE ABSOLUTE ERROR.
C
      FC = F(CENTR,FARR,N,FEXT,GEXT,NPR,PARAM)

      RESG = FC*WG(4)
      RESK = FC*WGK(8)
      RESABS = DABS(RESK)
      DO 10 J=1,3
        JTW = J*2
        ABSC = HLGTH*XGK(JTW)
        FVAL1 = F(CENTR-ABSC,FARR,N,FEXT,GEXT,NPR,PARAM)

        FVAL2 = F(CENTR+ABSC,FARR,N,FEXT,GEXT,NPR,PARAM)

        FV1(JTW) = FVAL1
        FV2(JTW) = FVAL2
        FSUM = FVAL1+FVAL2
        RESG = RESG+WG(J)*FSUM
        RESK = RESK+WGK(JTW)*FSUM
        RESABS = RESABS+WGK(JTW)*(DABS(FVAL1)+DABS(FVAL2))
   10 CONTINUE
      DO 15 J = 1,4
        JTWM1 = J*2-1
        ABSC = HLGTH*XGK(JTWM1)
        FVAL1 = F(CENTR-ABSC,FARR,N,FEXT,GEXT,NPR,PARAM)
        FVAL2 = F(CENTR+ABSC,FARR,N,FEXT,GEXT,NPR,PARAM)
        FV1(JTWM1) = FVAL1
        FV2(JTWM1) = FVAL2
        FSUM = FVAL1+FVAL2
        RESK = RESK+WGK(JTWM1)*FSUM
        RESABS = RESABS+WGK(JTWM1)*(DABS(FVAL1)+DABS(FVAL2))
   15 CONTINUE
      RESKH = RESK*5.0D-01
      RESASC = WGK(8)*DABS(FC-RESKH)
      DO 20 J=1,7
        RESASC = RESASC+WGK(J)*(DABS(FV1(J)-RESKH)+DABS(FV2(J)-RESKH))
   20 CONTINUE
      RESULT = RESK*HLGTH
      RESABS = RESABS*DHLGTH
      RESASC = RESASC*DHLGTH
      ABSERR = DABS((RESK-RESG)*HLGTH)
      IF(RESASC.NE.0.0D+00.AND.ABSERR.NE.0.0D+00)
     *  ABSERR = RESASC*DMIN1(1.0D+00,(2.0D+02*ABSERR/RESASC)**1.5D+00)
      IF(RESABS.GT.UFLOW/(5.0D+01*EPMACH)) ABSERR = DMAX1
     *  ((EPMACH*5.0D+01)*RESABS,ABSERR)
      RETURN
      END
c
      SUBROUTINE SRQSORTD(LIMIT,LAST,MAXERR,ERMAX,ELIST,IORD,NRMAX)
      DOUBLE PRECISION ELIST,ERMAX,ERRMAX,ERRMIN
      INTEGER I,IBEG,IDO,IORD,ISUCC,J,JBND,JUPBN,K,LAST,LIMIT,MAXERR,
     *  NRMAX
      DIMENSION ELIST(LAST),IORD(LAST)
C
C           CHECK WHETHER THE LIST CONTAINS MORE THAN
C           TWO ERROR ESTIMATES.
C
C***FIRST EXECUTABLE STATEMENT
      IF(LAST.GT.2) GO TO 10
      IORD(1) = 1
      IORD(2) = 2
      GO TO 90
C
C           THIS PART OF THE ROUTINE IS ONLY EXECUTED IF, DUE TO A
C           DIFFICULT INTEGRAND, SUBDIVISION INCREASED THE ERROR
C           ESTIMATE. IN THE NORMAL CASE THE INSERT PROCEDURE SHOULD
C           START AFTER THE NRMAX-TH LARGEST ERROR ESTIMATE.
C
   10 ERRMAX = ELIST(MAXERR)
      IF(NRMAX.EQ.1) GO TO 30
      IDO = NRMAX-1
      DO 20 I = 1,IDO
        ISUCC = IORD(NRMAX-1)
C***JUMP OUT OF DO-LOOP
        IF(ERRMAX.LE.ELIST(ISUCC)) GO TO 30
        IORD(NRMAX) = ISUCC
        NRMAX = NRMAX-1
   20    CONTINUE
C
C           COMPUTE THE NUMBER OF ELEMENTS IN THE LIST TO BE
C           MAINTAINED IN DESCENDING ORDER. THIS NUMBER
C           DEPENDS ON THE NUMBER OF SUBDIVISIONS STILL ALLOWED.
C
   30 JUPBN = LAST
      IF(LAST.GT.(LIMIT/2+2)) JUPBN = LIMIT+3-LAST
      ERRMIN = ELIST(LAST)
C
C           INSERT ERRMAX BY TRAVERSING THE LIST TOP-DOWN,
C           STARTING COMPARISON FROM THE ELEMENT ELIST(IORD(NRMAX+1)).
C
      JBND = JUPBN-1
      IBEG = NRMAX+1
      IF(IBEG.GT.JBND) GO TO 50
      DO 40 I=IBEG,JBND
        ISUCC = IORD(I)
C***JUMP OUT OF DO-LOOP
        IF(ERRMAX.GE.ELIST(ISUCC)) GO TO 60
        IORD(I-1) = ISUCC
   40 CONTINUE
   50 IORD(JBND) = MAXERR
      IORD(JUPBN) = LAST
      GO TO 90
C
C           INSERT ERRMIN BY TRAVERSING THE LIST BOTTOM-UP.
C
   60 IORD(I-1) = MAXERR
      K = JBND
      DO 70 J=I,JBND
        ISUCC = IORD(K)
C***JUMP OUT OF DO-LOOP
        IF(ERRMIN.LT.ELIST(ISUCC)) GO TO 80
        IORD(K+1) = ISUCC
        K = K-1
   70 CONTINUE
      IORD(I) = LAST
      GO TO 90
   80 IORD(K+1) = LAST
C
C           SET MAXERR AND ERMAX.
C
   90 MAXERR = IORD(NRMAX)
      ERMAX = ELIST(MAXERR)
      RETURN
      END
C
C-----------------------------------------------------------------------
C
      SUBROUTINE SRXERF(KODE,X,P)
C.......................................................................
C
C   COPYRIGHT 1992 Alfio Marazzi
C
C   AUTHOR : A. MARAZZI
C.......................................................................
C
      implicit double precision(a-h,o-z)
      EXTERNAL SRXEXPD
      DATA SPI/2.506628274631D0/
c      IF (KODE.NE.1.AND.KODE.NE.2) CALL MESSGE(500,'XERF  ',1)
      X2=-X*X/2.D0
      P=SRXEXPD(X2)
      IF (KODE.EQ.2) P=P/SPI
      RETURN
      END
C
C-----------------------------------------------------------------------
C
      SUBROUTINE SRGAUSSD (KODE,X,P)
C.......................................................................
C
C   COPYRIGHT 1992 Alfio Marazzi
C
C   AUTHOR : A. RANDRIAMIHARISOA
C
C.......................................................................
C
      DOUBLE PRECISION   P,X,SQR1D2,CD
      DATA               SQR1D2/.7071067811865475D0/
C
c      IF (KODE.NE.1.AND.KODE.NE.2) CALL MESSGE(500,'GAUSSD',1)
      CALL SRCERFD(-X*SQR1D2,CD)
      P = .5D0 * CD
      IF (KODE.EQ.2) P=1.D0-P
      RETURN
      END
C
C-----------------------------------------------------------------------
C
      SUBROUTINE SRCERFD(X,F)
C.......................................................................
C
C   COPYRIGHT 1992 Alfio Marazzi
C
C   AUTHOR : A. RANDRIAMIHARISOA
C
C.......................................................................
C
      DOUBLE PRECISION   F,X,SRXEXPD
      DIMENSION          P(5),Q(4),P1(9),Q1(8),P2(6),Q2(5)
      DOUBLE PRECISION   P,Q,P1,Q1,P2,Q2,XMIN,XLARGE,SQRPI,XX,
     *                   RES,XSQ,XNUM,XDEN,XI,XBIG
      INTEGER            ISW,I
      EXTERNAL SRXEXPD
C                                  COEFFICIENTS FOR 0.0 .LE. Y .LT.
C                                  .477
      DATA               P(1)/113.8641541510502D0/,
     *                   P(2)/377.4852376853020D0/,
     *                   P(3)/3209.377589138469D0/,
     *                   P(4)/.1857777061846032D0/,
     *                   P(5)/3.161123743870566D0/
      DATA               Q(1)/244.0246379344442D0/,
     *                   Q(2)/1282.616526077372D0/,
     *                   Q(3)/2844.236833439171D0/,
     *                   Q(4)/23.60129095234412D0/
C                                  COEFFICIENTS FOR .477 .LE. Y
C                                  .LE. 4.0
      DATA               P1(1)/8.883149794388376D0/,
     *                   P1(2)/66.11919063714163D0/,
     *                   P1(3)/298.6351381974001D0/,
     *                   P1(4)/881.9522212417691D0/,
     *                   P1(5)/1712.047612634071D0/,
     *                   P1(6)/2051.078377826071D0/,
     *                   P1(7)/1230.339354797997D0/,
     *                   P1(8)/2.153115354744038D-8/,
     *                   P1(9)/.5641884969886701D0/
      DATA               Q1(1)/117.6939508913125D0/,
     *                   Q1(2)/537.1811018620099D0/,
     *                   Q1(3)/1621.389574566690D0/,
     *                   Q1(4)/3290.799235733460D0/,
     *                   Q1(5)/4362.619090143247D0/,
     *                   Q1(6)/3439.367674143722D0/,
     *                   Q1(7)/1230.339354803749D0/,
     *                   Q1(8)/15.74492611070983D0/
C                                  COEFFICIENTS FOR 4.0 .LT. Y
      DATA               P2(1)/-3.603448999498044D-01/,
     *                   P2(2)/-1.257817261112292D-01/,
     *                   P2(3)/-1.608378514874228D-02/,
     *                   P2(4)/-6.587491615298378D-04/,
     *                   P2(5)/-1.631538713730210D-02/,
     *                   P2(6)/-3.053266349612323D-01/
      DATA               Q2(1)/1.872952849923460D0/,
     *                   Q2(2)/5.279051029514284D-01/,
     *                   Q2(3)/6.051834131244132D-02/,
     *                   Q2(4)/2.335204976268692D-03/,
     *                   Q2(5)/2.568520192289822D0/
C                                  CONSTANTS
      DATA               XMIN/1.0D-10/,XLARGE/6.375D0/
C                                  CERFD(XBIG) .APPROX. DETAP
      DATA               XBIG/13.3D0/
      DATA               SQRPI/.5641895835477563D0/
C
      Y=SNGL(X)
      XX = Y
      ISW = 1
      IF (XX.GE.0.0D0) GO TO 5
      ISW = -1
      XX = -XX
    5 IF (XX.LT..477D0) GO TO 10
      IF (XX.LE.4.0D0) GO TO 30
      IF (ISW .GT. 0) GO TO 40
      IF (XX.LT.XLARGE) GO TO 45
      RES = 2.0D0
      GO TO 70
C                                  ABS(Y) .LT. .477, EVALUATE
C                                  APPROXIMATION FOR CERFD
   10 IF (XX.LT.XMIN) GO TO 20
      XSQ = XX*XX
      XNUM = P(4)*XSQ+P(5)
      XDEN = XSQ+Q(4)
      DO 15 I = 1,3
         XNUM = XNUM*XSQ+P(I)
         XDEN = XDEN*XSQ+Q(I)
   15 CONTINUE
      RES = XX*XNUM/XDEN
      GO TO 25
   20 RES = XX*P(3)/Q(3)
   25 IF (ISW.EQ.-1) RES = -RES
      RES = 1.0D0-RES
      GO TO 70
C                                  .477 .LE. ABS(Y) .LE. 4.0
C                                  EVALUATE APPROXIMATION FOR CERFD
   30 XSQ = XX*XX
      XNUM = P1(8)*XX+P1(9)
      XDEN = XX+Q1(8)
      DO 35 I=1,7
         XNUM = XNUM*XX+P1(I)
         XDEN = XDEN*XX+Q1(I)
   35 CONTINUE
      RES = XNUM/XDEN
      GO TO 60
C                                  4.0 .LT. ABS(Y), EVALUATE
C                                  MINIMAX APPROXIMATION FOR CERFD
   40 IF (XX.GT.XBIG) GO TO 65
   45 XSQ = XX*XX
      XI = 1.0D0/XSQ
      XNUM= P2(5)*XI+P2(6)
      XDEN = XI+Q2(5)
      DO 50 I = 1,4
         XNUM = XNUM*XI+P2(I)
         XDEN = XDEN*XI+Q2(I)
   50 CONTINUE
      RES = (SQRPI+XI*XNUM/XDEN)/XX
   60 RES = RES*SRXEXPD(-XSQ)
      IF (ISW.EQ.-1) RES = 2.0D0-RES
      GO TO 70
   65 RES = 0.0D0
   70 F = RES
      RETURN
      END
C
C======================================================================
C
      DOUBLE PRECISION FUNCTION DNORM0(X)
      DOUBLE PRECISION TMP,SPI,X2,X,EXMIN
      DATA NCALL,EXMIN,SPI/0,0.D0,2.506628274631D0/
      IF (NCALL.EQ.0) THEN
        NCALL=1
        CALL SRMACHD(3,EXMIN)
      ENDIF
      X2=-X*X/2.D0
      TMP=0.D0
      IF (X2.GT.EXMIN) TMP=DEXP(X2)
      TMP=TMP/SPI
      DNORM0=TMP
      RETURN
      END
c
      DOUBLE PRECISION FUNCTION PNORM0(Z)
      DOUBLE PRECISION Z,TMP 
      CALL SRGAUSSD(1,Z,TMP)
      PNORM0=TMP
      RETURN
      END
c
      DOUBLE PRECISION FUNCTION SRChisk(S,K0)
      DOUBLE PRECISION S,S2,ABST,TMP,K0 
c     COMMON/PSIPR/IPSI,C,H1,H2,H3,XK,D
c      
c     K0=DBLE(XK)
      TMP=1.D0
      ABST=DABS(S)
      IF (ABST.GE.K0) GOTO 400
      S2=(S/K0)**2
      TMP=(S2*(S2-3.D0)+3.D0)*S2
  400 SRCHISK=TMP-0.5D0
      RETURN
      END
c
      DOUBLE PRECISION FUNCTION SRPSI1N(Z,U)
      DOUBLE PRECISION Z,U 
      SRPSI1N=0.D0
      IF (Z.LT.-U.OR.Z.GT.U) RETURN 
      SRPSI1N=Z
      RETURN
      END
c
      DOUBLE PRECISION FUNCTION SRPSI2N(Z,U)
      DOUBLE PRECISION Z,U
      SRPSI2N=0.D0
      IF (Z.LT.-U.OR.Z.GT.U) RETURN 
      SRPSI2N=Z*Z
      RETURN
      END
c
      DOUBLE PRECISION FUNCTION SRBETAN(u)
      double precision Alfa,u,sum,dnorm0,pnorm0,pnrm0
      external dnorm0,pnorm0
      pnrm0=pnorm0(u)
      Alfa=2.d0*pnrm0-1.d0
      sum=2.d0*(-u*dnorm0(u)+pnrm0-0.5d0)
      SRBETAN=SUM/Alfa
      RETURN
      END
c
      DOUBLE PRECISION FUNCTION IALPHAN(Z0,U,SIGMA,IS0)
      DOUBLE PRECISION Z0,U,SIGMA,IS0,ETA,RHO,TMP,DNORM0,pnorm0,
     *       XLGMN,YLGMN
      EXTERNAL DNORM0,pnorm0
      DATA NCALL,XLGMN,YLGMN/0,0.D0,0.D0/
      IF (NCALL.EQ.0) THEN
        CALL SRMACHD(4,XLGMN)
        CALL SRMACHD(5,YLGMN)
        NCALL=1
      ENDIF
      ETA=dnorm0(U)
      TMP=-YLGMN
      IF (ETA.GT.XLGMN) TMP=-dlog(ETA)      
      ETA=TMP  
      RHO=dnorm0(z0)
      TMP=-YLGMN
      IF (RHO.GT.XLGMN) TMP=-dlog(RHO)
      RHO=TMP
      TMP=2.d0*U*DNORM0(U)*IS0/SIGMA-(2.d0*pnorm0(U)-1.d0)
      IF (ETA.GE.RHO) TMP=TMP+1.D0
      IALPHAN=TMP
      RETURN
      END
c
      SUBROUTINE SRD1N(U,SIGMA,IT0,XtX,NP,VAL)
      DOUBLE PRECISION U,L,SIGMA,IT0(NP),XtX(NP,NP),TMP1,
     +       TMP,DNORM0,EZU,VAL(NP)
      EXTERNAL DNORM0
      L=-U
c     TMP2=(U*U-L*L)*IS0=0.D0
      TMP1=U-L
      EZU=DNORM0(U) 
      DO 200 I=1,NP
      TMP=0.D0
      DO 100 J=1,NP
      TMP=TMP+XtX(I,J)*IT0(J)
  100 CONTINUE
      VAL(I)=TMP1*TMP
      VAL(I)=EZU*(VAL(I))/SIGMA  
  200 CONTINUE
      RETURN
      END
c
      SUBROUTINE SRD2N(U,SIGMA,IS0,VAL)
      DOUBLE PRECISION L,U,SIGMA,IS0,TMP2,DNORM0,EZU,U2,VAL
      EXTERNAL DNORM0
      L=-U
      U2=U*U
      TMP2=(U*U2-L*U2)*IS0
c      TMP1=0.D0
      EZU=DNORM0(U) 
c      TMP=0.D0
c      DO 100 J=1,NP
c      TMP=TMP+XBAR(J)*IT0(J)
c  100 CONTINUE
c      VAL=TMP*TMP1
      VAL=EZU*TMP2/SIGMA  
      RETURN
      END
c
      SUBROUTINE if_tmlnf(y,n,k0,theta,sigma,invm0,its0)
      implicit double precision(A-H,O-Z)
      double precision y(n),k0,theta,invm0(2,2),its0(n,2),sc1(2)
      external SRpsimm, SRchisk
c     COMMON/PSIPR/IPSI,C,H1,H2,H3,XK,D
c     IPSI=4
c     K0=1.548D0
      do 500 k=1,n
        y0=y(k)
        z0=(y0-theta)/sigma
        sc1(1)=SRpsimm(z0,2,k0)
        sc1(2)=SRchisk(z0,k0)
        do 240 i=1,2     
        its0(k,i)=0.d0
        do 230 j=1,2
        its0(k,i)=its0(k,i)+invm0(i,j)*sc1(j)
  230   continue
  240   continue
  500 continue
      return
      end

c
      SUBROUTINE av_tmlnf(X,y,n,np,ncov,u,k0,theta,sigma,invm0,
     +           invm1,avts0,avts,xbar,XtX,sa,sc1,x0,its0,its)
      implicit double precision(A-H,O-Z)
      double precision X(n,np),k0,XtX(np,np),xbar(np),y(n),theta(np),
     +       x0(np),avts0(np+1,np+1),avts(np+1,np+1),invm0(np+1,np+1),
     +       invm1(np+1,np+1),its0(np+1),its(np+1),is0,ialf,IALPHAN,
     +       sa(ncov),sc1(ncov)
      external pnorm0,IALPHAN,SRPSI1N,SRPSI2N, SRpsimm,SRchisk,SRBETAN,
     +         SRD1N,SRD2N
c     COMMON/PSIPR/IPSI,C,H1,H2,H3,XK,D
c
c     IPSI=4
c     K0=1.548D0
c
      do 20 i=1,np+1
      do 10 j=1,np+1
        avts0(i,j)=0.d0
        avts(i,j)=0.d0
   10 continue
   20 CONTINUE
c      avs0=0.d0
c      avs=0.d0
      tmp1=xbar(1)
      en2=DBLE(n)*DBLE(n-np)
      pnrm0=pnorm0(u)
      alfa=2.d0*pnrm0-1.d0
      beta=SRBETAN(u)

      do 500 k=1,n
        y0=y(k)
        z0=y0
        do 220 j=1,np
        x0(j)=X(k,j)
        z0=z0-x0(j)*theta(j)
  220   continue
        z0=z0/sigma
        tmp1=SRpsimm(z0,2,k0)
        do 235 i=1,np
        sc1(i)=tmp1*x0(i)
  235   continue
        sc1(np+1)=SRchisk(z0,k0)
        do 240 i=1,np+1     
        its0(i)=0.d0
        do 230 j=1,np+1
        its0(i)=its0(i)+invm0(i,j)*sc1(j)
  230   continue
  240   continue
        is0=its0(np+1)
c
c        z0=y0
c        do 260 j=1,np
c        x0(j)=X(k,j)
c        z0=z0-x0(j)*theta1(j)
c  260   continue
c        z0=z0/sigma1
        ialf=IALPHAN(Z0,U,SIGMA,IS0)
        tmp1=SRPSI1N(z0,u)
        call SRD1N(U,SIGMA,ITS0,XtX,NP,SA)   
        call SRD2N(U,SIGMA,IS0,TMP2)
        tmp2=tmp2+SRPSI2N(z0,u)-alfa*beta-beta*ialf
        do 265 i=1,np
        sc1(i)=tmp1*x0(i)+sa(i)
  265   continue
        sc1(np+1)=tmp2
        do 280 i=1,np+1 
        its(i)=0.d0
        do 270 j=1,np+1
        its(i)=its(i)+invm1(i,j)*sc1(j)
  270   continue    
  280   continue
c
        do 320 i=1,np+1
        do 300 j=1,i
        avts0(i,j)=avts0(i,j)+its0(i)*its0(j)/en2
        if (i.ne.j) avts0(j,i)=avts0(i,j)
        avts(i,j)=avts(i,j)+its(i)*its(j)/en2
        if (i.ne.j) avts(j,i)=avts(i,j)
  300   continue
  320   continue
  500 continue
      return
      end
C
C********************************************************************************
C
      DOUBLE PRECISION FUNCTION SRRHOW(Z,CONST)
      implicit double precision(a-h,o-z)
      EXTERNAL SRXEXPD
c     COMMON/ZEZCOM/CONST
      SRRHOW=SRXEXPD(Z)-CONST-Z
      RETURN
      END
C
C-----------------------------------------------------------------------
C
      SUBROUTINE SRRGFL2(F,CONST,Y,A,B,TOL,MAXIT,X,ITERM)
      implicit double precision(a-h,o-z)
      EXTERNAL F
c     COMMON/ZEZCOM/CONST
c     DATA TL/1.D-16/
C
C  INITIALIZE
C
      TL=DMIN1(1.D-10,0.1D0*TOL)
      ITR=1
      MESS=0
   10 FA=F(A,CONST)-Y
      FB=F(B,CONST)-Y
C
C  REGULA FALSI ITERATION
C
   20 IF (DABS(FA-FB).GT.TL) GOTO 30
      MESS=MESS+1
      IF (MESS.LE.2) THEN
        A=A/10.D0
        GOTO 10
      ENDIF
C     CALL MESSGE(401,'RGFL2 ',0)
      RETURN
   30 XN=(A*FB-B*FA)/(FB-FA)
      FN=F(XN,CONST)-Y
C
C  TEST TO SEE IF MAXIMUM NUMBER OF ITERATIONS HAS BEEN EXECUTED
C
      IF (ITR.GE.MAXIT) GOTO 60
C
C  TEST TO SEE IF ROOT HAS BEEN FOUND
C
      IF (DABS(FN).LT.TOL) GOTO 70
      IF (FA*FN.LE.0.D0) GOTO 40
      A=XN
      FA=FN
      GOTO 50
   40 B=XN
      FB=FN
C
C  INCREMENT ITERATION COUNTER
C
   50 ITR=ITR+1
      GOTO 20
C
   60 ITERM=2
      X=XN
      RETURN
   70 ITERM=1
      X=XN
      RETURN
      END
C
      SUBROUTINE SRF0W(U,TOL,MAXIT,P)
      implicit double precision(a-h,o-z)
      DOUBLE PRECISION LOW
c     COMMON/ZEZCOM/CONST
      EXTERNAL SRXEXPD,SRRHOW
      P=0.D0
      IF (U.LE.1.D0) RETURN
      P=1.D0
      IF (U.GT.16.D0) RETURN
      CONST=U
      IF (U.GT.1.5D0) THEN
        LOW=-U
        UP=-U+1.5D0
        CALL SRRGFL2(SRRHOW,CONST,0.D0,LOW,UP,TOL,MAXIT,TL,ITRM)
      ELSE
        TLO=TOL
        IF (U-1.D0.LT.1.D-3) TLO=DMIN1(TOL,1.D-8)
        LOW=-U
        UP=0.D0
        CALL SRRGFL2(SRRHOW,CONST,0.D0,LOW,UP,TLO,MAXIT,TL,ITRM)
      ENDIF
      ALOGU=DLOG(U)
      CALL SRRGFL2(SRRHOW,CONST,0.D0,ALOGU,U,TOL,MAXIT,TU,ITRM)
      XU=SRXEXPD(TU)
      CALL SRPWEIBL(1.D0,1.D0,XU,PU)
      XL=SRXEXPD(TL)
      CALL SRPWEIBL(1.D0,1.D0,XL,PL)  
      P=PU-PL
      RETURN
      END
C
      SUBROUTINE SRPWEIBL(ALPHA,SIGMA,X,P)
      implicit double precision(a-h,o-z)
      DATA NCALL,EXMIN,XLGMN,YLGMN/0,0.D0,0.D0,0.D0/
      IF (NCALL.EQ.0) THEN
        NCALL=1
        CALL SRMACHD(3,EXMIN)
        CALL SRMACHD(4,XLGMN)
        CALL SRMACHD(5,YLGMN)
      ENDIF
c     IF (ALPHA.LE.0..OR.SIGMA.LE.0.) CALL MESSGE(500,'PWEIBL',1)
      P=0.D0
      IF (X.LE.0.D0) RETURN
      ALGXS=YLGMN
      XS=X/SIGMA
      IF (XS.GT.XLGMN) ALGXS=DLOG(XS)
      T=ALPHA*ALGXS
      EXPT=0.D0
      IF (T.GT.EXMIN) EXPT=DEXP(T)
      XXPT=0.D0
      IF (-EXPT.GT.EXMIN) XXPT=DEXP(-EXPT)
      P=1.D0-XXPT
      RETURN
      END
C
C======================================================================
C
      DOUBLE PRECISION FUNCTION CHIS1WP(DX,WGT,N,EXU,EXV,NP,PAR)
      implicit double precision(a-h,o-z)
      DIMENSION WGT(N),PAR(NP)
      EXTERNAL EXU,EXV
      CHIS1WP=0.D0 * PAR(1)
      I=IDINT(WGT(1))
      B1=WGT(2)
      ANS=EXU(DX)
      XZV=0.D0
      IF (I.GE.4) THEN
        VV=B1
        ZV=DX/VV
        XZV=EXV(ZV)
      ENDIF
c     GOTO (10,20,30,40,50) I (obsolescent)
      IF (I.EQ.2) GOTO 20
      IF (I.EQ.3) GOTO 30
      IF (I.EQ.4) GOTO 40
      IF (I.EQ.5) GOTO 50
      CHIS1WP=(EXV(DX-B1)-1.D0)*ANS
      RETURN
   20 CHIS1WP=EXV(DX-B1)*ANS
      RETURN
   30 CHIS1WP=DX*(EXV(DX)-1.D0)*ANS
      RETURN
   40 CHIS1WP=ZV*(XZV-1.D0)*ANS
      RETURN
   50 CHIS1WP=-ZV*(XZV-1.D0+XZV*ZV)*ANS/VV
      RETURN
      END
C
C----------------------------------------------------------------------
C
      SUBROUTINE SRINTMW(IWGT,TL,TU,B1,TIL,SUM)
      implicit double precision(a-h,o-z)
      DOUBLE PRECISION SRezez,ERRSTD,WORK,CHIS1WP,LO,HI,
     +       TL,TU,TIL,SUM
      DIMENSION WGT(2),IWORK(80),WORK(320)
      EXTERNAL SRezez,CHIS1WP,SRXEXPD
C
C     INITIALISATION
C
      LIMIT=80
      KEY=1
      LO=TL
      HI=TU
      WGT(1)=DBLE(IWGT)
      WGT(2)=B1
      CALL SRINTGRT(CHIS1WP,WGT,2,SRezez,SRXEXPD,LO,HI,TIL,TIL,
     1            KEY,LIMIT,SUM,ERRSTD,NEVAL,IER,WORK,IWORK,2,WGT)
      RETURN
      END

C======================================================================

      DOUBLE PRECISION FUNCTION SRezez(Z)
      implicit double precision(a-h,o-z)
      DOUBLE PRECISION Z,EXMIN,TMP,VAL 
      DATA NCALL,EXMIN/0,0.D0/
      IF (NCALL.EQ.0) THEN
        NCALL=1
        CALL SRMACHD(3,EXMIN)
      ENDIF
      TMP=Z
      IF (Z.GE.EXMIN) TMP=Z-DEXP(Z)        
      VAL=0.D0
      IF (TMP.GT.EXMIN) VAL=DEXP(TMP) 
      SRezez=VAL
      RETURN
      END
c
      DOUBLE PRECISION FUNCTION SRpezez(Z)
      implicit double precision(a-h,o-z)
      DOUBLE PRECISION Z,EXMIN,TMP,VAL 
      DATA NCALL,EXMIN/0,0.D0/
      IF (NCALL.EQ.0) THEN
        NCALL=1
        CALL SRMACHD(3,EXMIN)
      ENDIF
      TMP=0.D0
      IF (Z.GT.EXMIN) TMP=-DEXP(Z)        
      VAL=0.D0
      IF (TMP.GT.EXMIN) VAL=DEXP(TMP) 
      SRpezez=1.D0-VAL
      RETURN
      END
c
      DOUBLE PRECISION FUNCTION SRPSI1W(Z,L,U)
      DOUBLE PRECISION Z,L,U,EXMIN,TMP 
      DATA NCALL,EXMIN/0,0.D0/
      IF (NCALL.EQ.0) THEN
        NCALL=1
        CALL SRMACHD(3,EXMIN)
      ENDIF
      SRPSI1W=0.D0
      IF (Z.LT.L.OR.Z.GT.U) RETURN 
      TMP=-1.D0
      IF (Z.GT.EXMIN) TMP=DEXP(Z)-1.D0        
      SRPSI1W=TMP
      RETURN
      END
c
      DOUBLE PRECISION FUNCTION SRPSI2W(Z,L,U)
      DOUBLE PRECISION Z,L,U,EXMIN,TMP
      DATA NCALL,EXMIN/0,0.D0/
      IF (NCALL.EQ.0) THEN
        NCALL=1
        CALL SRMACHD(3,EXMIN)
      ENDIF
      SRPSI2W=0.D0
      IF (Z.LT.L.OR.Z.GT.U) RETURN 
      TMP=-Z
      IF (Z.GT.EXMIN) TMP=Z*(DEXP(Z)-1.D0)        
      SRPSI2W=TMP
      RETURN
      END
c
      DOUBLE PRECISION FUNCTION SRBetaw(l,u)
      implicit double precision(a-h,o-z)
      double precision Alfa,l,u,tild,sum,SRpezez
      external SRpezez
      Alfa=SRpezez(u)-SRpezez(l)
      tild=1.D-4
      CALL SRINTMW(3,L,U,0.D0,TILD,SUM)
      SRbetaw=SUM/Alfa
      RETURN
      END
c
      DOUBLE PRECISION FUNCTION SRIALFAW(Z0,L,U,SIGMA,IS0)
      implicit double precision(a-h,o-z)
      DOUBLE PRECISION Z0,L,U,SIGMA,IS0,ETA,RHO,TMP,SREZEZ,SRPEZEZ,EXMIN
      EXTERNAL SREZEZ,SRPEZEZ
      DATA NCALL,EXMIN/0,0.D0/
      IF (NCALL.EQ.0) THEN
        NCALL=1
        CALL SRMACHD(3,EXMIN)
      ENDIF
      ETA=DEXP(U)-U        
      RHO=-Z0
      IF (Z0.GT.EXMIN) RHO=DEXP(Z0)-Z0 
      TMP=(U*SREZEZ(U)-L*SREZEZ(L))*IS0/SIGMA-(SRPEZEZ(U)-SRPEZEZ(L))
      IF (ETA.GE.RHO) TMP=TMP+1.D0
      SRIALFAW=TMP
      RETURN
      END

c 
      SUBROUTINE SRD1W(L,U,SIGMA,IT0,IS0,XtX,XBAR,NP,VAL)
      implicit double precision(a-h,o-z)
      DOUBLE PRECISION L,U,SIGMA,IT0(NP),IS0,XtX(NP,NP),TMP,TMP1,TMP2,
     +       SREZEZ,EZU,VAL(NP),EXMIN,DXPL,XBAR(NP)
      EXTERNAL SREZEZ
      DATA NCALL,EXMIN/0,0.D0/
      IF (NCALL.EQ.0) THEN
        NCALL=1
        CALL SRMACHD(3,EXMIN)
      ENDIF
      DXPL=0.D0
      IF (L.GT.EXMIN) DXPL=DEXP(L)
      TMP1=(DEXP(U)-DXPL)
      TMP2=(U*DEXP(U)-U-L*DXPL+L)*IS0
      EZU=SREZEZ(U) 
      DO 200 I=1,NP
      TMP=0.D0
      DO 100 J=1,NP
      TMP=TMP+XtX(I,J)*IT0(J)
  100 CONTINUE
      VAL(I)=TMP1*TMP+TMP2*XBAR(I)
      VAL(I)=EZU*(VAL(I))/SIGMA  
  200 CONTINUE
      RETURN
      END
c 
      SUBROUTINE SRD2W(L,U,SIGMA,IT0,IS0,XBAR,NP,VAL)
      implicit double precision(a-h,o-z)
      DOUBLE PRECISION L,U,SIGMA,IT0(NP),IS0,XBAR(NP),TMP,TMP1,TMP2,
     +       SREZEZ,EZU,VAL,EXMIN,DXPL
      EXTERNAL SREZEZ
      DATA NCALL,EXMIN/0,0.D0/
      IF (NCALL.EQ.0) THEN
        NCALL=1
        CALL SRMACHD(3,EXMIN)
      ENDIF
      DXPL=0.D0
      IF (L.GT.EXMIN) DXPL=DEXP(L)
      TMP2=(U*U*(DEXP(U)-1.D0)-L*L*(DXPL-1.D0))*IS0
      TMP=U*(DEXP(U)-1.D0)-L*(DXPL-1.D0)
      EZU=SREZEZ(U) 
      TMP1=0.D0
      DO 100 J=1,NP
      TMP1=TMP1+XBAR(J)*IT0(J)
  100 CONTINUE
      VAL=TMP*TMP1
      VAL=EZU*(VAL+TMP2)/SIGMA  
      RETURN
      END
c
      SUBROUTINE av_tmlwf(X,y,n,np,ncov,l,u,xk,theta,sigma,invm0,
     +           invm1,avts0,avts,xbar,XtX,sa,sc1,x0,its0,its)          
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      double precision X(n,np),l,u,xbar(np),y(n),theta(np),XtX(np,np),
     +       invm0(np+1,np+1),invm1(np+1,np+1),x0(np),is0,sa(ncov),
     +       sc1(ncov),avts0(np+1,np+1),avts(np+1,np+1),ialf,
     +       its0(np+1),its(np+1)
      external SRpezez,SRIALFAW,SRPsi1w,SRPsi2w,SRchisk
c     COMMON/PSIPR/IPSI,C,H1,H2,H3,XK,D
c
c
      do 20 i=1,np+1
      do 10 j=1,np+1
        avts0(i,j)=0.d0
        avts(i,j)=0.d0
   10 continue
   20 continue
      en2=DBLE(n)*DBLE(n-np)
 
      alfa=SRpezez(u)-SRpezez(l)
      beta=SRBetaw(l,u)
      do 500 k=1,n
        y0=y(k)
        z0=y0
        do 220 j=1,np
        x0(j)=X(k,j)
        z0=z0-x0(j)*theta(j)
  220   continue
        z0=z0/sigma
        tmp1=SRpsimm(z0,2,xk)
        do 235 i=1,np
        sc1(i)=tmp1*x0(i)
  235   continue
        sc1(np+1)=SRchisk(z0,xk)
        do 240 i=1,np+1      
        its0(i)=0.d0
        do 230 j=1,np+1
        its0(i)=its0(i)+invm0(i,j)*sc1(j)
  230   continue    
  240   continue
        is0=its0(np+1)
        its0(1)=its0(1)+0.1352D0*is0
c
c        z0=y0
c        do 260 j=1,np
c        x0(j)=X(k,j)
c        z0=z0-x0(j)*theta1(j)
c  260   continue
c        z0=z0/sigma1
        ialf=SRIALFAW(Z0,L,U,SIGMA,IS0)
        tmp1=SRPsi1w(z0,l,u)
        call SRD1W(L,U,SIGMA,ITS0,IS0,XtX,XBAR,NP,SA)
        call SRD2W(L,U,SIGMA,ITS0,IS0,XBAR,NP,TMP2)
        tmp2=tmp2+SRPsi2w(z0,l,u)-alfa*beta-beta*ialf
        do 265 i=1,np
        sc1(i)=tmp1*x0(i)+sa(i)
  265   continue
        sc1(np+1)=tmp2
        do 280 i=1,np+1 
        its(i)=0.d0
        do 270 j=1,np+1
        its(i)=its(i)+invm1(i,j)*sc1(j)
  270   continue    
  280   continue
c
        do 320 i=1,np+1
        do 300 j=1,i
        avts0(i,j)=avts0(i,j)+its0(i)*its0(j)/en2
        if (i.ne.j) avts0(j,i)=avts0(i,j)
        avts(i,j)=avts(i,j)+its(i)*its(j)/en2
        if (i.ne.j) avts(j,i)=avts(i,j)
  300   continue
  320   continue
  500 continue
      return
      end

C================= ALREADY INCLUDED IN SPLUS 2000 ======================
      SUBROUTINE SRPSIAMM(N,SVALS,FVALS,IPS,XK)
C.......................................................................
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION SVALS(N),FVALS(N)
C-----------------------------------------------------------------------
C     COMPUTES THE VALUE OF PSI FUNCTION
C     IPS = 1: OPTIMAL FUNCTION
C     IPS = 2: RESCALED BISQUARE FUNCTION
C     IPS = 3: HUBER FUNCTION
C-----------------------------------------------------------------------
      DO 150 I=1,N
         FVALS(I)=SRPSIMM(SVALS(I),IPS,XK)
 150  CONTINUE
      RETURN
      END
C=======================================================================
      SUBROUTINE SRRHOAMM(N,SVALS,FVALS,IPS,XK)
C.......................................................................
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION SVALS(N),FVALS(N)
C-----------------------------------------------------------------------
C     COMPUTES THE VALUE OF RHO FUNCTION
C     IPS = 1: OPTIMAL FUNCTION
C     IPS = 2: RESCALED BISQUARE FUNCTION
C     IPS = 3: HUBER FUNCTION
C-----------------------------------------------------------------------
      DO 150 I=1,N
         FVALS(I)=SRRHOMM(SVALS(I),IPS,XK)
 150  CONTINUE
      RETURN
      END
C=======================================================================
      SUBROUTINE SRPSPAMM(N,SVALS,FVALS,IPS,XK)
C.......................................................................
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION SVALS(N),FVALS(N)
C-----------------------------------------------------------------------
C     COMPUTES THE VALUE OF PSP FUNCTION
C     IPS = 1: OPTIMAL FUNCTION
C     IPS = 2: RESCALED BISQUARE FUNCTION
C     IPS = 3: HUBER FUNCTION
C-----------------------------------------------------------------------
      DO 150 I=1,N
         FVALS(I)=SRPSPMM(SVALS(I),IPS,XK)
 150  CONTINUE
      RETURN
      END
C=======================================================================
      SUBROUTINE SRCHIAMM(N,SVALS,FVALS,IPS,XK)
C.......................................................................
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION SVALS(N),FVALS(N)
C-----------------------------------------------------------------------
C     COMPUTES THE VALUE OF PSP FUNCTION
C     IPS = 1: OPTIMAL FUNCTION
C     IPS = 2: RESCALED BISQUARE FUNCTION
C     IPS = 3: HUBER FUNCTION
C-----------------------------------------------------------------------
      DO 150 I=1,N
         FVALS(I)=SRCHIMM(SVALS(I),IPS,XK)
 150  CONTINUE
      RETURN
      END
C=======================================================================
      FUNCTION SRPSIMM(S,IPS,XK)
C.......................................................................
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
C     COMPUTES THE VALUE OF PSI FUNCTION
C     IPS = 1: OPTIMAL FUNCTION
C     IPS = 2: RESCALED BISQUARE FUNCTION
C     IPS = 3: HUBER FUNCTION
C     IPS = 4: SMOOTHED HUBER FUNCTION
C-----------------------------------------------------------------------
      ABST=DABS(S)
      IF (IPS .EQ. 1) GOTO 100
      IF (IPS .EQ. 2) GOTO 200
      IF (IPS .EQ. 3) GOTO 300
      IF (IPS .EQ. 4) GOTO 400
 100  R1= -1.944D0 
      R2=  1.728D0
      R3= -0.312D0
      R4=  0.016D0
      AX=ABST/XK
      IF (AX .GT. 3.D0) THEN 
         SRPSIMM=0.D0
      ELSE IF(AX .GT. 2.D0) THEN
         AX=S/XK 
         IF (AX .GT. 0.D0) THEN
            SRPSIMM=DMAX1(0.D0,XK*(R4*AX**7+R3*AX**5+R2*AX**3+R1*AX))
         ELSE
            SRPSIMM=-DABS(XK*(R4*AX**7+R3*AX**5+R2*AX**3+R1*AX))
         ENDIF
      ELSE
         SRPSIMM=S
      ENDIF 
      RETURN 
 200  SRPSIMM=0.D0
      IF (ABST .LT. XK) THEN
         SK=S/XK
         SRPSIMM=(6.D0*SK/XK)*(1.D0-SK*SK)*(1.D0-SK*SK)
      ENDIF
      RETURN
 300  SRPSIMM=DMIN1(XK,ABST)
      IF (S .LT. 0.D0) SRPSIMM=-SRPSIMM
      RETURN
 400  IF (ABST .LE. XK) THEN
         SRPSIMM=S
      ELSE
         SRPSIMM=S/ABST*XK*(1.D0+(1.D0-(ABST/XK)**(-3.D0))/3.D0)
      ENDIF
      RETURN
      END
C=======================================================================
      FUNCTION SRRHOMM(S,IPS,XK)
C.......................................................................
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
C     COMPUTES THE VALUE OF RHO FUNCTION
C     IPS = 1: OPTIMAL FUNCTION
C     IPS = 2: RESCALED BISQUARE FUNCTION
C     IPS = 3: HUBER FUNCTION 
C-----------------------------------------------------------------------
      ABST=DABS(S)
      S2=S*S
      IF (IPS .EQ. 1) GOTO 100
      IF (IPS .EQ. 2) GOTO 200
      IF (IPS .EQ. 3 .OR. IPS .EQ. 4) GOTO 300
 100  R1=-1.944D0/2.0D0
      R2= 1.728D0/4.0D0 
      R3=-0.312D0/6.0D0 
      R4= 0.016D0/8.0D0 
      AX=ABST/XK
      IF (AX .GT. 3.D0) THEN
         SRRHOMM=3.25D0*XK*XK
      ELSE IF (AX .GT. 2.D0) THEN
         SRRHOMM=XK*XK*(R1*AX**2+R2*AX**4+R3*AX**6+R4*AX**8+1.792D0)
      ELSE
         SRRHOMM=S2/2.D0
      ENDIF
      RETURN
 200  SRRHOMM=1.D0
      IF (ABST .LT. XK) THEN
         S2=S2/(XK**2)
         SRRHOMM=(S2*(S2-3.D0)+3.D0)*S2
      ENDIF
      RETURN
 300  SRRHOMM=S2/2.D0
      IF (ABST .GT. XK) SRRHOMM=XK*(ABST-XK/2.D0)
      RETURN
      END
C=======================================================================
      FUNCTION SRPSPMM(S,IPS,XK)
C.......................................................................
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
C     COMPUTES THE VALUE OF PSP FUNCTION
C     IPS = 1: OPTIMAL FUNCTION
C     IPS = 2: RESCALED BISQUARE FUNCTION
C     IPS = 3: HUBER FUNCTION
C     IPS = 4: SMOOTHED HUBER FUNCTION
C-----------------------------------------------------------------------
      ABST=DABS(S)
      IF (IPS .EQ. 1) GOTO 100
      IF (IPS .EQ. 2) GOTO 200
      IF (IPS .EQ. 3) GOTO 300
      IF (IPS .EQ. 4) GOTO 400
 100  R1= -1.944D0
      R2=  1.728D0
      R3= -0.312D0
      R4=  0.016D0
      AX=ABST/XK
      IF (AX .GT. 3.D0) THEN
         SRPSPMM=0.D0 
      ELSE IF (AX .GT. 2.D0) THEN 
         SRPSPMM=7.D0*R4*AX**6+5.D0*R3*AX**4+3.D0*R2*AX**2+R1
      ELSE 
         SRPSPMM=1.D0
      ENDIF
      RETURN
 200  SRPSPMM=0.D0
      IF (ABST .LT. XK) THEN
         S2=(S/XK)**2
         SRPSPMM=(6.D0/XK)*(1.D0-S2)*(1.D0-5.D0*S2)/XK
      ENDIF
      RETURN
 300  SRPSPMM=0.D0
      IF (ABST .LE. XK) SRPSPMM=1.D0
      RETURN
 400  SRPSPMM=1.D0
      IF (ABST .GT. XK) SRPSPMM=(ABST/XK)**(-3.D0)
      RETURN
      END
C=======================================================================
      FUNCTION SRCHIMM(S,IPS,XK)
C.......................................................................
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C-----------------------------------------------------------------------
C     COMPUTES THE VALUE OF CHI FUNCTION
C     IPS = 1: OPTIMAL FUNCTION
C     IPS = 2: RESCALED BISQUARE FUNCTION
C     IPS = 3: HUBER FUNCTION
C-----------------------------------------------------------------------
      IF (IPS .EQ. 1) GOTO 100
      IF (IPS .EQ. 2) GOTO 200
      IF (IPS .EQ. 3 .OR. IPS .EQ. 4) GOTO 300
 100  R1=-1.944D0/2.0D0
      R2= 1.728D0/4.0D0 
      R3=-0.312D0/6.0D0 
      R4= 0.016D0/8.0D0 
      AX=DABS(S/XK)
      IF (AX .GT. 3.D0) THEN
         SRCHIMM=3.25D0*XK*XK
      ELSE IF (AX .GT. 2.D0) THEN
         SRCHIMM=XK*XK*(R1*AX**2+R2*AX**4+R3*AX**6+R4*AX**8+1.792D0)
      ELSE 
         SRCHIMM=S*S/2.D0
      ENDIF
      RETURN
 200  SRCHIMM=1.D0
      ABST=DABS(S)
      IF (ABST .LT. XK) THEN
         S2=(S/XK)**2
         SRCHIMM=(S2*(S2-3.D0)+3.D0)*S2
      ENDIF
      RETURN
 300  ABST=DABS(S)
      PS=DMIN1(XK,ABST)
      SRCHIMM=PS*PS/2.D0
      RETURN
      END
c
      function pulk(a,n,k,b)
cc  Finds the kth order statistic of an array a of length n<=1000
      dimension a(n),b(n)
      do 15 j=1,n
      b(j)=a(j)
15    continue
      l=1
      lr=n
20    if (l.ge.lr) goto 90
      ax=b(k)
      jnc=l
      j=lr
30    if(jnc.gt.j)goto 80
40    if (b(jnc).ge.ax)goto 50
      jnc=jnc+1
      goto 40
50    if (b(j).le.ax)goto 60
      j=j-1
      goto 50
60    if(jnc.gt.j)goto 70
      buffer=b(jnc)
      b(jnc)=b(j)
      b(j)=buffer
      jnc=jnc+1
      j=j-1
70    goto 30
80    if(j.lt.k)l=jnc
      if(k.lt.jnc)lr=j
      goto 20
90    pulk=b(k)
      return
      end

cc
cc  Time-efficient algorithm for the scale estimator:
cc
cc       Qn = dn * 2.2219 * {|x_i-x_j|; i<j}_(k)
cc
cc  Parameters of the function Qn :
cc       x  : real array containing the observations
cc       n  : number of observations (n >=2)
cc
cc  The function Qn uses the procedures:
cc     whimed(a,iw,n): finds the weighted high median of an array
cc                     a of length n, using the array iw (also of
cc                     length n) with positive integer weights.
cc     pulk(a,n,k) : finds the k-th order statistic of an
cc                   array a of length n
cc
c     function Qn(x,n)   scratch length ns=500
      subroutine Qn(y,n,scale,sv,siw,sw,work,left,right,weight,Q,P,ns)
      dimension y(n),sv(n),sw(n),work(ns)
      integer left(ns),right(ns),weight(ns),Q(ns),P(ns),siw(ns)
      integer h,k,knew,jhelp,nL,nR,sumQ,sumP
      logical found
      h=n/2+1
      k=h*(h-1)/2
c     call sort(x,n,y,Q,P)
      do 20 i=1,n
          left(i)=n-i+2
          right(i)=n
20    continue 
      jhelp=n*(n+1)/2
      knew=k+jhelp
      nL=jhelp
      nR=n*n
      found=.false.
200   continue
      if ( (nR-nL.gt.n).and.(.not.found) ) then
          j=1
          do 30 i=2,n
          if (left(i).le.right(i)) then
              weight(j)=right(i)-left(i)+1
              jhelp=left(i)+weight(j)/2
              work(j)=y(i)-y(n+1-jhelp)
              j=j+1
          endif
30        continue
          trial=whimed(work,weight,j-1,sv,siw,sw)
          j=0
          do 40 i=n,1,-1
45            continue
              if ((j.lt.n).and.((y(i)-y(n-j)).lt.trial)) then
                  j=j+1
                  goto 45
              endif
              P(i)=j
40        continue
          j=n+1
          do 50 i=1,n
55            continue
              if ((y(i)-y(n-j+2)).gt.trial) then
                  j=j-1
                  goto 55
              endif
              Q(i)=j
50        continue
          sumP=0
          sumQ=0
          do 60 i=1,n
              sumP=sumP+P(i)
              sumQ=sumQ+Q(i)-1
60        continue
          if (knew.le.sumP) then
              do 70 i=1,n
                  right(i)=P(i)
70            continue
              nR=sumP
          else
              if (knew.gt.sumQ) then
                  do 80 i=1,n
                      left(i)=Q(i)
80                continue
                  nL=sumQ
              else
                  scale=trial
                  found=.true.
              endif
          endif
      goto 200
      endif
      if (.not.found) then
          j=1
          do 90 i=2,n
          if (left(i).le.right(i)) then
              do 100 jj=left(i),right(i)
                  work(j)=y(i)-y(n-jj+1)
                  j=j+1
100           continue
          endif
90        continue
          scale=pulk(work,j-1,knew-nL,sv)
      endif
      return
      end


cc
cc  Algorithm to compute the weighted high median in O(n) time.
cc
cc  The whimed is defined as the smallest a(j) such that the sum
cc  of the weights of all a(i) <= a(j) is strictly greater than
cc  half of the total weight.
cc
cc  Parameters of this function:
cc        a: real array containing the observations
cc        n: number of observations
cc       iw: array of integer weights of the observations.
cc
cc  This function uses the function pulk.
cc
cc  The size of acand, iwcand must be at least n.
cc
      function whimed(a,iw,n,acand,iwcand,sw)
      dimension a(n),iw(n),sw(n)
      dimension acand(n),iwcand(n)
      integer wtotal,wrest,wleft,wmid,wright
      nn=n
      wtotal=0
      do 20 i=1,nn
          wtotal=wtotal+iw(i)
20    continue
      wrest=0
100   continue
      trial=pulk(a,nn,nn/2+1,sw)
      wleft=0
      wmid=0
      wright=0
      do 30 i=1,nn
          if (a(i).lt.trial) then
              wleft=wleft+iw(i)
          else
              if (a(i).gt.trial) then
                  wright=wright+iw(i)
              else
                  wmid=wmid+iw(i)
              endif
          endif
30    continue
      if ((2*wrest+2*wleft).gt.wtotal) then
          kcand=0
          do 40 i=1,nn
              if (a(i).lt.trial) then
                  kcand=kcand+1
                  acand(kcand)=a(i)
                  iwcand(kcand)=iw(i)
              endif
40        continue
          nn=kcand
      else
          if ((2*wrest+2*wleft+2*wmid).gt.wtotal) then
              whimed=trial
              goto 99            
          else
              kcand=0
              do 50 i=1,nn
                  if(a(i).gt.trial) then
                      kcand=kcand+1
                      acand(kcand)=a(i)
                      iwcand(kcand)=iw(i)
                  endif
50            continue
              nn=kcand
              wrest=wrest+wleft+wmid
          endif
      endif
      do 60 i=1,nn
          a(i)=acand(i)
          iw(i)=iwcand(i)
60    continue
      go to 100
99    return
      end


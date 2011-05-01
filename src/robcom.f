C-----------------------------------------------------------------------
C
C                 R O B E T H  FORTRAN Source 
C
C  File DFCOMN.F  Routines for default values
C
C-----------------------------------------------------------------------
C
      SUBROUTINE DFRPAR(X,N,NP,MDX,ETYPE,UPAR,PSIPAR,
     1                  ITYPW,ITYPE,ISIGMA)
C.......................................................................
C
C   COPYRIGHT  1992  Alfio Marazzi
C
C   AUTHOR : A. RANDRIAMIHARISOA
C.......................................................................
C
C  SET THE PARAMETER VALUES OF THE COMMON/UCVPR/, /WWWPR/
C  AND /PSIPR/ ACCORDING TO THE ESTIMATOR TYPE
C
      CHARACTER*(7) ETYPE,CP*7,CC*1
      REAL X(MDX,NP)
      DOUBLE PRECISION SUMNRM
      LOGICAL NPRCHK
      COMMON/UCVPR/IUCV,A2,B2,CHK,CKW,BB,BT,CW
      COMMON/UCV56/EM,CR,VK,NNP,ENU,V7
      COMMON/ESTIM/IEST
      COMMON/WWWPR/IWWW
      COMMON/PSIPR/IPSI,C,H1,H2,H3,XK,D
      DATA K0,L0,L1,CP/-1,0,0,'       '/
      IF (K0.NE.-1) GOTO 5
      K0=ICHAR('A')
      L0=ICHAR('a')
      L1=L0+25
      IUCV=0
      IWWW=0
      IPSI=0
      A2=0.
      B2=0.
      CHK=0.
      CKW=0.
      BB=0.
      BT=0.
      CW=0.
      EM=0.
      CR=0.
      VK=0.
      NNP=NP
      C=0.
      H1=0.
      H2=0.
      H3=0.
      XK=0.
      D=0.
C
    5 NPRCHK=N.GT.0.AND.NP.GT.0
      IF (.NOT.NPRCHK) CALL MESSGE(500,'DFRPAR',1)
      CP=ETYPE
      DO 7 I=1,7
       CC=CP(I:I)
       J=ICHAR(CC)
       IF ((L0.LE.J).AND.(J.LE.L1)) THEN
         K=J-L0+K0
         CC=CHAR(K)
         CP(I:I)=CC
       ENDIF
    7 CONTINUE
      P=FLOAT(NP)
      C=0.
      IEST=5
      IF (CP(1:3).EQ.'OLS') GOTO 10
      IEST=6
      IF (CP(1:3).EQ.'LAR') GOTO 10
      IEST=1
      IF (CP(1:5).EQ.'HUBER') GOTO 15
      IEST=2
      IF (CP.EQ.'MAL-STD') GOTO 20
      IEST=3
      IF (CP.EQ.'KRA-WEL') GOTO 25
      IEST=4
      IF (CP.EQ.'MAL-HAM') GOTO 30
      IEST=7
      IF (CP.EQ.'HAM-KRA') GOTO 35
      IEST=8
      IF (CP.EQ.'MAL-UNS') GOTO 40
      IEST=9
      IF (CP.EQ.'MAL-TAU') GOTO 45
      IEST=10
      IF (CP.EQ.'SCH-TAU') GOTO 50
      IEST=11
      IF (CP(1:3).EQ.'LMS') GOTO 55
      IEST=12
      IF (CP(1:3).EQ.'LTS') GOTO 60
      IEST=13
      IF (CP(1:1).EQ.'S') GOTO 65
      IEST=14
      IF (CP.EQ.'ROCKE1') GOTO 70
      IEST=15
      IF (CP.EQ.'ROCKE2') GOTO 75
      IEST=0
      CALL MESSGE(500,'DFRPAR',1)
C
C ORDINARY LEAST SQUARES
C
   10 IUCV=0
      IWWW=0
      RETURN
C
C LEAST ABSOLUTE RESIDUALS
C
C     GOTO 10
C
C HUBER CASE
C
   15 IPSI=1
      IF (PSIPAR.LT.0.) THEN
        C=1.345
        D=1.345
      ELSE
        C=PSIPAR
        D=C
      ENDIF
      ITYPE=1
      ISIGMA=1
      IUCV=0
      IWWW=0
      RETURN
C
C MALLOWS STANDARD
C
   20 IWWW=3
   21 IUCV=1
      A2=0.
      IF (UPAR.LE.P) B2=1.05*1.05*P
      IF (UPAR.GT.P) B2=UPAR
      IPSI=1
      IF (PSIPAR.LT.0.) C=1.345
      IF (PSIPAR.GE.0.) C=PSIPAR
      ITYPW=1
      ITYPE=2
      ISIGMA=2
      RETURN
C
C KRASKER-WELSH
C
   25 IUCV=3
      IF (UPAR.LT.0.) UPAR=PSIPAR
      IF (UPAR.LE.SQRT(P)) CKW=1.05*SQRT(P)
      IF (UPAR.GT.SQRT(P)) CKW=UPAR
      IWWW=1
      IPSI=1
      C=CKW
      ITYPW=1
      ITYPE=3
      ISIGMA=2
      RETURN
C
C MALLOWS-HAMPEL
C
   30 IWWW=2
      GOTO 21
C
C HAMPEL-KRASKER
C
   35 IUCV=2
      IF (UPAR.LT.0.) UPAR=PSIPAR
      IF (UPAR.GE.0.) THEN
        CHK=UPAR
      ELSE
        SUMNRM=0.D0
        DO 36 I=1,N
        CALL NRM2(X(I,1),N,MDX,MDX*(NP-1)+1,XNRM)
        SUMNRM=SUMNRM+DBLE(XNRM)
   36   CONTINUE
        CHK=1.05*FLOAT(NP)*SQRT(1.5707963)
        CHK=CHK/(SNGL(SUMNRM)/FLOAT(N))
      ENDIF
      IWWW=1
      C=CHK
      IPSI=1
      ITYPW=2
      ITYPE=3
      ISIGMA=2
      RETURN
C
C MALLOWS : UNSTANDARDIZED CASE
C
   40 IUCV=4
      IF (UPAR.GE.0.) THEN
        BB=UPAR
      ELSE
        SUMNRM=0.D0
        DO 41 I=1,N
        CALL NRM2(X(I,1),N,MDX,MDX*(NP-1)+1,XNRM)
        SUMNRM=SUMNRM+DBLE(XNRM)
   41   CONTINUE
        BB=1.05*FLOAT(NP)
        BB=BB/(SNGL(SUMNRM)/FLOAT(N))
      ENDIF
      IWWW=2
      IPSI=1
      IF (PSIPAR.LT.0.) C=1.345
      IF (PSIPAR.GE.0.) C=PSIPAR
      ITYPW=2
      ITYPE=3
      ISIGMA=2
      RETURN
C
C MALLOWS ESTIMATOR IN THE TAU-TEST
C
   45 IUCV=4
      IF (UPAR.LE.0.) BB=9.999
      IF (UPAR.GT.0.) BB=UPAR
      IWWW=2
      IPSI=1
      IF (PSIPAR.LT.0.) C=1.345
      IF (PSIPAR.GE.0.) C=PSIPAR
      ITYPW=1
      ITYPE=2
      ISIGMA=2
      RETURN
C
C SCHWEPPE ESTIMATOR IN THE TAU-TEST
C
   50 IUCV=2
      IF (UPAR.LT.0.) UPAR=PSIPAR
      IF (UPAR.LE.0.) CHK=9.999
      IF (UPAR.GT.0.) CHK=UPAR
      IWWW=1
      IPSI=1
      C=CHK
      ITYPW=1
      ITYPE=3
      ISIGMA=2
      RETURN
C
C LMS-ESTIMATOR
C
   55 RETURN
C
C LTS-ESTIMATOR
C
   60 RETURN
C
C S-ESTIMATOR
C
   65 IPSI=4
      ITYPE=1
      IF (PSIPAR.LT.0.) XK=1.548
      IF (PSIPAR.GT.0.) XK=PSIPAR
      RETURN
C
C ROCKE ESTIMATOR 
C
   70 IUCV=5
      GOTO 76
   75 IUCV=6
   76 EM=PSIPAR
      IF (PSIPAR.LE.0.) EM=1.345
      CR=2
      IF (UPAR.GT.0.) CR=UPAR
      IWWW=2
      RETURN
      END
C
C-----------------------------------------------------------------------
C
      SUBROUTINE RPARDF(X,N,NP,MDX,RTYPE,UPAR,PSIPAR,
     1                  ITYPW,ITYPE,ISIGMA)
C.......................................................................
C
C   COPYRIGHT  1992  Alfio Marazzi
C
C   AUTHOR : A. RANDRIAMIHARISOA
C.......................................................................
C
C  SET THE PARAMETER VALUES OF THE COMMON/UCVPR/, /WWWPR/
C  AND /PSIPR/ ACCORDING TO THE ESTIMATOR TYPE
C
      INTEGER RTYPE
      REAL X(MDX,NP)
      DOUBLE PRECISION SUMNRM
      LOGICAL NPRCHK
      COMMON/UCVPR/IUCV,A2,B2,CHK,CKW,BB,BT,CW
      COMMON/UCV56/EM,CR,VK,NNP,ENU,V7
      COMMON/ESTIM/IEST
      COMMON/WWWPR/IWWW
      COMMON/PSIPR/IPSI,C,H1,H2,H3,XK,D
      IUCV=0
      IWWW=0
      IPSI=0
      A2=0.
      B2=0.
      CHK=0.
      CKW=0.
      BB=0.
      BT=0.
      CW=0.
      EM=0.
      CR=0.
      VK=0.
      NNP=NP
      C=0.
      H1=0.
      H2=0.
      H3=0.
      XK=0.
      D=0.
C
      NPRCHK=N.GT.0.AND.NP.GT.0
      IF (.NOT.NPRCHK) CALL MESSGE(500,'DFRPAR',1)
      P=FLOAT(NP)
      C=0.
      IEST=5
      IF (RTYPE.EQ.5) GOTO 10
      IEST=6
      IF (RTYPE.EQ.6) GOTO 10
      IEST=1
      IF (RTYPE.EQ.1) GOTO 15
      IEST=2
      IF (RTYPE.EQ.2) GOTO 20
      IEST=3
      IF (RTYPE.EQ.3) GOTO 25
      IEST=4
      IF (RTYPE.EQ.4) GOTO 30
      IEST=7
      IF (RTYPE.EQ.7) GOTO 35
      IEST=8
      IF (RTYPE.EQ.8) GOTO 40
      IEST=9
      IF (RTYPE.EQ.9) GOTO 45
      IEST=10
      IF (RTYPE.EQ.10) GOTO 50
      IEST=11
      IF (RTYPE.EQ.11) GOTO 55
      IEST=12
      IF (RTYPE.EQ.12) GOTO 60
      IEST=13
      IF (RTYPE.EQ.13) GOTO 65
      IEST=14
      IF (RTYPE.EQ.14) GOTO 70
      IEST=15
      IF (RTYPE.EQ.15) GOTO 75
      IEST=0
      CALL MESSGE(500,'DFRPAR',1)
C
C ORDINARY LEAST SQUARES
C
   10 IUCV=0
      IWWW=0
      RETURN
C
C LEAST ABSOLUTE RESIDUALS
C
C     GOTO 10
C
C HUBER CASE
C
   15 IPSI=1
      IF (PSIPAR.LT.0.) THEN
        C=1.345
        D=1.345
      ELSE
        C=PSIPAR
        D=C
      ENDIF
      ITYPE=1
      ISIGMA=1
      IUCV=0
      IWWW=0
      RETURN
C
C MALLOWS STANDARD
C
   20 IWWW=3
   21 IUCV=1
      A2=0.
      IF (UPAR.LE.P) B2=1.05*1.05*P
      IF (UPAR.GT.P) B2=UPAR
      IPSI=1
      IF (PSIPAR.LT.0.) C=1.345
      IF (PSIPAR.GE.0.) C=PSIPAR
      ITYPW=1
      ITYPE=2
      ISIGMA=2
      RETURN
C
C KRASKER-WELSH
C
   25 IUCV=3
      IF (UPAR.LT.0.) UPAR=PSIPAR
      IF (UPAR.LE.SQRT(P)) CKW=1.05*SQRT(P)
      IF (UPAR.GT.SQRT(P)) CKW=UPAR
      IWWW=1
      IPSI=1
      C=CKW
      ITYPW=1
      ITYPE=3
      ISIGMA=2
      RETURN
C
C MALLOWS-HAMPEL
C
   30 IWWW=2
      GOTO 21
C
C HAMPEL-KRASKER
C
   35 IUCV=2
      IF (UPAR.LT.0.) UPAR=PSIPAR
      IF (UPAR.GE.0.) THEN
        CHK=UPAR
      ELSE
        SUMNRM=0.D0
        DO 36 I=1,N
        CALL NRM2(X(I,1),N,MDX,MDX*(NP-1)+1,XNRM)
        SUMNRM=SUMNRM+DBLE(XNRM)
   36   CONTINUE
        CHK=1.05*FLOAT(NP)*SQRT(1.5707963)
        CHK=CHK/(SNGL(SUMNRM)/FLOAT(N))
      ENDIF
      IWWW=1
      C=CHK
      IPSI=1
      ITYPW=2
      ITYPE=3
      ISIGMA=2
      RETURN
C
C MALLOWS : UNSTANDARDIZED CASE
C
   40 IUCV=4
      IF (UPAR.GE.0.) THEN
        BB=UPAR
      ELSE
        SUMNRM=0.D0
        DO 41 I=1,N
        CALL NRM2(X(I,1),N,MDX,MDX*(NP-1)+1,XNRM)
        SUMNRM=SUMNRM+DBLE(XNRM)
   41   CONTINUE
        BB=1.05*FLOAT(NP)
        BB=BB/(SNGL(SUMNRM)/FLOAT(N))
      ENDIF
      IWWW=2
      IPSI=1
      IF (PSIPAR.LT.0.) C=1.345
      IF (PSIPAR.GE.0.) C=PSIPAR
      ITYPW=2
      ITYPE=3
      ISIGMA=2
      RETURN
C
C MALLOWS ESTIMATOR IN THE TAU-TEST
C
   45 IUCV=4
      IF (UPAR.LE.0.) BB=9.999
      IF (UPAR.GT.0.) BB=UPAR
      IWWW=2
      IPSI=1
      IF (PSIPAR.LT.0.) C=1.345
      IF (PSIPAR.GE.0.) C=PSIPAR
      ITYPW=1
      ITYPE=2
      ISIGMA=2
      RETURN
C
C SCHWEPPE ESTIMATOR IN THE TAU-TEST
C
   50 IUCV=2
      IF (UPAR.LT.0.) UPAR=PSIPAR
      IF (UPAR.LE.0.) CHK=9.999
      IF (UPAR.GT.0.) CHK=UPAR
      IWWW=1
      IPSI=1
      C=CHK
      ITYPW=1
      ITYPE=3
      ISIGMA=2
      RETURN
C
C LMS-ESTIMATOR
C
   55 RETURN
C
C LTS-ESTIMATOR
C
   60 RETURN
C
C S-ESTIMATOR
C
   65 IPSI=4
      ITYPE=1
      IF (PSIPAR.LT.0.) XK=1.548
      IF (PSIPAR.GT.0.) XK=PSIPAR
      RETURN
C
C ROCKE ESTIMATOR 
C
   70 IUCV=5
      GOTO 76
   75 IUCV=6
   76 EM=PSIPAR
      IF (PSIPAR.LE.0.) EM=1.345
      CR=2
      IF (UPAR.GT.0.) CR=UPAR
      IWWW=2
      RETURN
      END
C
C-----------------------------------------------------------------------
C
      SUBROUTINE DFCOMN(IPSI,C,H1,H2,H3,XK,D,BTA,BT0,IUCV,A2,B2,CHK,CKW,
     +                  BB,BT,CW,EM,CR,VK,NP,ENU,V7,IWWW)
C.......................................................................
C
C   COPYRIGHT  1992  Alfio Marazzi
C
C   AUTHOR : A. RANDRIAMIHARISOA
C.......................................................................
C
      COMMON/UCVPR/JUCV,AA,AB,PHK,PKW,PBB,PBT,PW
      COMMON/UCV56/PM,PCR,PK,NNP,PNU,P7
      COMMON/PSIPR/JPSI,PC,PH1,PH2,PH3,PXK,PD
      COMMON/WWWPR/JWWW
      COMMON/BETA/BETA,BET0
      IF (IPSI.GE.-5) JPSI=IPSI
      IF (C .GE.0.)  PC=C
      IF (H1.GE.0.)  PH1=H1
      IF (IPSI.EQ.10) PH1=H1
      IF (H2.GE.0.)  PH2=H2
      IF (H3.GE.0.)  PH3=H3
      IF (XK.GE.0.)  PXK=XK
      IF (D .GE.0.)  PD=D
      IF (BTA.GE.0.) BETA=BTA
      IF (BT0.GE.0.) BET0=BT0
      IF (IUCV.GE.0) JUCV=IUCV
      IF (A2.GE.0.)  AA=A2
      IF (B2.GE.0.)  AB=B2
      IF (CHK.GE.0.) PHK=CHK
      IF (CKW.GE.0.) PKW=CKW
      IF (BB.GE.0.)  PBB=BB
      IF (BT.GE.0.)  PBT=BT
      IF (CW.GE.0.)  PW=CW
      IF (EM.GT.0.)  PM=EM
      IF (CR.GT.0.)  PCR=CR
      IF (VK.GT.0.)  PK=VK
      IF (NP.GT.0)   NNP=NP
      IF (ENU.GT.0.) PNU=ENU
      IF (V7.GT.0.)  P7=V7
      IF (IWWW.GE.0) JWWW=IWWW
      RETURN
      END
C
C-----------------------------------------------------------------------
C
      SUBROUTINE COMVAL(IPSI,C,H1,H2,H3,XK,D,BTA,BT0,IUCV,A2,B2,CHK,
     +                  CKW,BB,BT,CW,EM,CR,VK,NP,ENU,V7,IWWW)
C.......................................................................
C
C   COPYRIGHT  1992  Alfio Marazzi
C
C   AUTHOR : A. RANDRIAMIHARISOA
C.......................................................................
C
      COMMON/UCVPR/JUCV,AA,AB,PHK,PKW,PBB,PBT,PW
      COMMON/UCV56/PM,PCR,PK,NNP,PNU,P7
      COMMON/PSIPR/JPSI,PC,PH1,PH2,PH3,PXK,PD
      COMMON/WWWPR/JWWW
      COMMON/BETA/BETA,BET0
      IPSI=JPSI
      C=PC
      H1=PH1
      H2=PH2
      H3=PH3
      XK=PXK
      D=PD
      BTA=BETA
      BT0=BET0
      IUCV=JUCV
      A2=AA
      B2=AB
      CHK=PHK
      CKW=PKW
      BB=PBB
      BT=PBT
      CW=PW
      EM=PM
      CR=PCR
      VK=PK
      NP=NNP
      ENU=PNU
      V7=P7
      IWWW=JWWW
      RETURN
      END
C
C-----------------------------------------------------------------------
C
      SUBROUTINE CHISQ(KODE,IFN,X,P)
C.......................................................................
C
C   AUTHORS :     I.D. HILL AND M.C. PIKE (1967)
C                 ALGORITHM 299: CHI-SQUARED INTEGRAL
C                 COMMUNICATION OF THE ACM, VOL.10, PP.243-244.
C                 ADAPTED FOR ROBETH BY A. RANDRIAMIHARISOA
C.......................................................................
C
C
C REMARK:   IF (X.LE.0.OR.IFN.LT.1).AND.(KODE.EQ.1) P=0. + MESSAGE 400
C ------    IF (X.LE.0.OR.IFN.LT.1).AND.(KODE.EQ.2) P=1. + MESSAGE 400
C
      LOGICAL EVEN,BIGX,ODD,SMLX
      DATA XLSPI,YLSPI/0.572364942925,0.564189583548/
C
      IF (KODE.NE.1.AND.KODE.NE.2) CALL MESSGE(500,'CHISQ ',1)
      S=1.
      FN=FLOAT(IFN)
      IF (X.GT.0..AND.FN.GE.1.) GOTO 5
      CALL MESSGE(400,'CHISQ ',0)
      GOTO 99
    5 NU=IFIX(FN+.5)
      CALL MACH(3,EXMIN)
      A=0.5*X
      BIGX=.FALSE.
      IF (-A.LE.EXMIN) BIGX=.TRUE.
      SMLX=.NOT.BIGX
      EVEN=(2*(NU/2).EQ.NU)
      ODD=.NOT.EVEN
      IF ((EVEN.OR.NU.GT.2).AND.SMLX) S=EXP(-A)
      IF (BIGX) S=0.
      Y=S
      IF (EVEN) GOTO 10
      SX=-SQRT(X)
      CALL GAUSS(1,SX,ANS)
      S=2.0*ANS
C
C  NU.LE.2
C
   10 IF (NU.LE.2) GOTO 99
C
C  NU.GT.2
C
      X1=0.5*(FN-1.0)
      IF (EVEN) Z=1.0
      IF (ODD ) Z=0.5
      IF (SMLX) GOTO 30
      IF (EVEN) E=0.0
      IF (ODD ) E=XLSPI
      C=ALOG(A)
   20 E=ALOG(Z)+E
      IF (C*Z-A-E.GT.EXMIN) S=EXP(C*Z-A-E)+S
      Z=Z+1.0
      IF (Z.LE.X1) GOTO 20
      GOTO 99
   30 IF (EVEN) E=1.0
      IF (ODD ) E=YLSPI/SQRT(A)
      C=0.0
   40 E=E*A/Z
      C=C+E
      Z=Z+1.0
      IF (Z.LE.X1) GOTO 40
      S=C*Y+S
   99 P=S
      IF (KODE.EQ.1) P=1.0-P
      RETURN
      END
C
C-----------------------------------------------------------------------
C
      SUBROUTINE GAUSS  (KODE,X,P)
C.......................................................................
C
C   COPYRIGHT 1992 Alfio Marazzi
C
C   AUTHOR : A. RANDRIAMIHARISOA
C                
C.......................................................................
C
      REAL               P,X,SQR1D2
      DATA               SQR1D2/.7071068/
C
      IF (KODE.NE.1.AND.KODE.NE.2) CALL MESSGE(500,'GAUSS ',1)
      CALL CERF(-X*SQR1D2,C)
      P = .5 * C
      IF (KODE.EQ.2) P=1.-P
      RETURN
      END
C
C-----------------------------------------------------------------------
C
      SUBROUTINE LMDD(X,Y,N,ISORT,XME,XMD,XSD)
C.......................................................................
C
C   COPYRIGHT 1992 Alfio Marazzi
C
C   AUTHORS : W. STAHEL / A. MARAZZI
C.......................................................................
C
      REAL X(N),Y(N)
C
      KM=(N+1)/2
      DO 20 I=1,N
   20 Y(I)=X(I)
      IF (ISORT.NE.0) CALL SRT1(Y,N,1,N)
      XME=Y(KM)
      IF (KM*2.EQ.N) XME=(XME+Y(KM+1))/2.
      K=0
      K1=KM
      K2=KM
      X1=0.
      X2=0.
   30 IF (K.GE.KM) GOTO 50
      K=K+1
      IF (X1.GT.X2) GOTO 40
      K1=K1-1
      IF (K1.EQ.0) GOTO 50
      X1=XME-Y(K1)
      GOTO 30
   40 K2=K2+1
      IF (K2.GT.N) GOTO 50
      X2=Y(K2)-XME
      GOTO 30
   50 XMD=AMIN1(X1,X2)
      XSD=XMD/.6745
      RETURN
      END
C
C-----------------------------------------------------------------------
C
      SUBROUTINE QRSS(RS,WGT,WGT2,EXRHO,N,ITYPE,SIGMA,CONST,QR)
C.......................................................................
C
C   COPYRIGHT 1992 Alfio Marazzi
C
C   AUTHOR : A. MARAZZI
C.......................................................................
C
      DIMENSION RS(N),WGT(N),WGT2(N)
      DOUBLE PRECISION TMP
      EXTERNAL EXRHO
      TMP=0.D0
      IF (ITYPE.NE.1) GOTO 15
C
C  HUBER-TYPE
C
      DO 10 I=1,N
        S=RS(I)/SIGMA
        TMP=TMP+DBLE(EXRHO(S))
   10 CONTINUE
      GOTO 50
C
C MALLOWS-TYPE
C
   15 IF (ITYPE.NE.2) GOTO 30
      DO 20 I=1,N
      IF (WGT(I).EQ.0..OR.WGT(I).EQ.-1.) GOTO 20
      S=RS(I)/SIGMA
      TMP=TMP+EXRHO(S)*DBLE(WGT(I))
   20 CONTINUE
      GOTO 50
C
C SCHWEPPE-TYPE
C
   30 DO 40 I=1,N
      IF (WGT(I).EQ.0..OR.WGT(I).EQ.-1.) GOTO 40
      S=RS(I)/(SIGMA*WGT(I))
      TMP=TMP+EXRHO(S)*DBLE(WGT2(I))
   40 CONTINUE
   50 QR=(SNGL(TMP)+CONST)*SIGMA
      RETURN
      END
C
C-----------------------------------------------------------------------
C
      SUBROUTINE NEWSIG(RS,WGT,WGT2,SIGMA,SIGMB,N,ITYPE,EXCHI)
C.......................................................................
C
C   COPYRIGHT 1992 Alfio Marazzi
C
C   AUTHORS : A. MARAZZI / A. RANDRIAMIHARISOA
C.......................................................................
C
C  PURPOSE
C  -------
C  COMPUTES A NEW VALUE SIGMB FOR THE ROBUST ESTIMATE OF THE
C  ERROR STANDARD DEVIATION IN THE HUBER'S ALGORITHM FOR REGRESSION.
C  NEWSIG CALLS THE FUNCTION EXCHI.
C
      DIMENSION RS(N),WGT(N),WGT2(N)
      EXTERNAL EXCHI
      COMMON/CONST/CONST
      TMP=0.0
      IF (ITYPE.NE.1) GOTO 20
C
C  HUBER-TYPE
C
      DO 10 I=1,N
      S=RS(I)/SIGMA
   10 TMP=TMP+EXCHI(S)
      GOTO 90
C
C  MALLOWS-TYPE
   20 IF (ITYPE.NE.2) GOTO 40
      DO 30 I=1,N
      S=RS(I)/SIGMA
      IF (WGT(I).LE.0.) GOTO 30
      TMP=TMP+EXCHI(S)*WGT(I)
   30 CONTINUE
      GOTO 90
C
C  SCHWEPPE-TYPE
C
   40 DO 50 I=1,N
      SW=SIGMA*WGT(I)
      IF (SW.EQ.0..OR.WGT(I).LE.0.) GOTO 50
      S=RS(I)/SW
      TMP=TMP+EXCHI(S)*WGT2(I)
   50 CONTINUE
   90 SIGMB=SQRT(TMP/CONST)*SIGMA
      RETURN
      END
C
C----------------------------------------------------------------------
C
      FUNCTION ICSIGM(SIGMA,SIGMB,TOL)
C.......................................................................
C
C   COPYRIGHT 1992 Alfio Marazzi
C
C   AUTHOR : A. RANDRIAMIHARISOA
C.......................................................................
C
      DS=ABS(SIGMA-SIGMB)/AMAX1(1.,SIGMA)
      ICSIGM=0
      IF (TOL.GE.DS) ICSIGM=1
      RETURN
      END
C
C-----------------------------------------------------------------------
C
      FUNCTION ICTHET(NP,NCOV,DELTA,SIGMA,S,TOL,ICNV)
C.......................................................................
C
C   COPYRIGHT 1992 Alfio Marazzi
C
C   AUTHOR : A. RANDRIAMIHARISOA
C.......................................................................
C
      DIMENSION DELTA(NP),S(NCOV)
      ICTHET=0
      TOL1=TOL*SIGMA
      IF (ICNV.EQ.2) GOTO 200
      IF (ICNV.EQ.3) GOTO 300
      L=0
      DO 100 J=1,NP
      L=L+J
      TOL2=TOL1*SQRT(S(L))
  100 IF (TOL2.LT.ABS(DELTA(J))) RETURN
      GOTO 500
  200 CALL XSY(DELTA,DELTA,S,NP,NCOV,TOL2)
      TOL2=SQRT(TOL2)
      IF (TOL1.GE.TOL2) ICTHET=1
      RETURN
  300 L=0
      DO 350 J=1,NP
      L=L+J
      TOL2=ABS(DELTA(J))*SQRT(S(L))
  350 IF (TOL1.LT.TOL2) RETURN
  500 ICTHET=1
      RETURN
      END
C
C-----------------------------------------------------------------------
C
      SUBROUTINE RYSIGM(RS,WGT,EXCHI,SIGMAI,N,NP,TOL,ITYPE,ISIGMA,
     1                  MAXIS,NIT,SIGMAF,SW,SC)
C.......................................................................
C
C   COPYRIGHT 1992 Alfio Marazzi
C
C   AUTHORS: A. MARAZZI / A. RANDRIAMIHARISOA
C.......................................................................
C
      DIMENSION RS(N),WGT(N),SW(N),SC(N)
      LOGICAL NPRCHK
      EXTERNAL EXCHI,ICSIGM
      COMMON/BETA/BETA,BET0
      COMMON/CONST/CONST
      DATA TL/1.E-10/
C
C  PARAMETER CHECK AND INITIALIZATION
C
      N0=N
      SIGMB=SIGMAI
      IASG=IABS(ISIGMA)
      NPRCHK=NP.GT.0.AND.N.GT.0.AND.(ITYPE.GE.1.AND.ITYPE.LE.3)
     1       .AND.((IASG.EQ.1.AND.MAXIS.GT.0.AND.TOL.GT.0..AND.
     2       SIGMAI.GT.0.).OR.IASG.EQ.2)
      IF (.NOT.NPRCHK) CALL MESSGE(500,'RYSIGM',1)
      ITYP=ITYPE
      IF (ITYP.EQ.1) GOTO 20
      IF (SIGMAI.EQ.SIGMAF) GOTO 20
      E=2.0
      IF (ITYP.EQ.2) E=0.5
      DO 10 I=1,N
      IF (WGT(I).LE.0.) THEN
        SW(I)=-1.
        N0=N0-1
      ELSE
        SW(I)=WGT(I)**E
      ENDIF
   10 CONTINUE
      IF (N0.EQ.0) ITYP=1
   20 CONTINUE
      IF (IASG.EQ.2) GOTO 500
      CONST=BETA*FLOAT(N-NP)
C
C  STEP 1. SET NIT := 1
C  -------
      NIT=1
C
C  STEP 2. COMPUTE A NEW VALUE SIGMB FOR SIGMA
C  -------
  100 SIGMA=SIGMB
      CALL NEWSIG(RS,WGT,SW,SIGMA,SIGMB,N,ITYP,EXCHI)
      IF (SIGMB.GT.TL) GOTO 300
      CALL MESSGE(460,'RYSIGM',0)
      RETURN
C
C  STEP 3. STOP ITERATIONS IF DESIRED PRECISION IS REACHED
C  -------
  300 IF (ICSIGM(SIGMA,SIGMB,TOL).EQ.1.OR.NIT.EQ.MAXIS) GOTO 400
      NIT=NIT+1
      GOTO 100
  400 SIGMAF=SIGMB
      RETURN
C
C COMPUTE SIGMA USING MEDIAN
C --------------------------
  500 IF (ITYPE.NE.1) GOTO 650
C
C  HUBER-TYPE
C
      DO 600 I=1,N
        SC(I)=ABS(RS(I))
  600 CONTINUE
      N0=N
      GOTO 900
C
C  MALLOWS
C
  650 IF (ITYPE.NE.2) GOTO 750
      N0=0
      DO 700 I=1,N
        IF (SW(I).LE.0.) GOTO 700
        N0=N0+1
        SC(N0)=ABS(RS(I))*SW(I)
  700 CONTINUE
      GOTO 900
C
C  SCHWEPPE-TYPE
C
  750 N0=0
      DO 800 I=1,N
        IF (WGT(I).EQ.0.) GOTO 800
        N0=N0+1
        SC(N0)=ABS(RS(I))
  800 CONTINUE
  900 MED=(N0/2)+1
      CALL FSTORD(SC,N0,MED,SIGMAF)
      SIGMAF=SIGMAF/BET0
      RETURN
      END 
C
C----------------------------------------------------------------------
C
      SUBROUTINE RANDOW(ISEED,RN)
C.......................................................................
C
C   COPYRIGHT 1992 Alfio Marazzi
C
C   AUTHOR : J. JOSS 
C.......................................................................
C
C  RANDOM NUMBER GENERATOR ACCORDING TO THE LINEAR CONGRUENT SCHEME
C                  ISEED=ISEED*5761+999 MODULO 65536
C  IMPROVED AFTER MACLAREN-MARSAGLIA
C
      DIMENSION T(128)
      DATA INIT,T/0,128*0./
      IF (INIT.EQ.0.OR.INIT.NE.ISEED) THEN
        ISEED=MOD(ISEED,65536)
        DO 100 I=1,128
        ISEED=ISEED*5761+999
        ISEED=MOD(ISEED,65536)
  100   T(I)=FLOAT(ISEED)/65536.0
      ENDIF
      ISEED=ISEED*5761+999
      ISEED=MOD(ISEED,65536)
      I=128*ISEED/65536
      RN=T(I+1)
      ISEED=ISEED*5761+999
      ISEED=MOD(ISEED,65536)
      T(I+1)=FLOAT(ISEED)/65536.0
      INIT=ISEED
      RETURN
      END
C
C----------------------------------------------------------------------
C
      SUBROUTINE FSTORD(Y,N,J,YJ)
C.......................................................................
C
C   AUTHOR :     P.J. ROUSSEEUW & A.M. LEROY
C                PROGRESS PACKAGE (SUBROUTINE PULL)
C                ADAPTED FOR ROBETH BY J. JOSS / A. RANDRIAMIHARISOA
C.......................................................................
C
C  FSTORD SEARCHES THE J-TH VALUE IN ORDER OF MAGNITUDE IN A VECTOR
C  OF LENGTH N.
C
      DIMENSION Y(N)
      IF (J.LE.0.OR.J.GT.N) CALL MESSGE(500,'FSTORD',1)
      L=1
      LR=N
   20 IF (L.GE.LR) GOTO 90
      AX=Y(J)
      JNC=L
      JJ=LR
   30 IF(JNC.GT.JJ) GOTO 80
   40 IF (Y(JNC).GE.AX) GOTO 50
      JNC=JNC+1
      GOTO 40
   50 IF(Y(JJ).LE.AX) GOTO 60
      JJ=JJ-1
      GOTO 50
   60 IF(JNC.GT.JJ) GOTO 70
      WA=Y(JNC)
      Y(JNC)=Y(JJ)
      Y(JJ)=WA
      JNC=JNC+1
      JJ=JJ-1
   70 GOTO 30
   80 IF(JJ.LT.J) L=JNC
      IF(J.LT.JNC) LR=JJ
      GOTO 20
   90 YJ=Y(J)
      RETURN
      END
C***********************************************************************
C**************************** H B A U X I ******************************
C
      FUNCTION ICNREP(N,NP,IOPT,IMODE)
C.......................................................................
C
C   R O B E T H  -  R O B S Y S   RELEASE 3.0 (COPYRIGHT) 1985, 1990
C
C   PROGRAMMER : J. JOSS
C.......................................................................
C
C  COMPUTE NUMBER OF REPETITIONS IN RYLMSR
C  M  NUMBER OF OBSERVATIONS
C  NP NUMBER OF PARAMETERS
C  IOPT  0  QUICK VERSION
C        1  EXTENDED VERSION
C        2  NOT USED
C        3  ALL COMBINATIONS
C
      DIMENSION NREPQ(8),NREPE(5)
      DATA NREPQ/150,300,400,500,600,700,850,1250/
      DATA NREPE/500,1000,1500,2000,2500/
      GOTO (1,2,3,4) IOPT+1
    1 IF(NP .GE. 9) THEN
         ICNREP=1500
      ELSE
         ICNREP=NREPQ(NP)
      ENDIF
      RETURN
    2 IF(NP .GE. 6) THEN
         ICNREP=3000
      ELSE
         ICNREP=NREPE(NP)
      ENDIF
    3 RETURN
    4 NN=N
      NR=1
      DO 10 I=1,NP
         NR=(NR*NN)/I
   10    NN=NN-1
      IF (IMODE.GE.3) NR=NR*2**(NP-1)
      ICNREP=NR
      RETURN
      END
C
C-----------------------------------------------------------------------
C
      SUBROUTINE NCOMB(N,NP,IT)
C.......................................................................
C
C   R O B E T H  -  R O B S Y S   RELEASE 3.0 (COPYRIGHT) 1985, 1990
C
C   PROGRAMMER : J .JOSS
C.......................................................................
C
C  COMPUTE ALL COMBINATIONS FOR RESAMPLING ALGORITHM
C
      DIMENSION IT(NP)
      IN=NP
   10 IT(IN)=IT(IN)+1
      IF(IT(IN).GT.N-NP+IN) THEN
         IN=IN-1
         GOTO 10
      ENDIF
      IF(IN.NE.NP) THEN
         DO 20 I=IN+1,NP
   20    IT(I)=IT(I-1)+1
      ENDIF
      RETURN
      END
C
C----------------------------------------------------------------------
C
      SUBROUTINE RICLL1(XT,Y,N,NP,MDXT,THETA,SH,SP)
      DIMENSION XT(MDXT,NP),Y(N),THETA(MDXT),SH(NP)
      INTEGER SP(NP)
C.......................................................................
C
C   R O B E T H  -  R O B S Y S   RELEASE 3.0 (COPYRIGHT) 1985, 1990
C
C   PROGRAMMER : J. JOSS
C.......................................................................
C
C  HOUSHOLDER TRANSFORMATION OF THE RIGHT SIDE
C
      DO 20 JJ=1,NP
      J=JJ
   20 CALL H12(2,J,J+1,N,XT(1,J),1,SH(J),Y,1,N,1,N)
C
C  SOLVE THE SYSTEM
C
      DO 30 I=1,N
   30 THETA(I)=Y(I)
      CALL SOLV(XT,THETA,NP,NP,MDXT,N)
C
C  TRANSFORM THE SOLUTION VECTOR FOR OUTPUT
C
      CALL PERM(THETA,SP,NP,NP)
      RETURN
      END
C
C-----------------------------------------------------------------------
C
      SUBROUTINE QRSSH(RS,EXRHO,N,NP,SIGMA,QR)
C.......................................................................
C
C   R O B E T H  -  R O B S Y S   RELEASE 3.0 (COPYRIGHT) 1985, 1990
C
C   PROGRAMMER : A. RANDRIAMIHARISOA
C.......................................................................
C
      DIMENSION RS(N)
      EXTERNAL EXRHO
      TMP=0.
      DO 10 I=1,N
        S=RS(I)/SIGMA
        TMP=TMP+EXRHO(S)
   10 CONTINUE
      QR=TMP/FLOAT(N-NP)
      RETURN
      END

c
c***********************************************************************
c
      subroutine int21(x,y,n,np,nq,ncov,mdx,mdw,mdi,iopt,intch,nrep,
     x           tols,tolr,tau,gam,maxit,maxs1,maxs2,expsi,expsp,
     x           exchi,iseed,ierr,smin,theta,rs,it1,cov,work,iwork)
C.......................................................................
C
C   COPYRIGHT 1992 Alfio Marazzi
C
C   AUTHOR : J. JOSS
C.......................................................................
C
       integer n,np,nq,ncov,mdx,mdw,mdi,iopt,intch,nrep,maxit,maxs1
       integer maxs2,expsi,expsp,exchi,iseed,ierr,it1(nq),iwork(mdi)
       real x(mdx,np),y(n),tols,tolr,tau,gam,smin,theta(n),rs(n),
     x      cov(ncov),work(mdw)
       external psy,userfs
       if (expsi.eq.1) then
         call int22(x,y,n,np,nq,ncov,mdx,mdw,mdi,iopt,intch,nrep,
     x              tols,tolr,tau,gam,maxit,maxs1,maxs2,psy,expsp,
     x              exchi,iseed,ierr,smin,theta,rs,it1,cov,work,iwork)
       else
         call int22(x,y,n,np,nq,ncov,mdx,mdw,mdi,iopt,intch,nrep,
     x              tols,tolr,tau,gam,maxit,maxs1,maxs2,userfs,expsp,
     x              exchi,iseed,ierr,smin,theta,rs,it1,cov,work,iwork)
       endif
       return
      end
      subroutine int22(x,y,n,np,nq,ncov,mdx,mdw,mdi,iopt,intch,nrep,
     x           tols,tolr,tau,gam,maxit,maxs1,maxs2,expsi,expsp,
     x           exchi,iseed,ierr,smin,theta,rs,it1,cov,work,iwork)
C.......................................................................
C
C   COPYRIGHT 1992 Alfio Marazzi
C
C   AUTHOR : J. JOSS
C.......................................................................
C
       integer n,np,nq,ncov,mdx,mdw,mdi,iopt,intch,nrep,maxit,maxs1
       integer maxs2,expsp,exchi,iseed,ierr,it1(nq),iwork(mdi)
       real x(mdx,np),y(n),tols,tolr,tau,gam,smin,theta(n),rs(n),
     x      cov(ncov),work(mdw)
       external expsi,psp,userfs
       if (expsp.eq.3) then
         call int23(x,y,n,np,nq,ncov,mdx,mdw,mdi,iopt,intch,nrep,
     x              tols,tolr,tau,gam,maxit,maxs1,maxs2,expsi,psp,
     x              exchi,iseed,ierr,smin,theta,rs,it1,cov,work,iwork)
       else
         call int23(x,y,n,np,nq,ncov,mdx,mdw,mdi,iopt,intch,nrep,
     x              tols,tolr,tau,gam,maxit,maxs1,maxs2,expsi,userfs,
     x              exchi,iseed,ierr,smin,theta,rs,it1,cov,work,iwork)
       endif
       return
      end
      subroutine int23(x,y,n,np,nq,ncov,mdx,mdw,mdi,iopt,intch,nrep,
     x           tols,tolr,tau,gam,maxit,maxs1,maxs2,expsi,expsp,
     x           exchi,iseed,ierr,smin,theta,rs,it1,cov,work,iwork)
C.......................................................................
C
C   COPYRIGHT 1992 Alfio Marazzi
C
C   AUTHOR : J. JOSS
C.......................................................................
C
       integer n,np,nq,ncov,mdx,mdw,mdi,iopt,intch,nrep,maxit,maxs1
       integer maxs2,exchi,iseed,ierr,it1(nq),iwork(mdi)
       real x(mdx,np),y(n),tols,tolr,tau,gam,smin,theta(n),rs(n),
     x      cov(ncov),work(mdw)
       external expsi,expsp,chi,userfs
       if (exchi.eq.4) then
         call hysest(x,y,n,np,nq,ncov,mdx,mdw,mdi,iopt,intch,nrep,
     x               tols,tolr,tau,gam,maxit,maxs1,maxs2,expsi,expsp,
     x               chi,iseed,ierr,smin,theta,rs,it1,cov,work,iwork)
       else
         call hysest(x,y,n,np,nq,ncov,mdx,mdw,mdi,iopt,intch,nrep,
     x               tols,tolr,tau,gam,maxit,maxs1,maxs2,expsi,expsp,
     x               userfs,iseed,ierr,smin,theta,rs,it1,cov,work,iwork)
       endif
       return
      end

C
C-----------------------------------------------------------------------
C
      SUBROUTINE HYSEST(X,Y,N,NP,NQ,NCOV,MDX,MDW,MDI,IOPT,INTCH,NREP,
     *           TOLS,TOLR,TAU,GAM,MAXIT,MAXS1,MAXS2,EXPSI,EXPSP,EXCHI,
     *           ISEED,IERR,SMIN,THETA,RS,IT1,COV,WORK,IWORK)
C.......................................................................
C
C   R O B E T H  -  R O B S Y S   RELEASE 3.0 (COPYRIGHT) 1985, 1990
C
C   PROGRAMMERS : A. MARAZZI / A. RANDRIAMIHARISOA
C.......................................................................
C
      DIMENSION X(MDX,NP),Y(N),THETA(N),RS(N),IT1(NQ),COV(NCOV)
      DIMENSION WORK(MDW),IWORK(MDI)
      EXTERNAL EXPSI,EXPSP,EXCHI
      COMMON/BETA/BETA,BET0
C
C   Resampling algorithm for the computation of S-estimates
C
      IMIN=NP+NQ
      NP1=NP+1
      MINW=NQ*(NP+2)+(MDX+3)*NP+N
      NN=NP*(NP+1)/2
      IF (N.LE.0 .OR. MDX.LT.N .OR. NP.LE.0 
     * .OR. NQ.LT.NP .OR. NCOV.NE.NN .OR. MDW.LT.MINW .OR. MDI.LT.IMIN
     * .OR. IOPT.LT.0 .OR. IOPT.GT.3 .OR. (IOPT.EQ.2 .AND. NREP.LE.0)
     * .OR. (INTCH.NE.0.AND.INTCH.NE.1) .OR. TOLS.LE.0. OR.
     *  TOLR.LE.0. .OR. TAU.LT.0. .OR. GAM.LE.0. .OR. GAM.GT.2. .OR.
     * MAXIT.LE.0 .OR. MAXS1.LE.0 .OR. MAXS2.LE.0)
     * CALL MESSGE(500,'HYSEST',1)
      N0=NP*NQ+1
      N1=N0+NQ
      N2=N1+NQ
      N3=N2+NP
      N4=N3+NP
      N5=N4+NP
      N6=N5+MDX*NP
      CALL HSEST2(X,Y,N,NP,NQ,NCOV,MDX,IOPT,INTCH,NREP,TOLS,TOLR,
     *            TAU,GAM,MAXIT,MAXS1,MAXS2,EXPSI,EXPSP,EXCHI,
     *            ISEED,IERR,SMIN,THETA,RS,IT1,COV,
     *            WORK(1),WORK(N0),WORK(N1),WORK(N2),WORK(N3),WORK(N4),
     *            WORK(N5),WORK(N6),IWORK(1),IWORK(NP1))
      RETURN
      END
C
C-----------------------------------------------------------------------
C
      SUBROUTINE HSEST2(X,Y,N,NP,NQ,NCOV,MDX,IOPT,INTCH,NREP,TOLS,TOLR,
     *           TAU,GAM,MAXIT,MAXS1,MAXS2,EXPSI,EXPSP,EXCHI,ISEED,
     *           IERR,SMIN,THETA,RS,IT1,COV,XX,YY,XTHETA,
     *           SF,SG,SH,SX,SZ,SP,IT)
C.......................................................................
C
C   R O B E T H  -  R O B S Y S   RELEASE 3.0 (COPYRIGHT) 1985, 1990
C
C   PROGRAMMERS : A. MARAZZI / A. RANDRIAMIHARISOA
C.......................................................................
C
      DIMENSION X(MDX,NP),Y(N),THETA(N),RS(N),COV(NCOV),XX(NQ,NP),YY(NQ)
      DIMENSION XTHETA(NQ),SF(NP),SG(NP),SH(NP),SX(MDX,NP),SZ(N)
      INTEGER IT1(NQ),SP(NP),IT(NQ)
      EXTERNAL EXPSI,EXPSP,EXCHI,ICNREP
      COMMON/BETA/BETA,BET0,/CONST/CONST
C
C   Resampling algorithm for the computation of S-estimates
C
      IF (N.LE.0 .OR. MDX.LT.N .OR. NP.LE.0
     * .OR. NQ.LT.NP .OR. NCOV.NE.NP*(NP+1)/2 .OR.
     *  IOPT.LT.0 .OR. IOPT.GT.3 .OR. (IOPT.EQ.2 .AND. NREP.LE.0)
     * .OR. (INTCH.NE.0.AND.INTCH.NE.1) .OR. TOLS.LE.0. OR.
     *  TOLR.LE.0. .OR. TAU.LT.0. .OR. GAM.LE.0. .OR. GAM.GT.2. .OR.
     * MAXIT.LE.0 .OR. MAXS1.LE.0 .OR. MAXS2.LE.0)
     * CALL MESSGE(500,'HSEST2',1)
C
C STEP 0: INITIALIZATIONS
C ------
      N2=N/2
      N2P=N-N2
      K1=N2+1
      NK1=N-K1+1
      CONST=BETA*FLOAT(N-NP)
      IF (IOPT.NE.2) NREP=ICNREP(N,NQ,IOPT,0)
      NIT=1
      IERR=2
      SMIN=0.
      ITYPE=1
      NITMON=0
      PSP0=EXPSP(0.)
C
C STEP 1: DRAW A SUBSAMPLE
C ------ 
c     irep=0 
  100 IF (IOPT.NE.3) THEN
        DO 130 K=1,NQ
  110     CALL RANDOW(ISEED,RND)
c         irep=irep+1
c         ver(irep)=rnd 
c         IF (irep.eq.5) then
c           call realpr('rnd',3,ver,5)
c           irep=0
c         endif
          ITK=RND*N+1
          DO 120 KK=1,K-1
          IF (ITK.EQ.IT(KK)) GOTO 110
  120     CONTINUE
          IT(K)=ITK
  130   CONTINUE
      ELSE
        IF (NIT.EQ.1) THEN
          DO 140 K=1,NQ
  140     IT(K)=K
        ELSE
          CALL  NCOMB(N,NQ,IT)
        ENDIF
      ENDIF
      DO 160 K=1,NQ
      ITK=IT(K)
      DO 150 J=1,NP
  150 XX(K,J)=X(ITK,J)
  160 YY(K)=Y(ITK)
C
C STEP 2: DECOMPOSE SAMPLE MATRIX
C -------
      CALL RIMTRF(XX,NQ,NP,NQ,INTCH,TAU,KK,SF,SG,SH,SP)
      IF(KK.NE.NP) GOTO 700
C
C STEP 3: SOLVE SYSTEM OF LINEAR EQUATIONS
C -------
      CALL RICLL1(XX,YY,NQ,NP,NQ,XTHETA,SH,SP)
C
C STEP 4: COMPUTE RESIDUALS
C -------
      DO 420 I=1,N
      S=Y(I)
      DO 410 J=1,NP
  410 S=S-XTHETA(J)*X(I,J)
  420 RS(I)=S
      IF (SMIN.EQ.0.) THEN
        S=1.0E7
        DO 430 I=1,N
        ARI=ABS(RS(I))
        SZ(I)=ARI
        IF (ARI.NE.0.) S=AMIN1(S,ARI)
  430   CONTINUE
        IF (S.EQ.1.0E7) GOTO 915
        CALL FSTORD(SZ,N,K1,S0)
        S0=2.*S0
        IF (S0.EQ.0.) S0=S
        SRES=S0
      ENDIF
  435 D=0.
      DO 440 I=1,N
  440 D=D+EXCHI(RS(I)/SRES)
      IF (SMIN.NE.0..AND.D.GT.CONST) GOTO 700
      IF (D.LE.CONST) GOTO 500
      S0=1.5*S0
      SRES=S0
      GOTO 435
C
C STEP 5: SOLVE FOR SRES
C -------
  500 CALL RYSIGM(RS,SZ,EXCHI,S0,N,NP,TOLR,ITYPE,1,MAXS1,NIS,SRES,SZ,SZ)
      IF (NIS.EQ.MAXS1) CALL MESSGE(110,'HSEST2',0)
C
C STEP 6: UPDATE BEST FIT
C ------
      IERR=0
      SMIN=SRES
      S0=SMIN
      DO 610 K=1,NP
      THETA(K)=XTHETA(K)
  610 continue
      DO 620 K=1,NQ
      IT1(K)=IT(K)
  620 CONTINUE
      IF (SRES .LE. TOLS) THEN
        IERR=1
        GOTO 800
      ENDIF
C
C STEP 7: END OF MAIN LOOP
C -------
  700 IF (NIT.EQ.NREP) GOTO 800
      NIT=NIT+1
      GOTO 100
C
C STEP 8: SOLVE SYSTEM OF EQUATIONS FOR THETA AND SRES
C --------
  800 IF (IERR.EQ.2) RETURN
      DO 820 I=1,N
      S=Y(I)
      DO 810 J=1,NP
  810 S=S-THETA(J)*X(I,J)
  820 RS(I)=S
      K=1
      MAXIW=1
      ISIGMA=-1
  830 SWI=0.
      DO 860 I=1,N
      WI=0.
      IF (RS(I).EQ.0.) GOTO 840
      T=RS(I)/SMIN
      WI=EXPSI(T)/T
      SWI=SWI+WI
      WI=SQRT(WI)
  840 DO 850 J=1,NP
  850 SX(I,J)=WI*X(I,J)
  860 CONTINUE
      CALL KFFACV(RS,EXPSI,EXPSP,N,NP,SMIN,FH)
      FACT=FH*SWI/FLOAT(N)
      IF (K.EQ.0) FACT=FACT*SMIN*SMIN
      CALL KTASKV(SX,N,NP,MDX,NCOV,TAU,FACT,XX,COV)
      IF (K.EQ.0) RETURN
      SRES=SMIN
      ICNV=1
      DO 870 J=1,NP
  870 XTHETA(J)=THETA(J)
      IF (MAXIW.EQ.1) CALL QRSSH(RS,EXCHI,N,NP,SRES,QR0)
  880 CONTINUE
       CALL RYWALG(X,Y,THETA,SZ,COV,PSP0,EXPSI,EXCHI,EXCHI,SRES,
     * N,NP,MDX,MDX,NCOV,TOLR,GAM,TAU,ITYPE,ISIGMA,ICNV,MAXIW,
     * MAXS2,NITMON,NIT8,SMIN,RS,YY,SZ,SF,SG,SH,SP,SZ,SX)
C
C STEP 9: EXIT
C -------
      IF (MAXIT.EQ.1) GOTO 900
      IF (MAXIW.EQ.1) THEN
        CALL QRSSH(RS,EXCHI,N,NP,SRES,QR1)
        IF (QR0.LE.QR1) GOTO 910
        ISIGMA=1
        MAXIW=MAXIW+MAXIT
        GOTO 880
      ENDIF
      IF (NIT8.EQ.MAXIW) CALL MESSGE(111,'HSEST2',0)
      IF (SMIN.GE.SRES) GOTO 910
  900 K=0
      GOTO 830
  910 CALL MESSGE(112,'HSEST2',0)
      SMIN=SRES
      FACT=SMIN*SMIN
      CALL SCAL(COV,FACT,NCOV,1,NCOV)
  915 DO 920 J=1,NP
  920 THETA(J)=XTHETA(J)
      DO 940 I=1,N
      S=Y(I)
      DO 930 J=1,NP
  930 S=S-THETA(J)*X(I,J)
  940 RS(I)=S
      RETURN
      END
C
C-----------------------------------------------------------------------
C
      SUBROUTINE XERF(KODE,X,P)
C.......................................................................
C
C   COPYRIGHT 1992 Alfio Marazzi
C
C   AUTHOR : A. MARAZZI
C.......................................................................
C
      EXTERNAL XEXP
      DATA SPI/2.506628274631/
C
C  EXMIN IS A MACHINE DEPENDENT PARAMETER SPECIFYING THE LARGEST NEGATIVE
C  REAL VALUE SUCH THAT EXP(EXMIN) CAN BE SUCCESSFULLY EVALUATED WITHOUT
C  UNDERFLOW.
C
      IF (KODE.NE.1.AND.KODE.NE.2) CALL MESSGE(500,'XERF  ',1)
      X2=-X*X/2.
      P=XEXP(X2)
      IF (KODE.EQ.2) P=P/SPI
      RETURN
      END
C
C-----------------------------------------------------------------------
C
      SUBROUTINE DOTP(X,Y,N,INCX,INCY,NX,NY,RESULT)
C.......................................................................
C
C   COPYRIGHT 1979 SOCIETY FOR INDUSTRIAL AND APPLIED MATHEMATICS.
C   ALL RIGHTS RESERVED.
C
C   AUTHOR :     LINPACK (SUBROUTINE SDOT)
C                REPRINTED WITH PERMISSION FROM 
C                LINPACK USER'S GUIDE.
C                ADAPTED FOR ROBETH BY A. MARAZZI
C.......................................................................
C
      REAL X(NX),Y(NY)
      DOUBLE PRECISION DTEMP
      LOGICAL NPRCHK
C
C  PARAMETER CHECK
C
      NPRCHK=INCX.NE.0.AND.IABS(INCX)*(N-1)+1.LE.NX
     1       .AND.INCY.NE.0.AND.IABS(INCY)*(N-1)+1.LE.NY
      IF (.NOT.NPRCHK) CALL MESSGE(500,'DOTP  ',1)
C
      DTEMP=0.D0
      RESULT=0.
      IF (N.LE.0) RETURN
      IF (INCX.EQ.1.AND.INCY.EQ.1) GOTO 20
C
C  CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS NOT EQUAL TO 1
C
      IX=1
      IY=1
      IF (INCX.LT.0) IX=(-N+1)*INCX+1
      IF (INCY.LT.0) IY=(-N+1)*INCY+1
      DO 10 I=1,N
      DTEMP=DTEMP+X(IX)*DBLE(Y(IY))
      IX=IX+INCX
      IY=IY+INCY
   10 CONTINUE
      RESULT=DTEMP
      RETURN
C
C  CODE FOR BOTH INCREMENTS EQUAL TO 1
C
   20 M=MOD(N,5)
      IF (M.EQ.0) GOTO 40
      DO 30 I=1,M
      DTEMP=DTEMP+X(I)*DBLE(Y(I))
   30 CONTINUE
      IF (N.LT.5) GOTO 60
   40 MP1=M+1
      DO 50 I=MP1,N,5
      DTEMP=DTEMP+X(I)*DBLE(Y(I))+X(I+1)*DBLE(Y(I+1))+
     1      X(I+2)*DBLE(Y(I+2))+X(I+3)*DBLE(Y(I+3))+
     1      X(I+4)*DBLE(Y(I+4))
   50 CONTINUE
   60 RESULT=DTEMP
      RETURN
      END
C
C-----------------------------------------------------------------------
C
      SUBROUTINE SCAL(X,SA,N,INCX,MDX)
C.......................................................................
C
C   COPYRIGHT 1979 SOCIETY FOR INDUSTRIAL AND APPLIED MATHEMATICS.
C   ALL RIGHTS RESERVED.
C
C   AUTHOR :     LINPACK (SUBROUTINE SSCAL)
C                REPRINTED WITH PERMISSION FROM 
C                LINPACK USER'S GUIDE.
C                ADAPTED FOR ROBETH BY A. MARAZZI
C.......................................................................
C
      REAL X(MDX)
      LOGICAL NPRCHK
C
C  PARAMETER CHECK
C
      NPRCHK=INCX.GT.0.AND.N.GE.0.AND.INCX*(N-1)+1.LE.MDX
      IF (.NOT.NPRCHK) CALL MESSGE(500,'SCAL  ',1)
C
      IF (N.LE.0) RETURN
      IF (INCX.EQ.1) GOTO 20
C
C  CODE FOR INCREMENT NOT EQUAL TO 1
C
      NINCX=N*INCX
      DO 10 I=1,NINCX,INCX
      X(I)=SA*X(I)
   10 CONTINUE
      RETURN
C
C  CODE FOR INCREMENT EQUAL TO 1
C
   20 M=MOD(N,5)
      IF (M.EQ.0) GOTO 40
      DO 30 I=1,M
      X(I)=SA*X(I)
   30 CONTINUE
      IF (N.LT.5) RETURN
   40 MP1=M+1
      DO 50 I=MP1,N,5
      X(I)=SA*X(I)
      X(I+1)=SA*X(I+1)
      X(I+2)=SA*X(I+2)
      X(I+3)=SA*X(I+3)
      X(I+4)=SA*X(I+4)
   50 CONTINUE
      RETURN
      END
C
C----------------------------------------------------------------------
C
      SUBROUTINE H12(MODE,LPIVOT,L1,M,U,IUE,UP,C,ICE,ICV,NCV,
     1               MDC)
C.......................................................................
C
C   AUTHORS :     CH.L. LAWSON & R.J. HANSON (1974)
C                 SOLVING LEAST SQUARES PROBLEMS 
C                 REPRINT FROM PP.290-291,308 BY PERMISSION OF 
C                 PRENTICE HALL, ENGLEWOOD CLIFFS, NEW JERSEY.
C                 ADAPTED FOR ROBETH BY A. MARAZZI
C.......................................................................
C
      REAL U(IUE,M),C(MDC)
      DOUBLE PRECISION SM,B
      ONE=1.
C
      IF (0.GE.LPIVOT.OR.LPIVOT.GE.L1.OR.L1.GT.M) RETURN
      CL=ABS(U(1,LPIVOT))
      IF (MODE.EQ.2) GOTO 60
C
C  CONSTRUCT THE TRANSFORMATION
C
      DO 10 J=L1,M
   10 CL=AMAX1(ABS(U(1,J)),CL)
      IF (CL) 130,130,20
   20 CLINV=ONE/CL
      SM=(DBLE(U(1,LPIVOT))*CLINV)**2
      DO 30 J=L1,M
   30 SM=SM+(DBLE(U(1,J))*CLINV)**2
C
C  CONVERT DBLE. PRE. SM TO SNGL. PREC. SM1
C
      SM1=SM
      CL=CL*SQRT(SM1)
      IF (U(1,LPIVOT)) 50,50,40
   40 CL=-CL
   50 UP=U(1,LPIVOT)-CL
      U(1,LPIVOT)=CL
      GOTO 70
C
C  APPLY THE TRANSFORMATION I+U*(U**T)/B TO C
C
   60 IF (CL) 130,130,70
   70 IF (NCV.LE.0) RETURN
      B=DBLE(UP)*U(1,LPIVOT)
C
C  B MUST BE NONPOSITIVE HERE. IF B=0., RETURN.
C
      IF (B) 80,130,130
   80 B=ONE/B
      I2=1-ICV+ICE*(LPIVOT-1)
      INCR=ICE*(L1-LPIVOT)
      DO 120 J=1,NCV
      I2=I2+ICV
      I3=I2+INCR
      I4=I3
      SM=C(I2)*DBLE(UP)
      DO 90 I=L1,M
      SM=SM+C(I3)*DBLE(U(1,I))
   90 I3=I3+ICE
      IF (SM) 100,120,100
  100 SM=SM*B
      C(I2)=C(I2)+SM*DBLE(UP)
      DO 110 I=L1,M
      C(I4)=C(I4)+SM*DBLE(U(1,I))
  110 I4=I4+ICE
  120 CONTINUE
  130 RETURN
      END
C
C-----------------------------------------------------------------------
C
      SUBROUTINE SOLV(X,THETA,NP,K,MDX,MDT)
C.......................................................................
C
C   COPYRIGHT 1992 Alfio Marazzi
C
C   AUTHOR : A. MARAZZI
C.......................................................................
C
C  PURPOSE
C  -------
C  LET U BE THE K BY K UPPER TRIANGULAR MATRIX WITH ELEMENTS
C  X(1,1)...X(1,K),X(2,2)...X(2,K),...X(K,K).
C  SOLV SOLVES THE TRIANGULAR SYSTEM U*THETA=Y (BY BACK SUBSTI-
C  TUTION). ON INPUT Y IS CONTAINED IN THETA.  ON OUTPUT
C  THETA(1)...THETA(K) CONTAIN THE DESIRED SOLUTION.
C
C  ERRORS
C  ------
C  1   AN ELEMENT OF THE PRINCIPAL DIAGONAL OF X IS =0.
C
      DIMENSION X(MDX,NP),THETA(MDT)
      DOUBLE PRECISION SM,DZERO
      DZERO=0.D0
      KP1=K+1
      DO 80 L=1,K
      SM=DZERO
      I=KP1-L
      IF (I.EQ.K) GOTO 60
      IP1=I+1
      DO 50 J=IP1,K
   50 SM=SM+X(I,J)*DBLE(THETA(J))
   60 SM1=SM
      IF (X(I,I)) 80,70,80
   70 CALL MESSGE(501,'SOLV  ',1)
   80 THETA(I)=(THETA(I)-SM1)/X(I,I)
      RETURN
      END
C
C-----------------------------------------------------------------------
C
      SUBROUTINE PERM(X,SP,N,NDIM)
C.......................................................................
C
C   COPYRIGHT 1992 Alfio Marazzi
C
C   AUTHOR : A. MARAZZI
C.......................................................................
C
C  PURPOSE
C  -------
C  PERMUTE COMPONENTS OF X TO COMPENSATE COLUMN INTERCH.
C
      DIMENSION X(NDIM)
      INTEGER SP(NDIM)
      DO 10 JJ=1,N
      J=N-JJ+1
      IF (SP(J).EQ.J) GOTO 10
      L=SP(J)
      TMP=X(L)
      X(L)=X(J)
      X(J)=TMP
   10 CONTINUE
      RETURN
      END
C
C-----------------------------------------------------------------------
C
      SUBROUTINE XSYD(X,Y,S,N,NN,RESULT)
C.......................................................................
C
C   COPYRIGHT 1992 Alfio Marazzi
C
C   AUTHORS : A. MARAZZI / A. RANDRIAMIHARISOA
C.......................................................................
C
      DOUBLE PRECISION X(N),Y(N),S(NN),RESULT,SM
C
C  PARAMETER CHECK
C
      NS=N*(N+1)/2
      IF (N.LE.0.OR.NN.NE.NS) CALL MESSGE(500,'XSYD  ',1)
C
      SM=0.D0
      L=0
      DO 20 I=1,N
      L=L+I
      L1=L-I+1
      K=0
      DO 20 J=L1,L
      K=K+1
      IF (J.EQ.L) GOTO 10
      SM=SM+S(J)*(X(I)*Y(K)+X(K)*Y(I))
      GOTO 20
   10 SM=SM+S(J)*X(I)*Y(I)
   20 CONTINUE
      RESULT=SM
      RETURN
      END
C
C-----------------------------------------------------------------------
C
      SUBROUTINE MINV(R,N,NN,TAU,ISING)
C.......................................................................
C
C   COPYRIGHT 1992 Alfio Marazzi
C
C   AUTHOR : A. MARAZZI
C.......................................................................
C
      REAL R(NN)
      DOUBLE PRECISION SM,DZERO
      LOGICAL NPRCHK
C
C  PARAMETER CHECK
C
      NPRCHK=N.GT.0.AND.NN.EQ.(N*(N+1)/2).AND.TAU.GE.0.
      IF (.NOT.NPRCHK) CALL MESSGE(500,'MINV  ',1)
C
      DZERO=0.D0
      ISING=0
      I1=0
      DO 10 I=1,N
      I1=I1+I
      IF (ABS(R(I1)).LE.TAU) GOTO 900
   10 R(I1)=1./R(I1)
      IF (N.EQ.1) RETURN
      I1=0
      NM1=N-1
      DO 40 I=1,NM1
      I1=I1+I
      J1=I1+I
      IP1=I+1
      DO 30 J=IP1,N
      SM=DZERO
      IL=I1
      LJ=J1
      JM1=J-1
      DO 20 L=I,JM1
      SM=SM+R(IL)*DBLE(R(LJ))
      LJ=LJ+1
   20 IL=IL+L
      R(J1)=-R(LJ)*SM
   30 J1=J1+J
   40 CONTINUE
      RETURN
  900 ISING=1
      RETURN
      END
C
C-----------------------------------------------------------------------
C
      SUBROUTINE MCHL(A,N,NN,INFO)
C.......................................................................
C
C   COPYRIGHT 1979 SOCIETY FOR INDUSTRIAL AND APPLIED MATHEMATICS.
C   ALL RIGHTS RESERVED.
C
C   AUTHOR :     LINPACK (SUBROUTINE SPPFA)
C                REPRINTED WITH PERMISSION FROM 
C                LINPACK USER'S GUIDE.
C                ADAPTED FOR ROBETH BY A. MARAZZI
C.......................................................................
C
      REAL A(NN)
      DOUBLE PRECISION S
      LOGICAL NPRCHK
C
C  PARAMETER CHECK
C
      NPRCHK=N.GT.0.AND.NN.EQ.(N*(N+1)/2)
      IF (.NOT.NPRCHK) CALL MESSGE(500,'MCHL  ',1)
C
      JJ=0
      DO 30 J=1,N
      INFO=J
      S=0.D0
      JM1=J-1
      KJ=JJ
      KK=0
      IF (JM1.LT.1) GOTO 20
      DO 10 K=1,JM1
      KJ=KJ+1
      CALL DOTP(A(KK+1),A(JJ+1),K-1,1,1,NN-KK,NN-JJ,DTP)
      T=A(KJ)-DTP
      KK=KK+K
      T=T/A(KK)
      A(KJ)=T
      S=S+T*DBLE(T)
   10 CONTINUE
   20 CONTINUE
      JJ=JJ+J
      S=DBLE(A(JJ))-S
      IF (S.LE.0.D0) GOTO 40
      A(JJ)=DSQRT(S)
   30 CONTINUE
      INFO=0
   40 CONTINUE
      RETURN
      END
C
C-----------------------------------------------------------------------
C
      SUBROUTINE MTT1(A,B,N,NN)
C.......................................................................
C
C   COPYRIGHT 1992 Alfio Marazzi
C
C   AUTHORS : A. MARAZZI / R. DUTTER
C.......................................................................
C
      REAL A(NN),B(NN)
      DOUBLE PRECISION SM,DZERO
      LOGICAL NPRCHK
C
C  PARAMETER CHECK
C
      NPRCHK=N.GT.0.AND.NN.EQ.(N*(N+1)/2)
      IF (.NOT.NPRCHK) CALL MESSGE(500,'MTT1  ',1)
C
      DZERO=0.D0
      IJ=0
      JJ=0
      DO 30 J=1,N
      DO 20 I=1,J
      IJ=IJ+1
      SM=DZERO
      IL=JJ+I
      JL=JJ+J
      DO 10 L=J,N
      SM=SM+A(IL)*DBLE(A(JL))
      IL=IL+L
   10 JL=JL+L
      B(IJ)=SNGL(SM)
   20 CONTINUE
   30 JJ=JJ+J
      RETURN
      END
C
C-----------------------------------------------------------------------
C
      SUBROUTINE RESIDU(X,Y,THETA,N,NP,MDX,RS)
C.......................................................................
C
C   COPYRIGHT 1992 Alfio Marazzi
C
C   AUTHOR : A. RANDRIAMIHARISOA
C.......................................................................
C
C  COMPUTES RESIDUALS
C  RS(I)=Y(I)-SUM X(I,J)*THETA(J)
C
      DIMENSION X(MDX,NP),Y(N),THETA(NP),RS(N)
      DOUBLE PRECISION SUM
      DO 200 I=1,N
        SUM=0.D0
        DO 100 J=1,NP
          SUM=SUM+X(I,J)*DBLE(THETA(J))
  100   CONTINUE
        RS(I)=Y(I)-SUM
  200 CONTINUE
      RETURN
      END
C
C-----------------------------------------------------------------------
C
      SUBROUTINE SWAP(X,Y,N,INCX,INCY,MDX,MDY)
C.......................................................................
C
C   COPYRIGHT 1979 SOCIETY FOR INDUSTRIAL AND APPLIED MATHEMATICS.
C   ALL RIGHTS RESERVED.
C
C   AUTHOR :     LINPACK (SUBROUTINE SSWAP)
C                REPRINTED WITH PERMISSION FROM 
C                LINPACK USER'S GUIDE.
C                ADAPTED FOR ROBETH BY A. MARAZZI
C.......................................................................
C
      REAL X(MDX),Y(MDY)
      LOGICAL NPRCHK
C
C  PARAMETER CHECK
C
      NPRCHK=N.GE.0.AND.INCX.NE.0.AND.IABS(INCX)*(N-1)+1.LE.MDX
     1       .AND.INCY.NE.0.AND.IABS(INCY)*(N-1)+1.LE.MDY
      IF (.NOT.NPRCHK) CALL MESSGE(500,'SWAP  ',1)
C
      IF (N.EQ.0) RETURN
      IF (INCX.EQ.1.AND.INCY.EQ.1) GOTO 20
C
C  CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS NOT
C  EQUAL TO 1
C
      IX=1
      IY=1
      IF (INCX.LT.0) IX=(-N+1)*INCX+1
      IF (INCY.LT.0) IY=(-N+1)*INCY+1
      DO 10 I=1,N
      TEMP=X(IX)
      X(IX)=Y(IY)
      Y(IY)=TEMP
      IX=IX+INCX
      IY=IY+INCY
   10 CONTINUE
      RETURN
C
C  CODE FOR BOTH INCREMENTS EQUAL TO 1
C
   20 M=MOD(N,3)
      IF (M.EQ.0) GOTO 40
      DO 30 I=1,M
      TEMP=X(I)
      X(I)=Y(I)
      Y(I)=TEMP
   30 CONTINUE
      IF (N.LT.3) RETURN
   40 MP1=M+1
      DO 50 I=MP1,N,3
      TEMP=X(I)
      X(I)=Y(I)
      Y(I)=TEMP
      TEMP=X(I+1)
      X(I+1)=Y(I+1)
      Y(I+1)=TEMP
      TEMP=X(I+2)
      X(I+2)=Y(I+2)
      Y(I+2)=TEMP
   50 CONTINUE
      RETURN
      END
C
C-----------------------------------------------------------------------
C
      SUBROUTINE XSY(X,Y,S,N,NN,RESULT)
C.......................................................................
C
C   COPYRIGHT 1992 Alfio Marazzi
C
C   AUTHOR : A. MARAZZI
C.......................................................................
C
      REAL X(N),Y(N),S(NN)
      DOUBLE PRECISION SM
C
C  PARAMETER CHECK
C
      NS=N*(N+1)/2
      IF (N.LE.0.OR.NN.NE.NS) CALL MESSGE(500,'XSY   ',1)
C
      SM=0.D0
      L=0
      DO 20 I=1,N
      L=L+I
      L1=L-I+1
      K=0
      DO 20 J=L1,L
      K=K+1
      IF (J.EQ.L) GOTO 10
      SM=SM+DBLE(S(J))*(X(I)*Y(K)+X(K)*Y(I))
      GOTO 20
   10 SM=SM+DBLE(S(J))*X(I)*Y(I)
   20 CONTINUE
      RESULT=SM
      RETURN
      END

C
C-----------------------------------------------------------------------
C
      SUBROUTINE FACS(RS,N,K,SIGMA,TL,XKAPPA,SUM2,PSY,PSP)
C.......................................................................
C
C   COPYRIGHT 1992 Alfio Marazzi
C
C   AUTHORS : A. MARAZZI / A. RANDRIAMIHARISOA
C.......................................................................
C
C  PURPOSE
C  -------
C  COMPUTES CORRECTION FACTORS XKAPPA AND SUM2 FOR
C  THE COVARIANCE MATRIX.
C  FACS CALLS THE FUNCTIONS PSI AND PSP
C
      DIMENSION RS(N)
      EXTERNAL PSY,PSP
      TMP1=0.
      TMP2=0.
      DO 10 J=1,N
      S=RS(J)/SIGMA
      TMP1=TMP1+PSP(S)
      PS=PSY(S)
   10 TMP2=TMP2+PS*PS
      XMU=TMP1/FLOAT(N)
      SUM2=TMP2
      VAR=0.
      DO 20 J=1,N
      S=RS(J)/SIGMA
   20 VAR=VAR+(PSP(S)-XMU)**2
      VAR=VAR/FLOAT(N)
      XKAPPA=0.
      IF (XMU.LE.TL) RETURN
      XMU2=XMU*XMU
      XKAPPA=1.+FLOAT(K)*VAR/FLOAT(N)/XMU2
      SUM2=SUM2/XMU2/FLOAT(N-K)
      RETURN
      END
C
C-----------------------------------------------------------------------
C
      FUNCTION DIFF(X,Y)
C.......................................................................
C
C   COPYRIGHT 1992 Alfio Marazzi
C
C   AUTHOR : A. MARAZZI
C.......................................................................
C
C  PURPOSE
C  -------
C  SEE (L6) IN H12, P.278.
C
      DIFF=X-Y
      RETURN
      END
C
C-----------------------------------------------------------------------
C
      FUNCTION CHIPHI(S,WGT,N,FCHI,FEXT)
C.......................................................................
C
C   COPYRIGHT 1992 Alfio Marazzi
C
C   AUTHORS : A. MARAZZI / A. RANDRIAMIHARISOA
C.......................................................................
C
C  PURPOSE
C  -------
C  AUXILIARY ROUTINE FOR RIBETU.
C
      DIMENSION WGT(N)
      EXTERNAL FCHI,FEXT
      COMMON/INTPAR/ITYPE,INTPAR,NEVAL,LIMIT,KEY
      CALL XERF(2,S,PHI)
      IF (ITYPE.EQ.3) GOTO 30
C
C  HUBER & MALLOWS CASE
C
      CHIPHI=FCHI(S)*PHI
      RETURN
C
C  SCHWEPPE CASE
C
   30 SM=0.
      DO 40 J=1,N
   40 IF (WGT(J).GT.0.) SM=SM+WGT(J)*WGT(J)*FCHI(S/WGT(J))
      CHIPHI=SM*PHI
      RETURN
      END
C
C-----------------------------------------------------------------------
C
      FUNCTION PSPPHI(S,WGT,N,FPSI,FEXT)
C.......................................................................
C
C   COPYRIGHT 1992 Alfio Marazzi
C
C   AUTHORS : A. MARAZZI / A. RANDRIAMIHARISOA
C.......................................................................
C
C  PURPOSE
C  -------
C  AUXILIARY ROUTINE FOR KIEDCU.
C
      DIMENSION WGT(N)
      EXTERNAL FPSI,FEXT
      COMMON/INTPAR/ITYPE,I,NEVAL,LIMIT,KEY
      R=S
      CALL XERF(2,R,PHI)
      PHI=R*PHI
      IF (ITYPE.EQ.3) R=R/WGT(I)
      PSPPHI=FPSI(R)*PHI
      RETURN
      END
C
C-----------------------------------------------------------------------
C
      FUNCTION PS2PHI(S,WGT,N,FPSI,FEXT)
C.......................................................................
C
C   COPYRIGHT 1992 Alfio Marazzi
C
C   AUTHORS : A. MARAZZI / A. RANDRIAMIHARISOA
C.......................................................................
C
C  PURPOSE
C  -------
C  AUXILIARY ROUTINE FOR KIEDCU
C
      DIMENSION WGT(N)
      EXTERNAL FPSI,FEXT
      COMMON/INTPAR/ITYPE,I,NEVAL,LIMIT,KEY
      R=S
      CALL XERF(2,R,PHI)
      IF (ITYPE.EQ.3) R=R/WGT(I)
      PS2PHI=FPSI(R)*FPSI(R)*PHI
      RETURN
      END
C
C-----------------------------------------------------------------------
C
      SUBROUTINE INTCHG(A,B)
C.......................................................................
C
C   COPYRIGHT 1992 Alfio Marazzi
C
C   AUTHOR : A. RANDRIAMIHARISOA
C.......................................................................
C
C  PURPOSE
C  -------
C  AUXILIARY ROUTINE FOR RILARS
C
      C=A
      A=B
      B=C
      RETURN
      END
C
C-----------------------------------------------------------------------
C
      SUBROUTINE COL(V1,V2,MLT,M,IOUT)
C.......................................................................
C
C   COPYRIGHT 1992 Alfio Marazzi
C
C   AUTHOR : A. RANDRIAMIHARISOA
C.......................................................................
C
C     PURPOSE
C     -------
C     AUXILIARY ROUTINE FOR RILARS
C
      REAL V1(M),V2(M),MLT
      DO 220 I=1,M
        IF (I.EQ.IOUT) GOTO 220
        V1(I)=V1(I)-V2(I)*MLT
  220 CONTINUE
      RETURN
      END

C
C-----------------------------------------------------------------------
C
      SUBROUTINE RIMTRF(X,N,NP,MDX,INTCH,TAU,K,SF,SG,SH,IP)
C.......................................................................
C
C   COPYRIGHT 1992 Alfio Marazzi
C
C   AUTHOR: A. MARAZZI
C.......................................................................
C
      DIMENSION X(MDX,NP),SF(NP),SG(NP),SH(NP)
      INTEGER IP(NP)
      LOGICAL NPRCHK
      EXTERNAL DIFF
C
C  PARAMETER CHECK AND INITIALIZATION
C
      FACTOR=0.001
      LDIAG=MIN0(N,NP)
      NPRCHK=LDIAG.GT.0.AND.N.LE.MDX.AND.TAU.GE.0.
      IF (.NOT.NPRCHK) CALL MESSGE(500,'RIMTRF',1)
C
      DO 80 JJ=1,LDIAG
      J=JJ
      IF (INTCH.EQ.0) IP(J)=J
      IF (INTCH.EQ.0) GOTO 70
      IF (J.EQ.1) GOTO 20
C
C  UPDATE SQUARED COLUMN LENGTHS AND FIND LMAX
C
      LMAX=J
      DO 10 L=J,NP
      SH(L)=SH(L)-X(J-1,L)**2
      IF(SH(L).GT.SH(LMAX)) LMAX=L
   10 CONTINUE
      IF (DIFF(HMAX+FACTOR*SH(LMAX),HMAX)) 20,20,50
C
C  COMPUTE SQUARED COLUMN LENGTHS AND FIND LMAX
C
   20 LMAX=J
      DO 40 L=J,NP
      SH(L)=0.
      DO 30 I=J,N
   30 SH(L)=SH(L)+X(I,L)**2
      IF (SH(L).GT.SH(LMAX)) LMAX=L
   40 CONTINUE
      HMAX=SH(LMAX)
C
C  LMAX HAS BEEN DETERMINED: INTERCHANGE COLUMNS IF NEEDED
C
   50 CONTINUE
      IP(J)=LMAX
      IF (IP(J).EQ.J) GOTO 70
      DO 60 I=1,N
      TMP=X(I,J)
      X(I,J)=X(I,LMAX)
   60 X(I,LMAX)=TMP
      SH(LMAX)=SH(J)
C
C  COMPUTE THE HOUSEHOLDER TRANSF. Q AND APPLY IT TO X
C
   70 MDC=NP-J
      IF (MDC.GT.0)
     1CALL H12(1,J,J+1,N,X(1,J),1,SH(J),X(1,J+1),1,MDX,MDC,MDX*MDC)
      IF (MDC.EQ.0)
     1CALL H12(1,J,J+1,N,X(1,J),1,SH(J),SF,1,1,0,1)
   80 CONTINUE
C
C  X CONTAINS NOW THE TRANSFORMED DESIGN MATRIX Q*X.
C  DETERMINE THE PSEUDORANK K USING THE TOLERANCE TAU
C
      DO 100 J=1,LDIAG
      IF (ABS(X(J,J)).LE.TAU) GOTO 110
  100 CONTINUE
      K=LDIAG
      GOTO 120
  110 K=J-1
  120 KP1=K+1
C
C  IF THE PSEUDORANK IS LESS THAN NP STORE THE FIRST K
C  DIAG.ELEMENTS OF X FOR FURTHER APPLICATIONS OF Q
C
      IF (K.EQ.NP) GOTO 130
      DO 125 I=1,K
  125 SF(I)=X(I,I)
  130 CONTINUE
C
C  SPECIAL FOR PSEUDORANK=0
C
      IF (K.GT.0) GOTO 140
      CALL MESSGE(401,'RIMTRF',0)
      RETURN
C
C  IF THE PSEUDORANK IS LESS THAN NP COMPUTE Q*X*V
C
  140 IF (K.EQ.NP) GOTO 160
      MDC=MDX*(NP-1)
      DO 150 II=1,K
      I=KP1-II
  150 CALL H12(1,I,KP1,NP,X(I,1),MDX,SG(I),X,MDX,1,I-1,MDC+I-1)
  160 CONTINUE
      RETURN
      END
C
C-----------------------------------------------------------------------
C
      SUBROUTINE RYWALG(X,Y,THETA,WGT,COV,PSP0,EXPSI,EXCHI,EXRHO,
     1                  SIGMAI,N,NP,MDX,MDT,NCOV,TOL,GAM,TAU,ITYPE,
     1                  ISIGMA,ICNV,MAXIT,MAXIS,NITMON,NIT,SIGMAF,RS,
     1                  DELTA,SC,SF,SG,SH,IP,SW,SX)
C.......................................................................
C
C   COPYRIGHT 1992 Alfio Marazzi
C
C   AUTHOR: A. RANDRIAMIHARISOA
C.......................................................................
C
C   PURPOSE
C   -------
C   W-ALGORITHM FOR ROBUST AND BOUNDED INFLUENCE LINEAR REGRESSION
C
      DIMENSION X(MDX,NP),Y(N),THETA(MDT),WGT(N),COV(NCOV),RS(N),
     1 DELTA(NP),SC(N),SF(NP),SG(NP),SH(NP),SW(N),SX(MDX,NP)
      INTEGER IP(NP)
      LOGICAL NPRCHK
      EXTERNAL EXPSI,EXCHI,EXRHO,ICTHET,ICSIGM
      COMMON/CONST/CONST
      COMMON/BETA/BETA,BET0
      DATA TL/1.E-10/
C     SAVE /CONST/
C
C   PARAMETER CHECK AND INITIALIZATION
C
      MDXP1=MDX+1
      LDIAG=MIN0(N,NP)
      MX=MAX0(N,NP)
      NN=NP*(NP+1)/2
      SIGMA=SIGMAI
      SIGMB=SIGMA
      IASG=IABS(ISIGMA)
      INTCH=1
      NPRCHK=N.LE.MDX.AND.MDT.GE.MX.AND.NCOV.EQ.NN.AND.GAM.GT.0.
     1       .AND.GAM.LT.2..AND.MAXIT.GT.0.AND.(MAXIS.GT.0.OR.
     1        IASG.NE.1).AND.LDIAG.GT.0.AND.SIGMA.GT.0..AND.TOL.GT.0.
     1       .AND.TAU.GE.0..AND.(ITYPE.EQ.1.OR.ITYPE.EQ.2.OR.
     1        ITYPE.EQ.3).AND.(ISIGMA.GE.-2.AND.ISIGMA.LE.2)
     1       .AND.(ICNV.EQ.1.OR.ICNV.EQ.2.OR.ICNV.EQ.3)
      IF (.NOT.NPRCHK) CALL MESSGE(500,'RYWALG',1)
      ITYP=ITYPE
      IF (ITYP.EQ.1) GOTO 15
      N0=N
      E=2.0
      IF (ITYP.EQ.2) E=0.5
      DO 10 I=1,N
        IF (WGT(I).LE.0.) THEN
          SW(I)=-1.
          N0=N0-1
        ELSE
          SW(I)=WGT(I)**E
        ENDIF
   10 CONTINUE
      IF (N0.EQ.0) ITYP=1
   15 IF (IASG.EQ.0) CONST=0.
      IF (IASG.EQ.1) CONST=BETA*FLOAT(N-NP)
      IF (IASG.EQ.2) CONST=BET0*FLOAT(N-NP)
C
C   STEP 1. SET NIT := 1
C   -------
      NIT=1
C
C   STEP 2. COMPUTE RESIDUALS R=Y-X*THETA
C   -------
  200 CALL RESIDU(X,Y,THETA,N,NP,MDX,RS)
C
C   STEP 3. COMPUTE A NEW VALUE SIGMB FOR SIGMA.
C   -------
      IF (ISIGMA.LT.0.AND.NIT.EQ.1) GOTO 300
      IF (ISIGMA.EQ.0) GOTO 300
      SIGMA=SIGMB
      CALL RYSIGM(RS,WGT,EXCHI,SIGMA,N,NP,TOL,ITYP,ISIGMA,MAXIS,
     1            NIS,SIGMB,SW,SC)
      IF (SIGMB.LE.TL) CALL MESSGE(460,'RYWALG',0)
      IF (SIGMB.LE.TL) RETURN
C
C  ITERATION MONITORING
C
      IF (NITMON.LE.0) GOTO 300
      IF (MOD(NIT,NITMON).NE.0) GOTO 300
      CALL QRSS(RS,WGT,SW,EXRHO,N,ITYP,SIGMB,CONST,QS)
      CALL MONITR(NIT,NP,GAM,QS/FLOAT(N),SIGMB,THETA,DELTA)
C
C   STEP 4. COMPUTE WEIGHTS AND APPLY THEM TO X; STORE RESULT IN SX
C   -------
  300 DO 430 I=1,N
        SC(I)=PSP0
        IF (RS(I).EQ.0.) GOTO 410
        T=RS(I)/SIGMB
        IF (ITYP.EQ.1) GOTO 400
        SC(I)=0.
        IF (WGT(I).LE.0.) GOTO 410
        IF (ITYP.EQ.2) GOTO 400
        T=T/WGT(I)
  400   SC(I)=EXPSI(T)/T
  410   PI=SQRT(SC(I))
        IF (ITYP.EQ.2) PI=PI*SW(I)
        RS(I)=PI*RS(I)
        DO 420 J=1,NP
          SX(I,J)=PI*X(I,J)
  420   CONTINUE
  430 CONTINUE
C   STEP 5. SOLVE FOR DELTA
C   -------
C
C   TRIANGULARIZATION OF SX
C
      CALL RIMTRF(SX,N,NP,MDX,INTCH,TAU,K,SF,SG,SH,IP)
      IF (K.EQ.0) CALL MESSGE(461,'RYWALG',0)
      IF (K.EQ.0) RETURN 
C
C   HOUSEHOLDER TRANSFORMATIONS OF THE RIGHT HAND SIDE
C
      KK=MDX*(K-1)+K
      IF (K.NE.NP) CALL SWAP(SX,SF,K,MDXP1,1,KK,K)
      DO 500 JJ=1,LDIAG
        J=JJ
        CALL H12(2,J,J+1,N,SX(1,J),1,SH(J),RS,1,N,1,N)
  500 CONTINUE
      IF (K.NE.NP) CALL SWAP(SX,SF,K,MDXP1,1,KK,K)
C
C   SOLVE FOR DELTA
C
      CALL SOLV(SX,RS,NP,K,MDX,N)
      IF (K.EQ.NP) GOTO 530
      KP1=K+1
      DO 510 J=KP1,NP
        RS(J)=0.0
  510 CONTINUE
      DO 520 J=1,K
        I=J
        CALL H12(2,I,KP1,NP,SX(I,1),MDX,SG(I),RS,1,N,1,NP)
  520 CONTINUE
  530 DO 540 J=1,NP
        DELTA(J)=GAM*RS(J)
  540 CONTINUE
      CALL PERM(DELTA,IP,LDIAG,NP)
C
C   STEP 6. COMPUTE NEW SOLUTION
C   -------
      DO 600 J=1,NP
        THETA(J)=THETA(J)+DELTA(J)
  600 CONTINUE
C
C   STEP 7. STOP ITERATIONS IF DESIRED PRECISION IS REACHED
C   -------
      IF (NIT.EQ.MAXIT) GOTO 800
      IF (ISIGMA.LT.0.AND.NIT.EQ.1) GOTO 700
      IF(ICTHET(NP,NCOV,DELTA,SIGMA,COV,TOL,ICNV).EQ.1
     +   .AND.ICSIGM(SIGMA,SIGMB,TOL).EQ.1) GOTO 800
  700 NIT=NIT+1
      GOTO 200
  800 SIGMAF=SIGMB
      CALL RESIDU(X,Y,THETA,N,NP,MDX,RS)
      RETURN
      END
C
C-----------------------------------------------------------------------
C
      SUBROUTINE KTASKV(X,N,NP,MDX,NCOV,TAU,F,A,COV)
C.......................................................................
C
C   COPYRIGHT 1992 Alfio Marazzi
C
C   AUTHOR : A. MARAZZI
C.......................................................................
C
      DIMENSION X(MDX,NP),A(NCOV),COV(NCOV)
      LOGICAL NPRCHK
      DOUBLE PRECISION SM1,DZERO
C
C  PARAMETER CHECK
C
      NN=NP*(NP+1)/2
      NPRCHK=NP.GT.0.AND.N.GE.NP.AND.N.LE.MDX.AND.NCOV.EQ.NN
     1       .AND.TAU.GE.0.
      IF (.NOT.NPRCHK) CALL MESSGE(500,'KTASKV',1)
 
      DZERO=0.D0
C
C  COMPUTE X**T*X AND STORE IT TEMPORARILY IN COV
C
      L=0
      DO 60 I=1,NP
      DO 60 J=1,I
      L=L+1
      SM1=DZERO
      DO 50 K=1,N
   50 SM1=SM1+DBLE(X(K,I))*X(K,J)
      COV(L)=SM1
   60 CONTINUE
C
C  COMPUTE A LOWER TRIANGULAR MATRIX A SUCH THAT
C  (X**T*X)**(-1)=A**T*A; SET COV=A**T*A.
C
      CALL MCHL(COV,NP,NN,INFO)
      IF (INFO.EQ.0) GOTO 65
      CALL MESSGE(400+INFO,'KTASKV',0)
      RETURN
   65 DO 70 L=1,NN
   70 A(L)=COV(L)
      CALL MINV(A,NP,NN,TAU,ISING)
      IF (ISING.EQ.0) GOTO 75
      CALL MESSGE(450,'KTASKV',0)
      RETURN
   75 CALL MTT1(A,COV,NP,NN)
      IF (F.GT.0.) CALL SCAL(COV,F,NCOV,1,NCOV)
      RETURN
      END
C
C-----------------------------------------------------------------------
C
      SUBROUTINE KFFACV(RS,EXPSI,EXPSP,N,NP,SIGMA,FH)
C.......................................................................
C
C   COPYRIGHT 1992 Alfio Marazzi
C
C   AUTHOR : A. MARAZZI
C.......................................................................
C
C  HUBER'S CORRECTION FACTOR FOR AS. COVARIANCE
C  MATRIX OF PARAMETER ESTIMATES
C
      DIMENSION RS(N)
      EXTERNAL EXPSI,EXPSP
      DATA TL/1.E-10/
C
C  PARAMETER CHECK
C
      IF (N.LT.NP.OR.NP.LE.0) CALL MESSGE(500,'KFFACV',1)
C
      FH=1.
      IF (NP.EQ.N) RETURN
      CALL FACS(RS,N,NP,SIGMA,TL,XKAPPA,SUM2,EXPSI,EXPSP)
      IF (XKAPPA.EQ.0.) CALL MESSGE(301,'KFFACV',0)
      IF (XKAPPA.EQ.0.) RETURN
      FH=(XKAPPA*XKAPPA)*SUM2
      RETURN
      END


C-----------------------------------------------------------------------
C
C                 R O B E T H  Interface to  S - P L U S
C
C  File USREXT.F  Interface for vector and user-defined weight
C                 functions (FORTRAN CODE)
C                 Note: This file also contains a version of SRDPRT.F
C                 used by the Interface. 
C
C-----------------------------------------------------------------------
C
      SUBROUTINE PSIA(N,SVALS,FVALS)
C.......................................................................
C
C   COPYRIGHT 1992 Alfio Marazzi  
C
C   AUTHOR : A. RANDRIAMIHARISOA
C.......................................................................
C
C  PURPOSE
C  -------
C  GIVES THE VALUE OF THE FUNCTION PSI FOR SVALS(I), I=1,...,N
C  This subroutine is mainly used by the S-interface
C
      DIMENSION SVALS(N),FVALS(N)
      COMMON/PSIPR/IPSI,C,H1,H2,H3,XK,D
      IPS=IABS(IPSI)
      IF (IPS.EQ.0) GOTO 100
      IF (IPS.EQ.1) GOTO 200
      IF (IPS.EQ.2) GOTO 300
      IF (IPS.EQ.3) GOTO 400
      IF (IPS.EQ.4) GOTO 500
      IF (IPS.EQ.10) GOTO 700
C
C  PSI(S)=USER PSI FUNCTION
C
C      DO 50 I=1,N
C      S=SVALS(I)
C   50 FVALS(I)=UPSI(S)
C      RETURN
C
C  PSI(S)=S
C
  100 DO 150 I=1,N
  150 FVALS(I)=SVALS(I)
      RETURN
C
C  PSI(S,C)=MAX(-C,MIN(C,S))
C
  200 DO 250 I=1,N
      S=SVALS(I)
      ABST=ABS(S)
      TMP=AMIN1(C,ABST)
      IF (S.LT.0.) TMP=-TMP
  250 FVALS(I)=TMP
      RETURN
C
C  PSI(S,H1,H2,H3)=-PSI(-S,H1,H2,H3)
C                 =S FOR 0 .LE. S .LE. H1
C                 =H1 FOR H1 .LE. S .LE. H2
C                 =H1*(H3-S)/(H3-H2) FOR H2 .LE. S .LE. H3
C                 =0 FOR S .GT. H3
C
  300 DO 350 I=1,N
      S=SVALS(I)
      ABST=ABS(S)
      TMP=0
      IF (ABST.GE.H3) GOTO 350
      IF (ABST.LE.H2) TMP=AMIN1(H1,ABST)
      IF (ABST.GT.H2) TMP=H1*(H3-ABST)/(H3-H2)
      IF (S.LT.0.) TMP=-TMP
  350 FVALS(I)=TMP
      RETURN
C
C  PSI(S)=S*[MAX(1-S**2,0)]**2
C
  400 DO 450 I=1,N
      S=SVALS(I)
      ABST=ABS(S)
      TMP=0.
      IF (ABST.GE.1.) GOTO 450
      TMP=S*(1.-S*S)*(1.-S*S)
  450 FVALS(I)=TMP
      RETURN
C
C  PSI(S)=(6/K)*(S/K)*[MAX(1-(S/K)**2,0)]**2
C
  500 DO 550 I=1,N
      S=SVALS(I)
      ABST=ABS(S)
      TMP=0.
      IF (ABST.GE.XK) GOTO 550
      SK=S/XK
      TMP=(6.*SK/XK)*(1.-SK*SK)*(1.-SK*SK)
  550 FVALS(I)=TMP
      RETURN
C
C  PSI(S,C)=MAX(-C,MIN(C,S))
C

  700 DO 750 I=1,N
      S=SVALS(I)
      TMP=AMIN1(H2,S)
      IF (TMP.LE.H1) TMP=H1
  750 FVALS(I)=TMP
      RETURN
      END
C
C-----------------------------------------------------------------------
C
      SUBROUTINE RHOA(N,SVALS,FVALS)
C.......................................................................
C
C   COPYRIGHT 1992 Alfio Marazzi
C
C   AUTHOR : A. RANDRIAMIHARISOA
C.......................................................................
C
C  PURPOSE
C  -------
C  GIVES THE VALUES OF THE INTEGRAL FROM 0 TO SVALS(I) OF PSI, I=1,..,N
C  This subroutine is mainly used by the S-interface
C
      DIMENSION SVALS(N),FVALS(N)
      COMMON/PSIPR/IPSI,C,H1,H2,H3,XK,D
      IPS=IABS(IPSI)
      IF (IPS.EQ.0) GOTO 100
      IF (IPS.EQ.1) GOTO 200
      IF (IPS.EQ.2) GOTO 300
      IF (IPS.EQ.3) GOTO 400
      IF (IPS.EQ.4) GOTO 500
      IF (IPS.EQ.10) GOTO 700
C      DO 50 I=1,N
C      S=SVALS(I)
C      RHO=URHO(S)
C   50 FVALS(I)=TMP
C      RETURN
  100 DO 150 I=1,N
      S=SVALS(I)
  150 FVALS(I)=S*S/2.
      RETURN
  200 DO 250 I=1,N
      S=SVALS(I)
      ABST=ABS(S)
      TMP=S*S/2.
      IF (ABST.GT.C) TMP=C*(ABST-C/2.)
  250 FVALS(I)=TMP
      RETURN
  300 DO 350 I=1,N
      S=SVALS(I)
      ABST=ABS(S)
      IF (ABST.GT.H2) GOTO 325
      TMP=S*S/2.
      IF (ABST.GT.H1) TMP=H1*(ABST-H1/2.)
  325 TMP=0.5*H1*(H3+H2-H1)
      IF (ABST.LT.H3) TMP=TMP-.5*H1*(H3-ABST)**2/(H3-H2)
  350 FVALS(I)=TMP
      RETURN
  400 DO 450 I=1,N
      S=SVALS(I)
      ABST=ABS(S)
      TMP=1./6.
      IF (ABST.GE.1.) GOTO 450
      S2=S*S
      TMP=(S2*(S2-3)+3)*S2/6.
  450 FVALS(I)=TMP
      RETURN
  500 DO 550 I=1,N
      S=SVALS(I)
      ABST=ABS(S)
      TMP=1.
      IF (ABST.GE.XK) GOTO 550
      S2=(S/XK)**2
      TMP=(S2*(S2-3)+3)*S2
  550 FVALS(I)=TMP
      RETURN
  700 DO 750 I=1,N
      S=SVALS(I)
      TMP=S*S/2.
      IF (S.LT.H1) TMP=H1*(S-H1/2.)
      IF (S.GT.H2) TMP=H2*(S-H2/2.)
  750 FVALS(I)=TMP
      RETURN
      END
C
C-----------------------------------------------------------------------
C
      SUBROUTINE PSPA(N,SVALS,FVALS)
C.......................................................................
C
C   COPYRIGHT 1992 Alfio Marazzi
C
C   AUTHOR : A. RANDRIAMIHARISOA
C.......................................................................
C
C  PURPOSE
C  -------
C  GIVES THE VALUE AT THE POINTS S=SVALS(I), I=1,...,N
C  OF THE FIRST DERIVATIVE OF THE FUNCTION PSI.
C  This subroutine is mainly used by the S-interface
C
      DIMENSION SVALS(N),FVALS(N)
      COMMON/PSIPR/IPSI,C,H1,H2,H3,XK,D
      IPS=IABS(IPSI)
      IF (IPS.EQ.0) GOTO 100
      IF (IPS.EQ.1) GOTO 200
      IF (IPS.EQ.2) GOTO 300
      IF (IPS.EQ.3) GOTO 400
      IF (IPS.EQ.4) GOTO 500
      IF (IPS.EQ.10) GOTO 700
C      DO 50 I=1,N
C      S=SVALS(I)
C   50 FVALS(I)=UPSP(S)
C      RETURN
  100 DO 150 I=1,N
  150 FVALS(I)=1.
      RETURN
  200 DO 250 I=1,N
      S=SVALS(I)
      ABST=ABS(S)
      TMP=0.
      IF (ABST.LE.C) TMP=1.
  250 FVALS(I)=TMP
      RETURN
  300 DO 350 I=1,N
      S=SVALS(I)
      ABST=ABS(S)
      TMP=1.
      IF (ABST.LT.H1) GOTO 350
      TMP=0.
      IF ((ABST.LE.H2).OR.(ABST.GE.H3)) GOTO 350
      TMP=H1/(H2-H3)
  350 FVALS(I)=TMP
      RETURN
  400 DO 450 I=1,N
      S=SVALS(I)
      ABST=ABS(S)
      TMP=0.
      IF (ABST.GE.1.) GOTO 450
      S2=S*S
      TMP=(1.-S2)*(1.-5.*S2)
  450 FVALS(I)=TMP
      RETURN
  500 DO 550 I=1,N
      S=SVALS(I)
      ABST=ABS(S)
      TMP=0.
      IF (ABST.GE.XK) GOTO 550
      S2=(S/XK)**2
      TMP=(6./XK)*(1.-S2)*(1.-5.*S2)/XK
  550 FVALS(I)=TMP
      RETURN
  700 DO 750 I=1,N
      S=SVALS(I)
      TMP=0.
      IF (S.GE.H1.AND.S.LE.H2) TMP=1.
  750 FVALS(I)=TMP
      RETURN
      END
C
C-----------------------------------------------------------------------
C
      SUBROUTINE CHIA(N,SVALS,FVALS)
C.......................................................................
C
C   COPYRIGHT 1992 Alfio Marazzi
C
C   AUTHOR : A. RANDRIAMIHARISOA
C.......................................................................
C
C  PURPOSE
C  -------
C  FOR S=SVALS(I), I=1,...,N, CHIA GIVES THE VALUE OF :
C  THE FUNCTION CHI(S)=S*S/2 IF IPSI=0,
C  THE HUBER'S CHI FUNCTION IF ABS(IPSI)<4, AND
C  THE FUNCTION CHI(S)=CHIK(S) FOR S-ESTIMATES IF IPS=4.
C  This subroutine is mainly used by the S-interface
C
      DIMENSION SVALS(N),FVALS(N)
      COMMON/PSIPR/IPSI,C,H1,H2,H3,XK,D
      IF (IPSI.EQ.0) GOTO 100
      IPS=IABS(IPSI)
      IF (IPS.LE.3) GOTO 200
      IF (IPS.EQ.4) GOTO 300
      IF (IPS.EQ.10) GOTO 400
C      DO 50 I=1,N
C      S=SVALS(I)
C   50 FVALS(I)=UCHI(S)
C      RETURN
  100 DO 150 I=1,N
      S=SVALS(I)
  150 FVALS(I)=S*S/2.
      RETURN
  200 DO 250 I=1,N
      S=SVALS(I)
      ABST=ABS(S)
      PS=AMIN1(D,ABST)
  250 FVALS(I)=PS*PS/2.
      RETURN
  300 DO 350 I=1,N
      S=SVALS(I)
      ABST=ABS(S)
      TMP=1.
      IF (ABST.GE.XK) GOTO 350
      S2=(S/XK)**2
      TMP=(S2*(S2-3)+3)*S2
  350 FVALS(I)=TMP
      RETURN
  400 DO 450 I=1,N
      S=SVALS(I)
      PS=AMIN1(H2,S)
      IF (PS.LT.H1) PS=H1
  450 FVALS(I)=PS*PS/2.
      RETURN
      END
C
      SUBROUTINE MONITR(NIT,NP,GAM,Q,SIGMA,THETA,DELTA)
C.......................................................................
C
C   COPYRIGHT 1992 Alfio Marazzi
C
C   AUTHOR : A. RANDRIAMIHARISOA
C.......................................................................
C
      DIMENSION THETA(NP),DELTA(NP),TMP(2)
      DATA NEXT,INIT/0,0/
      IF (NEXT.NE.NIT) NEXT=0
      IF (NEXT.EQ.0) INIT=NIT
      IF (NEXT.EQ.0) CALL INTPR
     + ('* * * I T E R A T I O N   M O N I T O R I N G * * *',51,0,0)
      NEXT=NIT+INIT
      TMP(1)=Q
      TMP(2)=GAM
      CALL INTPR('Nb of iterations',16,NIT,1)
      CALL REALPR('Qs, Gamma',9,TMP,2)
      CALL REALPR('Theta',5,THETA,NP)
      CALL REALPR('Sigma',5,SIGMA,1)
      CALL REALPR('Delta',5,DELTA,NP)
      RETURN
      END
C
      SUBROUTINE MONITW(NIT,NP,NCOV,A,TOLA)
C.......................................................................
C
C   COPYRIGHT 1992 Alfio Marazzi
C
C   AUTHOR : A. RANDRIAMIHARISOA
C.......................................................................
C
      DOUBLE PRECISION A(NCOV)
      DATA NEXT,INIT/0,0/
      TMP=NP
      IF (NEXT.NE.NIT) NEXT=0
      IF (NEXT.EQ.0) INIT=NIT
      IF (NEXT.EQ.0) CALL INTPR
     + ('* * * I T E R A T I O N   M O N I T O R I N G * * *',51,0,0)
      NEXT=NIT+INIT
      CALL INTPR('Nb of iterations',16,NIT,1)
      CALL REALPR('TOLA',4,TOLA,1)
      CALL DBLEPR('A matrix',8,A,NCOV)
      RETURN
      END
C
      SUBROUTINE MONITC(NIT,NVAR,NCOV,B,A,TOLB,TOLA)
C.......................................................................
C
C   COPYRIGHT 1992 Alfio Marazzi
C
C   AUTHOR : A. RANDRIAMIHARISOA
C.......................................................................
C
      DOUBLE PRECISION A(NCOV)
      REAL B(NVAR),TOL(2)
      DATA NEXT,INIT/0,0/
      TOL(1)=TOLA
      TOL(2)=TOLB
      IF (NEXT.NE.NIT) NEXT=0
      IF (NEXT.EQ.0) INIT=NIT
      IF (NEXT.EQ.0) CALL INTPR
     + ('* * * I T E R A T I O N   M O N I T O R I N G * * *',51,0,0)
      NEXT=NIT+INIT
      CALL INTPR('Nb of iterations',16,NIT,1)
      CALL DBLEPR('A matrix',8,A,NCOV)
      CALL REALPR('B vector',8,B,NVAR)
      RETURN
      END
C
      SUBROUTINE MONITA(NIT,NVAR,NCOV,B,A,TOLB,TOLA)
C.......................................................................
C
C   COPYRIGHT 1992 Alfio Marazzi
C
C   AUTHOR : A. RANDRIAMIHARISOA
C.......................................................................
C
      DOUBLE PRECISION A(NCOV)
      REAL B,TMP(3)
      DATA NEXT,INIT/0,0/
      TMP(1)=TOLA
      TMP(2)=TOLB
      TMP(3)=NVAR
      IF (NEXT.NE.NIT) NEXT=0
      IF (NEXT.EQ.0) INIT=NIT
      IF (NEXT.EQ.0) CALL INTPR
     + ('* * * I T E R A T I O N   M O N I T O R I N G * * *',51,0,0)
      NEXT=NIT+INIT
      CALL INTPR('Nb of iterations',16,NIT,1)
      CALL REALPR('B',1,B,1)
      CALL DBLEPR('A matrix',8,A,NCOV)
      RETURN
      END
C
      SUBROUTINE MONITG(NIT,NP,GAM,Q,THETA,DELTA)
C.......................................................................
C
C   COPYRIGHT 1992 Alfio Marazzi
C
C   AUTHOR : A. RANDRIAMIHARISOA
C.......................................................................
C
      DIMENSION THETA(NP),DELTA(NP),TMP(2)
      COMMON/OUT/IOUT
      DATA NEXT,INIT/0,0/
      IF (NEXT.NE.NIT) NEXT=0
      IF (NEXT.EQ.0) INIT=NIT
      IF (NEXT.EQ.0) CALL INTPR
     + ('* * * I T E R A T I O N   M O N I T O R I N G * * *',51,0,0)
      NEXT=NIT+INIT
      TMP(1)=Q
      TMP(2)=GAM
      CALL INTPR('Nb of iterations',16,NIT,1)
      CALL REALPR('Q0, Gamma',9,TMP,2)
      CALL REALPR('Theta',5,THETA,NP)
      CALL REALPR('Delta',5,DELTA,NP)
      RETURN
      END
C
C-----------------------------------------------------------------------
C
      FUNCTION USERFS(S)
C
C  PURPOSE
C  -------
C  Dummy function (Single precision)
C
      USERFS=S
      RETURN
      END
c
c***********************************************************************
c
      subroutine int51(rs,wgt,exchi,sigmai,n,np,tol,itype,isigma,maxis,
     x                 nit,sigmaf,sw,sc)
C.......................................................................
C
C   COPYRIGHT 1992 Alfio Marazzi
C
C   AUTHOR : J. JOSS
C.......................................................................
C
       integer n,np,exchi,itype,isigma,maxis,nit
       real rs(n),wgt(n),sigmai,tol,sigmaf,sw(n),sc(n)
       external chi,userfs
       if (exchi.eq.4) then
         call rysigm(rs,wgt,chi,sigmai,n,np,tol,itype,isigma,maxis,
     x               nit,sigmaf,sw,sc)
       else
         call rysigm(rs,wgt,userfs,sigmai,n,np,tol,itype,isigma,maxis,
     x               nit,sigmaf,sw,sc)
       endif
       return
      end
c
c***********************************************************************
c
      subroutine int59(s,result)
C.......................................................................
C
C   COPYRIGHT 1992 Alfio Marazzi
C
C   AUTHOR : J. JOSS
C.......................................................................
C
       real s,result,psy
       external psy
       result=psy(s)
      end
c
c***********************************************************************
c
      subroutine int60(s,result)
C.......................................................................
C
C   COPYRIGHT 1992 Alfio Marazzi
C
C   AUTHOR : J. JOSS
C.......................................................................
C
       real s,result,rho
       external rho
       result=rho(s)
      end
c
c***********************************************************************
c
      subroutine int61(s,result)
C.......................................................................
C
C   COPYRIGHT 1992 Alfio Marazzi
C
C   AUTHOR : J. JOSS
C.......................................................................
C
       real s,result,psp
       external psp
       result=psp(s)
      end
c
c***********************************************************************
c
      subroutine int62(s,result)
C.......................................................................
C
C   COPYRIGHT 1992 Alfio Marazzi
C
C   AUTHOR : J. JOSS
C.......................................................................
C
       real s,result,chi
       external chi
       result=chi(s)
      end
cc
cc======================================================================
cc
      subroutine int92(y,theta,psp0,expsi,exchi,exrho,sigmai,
     x                 n,tol,gam,isigma,maxit,maxis,nitmon,
     x                 nit,sigmaf,rs,sc)
C.......................................................................
C
C   COPYRIGHT 1992 Alfio Marazzi
C
C   AUTHOR : J. JOSS
C.......................................................................
C
       integer expsi,exchi,exrho,n,isigma,maxit,maxis,nitmon,nit
       real y(n),theta,psp0,sigmai,tol,gam,sigmaf,rs(n),sc(n)
       external psy,userfs
       if (expsi.eq.1) then
         call int93(y,theta,psp0,psy,exchi,exrho,sigmai,
     x              n,tol,gam,isigma,maxit,maxis,nitmon,
     x              nit,sigmaf,rs,sc)
       else
         call int93(y,theta,psp0,userfs,exchi,exrho,sigmai,
     x              n,tol,gam,isigma,maxit,maxis,nitmon,
     x              nit,sigmaf,rs,sc)
       endif
       return
      end
      subroutine int93(y,theta,psp0,expsi,exchi,exrho,sigmai,
     x                 n,tol,gam,isigma,maxit,maxis,nitmon,
     x                 nit,sigmaf,rs,sc)
C.......................................................................
C
C   COPYRIGHT 1992 Alfio Marazzi
C
C   AUTHOR : J. JOSS
C.......................................................................
C
       integer exchi,exrho,n,isigma,maxit,maxis,nitmon,nit
       real y(n),theta,psp0,sigmai,tol,gam,sigmaf,rs(n),sc(n)
       external expsi,chi,userfs
       if (exchi.eq.4) then
         call int94(y,theta,psp0,expsi,chi,exrho,sigmai,
     x              n,tol,gam,isigma,maxit,maxis,nitmon,
     x              nit,sigmaf,rs,sc)
       else
         call int94(y,theta,psp0,expsi,userfs,exrho,sigmai,
     x              n,tol,gam,isigma,maxit,maxis,nitmon,
     x              nit,sigmaf,rs,sc)
       endif
       return
      end
      subroutine int94(y,theta,psp0,expsi,exchi,exrho,sigmai,
     x                 n,tol,gam,isigma,maxit,maxis,nitmon,
     x                 nit,sigmaf,rs,sc)
C.......................................................................
C
C   COPYRIGHT 1992 Alfio Marazzi
C
C   AUTHOR : J. JOSS
C.......................................................................
C
       integer n,isigma,maxit,maxis,nitmon,nit,exrho 
       real y(n),theta,psp0,sigmai,tol,gam,sigmaf,rs(n),sc(n)
       external expsi,exchi,rho,userfs
       if (exrho.eq.2) then
         call lywalg(y,theta,psp0,expsi,exchi,rho,sigmai,
     x               n,tol,gam,isigma,maxit,maxis,nitmon,
     x               nit,sigmaf,rs,sc)
       else
         call lywalg(y,theta,psp0,expsi,exchi,userfs,sigmai,
     x               n,tol,gam,isigma,maxit,maxis,nitmon,
     x               nit,sigmaf,rs,sc)
       endif
       return
      end
C
C-----------------------------------------------------------------------
C
      SUBROUTINE LYWALG(Y,THETA,PSP0,EXPSI,EXCHI,EXRHO,SIGMAI,
     1                  N,TOL,GAM,ISIGMA,MAXIT,MAXIS,NITMON,
     1                  NIT,SIGMAF,RS,SC)
C.......................................................................
C
C   COPYRIGHT 1996 Alfio Marazzi
C
C   AUTHOR: A. RANDRIAMIHARISOA
C.......................................................................
C
C   PURPOSE
C   -------
C   W-ALGORITHM FOR ROBUST LOCATION (Huber case)
C
      DIMENSION Y(N),RS(N),SC(N),TETA(1),DLTA(1)
      LOGICAL NPRCHK
      EXTERNAL EXPSI,EXCHI,EXRHO
      COMMON/CONST/CONST
      COMMON/BETA/BETA,BET0
      DATA TL/1.E-10/
C     SAVE /CONST/
C
C   PARAMETER CHECK AND INITIALIZATION
C
      SIGMA=SIGMAI
      SIGMB=SIGMA
      IASG=IABS(ISIGMA)
      NPRCHK=GAM.GT.0..AND.GAM.LT.2..AND.MAXIT.GT.0.AND.(MAXIS.GT.0
     1       .OR.IASG.NE.1).AND.SIGMA.GT.0..AND.TOL.GT.0..AND.
     1       (ISIGMA.GE.-2.AND.ISIGMA.LE.2)
      IF (.NOT.NPRCHK) CALL MESSGE(500,'LYWALG',1)
      ITYP=1
      NP=1
      IF (IASG.EQ.0) CONST=0.
      IF (IASG.EQ.1) CONST=BETA*FLOAT(N-1)
      IF (IASG.EQ.2) CONST=BET0*FLOAT(N-1)
C
C   STEP 1. SET NIT := 1
C   -------
      NIT=1
C
C   STEP 2. COMPUTE RESIDUALS R=Y-X*THETA
C   -------
  200 DO 210 I=1,N
  210 RS(I)=Y(I)-THETA
c      call intpr('Nit',3,nit,1)
C
C   STEP 3. COMPUTE A NEW VALUE SIGMB FOR SIGMA.
C   -------
      IF (ISIGMA.LT.0.AND.NIT.EQ.1) GOTO 400
      IF (ISIGMA.EQ.0) GOTO 400
      SIGMA=SIGMB
      CALL RYSIGM(RS,SC,EXCHI,SIGMA,N,1,TOL,ITYP,ISIGMA,MAXIS,
     1            NIS,SIGMB,SC,SC)
      IF (SIGMB.LE.TL) CALL MESSGE(460,'LYWALG',0)
      IF (SIGMB.LE.TL) RETURN
C
C  ITERATION MONITORING
C
      IF (NITMON.LE.0) GOTO 400
      IF (MOD(NIT,NITMON).NE.0) GOTO 400
      CALL QRSS(RS,SC,SC,EXRHO,N,ITYP,SIGMB,CONST,QS)
      TETA(1)=THETA
      DLTA(1)=DELTA
      CALL MONITR(NIT,NP,GAM,QS/FLOAT(N),SIGMB,TETA,DLTA)
C
C   STEP 4. COMPUTE WEIGHTS AND APPLY THEM TO X; STORE RESULT IN SX
C   -------
  400 SUM1=0.
      SUM2=0.
      DO 430 I=1,N
        SC(I)=PSP0
        IF (RS(I).EQ.0.) GOTO 410
        T=RS(I)/SIGMB
        SC(I)=EXPSI(T)/T
  410   S1=SC(I)
        S2=RS(I)*SC(I)                               
        SUM1=SUM1+S1 
        SUM2=SUM2+S2 
  430 CONTINUE
C
C   STEP 5. SOLVE FOR DELTA
C   -------
      IF (ABS(SUM1).LT.1E-6) THEN
        CALL MESSGE(450,'LYWALG',0)
        SUM1=SIGN(1E-6,SUM1)
      ENDIF
      DELTA=GAM*SUM2/SUM1
C
C   STEP 6. COMPUTE NEW SOLUTION
C   -------
      THETA=THETA+DELTA
C
C   STEP 7. STOP ITERATIONS IF DESIRED PRECISION IS REACHED
C   -------
      IF (NIT.EQ.MAXIT) GOTO 800
      IF (ISIGMA.LT.0.AND.NIT.EQ.1) GOTO 700
      IF(ABS(DELTA).LT.TOL*AMAX1(1.,ABS(THETA)).AND.
     +   ABS(SIGMA-SIGMB).LT.TOL*AMAX1(1.,ABS(SIGMB))) GOTO 800
  700 NIT=NIT+1
      GOTO 200
  800 SIGMAF=SIGMB
      DO 900 I=1,N
        RS(I)=Y(I)-THETA
  900 CONTINUE
      RETURN
      END



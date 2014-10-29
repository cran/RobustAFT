      SUBROUTINE SMINCC(K,I,X,Y,DELTA,SIGMA0,MU0,S0,IPSI,XK,B,BETA,
     *           GAMA,CNST,N,NP,NS,MDX,LINT,METH,IALG,MAXIT,TOL,
     *           SIGMA,THETA,RS,YY,DD,SBETA,SGAMA,SX,SZ,SW,IT,MES2) 
C
C.......................................................................
C   PROGRAMMER : A. RANDRIAMIHARISOA
C.......................................................................
C
      DIMENSION X(MDX,NP),Y(N),DELTA(N),BETA(NS,NP),THETA(N),RS(N),
     +          CNST(2),YY(N),DD(N),SBETA(NP),SGAMA(NP),MES2(4),
     +          GAMA(NS,NP),SX(N,NP),SZ(N),SW(N),SIGMA(1),SIGI(1)
      REAL MU0
      INTEGER LINT,IT(N) 
      COMMON/PSIPR/JPSI,C,H1,H2,H3,YK,D
C
      IF (N.LE.0.OR.MDX.LT.N.OR.SIGMA0.LT.0.0.OR.XK.LE.0.0) 
     * CALL MESSGE(500,'SMINCC',1)

          DO 40 J=1,NP
          SBETA(J)=BETA(K,J)
          IF (METH.EQ.5) SBETA(J)=SBETA(J)+GAMA(I,J)
          SGAMA(J)=0.0
   40     CONTINUE
          monit=0
          CNST(2)=SIGMA0/S0
          SIGI(1)=CNST(2)
          NIT=1
          SIGMA(1)=SIGMA0
          IF (LINT.NE.0) CALL SIGAMA(X,Y,DELTA,SIGMA0,MU0,S0,IPSI,XK,B,
     *          SBETA,SGAMA,CNST,N,NP,1,MDX,LINT,4,IALG,MAXIT,TOL,NIT,
     *          SIGMA,THETA,RS,YY,DD,SBETA,SGAMA,SX,SZ,SW,SIGI,IT,MES2)
         DO 50 J=1,NP
          IF (METH.EQ.5) SBETA(J)=BETA(K,J)
          SGAMA(J)=GAMA(I,J)
   50     CONTINUE
          SIGB=SIGMA(1)
          CNST(2)=SIGB
          SIGI(1)=CNST(2)
          NIT=1
          CALL SIGAMA(X,Y,DELTA,SIGMA0,MU0,S0,IPSI,XK,B,SBETA,SGAMA,
     *            CNST,N,NP,1,MDX,LINT,1,IALG,MAXIT,TOL,NIT,SIGMA,
     *            THETA,RS,YY,DD,SBETA,SGAMA,SX,SZ,SW,SIGI,IT,MES2)
          SGAMA(1)=SIGB
          RETURN
          END
C
C------------------------------------------------------------------------
C 
      SUBROUTINE SMINAC(JJJ,X,Y,DELTA,SIGMA0,MU0,S0,IPSI,XK,B,BETA,
     *           GAMA,N,NP,NS,MDX,LGAM,LINT,METH,IALG,MAXIT,TOL,
     *           BETMIN,SMIN,KAPPA,SIGMA,THETA,RS,YY,DD,SBETA,SGAMA,
     *           SX,SZ,SW,SIG5,GNRM,IT,MES2)     
C
C.......................................................................
C   PROGRAMMER : A. RANDRIAMIHARISOA
C.......................................................................
C
      DIMENSION X(MDX,NP),Y(N),DELTA(N),BETA(NS,NP),GAMA(NS,NP),RS(N),
     +    THETA(N),BETMIN(NP),SIGMA(NS),CNST(2),YY(N),DD(N),MES2(4),
     +    SBETA(NP),SGAMA(NP),SX(N,NP),SZ(N),SW(N),SIG5(NS),GNRM(NS),
     +    SCNST(1,2) 
      REAL MU0,KAPPA,SIGMS(1)
      INTEGER LINT,JJJ(NS),IT(N) 
      COMMON/PSIPR/JPSI,C,H1,H2,H3,YK,D
C
      IF (N.LE.0.OR.MDX.LT.N.OR.SIGMA0.LT.0.0.OR.XK.LE.0.0) 
     * CALL MESSGE(500,'SMINAC',1)
      KAPPA=9.E9
      SMINK=9.E9
      MSTOR=1
      KSTOR=1
      CNST(1)=0.0
      IF (METH.EQ.2) CNST(2)=SMIN
      JPSI=IPSI
      YK=XK
      SMAD=0.0
      NETH=METH
      DO 900 K=1,NS
        DO 10 I=1,NP
   10   SBETA(I)=BETA(K,I)
        IF (METH.EQ.1) THEN
c         L1INT=LINT
c         IF (LINT.EQ.0) L1INT=3
          CALL KMEDMAD(X,Y,IT,DELTA,SBETA,N,NP,MDX,1,1,LINT,SCNST,
     +    THETA,RS,DD,YY,SZ,SW)
c         CNST2=SCNST(1,2)
          CNST2=SIGMA0
        ENDIF
        IF (METH.LT.4) THEN
          DO 25 I=1,NS
   25     SIG5(I)=CNST2
        ENDIF
        MS=0
        DO 110 I=1,NS
        DO 30 J=1,NP
        IF (LGAM.EQ.0) THEN
          SGAMA(J)=BETA(I,J)-BETA(K,J)
          GAMA(I,J)=SGAMA(J)
        ELSE
          SGAMA(J)=GAMA(I,J)
        ENDIF
   30   CONTINUE
        CALL NRM2(SGAMA,NP,1,NP,TMP)
        GNRM(I)=TMP
        IF (TMP.LE.KAPPA) THEN
         MS=MS+1
         JJJ(MS)=I
        ELSE
         GNRM(I)=-TMP-1.0
         GOTO 110
        ENDIF
        IF (METH.GE.4) THEN
          CALL SMINCC(K,I,X,Y,DELTA,SIGMA0,MU0,S0,IPSI,XK,B,BETA,
     *         GAMA,CNST,N,NP,NS,MDX,LINT,METH,IALG,MAXIT,TOL,
     *         SIGMS,THETA,RS,YY,DD,SBETA,SGAMA,SX,SZ,SW,IT,MES2)
          SIG5(MS)=SGAMA(1)
          SIGMA(MS)=SIGMS(1) 
        ENDIF          
  110   CONTINUE
        IF (MS.EQ.0) GOTO 900
        IF (METH.LT.4) THEN
         STEST=0.0
         CALL SIGMAJL(JJJ,X,Y,DELTA,SIGMA0,MU0,S0,IPSI,XK,STEST,BETA,
     *      GAMA,CNST,K,B,MS,N,NP,NS,MDX,LINT,NETH,IALG,MAXIT,TOL,NIT,
     *      MES2,SIGMA,THETA,RS,YY,DD,SBETA,SGAMA,SX,SZ,SW,SIG5,IT)
        ENDIF  
      KSTAR=-1  
      OMEGA=1.E9 
      ONMIN=1.E9
      DO 120 J=1,MS
      IF (SIGMA(J).LE.OMEGA) OMEGA=SIGMA(J)
  120 CONTINUE
      DO 130 J=1,MS
      IF (SIGMA(J).GT.OMEGA) GOTO 130
      JM=JJJ(J)
      IF (GNRM(JM).LT.ONMIN) THEN
        MSTOR=J
        KSTOR=JM
        ONMIN=GNRM(JM)
      ENDIF
  130 CONTINUE
      IF (ABS(ONMIN-KAPPA).LT.1.E-6.AND.OMEGA.GT.SMINK) GOTO 900
      MG=0  
      IND=0
      IF (METH.GE.4) IND=MS
      DO 200 I=1,NS    
      IF (GNRM(I).GE.0.0) GOTO 200
      MG=MG+1
      JJJ(IND+MG)=I
  200 CONTINUE
      STEST=OMEGA
      IF (MG.GT.0) THEN
        IF (METH.GE.4) THEN
         DO 250 I=1,MG
         JI=JJJ(IND+I)
         CALL SMINCC(K,JI,X,Y,DELTA,SIGMA0,MU0,S0,IPSI,XK,B,BETA,
     *        GAMA,CNST,N,NP,NS,MDX,LINT,METH,IALG,MAXIT,TOL,
     *        SIGMS,THETA,RS,YY,DD,SBETA,SGAMA,SX,SZ,SW,IT,MES2)
          IF (SIGMS(1).LT.OMEGA) GOTO 900 
  250    CONTINUE
        ELSE
          CALL SIGMAJL(JJJ,X,Y,DELTA,SIGMA0,MU0,S0,IPSI,XK,STEST,BETA,
     *      GAMA,CNST,K,B,MG,N,NP,NS,MDX,LINT,NETH,IALG,MAXIT,TOL,NIT,
     *      MES2,SIGMA,THETA,RS,YY,DD,SBETA,SGAMA,SX,SZ,SW,SIG5,IT)
        ENDIF
      ENDIF
      IF (STEST.GT.0.0) THEN
        SMINK=OMEGA
        KAPPA=ONMIN
        MSTAR=MSTOR
        KSTAR=KSTOR       
        DO 300 J=1,NP
  300   BETMIN(J)=BETA(K,J)
        SMIN2=SMINK
        IF (METH.EQ.1) SMIN=CNST2
        IF (LINT.EQ.0) THEN
           SMIN=SMIN2
C          SMIN2=CNST2
         ENDIF
        IF (METH.GE.4) SMIN=SIG5(MSTAR)
      ENDIF
  900 CONTINUE
      SIGMA(1)=SMIN2
      RETURN
      END
C
C------------------------------------------------------------------------
C 
      SUBROUTINE KMEDMAD(X,Y,IT,DELTA,BETA,N,NP,MDX,NB,MDB,LINT,
     *           FMEDMAD,THETA,RS,DD,YY,SZ,SW)        
C
C.......................................................................
C   PROGRAMMER : A. RANDRIAMIHARISOA
C.......................................................................
C
      REAL X(MDX,NP),Y(N),DELTA(N),THETA(N),RS(N),DD(N),SZ(N),SW(N)
      REAL BETA(MDB,NP),FMEDMAD(NB,2),YY(N),MEDF0(3),MADF0(3)
      INTEGER IT(N)
      DATA MEDF0/0.3665129,-0.3665129,0/
      DATA MADF0/0.7670493,0.7670493,0.6745/
C
C   Compute MED and MAD for Kaplan Meir (Version A. Marazzi)
C
      IF (N.LE.0.OR.MDX.LT.N.OR.NP.LE.0) CALL MESSGE(500,'KMEDMAD',1)
C
C STEP 0: INITIALIZATIONS
C ------

      ND1=0
      NUP=0
      IIR=-1
      DO 10 K=1,N
      IF (DELTA(K).NE.0.0) THEN
        ND1=ND1+1
c       YY(ND1)=FLOAT(K)
      ENDIF
      THETA(K)=0.0 
      YY(K)=0.0
      SZ(K)=0.0
      DD(K)=0.0
   10 CONTINUE
      LI=LINT
      IF (LI.EQ.0) LI=3
      DO 777 II=1,NB 
      FMEDMAD(II,1)=0.0
      FMEDMAD(II,2)=1.0
      I0=ND1
      I1=0
      RS1MAX=-9.E9
      RS0MAX=-9.E9
      DO 220 I=1,N
      RI=Y(I)
      DO 210 J=1,NP
      RI=RI-BETA(II,J)*X(I,J)
  210 CONTINUE
      IF (DELTA(I).EQ.1.0) THEN
        I1=I1+1
        RS(I1)=RI
        SW(I1)=FLOAT(I)
        IF (RI.GE.RS1MAX) RS1MAX=RI
c       YY(I)=I
      ELSE
        I0=I0+1
        RS(I0)=RI
        IF (RI.LT.RS0MAX) GOTO 220 
        IR0=I0
        IIR=I
        RS0MAX=RI
      ENDIF
  220 CONTINUE
      IF (RS0MAX.GT.RS1MAX+1.e-5) THEN
         ND1=ND1+1
         I1=I1+1
         TMP=RS(ND1)
         RS(ND1)=RS0MAX
         RS(IR0)=TMP
       ENDIF
      CALL SRT2(RS,SW,ND1,1,ND1)
      RSMAX=RS(ND1)
c     CALL SRT2(RS,YY,ND1,1,ND1) 
C     KpMr(RS(1:I1),RS(I1+1,N))
      NU=1
      TI=RS(1)
      SZ(1)=TI
      IT(1)=INT(SW(1))
C     YY(NU)=YY(1)
      THETA(1)=1.0
      DO 225 I=2,ND1
      IF (ABS(RS(I)-TI).LT.1.E-5) THEN
        THETA(NU)=THETA(NU)+1.0
      ELSE
        NU=NU+1
        THETA(NU)=1.0
      ENDIF
      SZ(NU)=RS(I)
      IT(NU)=INT(SW(I))
      TI=RS(I)
  225 CONTINUE
      N1=N
      IF (ND1.LT.N) THEN
      DO 230 I=ND1+1,N
      IF (RS(I).LT.RS(1)) N1=N1-1
      IF (RS(I).GT.RSMAX) NUP=NUP+1
  230 CONTINUE
      ENDIF
      DO 240 I=1,NU-1
      SW(I)=FLOAT(IT(I))
      IT(I)=0
      IF (ND1.LT.N) THEN
      DO 235 J=ND1+1,N
      RSJ=RS(J)+1.E-5
      IF (SZ(I).LE.RSJ.AND.RSJ.LT.SZ(I+1)) IT(I)=IT(I)+1
  235 CONTINUE
      ENDIF
  240 CONTINUE   
      SW(NU)=FLOAT(IT(NU))
      IT(NU)=0 
      CS1=0.0
      CS2=0.0
      DO 250 I=1,NU
      YY(I)=FLOAT(N1)-CS1-CS2
c     IF (I.EQ.NU) GOTO 250
      CS1=CS1+THETA(I)
      CS2=CS2+FLOAT(IT(I))
  250 CONTINUE
      CP=1.0
      DO 260 I=1,NU     
      IF (ABS(CP).GT.1.E-6) THEN
        YY(I)=CP*(YY(I)-THETA(I))/YY(I)
      ELSE
        YY(I)=0.0
      ENDIF
      CP=YY(I)
  260 CONTINUE
      T0=0.0
      SUMDD=0.0
      DO 270 I=1,NU
      THETA(I)=1.0-YY(I)-T0
      DD(I)=THETA(I)
      T0=1.0-YY(I)
      IF (I.NE.NU) SUMDD=SUMDD+DD(I)
  270 CONTINUE
      IF (NU.LT.N) THETA(NU+1)=YY(NU)
      THETA(NU)=THETA(NU)+YY(NU)
      DD(NU)=1.0-SUMDD
      IT(1)=NU
c   k=nu=it(1)
c   tu=sz
c   Ft=YY
c   wkm=DD

      IF (LINT.EQ.0.AND.NB.EQ.1) THEN
C     SET Ft(n)=FT(n-1) and Ft(1)=1.0
        NU1=NU-1
        DO 280 I=1,NU1
        INU=NU-I+1
        YY(INU)=YY(INU-1)
  280   CONTINUE
        YY(1)=1.0
        DD(N)=DD(NU)
        SZ(N)=FLOAT(IIR)
        IF (IIR.GT.0) SW(NU)=FLOAT(IIR)
        RETURN
      ENDIF
      DO 290 I=1,NU
      SW(I)=1.0-YY(I)
  290 CONTINUE
      I=0
  300 I=I+1
      IF (SW(I).GE.0.5) GOTO 310 
      IF (I.LT.NU) GOTO 300
  310 FKMED=SZ(I)
      IMED=I
      DO 320 I=1,NU
      RS(I)=ABS(SZ(I)-FKMED)
  320 CONTINUE
      CALL SRT1(RS,NU,1,NU)
      DO 380 I=1,NU
      XF=RS(I)+FKMED
      SUM1=0.0
      DO 340 J=1,NU
        IF (SZ(J).GT.XF) GOTO 340
        SUM1=SUM1+THETA(J)
  340 CONTINUE
      XF=-RS(I)+FKMED
      SUM2=0.0
      DO 360 J=1,NU
        IF (SZ(J).GT.XF) GOTO 360
        SUM2=SUM2+THETA(J)
  360 CONTINUE   
      SW(I)=SUM1-SUM2
  380 CONTINUE    

      I=0
  400 I=I+1
      IF (SW(I).GE.0.5) GOTO 410 
      IF (I.LT.NU) GOTO 400
  410 FKMAD=RS(I)
      FKMAD=FKMAD/MADF0(LI)
      FKMED=FKMED-MEDF0(LI)*FKMAD
      FMEDMAD(II,1)=FKMED
      FMEDMAD(II,2)=FKMAD
c       if (ii.eq.1) call realpr('fkmad',4,fkmad,1)
  777 CONTINUE 
      RETURN
      END
C
C----------------------------------------------------------------------
C
      SUBROUTINE NTRP0L(YMBX,NU,TAB,IND)
      REAL TAB(NU)
C     SEARCH THE INDEX IN THE NON-CENSORED TAB VALUE
C
      NUP=NU+1
      IF (YMBX.LT.TAB(1)-1.E-5) THEN
        IND=0
      ELSEIF (YMBX.GT.TAB(NU)+1.E-5) THEN
        IND=NUP
        NUP=NUP+1
      ELSE 
        DO 100 I=2,NU
C       IF (YMBX.GE.TAB(I)) GOTO 100
C       IND=I-1
        IF (ABS(YMBX-TAB(I)).LE.1.E-5) THEN
          IND=I
          RETURN
        ENDIF
        IF (YMBX.GT.TAB(I)) GOTO 100
        IND=I-1
        RETURN 
  100   CONTINUE
        IND=NU
      ENDIF
      RETURN
      END
C
C------------------------------------------------------------------------
C 
      SUBROUTINE SIGSCENS(X,Y,DELTA,SIG,MU0,S0,BETA,GAMMA,B,N,NP,MDX,
     *           METH,NIT,SIGMA,SIGBET,THETA,RS,DD,YY,SX,SZ,SW,IT,EQB)        
C
C.......................................................................
C   PROGRAMMER : A. RANDRIAMIHARISOA
C.......................................................................
C
      DIMENSION X(MDX,NP),Y(N),DELTA(N),BETA(NP),GAMMA(NP),THETA(N),
     + RS(N),DD(N),YY(N),SX(N,NP),SZ(N),SW(N),SIGBET(2),DUMMY(1,2)
      DOUBLE PRECISION SUMI,SUMJ
      INTEGER IT(N)
      REAL MU0
      EXTERNAL CHI
      COMMON/PSIPR/IPSI,C,H1,H2,H3,XK,D
      DATA NU/0/
C
      IF (N.LE.0.OR.MDX.LT.N.OR.NP.LE.0) CALL MESSGE(500,'SIGSCENS',1)
      NB=1
      LINT=0
C     BET1=SIGBET(1)   30/10/07
      SIGMN=SIGBET(2)
      IF (NIT.GT.1) GOTO 300
        TMP=MU0*FLOAT(METH)+S0+SIGBET(1) ! not used 
        DO 10 J=1,NP
   10   SX(1,J)=BETA(J)
        CALL KMEDMAD(X,Y,IT,DELTA,SX,N,NP,MDX,NB,N,LINT,DUMMY,THETA,
     *               RS,DD,YY,SZ,SW)
        NU=IT(1)
      DO 160 I=1,N
      TMP=Y(I)
      DO 140 J=1,NP
      TMP=TMP-BETA(J)*X(I,J)
  140 CONTINUE
      SW(I)=TMP  !-SIGMN*MU0/S0
c-    RS(I)=FLOAT(I)
  160 CONTINUE
c-    CALL SRT2(SW,RS,N,1,N)
c     DO 180 I=1,N
c     CALL NTRP0L(SW(I),NU,SZ,IR)
c     IT(I)=IR
c 180 CONTINUE
      DO 220 I=1,N
c-    IO=INT(RS(I))
      GMX=0.0
      DO 200 J=1,NP
c-    XIO=X(IO,J)
c-    SX(I,J)=XIO
      GMX=GMX+GAMMA(J)*X(I,J) !XIO
  200 CONTINUE
C     RS(I)=DELTA(IO)
C     SW(I)=SW(I)-FKMED   
      THETA(I)=GMX
  220 CONTINUE
C
  300 CONTINUE
      SUMI=0.D0
      JJO=0
      SWI=-9.E9
      DO 370 I=1,N
      TMP=(SW(I)-THETA(I))/SIG
C     IF (METH.GE.4) SIGMN=SIG
c-    IS=INT(RS(I))
      IF (DELTA(I).NE.0.0) THEN
C       IF (METH.EQ.3) TMP=SW(I)-MU0*SIG/S0
        CI=CHI(TMP)
        SUMI=SUMI+DBLE(CI)
      ELSE  
        CALL NTRP0L(SW(I),NU,SZ,II)      
C       II=IT(I)
        IF (II.GE.NU-1.OR.TMP.GE.XK) THEN
         IF (II.GE.NU-1) TMP=(SZ(NU)-THETA(I))/SIG
         CI=CHI(TMP)
         SUMI=SUMI+DBLE(CI)
         GOTO 370
c       ELSE
c        IF (JJO.EQ.II.AND.ABS(SW(I)-SWI).LT.1.E-5) GOTO 370
c        JJO=II
c        SWI=SW(I)
        ENDIF
        AI=YY(II+1)  ! No more Tail censured obs. (delta(n)=-1)
        IF (ABS(AI).LT.0.00001) AI=1./FLOAT(N)
        SUMJ=0.D0
        DO 350 J=II+1,NU
        TMP=(SZ(J)-THETA(I))/SIG 
        IF (TMP.GE.XK.AND.SUMJ.EQ.0.D0) THEN
          SUMI=SUMI+1.D0
          GOTO 370
        ENDIF
        CI=CHI(TMP)
        PIJ=DD(J)
        SUMJ=SUMJ+DBLE(CI*PIJ)
  350   CONTINUE
        SUMI=SUMI+SUMJ/DBLE(AI)
      ENDIF 
C 360 continue 
  370 CONTINUE
  380 SUMI=SUMI/DFLOAT(N)
      EQB=SNGL(SUMI)
      SIGMA=SIG*SQRT(EQB/B)
      RETURN
      END
C
C------------------------------------------------------------------------
C 
      SUBROUTINE FSIGMA(X,Y,DELTA,SIGMI,MU0,S0,B,CNST,N,NP,MDX,
     *           LINT,METH,NIT,SIGMAF,THETA,RS,YY,DD,
     *           SBETA,SGAMA,SX,SZ,SW,IT,EQB)     
C
C.......................................................................
C   PROGRAMMER : A. RANDRIAMIHARISOA
C.......................................................................
C 
      DIMENSION X(MDX,NP),Y(N),DELTA(N),THETA(N),CNST(2),RS(N),
     +          YY(N),DD(N),SBETA(NP),SGAMA(NP),SX(N,NP),SZ(N),SW(N)
      REAL MU0
      INTEGER LINT,IT(N)
      COMMON/PSIPR/JPSI,C,H1,H2,H3,YK,D
      IF (METH.EQ.3) THEN
       CNST(1)=-SIGMI*MU0/S0  ! => vi = residuals
       CNST(2)=SIGMI/S0 
      ENDIF
      IF (METH.GE.4) CNST(2)=SIGMI
      IF (LINT.EQ.0) THEN 
c       XK=1.5477
        CALL SIGSCENS(X,Y,DELTA,SIGMI,MU0,S0,SBETA,SGAMA,B,N,NP,MDX,
     *       METH,NIT,SIGMAF,CNST,THETA,RS,DD,YY,SX,SZ,SW,IT,EQB)
      ELSEIF (LINT.EQ.3) THEN
c       XK=1.5477
        CALL SIGSNRM(X,Y,DELTA,SIGMI,MU0,S0,SBETA,SGAMA,B,N,NP,MDX,
     *               METH,NIT,SIGMAF,CNST,RS,SX,SZ,SW,EQB)
      ELSE
c       XK=1.717816
        CALL SIGSGMB(X,Y,DELTA,SIGMI,MU0,S0,SBETA,SGAMA,B,N,NP,MDX,LINT,
     *               METH,NIT,SIGMAF,CNST,RS,SX,SZ,SW,EQB)
      ENDIF
      RETURN
      END
C
C------------------------------------------------------------------------
C 
      SUBROUTINE BISIGAM(X,Y,DELTA,SIG,MU0,S0,BB,CNST,N,NP,MDX,
     *           LINT,METH,IALG,TOL,MAXIT,SIGMA,THETA,RS,YY,DD,
     *           SBETA,SGAMA,SX,SZ,SW,IT,ITR,ITERM)     
C
C.......................................................................
C   PROGRAMMER : A. RANDRIAMIHARISOA
C   IALG=1 (Fixpoint), IALG=2 (Regula falsi), IALG=3 (Bisection)
C.......................................................................
C
      DIMENSION X(MDX,NP),Y(N),DELTA(N),THETA(N),CNST(2),RS(N),
     +          YY(N),DD(N),SBETA(NP),SGAMA(NP),SX(N,NP),SZ(N),SW(N)

      REAL MU0
      INTEGER LINT,IT(N)
      COMMON/PSIPR/JPSI,C,H1,H2,H3,YK,D
C      
      IF (N.LE.0 .OR. MDX.LT.N .OR. NP.LE.0) 
     * CALL MESSGE(500,'BISIGAM',1)

      DATA TL/0.0001/
C
C  INITIALIZE
C
      ITR=1
      NIT=1
      A=SIG
      IT0=0
      INFA=0
      ITERM=1
   5  CONTINUE
      CALL FSIGMA(X,Y,DELTA,A,MU0,S0,BB,CNST,N,NP,MDX,LINT,METH,
     *     NIT,SIGMA,THETA,RS,YY,DD,SBETA,SGAMA,SX,SZ,SW,IT,EQB)           
      FA=EQB-BB
      IF (A.LE.1.E-5) CALL MESSGE(500,'BISIGAM',1)
      IF (FA.GE.1.E6) THEN
        A=A*0.5
        IF (INFA.EQ.0) INFA=1
        GOTO 5
      ENDIF
      SIGMA=A
      IF (ABS(FA).LT.TOL) RETURN
      IF (NIT.EQ.2) GOTO 15
      NIT=2
      B=2.0*SIG
      INFB=0
      IF (INFA.NE.0) B=A+2.0
   10 CONTINUE
      CALL FSIGMA(X,Y,DELTA,B,MU0,S0,BB,CNST,N,NP,MDX,LINT,METH,
     *     NIT,SIGMA,THETA,RS,YY,DD,SBETA,SGAMA,SX,SZ,SW,IT,EQB)      
      FB=EQB-BB
      SIGMA=B
      IF (ABS(FB).LT.TOL) RETURN
      IF (FB.GE.1.E6) THEN
        B=B-0.1
        IF (INFB.EQ.0) INFB=1
        IF (B.LE.0.0) CALL MESSGE(500,'BISIGAM',1)
        GOTO 10
      ENDIF
   15 CONTINUE
      IF (FA*FB.LT.0.0) GOTO 20
      IT0=IT0+1
      IF (FB.GT.FA) THEN
        IF (FB.LT.0.0) THEN
          A=B
          FA=FB
          IF (INFB.EQ.1) THEN
            B=AMIN1(B+1.0,2.0*B)
          ELSE  
            B=2.0*B
          ENDIF
          LBL=2
        ELSE
          B=A
          FB=FA
          A=A*0.5
          LBL=1
        ENDIF
      ELSE
        IF (FA.LT.0.0) THEN 
          B=A
          FB=FA
          A=0.5*A 
          IF (ABS(B-A).LT.TOL) A=0.5*A
          LBL=1
        ELSE
          A=B
          FA=FB
          IF (INFB.EQ.1) THEN
            B=AMIN1(B+1.0,2.0*B)
          ELSE  
            B=2.0*B
          ENDIF
          LBL=2
        ENDIF        
      ENDIF
      IF (IT0.EQ.5.AND.LINT.EQ.0) THEN
        NIT=1
        ITR=1
        A=A+5.0
        GOTO 5
      ENDIF
      IF (IT0.EQ.12.AND.LINT.NE.0) THEN
        ITERM=3
        RETURN
      ENDIF
      IF (IT0.EQ.17.AND.LINT.EQ.0) THEN
        ITERM=3
        RETURN
      ENDIF
      GOTO (5,10) LBL
C
C  REGULA FALSI ITERATION
C
   20 IF (ABS(FA-FB).GT.TL) GOTO 30
      SIGMA=0.5*(A+B)
      ITERM=4
      RETURN
   30 XN=0.5*(A+B)
      IF (IALG.EQ.2) XN=(A*FB-B*FA)/(FB-FA)
      CALL FSIGMA(X,Y,DELTA,XN,MU0,S0,BB,CNST,N,NP,MDX,LINT,METH,
     *     NIT,SIGMA,THETA,RS,YY,DD,SBETA,SGAMA,SX,SZ,SW,IT,EQB)     
      FN=EQB-BB
C
C  TEST TO SEE IF MAXIMUM NUMBER OF ITERATIONS HAS BEEN EXECUTED
C
      IF (ITR.GE.MAXIT) GOTO 60
C
C  TEST TO SEE IF ROOT HAS BEEN FOUND
C
      IF (ABS(FN).LT.TOL.AND.ABS(A-B).LE.TOL) GOTO 70
      IF (FA*FN.LE.0.0) GOTO 40
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
      SIGMA=XN
      RETURN
   70 ITERM=1
      SIGMA=XN
      RETURN
      END
C
C------------------------------------------------------------------------
C 
      SUBROUTINE SIGAMA(X,Y,DELTA,SIG,MU0,S0,IPSI,XK,B,BETA,GAMMA,CNST,
     *           N,NP,NS,MDX,LINT,METH,IALG,MAXIT,TOL,NIT,SIGMA,
     *           THETA,RS,YY,DD,SBETA,SGAMA,SX,SZ,SW,SIG5,IT,MES2)     
C
C.......................................................................
C   PROGRAMMER : A. RANDRIAMIHARISOA
C.......................................................................
C
C  BETA and GAMMA are INPUT VARIABLES
  
      DIMENSION X(MDX,NP),Y(N),DELTA(N),BETA(NS,NP),GAMMA(NS,NP),
     + THETA(N),SIGMA(NS),CNST(2),RS(N),YY(N),DD(N),SBETA(NP),
     + SGAMA(NP),SX(N,NP),SZ(N),SW(N),SIG5(NS),MES2(4)
      REAL MU0
      INTEGER LINT,IT(N)
      COMMON/PSIPR/JPSI,C,H1,H2,H3,YK,D
C      
      IF (N.LE.0 .OR. MDX.LT.N .OR. NP.LE.0 .OR. NS.LE.0) 
     * CALL MESSGE(500,'SIGAMA',1)
      CNST(1)=0.0
      JPSI=IPSI
      YK=XK

      DO 500 J=1,NS 
      NIT=0

      DO 100 I=1,NP
      SBETA(I)=BETA(J,I)
      SGAMA(I)=GAMMA(J,I)
 100  CONTINUE
      SIGMI=SIG
      CNST(2)=SIG5(J)
      IF (METH.EQ.3) THEN
       CNST(1)=-SIGMI*MU0/S0  ! => vi = residuals
       CNST(2)=SIGMI/S0 
      ENDIF
      IF (METH.GE.4) CNST(2)=SIGMI
      IF (IALG.NE.1) THEN
         CALL BISIGAM(X,Y,DELTA,SIGMI,MU0,S0,B,CNST,N,NP,MDX,LINT,
     +   METH,IALG,TOL,MAXIT,SIGMAF,THETA,RS,YY,DD,SBETA,SGAMA,
     +   SX,SZ,SW,IT,ITR,ITERM)
         MES2(ITERM)=MES2(ITERM)+1
         NIT=ITR 
         GOTO 400
      ENDIF
 300  NIT=NIT+1
c     B=0.5
      IF (LINT.EQ.0) THEN 
c       XK=1.5477
        CALL SIGSCENS(X,Y,DELTA,SIGMI,MU0,S0,SBETA,SGAMA,B,N,NP,MDX,
     *       METH,NIT,SIGMAF,CNST,THETA,RS,DD,YY,SX,SZ,SW,IT,EQB)
      ELSEIF (LINT.EQ.3) THEN
c       XK=1.5477
        CALL SIGSNRM(X,Y,DELTA,SIGMI,MU0,S0,SBETA,SGAMA,B,N,NP,MDX,
     *               METH,NIT,SIGMAF,CNST,RS,SX,SZ,SW,EQB)
      ELSE
c       XK=1.717816
        CALL SIGSGMB(X,Y,DELTA,SIGMI,MU0,S0,SBETA,SGAMA,B,N,NP,MDX,
     *               LINT,METH,NIT,SIGMAF,CNST,RS,SX,SZ,SW,EQB)
      ENDIF
      IF (ABS(EQB-B).LT.TOL.AND.ABS(SIGMI-SIGMAF).LE.TOL) THEN
        MES2(1)=MES2(1)+1
        GOTO 400
      ENDIF
      IF (NIT.EQ.MAXIT) THEN
        MES2(2)=MES2(2)+1
        GOTO 400
      ENDIF
      SIGMI=SIGMAF
      IF (METH.GE.3) CNST(2)=SIGMI/S0
      GOTO 300
 400  SIGMA(J)=SIGMAF
 500  CONTINUE
      RETURN
      END
C
C------------------------------------------------------------------------
C 
      SUBROUTINE SIGMAJL(JJJ,X,Y,DELTA,SIG,MU0,S0,IPSI,XK,STEST,BETA,
     *      GAMA,CNST,K,B,NJ,N,NP,NS,MDX,LINT,METH,IALG,MAXIT,TOL,NIT,
     *      MES2,SIGMA,THETA,RS,YY,DD,SBETA,SGAMA,SX,SZ,SW,SIG5,IT)     
C
C.......................................................................
C   PROGRAMMER : A. RANDRIAMIHARISOA
C.......................................................................
C
      DIMENSION X(MDX,NP),Y(N),DELTA(N),BETA(NS,NP),GAMA(NS,NP),RS(N),
     + THETA(N),SIGMA(NJ),CNST(2),YY(N),DD(N),SBETA(NP),SGAMA(NP),
     + SX(N,NP),SZ(N),SW(N),SIG5(NS),MES2(4) 
      REAL MU0
      INTEGER LINT,JJJ(NJ),IT(N)
      COMMON/PSIPR/JPSI,C,H1,H2,H3,YK,D
C      
      IF (N.LE.0 .OR. MDX.LT.N .OR. NP.LE.0 .OR. NJ.LE.0) 
     * CALL MESSGE(500,'SIGMAJL',1)
      CNST(1)=0.0
      JPSI=IPSI  !4
      YK=XK
      DO 500 JJ=1,NJ
      NIT=0
      J=JJJ(JJ)       
c      iver(jj)=j  
      DO 100 I=1,NP
      SBETA(I)=BETA(K,I)
 100  SGAMA(I)=GAMA(J,I)
      SIGMI=SIG
      CNST(2)=SIG5(J)
      IF (METH.EQ.3) THEN
       CNST(1)=-SIGMI*MU0/S0  ! => vi = residuals
       CNST(2)=SIGMI/S0 
      ENDIF
      IF (METH.GE.4) CNST(2)=SIGMI
      IF (IALG.NE.1) THEN
         CALL BISIGAM(X,Y,DELTA,SIGMI,MU0,S0,B,CNST,N,NP,MDX,LINT,
     +   METH,IALG,TOL,MAXIT,SIGMAF,THETA,RS,YY,DD,SBETA,SGAMA,
     +   SX,SZ,SW,IT,ITR,ITERM)
         MES2(ITERM)=MES2(ITERM)+1
         NIT=ITR
         GOTO 400
      ENDIF
 300  NIT=NIT+1
c     B=0.5
      IF (LINT.EQ.0) THEN 
c       XK=1.5477
        CALL SIGSCENS(X,Y,DELTA,SIGMI,MU0,S0,SBETA,SGAMA,B,N,NP,MDX,
     *       METH,NIT,SIGMAF,CNST,THETA,RS,DD,YY,SX,SZ,SW,IT,EQB)
      ELSEIF (LINT.EQ.3) THEN
c       XK=1.5477
        CALL SIGSNRM(X,Y,DELTA,SIGMI,MU0,S0,SBETA,SGAMA,B,N,NP,MDX,
     *               METH,NIT,SIGMAF,CNST,RS,SX,SZ,SW,EQB)
      ELSE
c       XK=1.717816
        CALL SIGSGMB(X,Y,DELTA,SIGMI,MU0,S0,SBETA,SGAMA,B,N,NP,MDX,LINT,
     *               METH,NIT,SIGMAF,CNST,RS,SX,SZ,SW,EQB)
      ENDIF
      IF (ABS(EQB-B).LT.TOL.AND.ABS(SIGMI-SIGMAF).LE.TOL) THEN
        MES2(1)=MES2(1)+1
        GOTO 400
      ENDIF
      IF (NIT.EQ.MAXIT) THEN
        MES2(2)=MES2(2)+1
        GOTO 400
      ENDIF
      SIGMI=SIGMAF
      IF (METH.GE.3) CNST(2)=SIGMI/S0
      GOTO 300
 400  SIGMA(JJ)=SIGMAF
      IF (STEST.EQ.0.0) GOTO 500
      IF (SIGMAF.GE.STEST) GOTO 500
      STEST=SIGMAF-STEST
      GOTO 900
 500  CONTINUE
 900  CONTINUE
      RETURN
      END
C
      SUBROUTINE SIGSNRM(X,Y,DELTA,SIG,MU0,S0,BETA,GAMMA,B,N,NP,MDX,
     *           METH,NIT,SIGMA,SIGBET,RS,SX,SZ,SW,EQB)        
C
C.......................................................................
C   PROGRAMMER : A. RANDRIAMIHARISOA
C.......................................................................
C
      DIMENSION X(MDX,NP),Y(N),DELTA(N),BETA(NP),GAMMA(NP),
     +          RS(N),SX(N,NP),SZ(N),SW(N),SIGBET(2),WGT(3)
      DOUBLE PRECISION SUMI,SUMJ
      REAL MU0
      EXTERNAL CHI,FGAUSS
      COMMON/PSIPR/IPSI,C,H1,H2,H3,XK,D
C
      IF (N.LE.0.OR.MDX.LT.N.OR.NP.LE.0) CALL MESSGE(500,'SIGSNRM',1)
      SUMI=0.D0
      BET1=SIGBET(1)
      SIGMN=SIGBET(2)
      IF (NIT.NE.1) GOTO 300
      MU0=0.0*METH/S0        
      DO 160 I=1,N
      TMP=Y(I)
      GMX=0.0
      DO 140 J=1,NP
      TMP=TMP-BETA(J)*X(I,J)
      GMX=GMX+GAMMA(J)*X(I,J)
  140 CONTINUE
      RS(I)=TMP-BET1
      SW(I)=GMX
  160 CONTINUE
      SX(1,1)=0.0
      SZ(1)=0.0
c     CALL SRT2(SW,RS,N,1,N)
c      DO 220 I=1,N
c     IO=RS(I)
c      DO 200 J=1,NP
c      SX(I,J)=X(IO,J)
c  200 CONTINUE
c      SZ(I)=DELTA(IO)
c  220 CONTINUE
C
  300 CONTINUE
      DO 370 I=1,N
      VI=RS(I)
c     IF (METH.EQ.3) VI=VI-SIG*MU0/S0 
      ci=0.0 
      ai=0.0
      sumj=0.D0
      ak=0.0
      TMP=(VI-SW(I))/SIG
      IF (METH.GE.4) SIGMN=SIG
      IF (DELTA(I).NE.0.0) THEN
        IF (METH.GE.4) TMP=TMP/S0
        CI=CHI(TMP)
        SUMI=SUMI+DBLE(CI)
      ELSE 
        IF (TMP.GE.XK) THEN
          SUMI=SUMI+1.D0
          GOTO 370
        ENDIF
        TT=VI/SIGMN
        AI=1.0001-FGAUSS(TT)
        WGT(1)=SW(I)
        WGT(2)=SIG
        IF (METH.GE.4) WGT(2)=S0*SIG
        WGT(3)=SIGMN
        CALL RHONRM(VI,WGT,SUMJ)
        SUMI=SUMI+SUMJ/DBLE(AI)
      ENDIF    
C 360 continue
  370 CONTINUE
      SUMJ=SUMI/DFLOAT(N)
      EQB=SNGL(SUMJ)
      SIGMA=SIG*SQRT(EQB/B)
      RETURN
      END
C
      SUBROUTINE SIGSGMB(X,Y,DELTA,SIG,MU0,S0,BETA,GAMMA,B,N,NP,MDX,
     *           LINT,METH,NIT,SIGMA,SIGBET,RS,SX,SZ,SW,EQB)        
C
C.......................................................................
C   PROGRAMMER : A. RANDRIAMIHARISOA
C   Gumbel Type I => f(x)=exp(-x) exp(-exp(-x)) => lint=1
C   Gumbel Type II => f(x)=exp(x) exp(-exp(x))  => lint=2
C.......................................................................
C
      DIMENSION X(MDX,NP),Y(N),DELTA(N),BETA(NP),GAMMA(NP),
     +          RS(N),SX(N,NP),SZ(N),SW(N),SIGBET(2),WGT(5)
c      dimension vervi(75),vera(75),verk(75),verak(75)  
c      logical iver
      DOUBLE PRECISION SUMI,SUMJ,TT,AA,FGUMBL
      REAL MU0
      EXTERNAL CHI,FGUMBL
      COMMON/PSIPR/IPSI,C,H1,H2,H3,XK,D
C
      IF (N.LE.0.OR.MDX.LT.N.OR.NP.LE.0) CALL MESSGE(500,'SIGSGMB',1)
      SUMI=0.D0
      BET1=SIGBET(1)
      SIGMN=SIGBET(2)
      IF (NIT.NE.1) GOTO 300
      DO 160 I=1,N
      TMP=Y(I)
      GMX=0.0
      DO 140 J=1,NP
      TMP=TMP-BETA(J)*X(I,J)
      GMX=GMX+GAMMA(J)*X(I,J)
  140 CONTINUE
      RS(I)=TMP-BET1
      SW(I)=GMX
  160 CONTINUE
      SX(1,1)=0.0
      SZ(1)=0.0
c      DO 220 I=1,N
c      IO=RS(I)
c      DO 200 J=1,NP
c      SX(I,J)=X(IO,J)
c  200 CONTINUE
c      SZ(I)=DELTA(IO)
c  220 CONTINUE
C
  300 CONTINUE
      DO 370 I=1,N
      SUMJ=0.D0
      VI=RS(I)
      IF (METH.EQ.3) SIGMN=SIG/S0
      IF (METH.GE.4) SIGMN=SIG
      TMP=(VI-SW(I))/SIG-MU0
      IF (METH.GE.4) TMP=TMP/S0
      IF (DELTA(I).NE.0.0) THEN
        CI=CHI(TMP)
        SUMI=SUMI+DBLE(CI)
      ELSE 
        IF (TMP.GE.XK) THEN
          SUMI=SUMI+1.D0
          GOTO 370 
        ENDIF
        TT=DBLE(VI)/DBLE(SIGMN)
        AA=1.0001D0-FGUMBL(TT,LINT)
        WGT(1)=MU0
        WGT(2)=SW(I)
        WGT(3)=SIG
        IF (METH.GE.4) WGT(3)=S0*SIG
        WGT(4)=SIGMN
        WGT(5)=FLOAT(LINT)
        CALL RHOGMB(VI,WGT,SUMJ)
        SUMI=SUMI+SUMJ/AA
      ENDIF    
  370 CONTINUE
      SUMJ=SUMI/DFLOAT(N)
      EQB=SNGL(SUMJ)
      SIGMA=SIG*SQRT(EQB/B)
      RETURN
      END
C
C==========================================================================
C
      DOUBLE PRECISION FUNCTION ROGMBL(DX,WGT,N,EXU,EXV)
      DIMENSION WGT(N)
      DOUBLE PRECISION ANS,EXU,DX,AA,G,MU0,SIG,SSN
      EXTERNAL EXU,EXV
      COMMON/PSIPR/IPSI,C,H1,H2,H3,XK,D
C
C  Initializations
C
      MU0=DBLE(WGT(1))
      G=DBLE(WGT(2))
      SIG=DBLE(WGT(3))
      SSN=DBLE(WGT(4))
      ITYP=WGT(5)
      IF (SSN.LT.1.D-4) SSN=1.D-4
      AA=(SIG*(DX+MU0)+G)/SSN
      ANS=EXU(AA,ITYP)
      ROGMBL=0.D0
      IF (ANS.EQ.0.D0) RETURN
      V=SNGL(DX)
      TMP=EXV(V)
      ROGMBL=DBLE(TMP)*ANS
      RETURN
      END
C
C-----------------------------------------------------------------------
C
      SUBROUTINE RHOGMB(VI,WGT,SUM)
c 
c  INT  Rho[...]*Gumbel(x) dx  for  x=vi to Inf
c
      DOUBLE PRECISION LOW,HI,SUM,ROGMBL,DGUMBL,FGUMBL,TILD,ERRSTD,
     +                 WORK,SUMA,SUMB,SUMC,TMP1,TMP2,S,G,MU0,SIGMA,XXK
      REAL WGT(5)
      EXTERNAL ROGMBL,DGUMBL,FGUMBL,CHI
      COMMON/INTEGN/AINTEG(4),IWORK(80),WORK(320)
      COMMON/PSIPR/IPSI,C,H1,H2,H3,XK,D
      DATA KEY,LIMIT,TILD/1,80,0.001D0/
      MU0=DBLE(WGT(1))
      G=DBLE(WGT(2))
      S=DBLE(WGT(3))
      SIGMA=DBLE(WGT(4))
      ITYP=WGT(5)
      XXK=DBLE(XK)
      LOW=(DBLE(VI)-G)/S-MU0
      SUMA=0.D0
      SUMB=0.D0
      SUMC=0.D0
      IF (LOW.LT.-XXK) THEN
        TMP1=(S*(MU0-XXK)+G)/SIGMA
        TMP2=(S*(MU0+LOW)+G)/SIGMA
        SUMA=FGUMBL(TMP1,ITYP)-FGUMBL(TMP2,ITYP)
        LOW=-XXK
      ENDIF
      IF (LOW.LT.XXK) THEN
        HI=XXK
        CALL INTGRD(ROGMBL,WGT,5,DGUMBL,CHI,LOW,HI,TILD,TILD,KEY,LIMIT,
     +  SUMB,ERRSTD,NEVAL,IER,WORK,IWORK)
        IF (IER.NE.0) CALL MESSGE(400+IER,'RHOGMB',0)
        LOW=XXK
      ENDIF
      TMP2=(S*(MU0+LOW)+G)/SIGMA
      SUMC=1.D0-FGUMBL(TMP2,ITYP)
      SUM=SUMA+(S/SIGMA)*SUMB+SUMC
      RETURN
      END
C
C==========================================================================
C
      DOUBLE PRECISION FUNCTION RONORM(DX,WGT,N,EXU,EXV)
      DIMENSION WGT(N)
      DOUBLE PRECISION ANS,EXU,DX,AA,G,SIG,SSN
      EXTERNAL EXU,EXV
      COMMON/PSIPR/IPSI,C,H1,H2,H3,XK,D
C
C  Initializations
C
      G=DBLE(WGT(1))
      SIG=DBLE(WGT(2))
      SSN=DBLE(WGT(3))
      IF (SSN.LT.1.D-4) SSN=1.D-4
      AA=(SIG*DX+G)/SSN
      ANS=EXU(AA)
      RONORM=0.D0
      IF (ANS.EQ.0.D0) RETURN
      TMP=SNGL(DX)
      TMP=EXV(TMP)
      RONORM=DBLE(TMP)*SIG*ANS/SSN
      RETURN
      END
C
C-----------------------------------------------------------------------
C
      SUBROUTINE RHONRM(VI,WGT,SUM)
c 
c  INT  Rho[...]*dnorm(x) dx  for  x=vi to Inf
c
      DOUBLE PRECISION LOW,HI,SUM,RONORM,DGAUSS,TILD,ERRSTD,WORK,
     +                 SUMA,SUMB,SUMC,TMP1,TMP2,S,G,SIGMA,XXK
      REAL WGT(3)
      EXTERNAL RONORM,DGAUSS,FGAUSS,CHI
      COMMON/INTEGN/AINTEG(4),IWORK(80),WORK(320)
      COMMON/PSIPR/IPSI,C,H1,H2,H3,XK,D
      DATA KEY,LIMIT,TILD/1,80,0.0001D0/
      G=DBLE(WGT(1))
      S=DBLE(WGT(2))
      SIGMA=DBLE(WGT(3))
      XXK=DBLE(XK)
      LOW=(DBLE(VI)-G)/S
      SUMA=0.D0
      SUMB=0.D0
      SUMC=0.D0
      IF (LOW.LT.-XXK) THEN
        TMP1=(-S*XXK+G)/SIGMA
        TMP2=(S*LOW+G)/SIGMA
        TMP=FGAUSS(SNGL(TMP1))-FGAUSS(SNGL(TMP2))
        SUMA=DBLE(TMP)
        LOW=-XXK
      ENDIF
      IF (LOW.LT.XXK) THEN
        HI=XXK
        CALL INTGRD(RONORM,WGT,3,DGAUSS,CHI,LOW,HI,TILD,TILD,KEY,LIMIT,
     +  SUMB,ERRSTD,NEVAL,IER,WORK,IWORK)
        IF (IER.NE.0) CALL MESSGE(400+IER,'RHONRM',0)
        LOW=XXK
      ENDIF
      TMP2=(S*LOW+G)/SIGMA
      TMP=1.0-FGAUSS(SNGL(TMP2))
      SUMC=DBLE(TMP)
      SUM=SUMA+SUMB+SUMC
      RETURN
      END
C
      DOUBLE PRECISION FUNCTION SIGMBL(DX,WGT,N,EXU,EXV)
      DIMENSION WGT(N)
      DOUBLE PRECISION ANS,ARG,EXU,DX
      EXTERNAL EXU,EXV
C
C  Initializations
C
      ITYP=WGT(5)
      ANS=EXU(DX,ITYP)
      SIGMBL=0.D0
      IF (ANS.EQ.0.D0) RETURN
      S=WGT(1)
      T=WGT(2)
      SM=WGT(3)
      SI=WGT(4)
      ARG=(DX-DBLE(T))/DBLE(S)
      XTS=SNGL(ARG)
      EMX=EXV(SM*XTS)
      IF (SI.LE.2.0) TMP=SM*(EMX-1.0)
      IF (SI.EQ.2.0) TMP=XTS*TMP-1.0
      IF (SI.EQ.3.0) TMP=EMX
      IF (SI.EQ.4.0) TMP=XTS*EMX
      IF (SI.GE.5.0) TMP=SM*(EMX-1.0)+XTS*EMX
      IF (SI.EQ.6.0) TMP=XTS*TMP
      SIGMBL=DBLE(TMP)*ANS
      RETURN
      END
C
C-----------------------------------------------------------------------
C
      SUBROUTINE DRGFLI(F,L,Y,A,B,TOL,MAXIT,X,ITERM)
C.......................................................................
C   COPYRIGHT 1992 Alfio Marazzi
C.......................................................................
C
      DOUBLE PRECISION F,A,B,X,Y,FA,FB,FN,XN,TL,TOL
      EXTERNAL F
      LOGICAL NPRCHK
      DATA TL/1.D-10/
C
C  PARAMETER CHECK
C
      NPRCHK=A.LE.B.AND.TOL.GT.0..AND.MAXIT.GT.1
      IF (.NOT.NPRCHK) CALL MESSGE(500,'DRGFLI',1)
C
C  INITIALIZE
C
      ITR=1
      FA=F(A,L)-Y
      FB=F(B,L)-Y
C
C  REGULA FALSI ITERATION
C
   20 IF (DABS(FA-FB).GT.TL) GOTO 30
      CALL MESSGE(401,'DRGFLI',0)
      RETURN
   30 XN=(A*FB-B*FA)/(FB-FA)
      FN=F(XN,L)-Y
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
      DOUBLE PRECISION FUNCTION DGAUSI(X,IOPT)
      DOUBLE PRECISION X,SPI,X2,XEXPD
      EXTERNAL XEXPD
      DATA SPI/2.506628274631D0/
      X2=DFLOAT(IOPT)
      X2=-X*X/2.D0
      DGAUSI=XEXPD(X2)/SPI
      RETURN
      END
C
      DOUBLE PRECISION FUNCTION DGAUSS(X)
      DOUBLE PRECISION X,SPI,X2,XEXPD
      EXTERNAL XEXPD
      DATA SPI/2.506628274631D0/
C
      X2=-X*X/2.D0
      DGAUSS=XEXPD(X2)/SPI
      RETURN
      END
C
      FUNCTION FGAUSS(X)
      REAL               X,SQR1D2
      DATA               SQR1D2/.7071068/
      CALL CERF(-X*SQR1D2,C)
      FGAUSS = 0.5 * C
      RETURN
      END
C
      SUBROUTINE CERF(X,F)
C.......................................................................
C
C   COPYRIGHT 1992 Alfio Marazzi
C
C   AUTHOR : A. RANDRIAMIHARISOA
C                
C.......................................................................
C
      REAL               X,F
      INTEGER            ISW,I
      DIMENSION          P(3),Q(2),P1(5),Q1(4),P2(3),Q2(2)
      REAL               P,Q,P1,Q1,P2,Q2,XMIN,XLARGE,SQRPI,XX,
     *                   RES,XSQ,XNUM,XDEN,XI,XBIG
      EXTERNAL XEXP
      DATA               P(1)/.3166529/,P(2)/1.722276/,
     *                   P(3)/21.38533/
      DATA               Q(1)/7.843746/,Q(2)/18.95226/
      DATA               P1(1)/.5631696/,P1(2)/3.031799/,
     *                   P1(3)/6.865018/,P1(4)/7.373888/,
     *                   P1(5)/4.318779E-5/
      DATA               Q1(1)/5.354217/,Q1(2)/12.79553/,
     *                   Q1(3)/15.18491/,Q1(4)/7.373961/
      DATA               P2(1)/-5.168823E-2/,P2(2)/-.1960690/,
     *                   P2(3)/-4.257996E-2/
      DATA               Q2(1)/.9214524/,Q2(2)/.1509421/
      DATA               XMIN/1.0E-5/,XLARGE/4.1875E0/
      DATA               XBIG/9.0/
      DATA               SQRPI/.5641896/
C
      Y=X
      XX = Y
      ISW = 1
      IF (XX.GE.0.0E0) GO TO 5
      ISW = -1
      XX = -XX
    5 IF (XX.LT..477E0) GO TO 10
      IF (XX.LE.4.0E0) GO TO 25
      IF (ISW .GT. 0) GO TO 35
      IF (XX.LT.XLARGE) GO TO 40
      RES = 2.0E0
      GO TO 55
   10 IF (XX.LT.XMIN) GO TO 15
      XSQ = XX*XX
      XNUM = (P(1)*XSQ+P(2))*XSQ+P(3)
      XDEN = (XSQ+Q(1))*XSQ+Q(2)
      RES = XX*XNUM/XDEN
      GO TO 20
   15 RES = XX*P(3)/Q(2)
   20 IF (ISW.EQ.-1) RES = -RES
      RES = 1.0E0-RES
      GO TO 55
   25 XSQ = XX*XX
      XNUM = P1(5)*XX+P1(1)
      XDEN = XX+Q1(1)
      DO 30 I=2,4
         XNUM = XNUM*XX+P1(I)
         XDEN = XDEN*XX+Q1(I)
   30 CONTINUE
      RES = XNUM/XDEN
      GO TO 45
   35 IF (XX.GT.XBIG) GO TO 50
   40 XSQ = XX*XX
      XI = 1.0E0/XSQ
      XNUM = (P2(1)*XI+P2(2))*XI+P2(3)
      XDEN = (XI+Q2(1))*XI+Q2(2)
      RES = (SQRPI+XI*XNUM/XDEN)/XX
   45 RES = RES*XEXP(-XSQ)
      IF (ISW.EQ.-1) RES = 2.0E0-RES
      GO TO 55
   50 RES = 0.0E0
   55 F = RES
      RETURN
      END
C
      DOUBLE PRECISION FUNCTION DGUMBL(X,ITYP)
      DOUBLE PRECISION X,TMP,LOWER,UPPER,XEXPD
      EXTERNAL XEXPD
      CALL GMBLIM(0.D0,1.D0,ITYP,LOWER,UPPER)
      DGUMBL=0.D0
      IF (X.LE.LOWER) RETURN
      IF (X.GE.UPPER) RETURN
      IF (ITYP.EQ.1) THEN
       TMP=-X-XEXPD(-X)
      ELSE
       TMP=X-XEXPD(X)
      ENDIF
      DGUMBL=XEXPD(TMP)
      RETURN
      END
c
      DOUBLE PRECISION FUNCTION FGUMBL(X,ITYP)
      DOUBLE PRECISION X,TMP,LOWER,UPPER,XEXPD
      EXTERNAL XEXPD
      CALL GMBLIM(0.D0,1.D0,ITYP,LOWER,UPPER)
      FGUMBL=0.D0
      IF (X.LE.LOWER) RETURN
      FGUMBL=1.D0
      IF (X.GE.UPPER) RETURN
      IF (ITYP.EQ.1) THEN
       TMP=XEXPD(-X)
       TMP=XEXPD(-TMP)
      ELSE
       TMP=XEXPD(X)
       TMP=1.D0-XEXPD(-TMP)
      ENDIF
      FGUMBL=TMP
      RETURN
      END
C
      SUBROUTINE GMBLIM(TAU,V,ITYP,LOWER,UPPER)
      DOUBLE PRECISION TAU,V,LOWER,UPPER,ZUP,ZLOW
      DATA ZUP,ZLOW/3.5D0,-2.8D1/
      IF (ITYP.LT.1.OR.ITYP.GT.2) CALL MESSGE(500,'GMBLIM',1)
      LOWER=ZLOW*V + TAU
      UPPER=ZUP*V + TAU
      IF (ITYP.EQ.1) THEN
        LOWER=-ZUP*V + TAU
        UPPER=-ZLOW*V + TAU
      ENDIF
      RETURN
      END
C
C-----------------------------------------------------------------------
C
      SUBROUTINE INTGAMD(X,P,G)
C.......................................................................
C
C   AUTHOR :   G. P. BHATTACHARJEE (1970)
C              ALGORITHM AS 32 "THE (INCOMPLETE) GAMA INTEGRAL" MODIFIED
C              APPLIED STATISTICS, VOL.19.
C              REPRINT FROM PP.285-287 WITH THE PERMISSION OF 
C              BLACKWELL PUBLISHERS.
C.......................................................................
C
      DOUBLE PRECISION X,PN(6),P,G,GP,GIN,OFLO,FACTOR,RN,TERM,TOL,
     1                 A,B,AN,DIF,XEXPD
      EXTERNAL XEXPD
      DATA TOL/1.0D-8/
C
      G=0.D0
      IF (X.EQ.0.D0) RETURN
      IF (X.LT.0.D0.OR.P.LE.0.D0) CALL MESSGE(500,'INTGAMD',1)
      CALL MACHD(6,OFLO)
      OFLO=OFLO*1.D-15
      CALL LGAMAD(P,GP)
      GIN=0.D0
      FACTOR=XEXPD(P*DLOG(X)-X-GP)
      IF (X.GT.1.D0.AND.X.GE.P) GOTO 30
C
C  CALCULATION BY SERIES EXPANSION
C
      GIN=1.D0
      TERM=1.D0
      RN=P
   20 RN=RN+1.D0
      TERM=TERM*X/RN
      GIN=GIN+TERM
      IF (TERM.GT.TOL) GOTO 20
      GIN=GIN*FACTOR/P
      GOTO 50
C
C  CALCULATION BY CONTINUED FRACTION
C
   30 A=1.D0-P
      B=A+X+1.D0
      TERM=0.D0
      PN(1)=1.D0
      PN(2)=X
      PN(3)=X+1.D0
      PN(4)=X*B
      GIN=PN(3)/PN(4)
   32 A=A+1.D0
      B=B+2.D0
      TERM=TERM+1.D0
      AN=A*TERM
      DO 33 I=1,2
   33 PN(I+4)=B*PN(I+2)-AN*PN(I)
      IF (PN(6).EQ.0.D0) GOTO 35
      RN=PN(5)/PN(6)
      DIF=DABS(GIN-RN)
      IF (DIF.GT.TOL) GOTO 34
      IF (DIF.LE.TOL*RN) GOTO 42
   34 GIN=RN
   35 DO 36 I=1,4
   36 PN(I)=PN(I+2)
      IF (DABS(PN(5)).LT.OFLO) GOTO 32
      DO 41 I=1,4
   41 PN(I)=PN(I)/OFLO
      GOTO 32
   42 GIN=1.D0-FACTOR*GIN
   50 G=GIN*XEXPD(GP)
      RETURN
      END
C
C-----------------------------------------------------------------------
C
      SUBROUTINE LGAMAD(X,GL)
C.......................................................................
C
C   AUTHORS :     M.C. PIKE AND I.D. HILL (1966)
C                 ALGORITHM 291: LOGARITHM OF GAMMA FUNCTION.
C                 COMMUNICATIONS OF THE ACM, VOL.9, P 684.
C                 ADAPTED FOR ROBETH BY A. RANDRIAMIHARISOA
C.......................................................................
C
      DOUBLE PRECISION X,GL,V,F,Z
      IF (X.LE.0.D0) CALL MESSGE(500,'LGAMAD',1)
      V=X
      F=0.D0
      IF (X.GE.7.D0) GOTO 300
      F=1.D0
      Z=X-1.D0
  100 Z=Z+1.D0
      IF (Z.GE.7.D0) GOTO 200
      V=Z
      F=F*Z
      GOTO 100
  200 V=V+1.D0
      F=-DLOG(F)
  300 Z=1.D0/V**2
      GL=F+(V-0.5D0)*DLOG(V)-V+.9189385332D0+(((-.000595238D0*Z+
     +   .0007936507D0)*Z - .0027777778D0)*Z+.0833333333D0)/V
      RETURN
      END
C
C-----------------------------------------------------------------------
C
      SUBROUTINE WHSKNRM(VI,WGT,SUM)
c 
c  INT  WH(x)*sk((x-xbet)/dsig))*dnorm(x) dx  for  x=vi to Inf
c
      DOUBLE PRECISION LOW,HI,SUM,WSKNORM,DGAUSS,TILD,ERRSTD,WORK,
     +                 SUMA,SUMB,SUMC,XXK
      DIMENSION WGT(4)
      EXTERNAL WSKNORM,DGAUSS,W0HMP,W0TUK
      COMMON/INTEGN/AINTEG(4),IWORK(80),WORK(320)
      COMMON/PSIPR/IPSI,CPSI,H1,H2,H3,XK,DCHI
      DATA KEY,LIMIT,TILD/1,80,0.00001D0/
      TU=WGT(1)
      CC=WGT(2)
      IWGT=INT(WGT(4))
C     ETA2=SQRT(ETA)
C     IF (IWGT.EQ.2) ETA2=XK
      LOW=DBLE(VI)
      IF (VI.LT.-TU) LOW=-DBLE(TU)
      SUMA=0.D0
      SUMB=0.D0
      SUMC=0.D0
      IF (IWGT.EQ.2) GOTO 10
      IF (CC.LE.0.2) GOTO 10
      XXK=DBLE(SQRT(TU*TU-2.0*CC))
      IF (LOW.LT.-XXK) THEN
        HI=-XXK
        CALL INTGRD(WSKNORM,WGT,4,DGAUSS,W0HMP,LOW,HI,TILD,TILD,KEY,
     +       LIMIT,SUMA,ERRSTD,NEVAL,IER,WORK,IWORK)
        IF (IER.NE.0) CALL MESSGE(400+IER,'WHSKNRM',0)
        LOW=-XXK
      ENDIF
      IF (LOW.LT.XXK) THEN
        HI=XXK
        CALL INTGRD(WSKNORM,WGT,4,DGAUSS,W0HMP,LOW,HI,TILD,TILD,KEY,
     +       LIMIT,SUMB,ERRSTD,NEVAL,IER,WORK,IWORK)
        IF (IER.NE.0) CALL MESSGE(400+IER,'WHSKNRM',0)
        LOW=XXK
      ENDIF
  10  HI=DBLE(TU)
      CALL INTGRD(WSKNORM,WGT,4,DGAUSS,W0HMP,LOW,HI,TILD,TILD,KEY,
     +     LIMIT,SUMC,ERRSTD,NEVAL,IER,WORK,IWORK)
      IF (IER.NE.0) CALL MESSGE(400+IER,'WHSKNRM',0)
      SUM=SUMA+SUMB+SUMC
      RETURN
      END
C
C==========================================================================
C
      DOUBLE PRECISION FUNCTION WSKNORM(DX,WGT,N,EXU,EXV)
      DIMENSION WGT(N)
      DOUBLE PRECISION ANS,EXU,DX,XBT,DS,SK
      EXTERNAL EXU,EXV,RHO,W0TUK
      COMMON/PSIPR/IPSI,CPSI,H1,H2,H3,XK,DCHI
      COMMON/WGTML/XBT,DS,T2
C
C  Initializations
C
      ANS=EXU(DX)
      WSKNORM=0.D0
      IF (ANS.EQ.0.D0) RETURN
      TU=WGT(1)
      CC=WGT(2)
      POW=WGT(3)
      IWGT=INT(WGT(4))
      IF (DS.LT.1.D-6) DS=1.D-6
      TMP=SNGL(DX)
      IF (IWGT.EQ.1) TMP=EXV(TMP,TU,CC)
      IF (IWGT.EQ.2) TMP=1.0-RHO(TMP)
      IF (IWGT.EQ.3) TMP=W0TUK(TMP,TU,CC)
      SK=1.D0
      IF (POW.GE.1.0) THEN 
        SK=(DX-XBT)/DS 
        IF (POW.EQ.2.0) SK=SK*SK
      ENDIF
      WSKNORM=DBLE(TMP)*SK*ANS
      RETURN
      END
C
C-----------------------------------------------------------------------
C
c      SUBROUTINE W3SKNRM(VI,WGT,SUM)
c c 
c c  INT  WH(dsig*x+xbet)*sk(x)*dnorm(x) dx  for  x=(vi-xbet)/dsig to Inf
c c
c       DOUBLE PRECISION LOW,HI,SUM,W3SKNORM,DGAUSS,TILD,ERRSTD,WORK,
c      +                 XBT,DS,SUMA,SUMB,SUMC,XXK,LOW2,HI2
c       DIMENSION WGT(4)
c       EXTERNAL W3SKNORM,DGAUSS,W0HMP,W0TUK
c       COMMON/INTEGN/AINTEG(4),IWORK(80),WORK(320)
c       COMMON/PSIPR/IPSI,CPSI,H1,H2,H3,XK,DCHI
c       COMMON/WGTML/XBT,DS,T2
c       DATA KEY,LIMIT,TILD/1,80,0.00001D0/
c       TU=WGT(1)
c       CC=WGT(2)
c       IWGT=INT(WGT(4))
c C     ETA2=SQRT(ETA)
c C     IF (IWGT.EQ.2) ETA2=XK
c       LOW=DBLE(VI)
c       IF (VI.LT.-TU) LOW=-DBLE(TU)
c       SUMA=0.D0
c       SUMB=0.D0
c       SUMC=0.D0
c       IF (IWGT.EQ.2) GOTO 10
c       IF (CC.LE.0.1) GOTO 10
c       XXK=DBLE(SQRT(TU*TU-2.0*CC))
c       IF (LOW.LT.-XXK) THEN
c         HI=-XXK
c         LOW2=(LOW-XBT)/DS
c         HI2=(HI-XBT)/DS
c         CALL INTGRD(W3SKNORM,WGT,4,DGAUSS,W0HMP,LOW2,HI2,TILD,TILD,KEY,
c      +       LIMIT,SUMA,ERRSTD,NEVAL,IER,WORK,IWORK)
c         IF (IER.NE.0) CALL MESSGE(400+IER,'WHSKNRM',0)
c         LOW=-XXK
c       ENDIF
c       IF (LOW.LT.XXK) THEN
c         HI=XXK
c         LOW2=(LOW-XBT)/DS
c         HI2=(HI-XBT)/DS
c         CALL INTGRD(W3SKNORM,WGT,4,DGAUSS,W0HMP,LOW2,HI2,TILD,TILD,KEY,
c      +       LIMIT,SUMB,ERRSTD,NEVAL,IER,WORK,IWORK)
c         IF (IER.NE.0) CALL MESSGE(400+IER,'WHSKNRM',0)
c         LOW=XXK
c       ENDIF
c   10  HI=DBLE(TU)
c       LOW2=(LOW-XBT)/DS
c       HI2=(HI-XBT)/DS
c       CALL INTGRD(W3SKNORM,WGT,4,DGAUSS,W0HMP,LOW2,HI2,TILD,TILD,KEY,
c      +     LIMIT,SUMC,ERRSTD,NEVAL,IER,WORK,IWORK)
c       IF (IER.NE.0) CALL MESSGE(400+IER,'WHSKNRM',0)
c       SUM=SUMA+SUMB+SUMC
c       RETURN
c       END
c C
c C==========================================================================
c C
c       DOUBLE PRECISION FUNCTION W3SKNORM(DX,WGT,N,EXU,EXV)
c       DIMENSION WGT(N)
c       DOUBLE PRECISION ANS,EXU,DX,XBT,DS,SK
c       EXTERNAL EXU,EXV,RHO,W0TUK
c       COMMON/PSIPR/IPSI,CPSI,H1,H2,H3,XK,DCHI
c       COMMON/WGTML/XBT,DS,T2
c C
c C  Initializations
c C
c       ANS=EXU(DX)
c       W3SKNORM=0.D0
c       IF (ANS.EQ.0.D0) RETURN
c       TU=WGT(1)
c       CC=WGT(2)
c       POW=WGT(3)
c       IWGT=INT(WGT(4))
c       IF (DS.LT.1.D-6) DS=1.D-6
c       TMP=SNGL(DS*DX+XBT)
c       IF (IWGT.EQ.1) TMP=EXV(TMP,TU,CC)
c       IF (IWGT.EQ.2) TMP=1.0-RHO(TMP)
c       IF (IWGT.EQ.3) TMP=W0TUK(TMP,TU,CC)
c       SK=1.D0
c       IF (POW.GE.1.0) THEN 
c C       SK=(DX-XBT)/DS 
c         SK=DX
c         IF (POW.EQ.2.0) SK=SK*SK
c       ENDIF
c       W3SKNORM=DBLE(TMP)*SK*ANS
c       RETURN
c       END
c
      FUNCTION W0HMP(X,TU,CC)
C     tu = sqrt(2*(eta-log(sqrt(2*pi)))) ; tl = -tu
      TMP=0.5*(TU*TU-X*X)/CC
      IF (TMP.GE.1.0) TMP=1.0
      IF (TMP.LE.0.0) TMP=0.0
      W0HMP=TMP
      RETURN
      END
C
      SUBROUTINE WHAMP(N,TU,TL,CC,U,WU)
C     tu = sqrt(2*(eta-log(sqrt(2*pi)))) ; tl = -tu
      DIMENSION U(N),WU(N)
      DO 100 I=1,N
      TMP=-0.5*(U(I)*U(I)+TU*TL)/CC
      IF (TMP.GE.1.0) TMP=1.0
      IF (TMP.LE.0.0) TMP=0.0
      WU(I)=TMP
  100 CONTINUE
      RETURN
      END 
C
      FUNCTION W0TUK(X,TU,CC) 
c  tu = sqrt(2*(eta-log(sqrt(2*pi)))) ; tl = -tu
      TMP=0.5*(X*X-TU*TU)
      CHC=1.0
      IF (ABS(TMP).GE.CC) GOTO 100
      S2=(TMP/CC)**2
      CHC=(S2*(S2-3.0)+3.0)*S2
  100 IF (ABS(X).GT.TU) CHC=0.0
      W0TUK=CHC
      RETURN
      END
C
      SUBROUTINE W1TUK(N,X,TU,CC,WX) 
c  tu = sqrt(2*(eta-log(sqrt(2*pi)))) ; tl = -tu
      DIMENSION X(N),WX(N)
      DO 200 I=1,N
      TMP=0.5*(X(I)*X(I)-TU*TU)
      CHC=1.0
      IF (ABS(TMP).GE.CC) GOTO 100
      S2=(TMP/CC)**2
      CHC=(S2*(S2-3.0)+3.0)*S2
  100 IF (ABS(X(I)).GT.TU) CHC=0.0
      WX(I)=CHC
  200 CONTINUE
      RETURN
      END
C
      FUNCTION W0GMB(X,T2,CC,ITYP)
      EXTERNAL XEXP
      ONE=1.0
      IF (ITYP.EQ.1) ONE=-1.0
      TMP=(-XEXP(ONE*X)+ONE*X + T2)/CC
      IF (TMP.GE.1.0) TMP=1.0
      IF (TMP.LE.0.0) TMP=0.0
      W0GMB=TMP
      RETURN
      END 
C
      SUBROUTINE WGMBL(N,T2,CC,ITYP,U,WU)
      DIMENSION U(N),WU(N)
      EXTERNAL XEXP
      ONE=1.0
      IF (ITYP.EQ.1) ONE=-1.0
      DO 100 I=1,N
      TMP=(-XEXP(ONE*U(I))+ONE*U(I) + T2)/CC
      IF (TMP.GE.1.0) TMP=1.0
      IF (TMP.LE.0.0) TMP=0.0
      WU(I)=TMP
  100 CONTINUE
      RETURN
      END 
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      SUBROUTINE WHSKGMB(VI,WGT,SUM)
c 
c  INT  WH(x)*sk((x-xbet)/dsig))*dgumbel(x) dx  for  x=vi to Inf
c
      DOUBLE PRECISION LOW,HI,SUM,WSKGMBL,DGUMBL,TILD,ERRSTD,WORK,
     +                 SUMA,SUMB,SUMC,XXK
      DIMENSION WGT(8)
      EXTERNAL WSKGMBL,DGUMBL,W0GMB
      COMMON/INTEGN/AINTEG(4),IWORK(80),WORK(320)
      COMMON/PSIPR/IPSI,CPSI,H1,H2,H3,XK,DCHI
      DATA KEY,LIMIT,TILD/1,80,0.00001D0/
      TU=WGT(1)
      CC=WGT(2)
      IWGT=INT(WGT(4))
      LINT=INT(WGT(6))
      TL=WGT(7)
      T2=WGT(8)
      LOW=DBLE(VI)
      IF (VI.LT.TL) LOW=DBLE(TL)
      SUMA=0.D0
      SUMB=0.D0
      SUMC=0.D0
      IF (IWGT.EQ.2) GOTO 10
      IF (CC.LE.0.2) GOTO 10
      CALL SOLVT2(T2-CC,ITYP,TILD,MAXIT,TCU,TCL)
      XXK=DBLE(TCL)
      IF (LOW.LT.XXK) THEN
        HI=XXK
        CALL INTGRD(WSKGMBL,WGT,8,DGUMBL,W0GMB,LOW,HI,TILD,TILD,KEY,
     +       LIMIT,SUMA,ERRSTD,NEVAL,IER,WORK,IWORK)
        IF (IER.NE.0) CALL MESSGE(400+IER,'WHSKGMB',0)
        LOW=XXK
      ENDIF
      XXK=DBLE(TCU)     
      IF (LOW.LT.XXK) THEN
        HI=XXK
        CALL INTGRD(WSKGMBL,WGT,8,DGUMBL,W0GMB,LOW,HI,TILD,TILD,KEY,
     +       LIMIT,SUMB,ERRSTD,NEVAL,IER,WORK,IWORK)
        IF (IER.NE.0) CALL MESSGE(400+IER,'WHSKGMB',0)
        LOW=XXK
      ENDIF
  10  HI=DBLE(TU)
      CALL INTGRD(WSKGMBL,WGT,8,DGUMBL,W0GMB,LOW,HI,TILD,TILD,KEY,
     +     LIMIT,SUMC,ERRSTD,NEVAL,IER,WORK,IWORK)
      IF (IER.NE.0) CALL MESSGE(400+IER,'WHSKGMB',0)
      SUM=SUMA+SUMB+SUMC
      RETURN
      END
C
C==========================================================================
C
      DOUBLE PRECISION FUNCTION WSKGMBL(DX,WGT,N,EXU,EXV)
      DIMENSION WGT(N)
      DOUBLE PRECISION ANS,EXU,DX,XBT,XTS,SS,DS,SK,XEXPD
      EXTERNAL EXU,EXV,RHO,XEXPD
      COMMON/PSIPR/IPSI,CPSI,H1,H2,H3,XK,DCHI
      COMMON/WGTML/XBT,DS,T2
C
C  Initializations
C
      LINT=INT(WGT(6))
      ANS=EXU(DX,LINT)
      WSKGMBL=0.D0
      IF (ANS.EQ.0.D0) RETURN
      TU=WGT(1)
      CC=WGT(2)
      POW=WGT(3)
      IWGT=INT(WGT(4))
      SS=DBLE(WGT(5))
      TL=WGT(7)
      T2=WGT(8)
      IF (DS.LT.1.D-6) DS=1.D-6
      TMP=SNGL(DX)
      IF (IWGT.EQ.1) TMP=EXV(TMP,T2,CC,LINT)
      IF (IWGT.EQ.2) TMP=1.0-RHO(TMP-0.5*(TU+TL))
      SK=1.D0
      IF (POW.GE.1.0) THEN 
        XTS=(DX-XBT)/DS
        SK=SS*(XEXPD(SS*XTS)-1.D0) 
        IF (POW.EQ.2.0) SK=XTS*SK
      ENDIF
      WSKGMBL=DBLE(TMP)*SK*ANS
      RETURN
      END
c C
c C********************************************************************************
c C
      DOUBLE PRECISION FUNCTION SRRHOG(Z,CONST,ITYP)
      implicit double precision(a-h,o-z)
      EXTERNAL XEXPD
      X=Z*DFLOAT(2*ITYP-3)  
      SRRHOG=XEXPD(X)-CONST-X
      RETURN
      END
C
C-----------------------------------------------------------------------
C
      SUBROUTINE SRRGFL3(F,CONST,ITYP,Y,A,B,TOL,MAXIT,X,ITERM)
      implicit double precision(a-h,o-z)
      EXTERNAL F
C
C  INITIALIZE
C
      TL=DMIN1(1.D-10,0.1D0*TOL)
      ITR=1
      MESS=0
   10 FA=F(A,CONST,ITYP)-Y
      FB=F(B,CONST,ITYP)-Y
C
C  REGULA FALSI ITERATION
C
   20 IF (DABS(FA-FB).GT.TL) GOTO 30
      MESS=MESS+1
      IF (MESS.LE.2) THEN
        A=A/10.D0
        GOTO 10
      ENDIF
C     CALL MESSGE(401,'RGFL3 ',0)
      RETURN
   30 XN=(A*FB-B*FA)/(FB-FA)
      FN=F(XN,CONST,ITYP)-Y
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
      SUBROUTINE SRF0G(U,TOL,MAXIT,P)
      implicit double precision(a-h,o-z)
      DOUBLE PRECISION LOW
c     Same result if ITYP=1 or ITYP=2
      EXTERNAL XEXPD,SRRHOG,FGUMBL
      P=0.D0
      IF (U.LE.1.D0) RETURN
      P=1.D0
      IF (U.GT.16.D0) RETURN
      CONST=U
      IF (U.GT.1.5D0) THEN
        LOW=-U
        UP=-U+1.5D0
C       IF (ITYP.EQ.1) THEN
C         LOW=U-1.5D0
C         UP=U
C       ENDIF
        CALL SRRGFL3(SRRHOG,CONST,2,0.D0,LOW,UP,TOL,MAXIT,TL,ITRM)
      ELSE
        TLO=TOL
        IF (U-1.D0.LT.1.D-3) TLO=DMIN1(TOL,1.D-8)
        LOW=-U
        UP=0.D0
C       IF (ITYP.EQ.1) THEN
C         LOW=0.D0
C         UP=U
C       ENDIF
        CALL SRRGFL3(SRRHOG,CONST,2,0.D0,LOW,UP,TOL,MAXIT,TL,ITRM)
      ENDIF
      LOW=DLOG(U)
      UP=U
C     IF (ITYP.EQ.1) THEN
C       UP=-LOW
C       LOW=-U
C     ENDIF
      CALL SRRGFL3(SRRHOG,CONST,2,0.D0,LOW,UP,TOL,MAXIT,TU,ITRM)
C     IF (ITYP.EQ.1) THEN
C       T=TU
C       TU=TL
C       TL=T
C     ENDIF
      P=FGUMBL(TU,2)-FGUMBL(TL,2)
      RETURN
      END
C
      SUBROUTINE SOLVT2(T2,ITYP,TOL,MAXIT,TU,TL)
      DOUBLE PRECISION LOW,LOGT2,TT2,UP2,UP,TOL,TTU,TTL,SRRHOG
      EXTERNAL SRRHOG
      TU=0.0
      TL=0.0
      IF (T2.LE.1.0) RETURN
      TT2=DBLE(T2)
      LOGT2=DLOG(TT2)
      UP2=DLOG(TT2+1.2D0*LOGT2)
      IF (T2.GT.1.D0.AND.T2.LE.1.5) THEN
        LOW=-TT2
        UP=0.D0
        CALL SRRGFL3(SRRHOG,TT2,2,0.D0,LOW,UP,TOL,MAXIT,TTL,ITRM)
        TL=SNGL(TTL)
      ENDIF
      IF (T2.GT.1.5D0.AND.T2.LE.16.0) THEN
        LOW=-TT2
        UP=-TT2+1.5D0
        CALL SRRGFL3(SRRHOG,TT2,2,0.D0,LOW,UP,TOL,MAXIT,TTL,ITRM)
        TL=SNGL(TTL)
      ENDIF
      IF (T2.GT.16.0) TL=-T2
      IF (T2.GT.1.D0.AND.T2.LE.50.0) THEN
        LOW=LOGT2
        UP=TT2
        CALL SRRGFL3(SRRHOG,TT2,2,0.D0,LOW,UP,TOL,MAXIT,TTU,ITRM)
        TU=SNGL(TTU)
      ENDIF
      IF (T2.GT.50.D0) THEN
        LOW=LOGT2
        UP=UP2
        CALL SRRGFL3(SRRHOG,TT2,2,0.D0,LOW,UP,TOL,MAXIT,TTU,ITRM)
        TU=SNGL(TTU)
      ENDIF
      IF (ITYP.EQ.1) THEN
        TMP=TU
        TU=-TL
        TL=-TMP
      ENDIF
      RETURN
      END
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      SUBROUTINE REFIRES(DPAR,X,Y,DELTA,N,NP,MDX,LINT,RES)  
C
C.......................................................................
C   PROGRAMMER : A. RANDRIAMIHARISOA
C.......................................................................
C
      DIMENSION DELTA(N)
      DOUBLE PRECISION X(MDX,NP),Y(N),DPAR(NP+1),RES(NP+1),DS
      COMMON/PSIPR/IPSI,CPSI,H1,H2,H3,XK,DCHI
C
      DS= DPAR(NP+1) 
      IF (DS.LT.1.D-6) DS=1.D-6
      IF (LINT.EQ.0.OR.LINT.EQ.3) THEN 
        CALL NRMRRES(DPAR,DS,X,DELTA,Y,N,NP,MDX,RES)
      ELSE
        CALL GMBRRES(DPAR,DS,X,DELTA,Y,LINT,N,NP,MDX,RES)
      ENDIF
      RETURN
      END
C
      SUBROUTINE REF1RES(DPAR,DS,X,Y,DELTA,N,NP,MDX,LINT,RES)  
C
C.......................................................................
C   PROGRAMMER : A. RANDRIAMIHARISOA
C.......................................................................
C
      DIMENSION DELTA(N)
      DOUBLE PRECISION X(MDX,NP),Y(N),DPAR(NP),RES(NP),DS
      COMMON/PSIPR/IPSI,CPSI,H1,H2,H3,XK,DCHI
C
      IF (DS.LT.1.D-6) DS=1.D-6
      IF (LINT.EQ.0.OR.LINT.EQ.3) THEN 
        CALL NRM1RES(DPAR,DS,X,DELTA,Y,N,NP,MDX,RES)
      ELSE
        CALL GMB1RES(DPAR,DS,X,DELTA,Y,LINT,N,NP,MDX,RES)
      ENDIF
      RETURN
      END
C
      SUBROUTINE REF2RES(DS,DPAR,X,Y,DELTA,N,NP,MDX,LINT,RES)  
C
C.......................................................................
C   PROGRAMMER : A. RANDRIAMIHARISOA
C.......................................................................
C
      DIMENSION DELTA(N)
      DOUBLE PRECISION X(MDX,NP),Y(N),DPAR(NP),RES,DS
      COMMON/PSIPR/IPSI,CPSI,H1,H2,H3,XK,DCHI
C
      IF (DS.LT.1.D-6) DS=1.D-6
      IF (LINT.EQ.0.OR.LINT.EQ.3) THEN 
        CALL NRM2RES(DS,DPAR,X,DELTA,Y,N,NP,MDX,RES)
      ELSE
        CALL GMB2RES(DS,DPAR,X,DELTA,Y,LINT,N,NP,MDX,RES)
      ENDIF
      RETURN
      END
C
C------------------------------------------------------------------------
C 
      SUBROUTINE SUMRRES(DBETA,X,Y,DELTA,N,NP,MDX,
     *           BETA,THETA,RS,DD,YY,SZ,SW,IT,RES)        
C
C.......................................................................
C   PROGRAMMER : A. RANDRIAMIHARISOA
C.......................................................................
C
      DIMENSION X(MDX,NP),Y(N),DELTA(N),BETA(1,NP),THETA(N),
     + RS(N),DD(N),YY(N),SZ(N),SW(N),DUMMY(1,2)
      DOUBLE PRECISION DBETA(NP+1),DS,SUMI,SUMP,SUMR,RES(NP+1)
      INTEGER IT(N)
      EXTERNAL PSY,RHO
      COMMON/PSIPR/IPSI,C,H1,H2,H3,XK,D
      DATA NU/0/
C
      IF (N.LE.0.OR.MDX.LT.N.OR.NP.LE.0) CALL MESSGE(500,'SUMRRES',1)
      NB=1
      B=0.5
      TU=1.E9
      IF (IPSI.EQ.2) TU=H3
      IF (IPSI.EQ.3) TU=1.0
      IF (IPSI.EQ.4) TU=XK
      TL=-TU
      LINT=0
      DS=DBETA(NP+1)
        DO 10 J=1,NP
        RES(J)=0.D0
        BETA(1,J)=SNGL(DBETA(J))
   10   CONTINUE
        CALL KMEDMAD(X,Y,IT,DELTA,BETA,N,NP,MDX,NB,1,LINT,DUMMY,THETA,
     *               RS,DD,YY,SZ,SW)
        NU=IT(1)
        RES(NP+1)=0.D0
        SIG=SNGL(DS)
      DO 160 I=1,N
      TMP=Y(I)
      DO 140 J=1,NP
      TMP=TMP-BETA(1,J)*X(I,J)
  140 CONTINUE
      SW(I)=TMP 
  160 CONTINUE
C
      DO 370 I=1,N
      TMP=SW(I)/SIG
      IF (DELTA(I).NE.0.0) THEN
        IF (TMP.GE.TU) GOTO 350
        IF (TMP.LE.TL) GOTO 350
        WRI=PSY(TMP)
        SUMI=DBLE(WRI)
        GOTO 330
      ELSE
        IF (TMP.GE.TU) GOTO 350        
        CALL NTRP0L(SW(I),NU,SZ,II)      
        IF (II.GE.NU-1) THEN
          TMP=SZ(NU)/SIG
          WRI=PSY(TMP)
          SUMI=DBLE(WRI)
          GOTO 330
        ENDIF
        AI=YY(II+1)  ! No more Tail censured obs. (delta(n)=-1)
        IF (ABS(AI).LT.0.00001) AI=1./FLOAT(N)
        SUMP=0.D0
        SUMR=0.D0
        DO 200 J=II+1,NU
        TMP=SZ(J)/SIG 
        IF (TMP.GE.TU.AND.SUMP.EQ.0.D0) GOTO 350
        WRI=PSY(TMP)
        PIJ=DD(J)
        SUMP=SUMP+DBLE(WRI*PIJ)
        WRI=RHO(TMP)
        SUMR=SUMR+DBLE(WRI*PIJ)
  200   CONTINUE
        SUMI=SUMP/DBLE(AI)
        DO 220 J=1,NP
  220   RES(J)=RES(J)+SUMI*DBLE(X(I,J))
        RES(NP+1)=RES(NP+1)+SUMR/DBLE(AI)
        GOTO 370
      ENDIF 
  330 DO 340 J=1,NP
  340 RES(J)=RES(J)+SUMI*DBLE(X(I,J))
  350 WRI=RHO(TMP)
      RES(NP+1)=RES(NP+1)+DBLE(WRI)
  370 CONTINUE
      DO 380 J=1,NP+1
  380 RES(J)=RES(J)/DFLOAT(N)
      RES(NP+1)=RES(NP+1)-0.5D0
      RETURN
      END
C
      SUBROUTINE SUM1RES(DBETA,DS,X,Y,DELTA,N,NP,MDX,
     *           BETA,THETA,RS,DD,YY,SZ,SW,IT,RES)        
C
C.......................................................................
C   PROGRAMMER : A. RANDRIAMIHARISOA
C.......................................................................
C
      DIMENSION X(MDX,NP),Y(N),DELTA(N),BETA(1,NP),THETA(N),
     + RS(N),DD(N),YY(N),SZ(N),SW(N),DUMMY(1,2)
      DOUBLE PRECISION DBETA(NP),DS,SUMI,SUMP,RES(NP)
      INTEGER IT(N)
      EXTERNAL PSY
      COMMON/PSIPR/IPSI,C,H1,H2,H3,XK,D
      DATA NU/0/
C
      IF (N.LE.0.OR.MDX.LT.N.OR.NP.LE.0) CALL MESSGE(500,'SUM1RES',1)
      NB=1
      B=0.5
      TU=1.E9
      IF (IPSI.EQ.2) TU=H3
      IF (IPSI.EQ.3) TU=1.0
      IF (IPSI.EQ.4) TU=XK
      TL=-TU
      LINT=0
c     DS=DBETA(NP+1)
        DO 10 J=1,NP
        RES(J)=0.D0
        BETA(1,J)=SNGL(DBETA(J))
   10   CONTINUE
        CALL KMEDMAD(X,Y,IT,DELTA,BETA,N,NP,MDX,NB,1,LINT,DUMMY,THETA,
     *               RS,DD,YY,SZ,SW)
        NU=IT(1)
        SIG=SNGL(DS)
      DO 160 I=1,N
      TMP=Y(I)
      DO 140 J=1,NP
      TMP=TMP-BETA(1,J)*X(I,J)
  140 CONTINUE
      SW(I)=TMP 
  160 CONTINUE
C
      DO 370 I=1,N
      TMP=SW(I)/SIG
      IF (DELTA(I).NE.0.0) THEN
        IF (TMP.GE.TU) GOTO 370
        IF (TMP.LE.TL) GOTO 370
        WRI=PSY(TMP)
        SUMI=DBLE(WRI)
        GOTO 330
      ELSE
        IF (TMP.GE.TU) GOTO 370        
        CALL NTRP0L(SW(I),NU,SZ,II)      
        IF (II.GE.NU-1) THEN
          TMP=SZ(NU)/SIG
          WRI=PSY(TMP)
          SUMI=DBLE(WRI)
          GOTO 330
        ENDIF
        AI=YY(II+1)  ! No more Tail censured obs. (delta(n)=-1)
        IF (ABS(AI).LT.0.00001) AI=1./FLOAT(N)
        SUMP=0.D0
c       SUMR=0.D0
        DO 200 J=II+1,NU
        TMP=SZ(J)/SIG 
        IF (TMP.GE.TU.AND.SUMP.EQ.0.D0) GOTO 370 
        WRI=PSY(TMP)
        PIJ=DD(J)
        SUMP=SUMP+DBLE(WRI*PIJ)
C       WRI=RHO(TMP)
C       SUMR=SUMR+DBLE(WRI*PIJ)
  200   CONTINUE
        SUMI=SUMP/DBLE(AI)
        DO 220 J=1,NP
  220   RES(J)=RES(J)+SUMI*DBLE(X(I,J))
C       RES(NP+1)=RES(NP+1)+SUMR/DBLE(AI)
        GOTO 370
      ENDIF 
  330 DO 340 J=1,NP
  340 RES(J)=RES(J)+SUMI*DBLE(X(I,J))
c 350 WRI=RHO(TMP)
c     RES(NP+1)=RES(NP+1)+DBLE(WRI)
  370 CONTINUE
      DO 380 J=1,NP
  380 RES(J)=RES(J)/DFLOAT(N)
c     RES(NP+1)=RES(NP+1)-0.5D0
      RETURN
      END
C
      SUBROUTINE SUM2RES(DS,DBETA,X,Y,DELTA,N,NP,MDX,
     *           BETA,THETA,RS,DD,YY,SZ,SW,IT,RES)        
C
C.......................................................................
C   PROGRAMMER : A. RANDRIAMIHARISOA
C.......................................................................
C
      DIMENSION X(MDX,NP),Y(N),DELTA(N),BETA(1,NP),THETA(N),
     + RS(N),DD(N),YY(N),SZ(N),SW(N),DUMMY(1,2)
      DOUBLE PRECISION DBETA(NP),DS,SUMR,RES
      INTEGER IT(N)
      EXTERNAL RHO
      COMMON/PSIPR/IPSI,C,H1,H2,H3,XK,D
      DATA NU/0/
C
      IF (N.LE.0.OR.MDX.LT.N.OR.NP.LE.0) CALL MESSGE(500,'SUMRRES',1)
      NB=1
      B=0.5
      TU=1.E9
      IF (IPSI.EQ.2) TU=H3
      IF (IPSI.EQ.3) TU=1.0
      IF (IPSI.EQ.4) TU=XK
      TL=-TU
      LINT=0
c     DS=DBETA(NP+1)
        RES=0.D0
        DO 10 J=1,NP
        BETA(1,J)=SNGL(DBETA(J))
   10   CONTINUE
        CALL KMEDMAD(X,Y,IT,DELTA,BETA,N,NP,MDX,NB,1,LINT,DUMMY,THETA,
     *               RS,DD,YY,SZ,SW)
        NU=IT(1)
        SIG=SNGL(DS)
      DO 160 I=1,N
      TMP=Y(I)
      DO 140 J=1,NP
      TMP=TMP-BETA(1,J)*X(I,J)
  140 CONTINUE
      SW(I)=TMP 
  160 CONTINUE
C
      DO 370 I=1,N
      TMP=SW(I)/SIG
      IF (DELTA(I).NE.0.0) THEN
        WRI=RHO(TMP)
        RES=RES+DBLE(WRI)
        GOTO 370
      ELSE
        IF (TMP.GE.TU) GOTO 350        
        CALL NTRP0L(SW(I),NU,SZ,II)      
        IF (II.GE.NU-1) THEN
          TMP=SZ(NU)/SIG
          GOTO 350
        ENDIF
        AI=YY(II+1)  ! No more Tail censured obs. (delta(n)=-1)
        IF (ABS(AI).LT.0.00001) AI=1./FLOAT(N)
c       SUMP=0.D0
        SUMR=0.D0
        DO 200 J=II+1,NU
        TMP=SZ(J)/SIG 
        IF (TMP.GE.TU.AND.SUMR.EQ.0.D0) GOTO 350
c       WRI=PSY(TMP)
        PIJ=DD(J)
c       SUMP=SUMP+DBLE(WRI*PIJ)
        WRI=RHO(TMP)
        SUMR=SUMR+DBLE(WRI*PIJ)
  200   CONTINUE
c       SUMI=SUMP/DBLE(AI)
c       DO 220 J=1,NP
c 220   RES(J)=RES(J)+SUMI*DBLE(X(I,J))
        RES=RES+SUMR/DBLE(AI)
        GOTO 370
      ENDIF 
c 330 DO 340 J=1,NP
c 340 RES(J)=RES(J)+SUMI*DBLE(X(I,J))
  350 WRI=RHO(TMP)
      RES=RES+DBLE(WRI)
  370 CONTINUE
c     DO 380 J=1,NP+1
c 380 RES(J)=RES(J)/DFLOAT(N)
      RES=RES/DFLOAT(N-NP)-0.5D0
      RETURN
      END
C
      SUBROUTINE SUMRJAC(DBETA,X,Y,DELTA,N,NP,MDX,
     *           BETA,THETA,RS,DD,YY,SZ,SW,IT,RES)        
C
C.......................................................................
C   PROGRAMMER : A. RANDRIAMIHARISOA
C.......................................................................
C
      DIMENSION X(MDX,NP),Y(N),DELTA(N),BETA(1,NP),THETA(N),
     + RS(N),DD(N),YY(N),SZ(N),SW(N),DUMMY(1,2)
      DOUBLE PRECISION DBETA(NP+1),DS,SUM1,SUM2,SUMP,SUMR,
     +       RES(NP+1,NP+1)
      INTEGER IT(N)
      EXTERNAL PSP,PSY,RHO
      COMMON/PSIPR/IPSI,CPSI,H1,H2,H3,XK,DCHI
      DATA NU/0/
C
      IF (N.LE.0.OR.MDX.LT.N.OR.NP.LE.0) CALL MESSGE(500,'SUMRJAC',1)
      NB=1
      B=0.5
      TU=1.E9
      IF (IPSI.EQ.2) TU=H3
      IF (IPSI.EQ.3) TU=1.0
      IF (IPSI.EQ.4) TU=XK
      TL=-TU
      LINT=0
      DS=DBETA(NP+1)
      NP1=NP+1
      DO 10 I=1,NP1
      DO 10 J=1,NP1
      RES(I,J)=0.D0
   10 CONTINUE
        DO 20 J=1,NP
        BETA(1,J)=SNGL(DBETA(J))
   20   CONTINUE
        CALL KMEDMAD(X,Y,IT,DELTA,BETA,N,NP,MDX,NB,1,LINT,DUMMY,THETA,
     *               RS,DD,YY,SZ,SW)
        NU=IT(1)
        SIG=SNGL(DS)
      DO 150 I=1,NU
      IT(I)=INT(SW(I))
  150 CONTINUE
      DO 170 I=1,N
      TMP=Y(I)
      DO 160 J=1,NP
      TMP=TMP-BETA(1,J)*X(I,J)
  160 CONTINUE
      SW(I)=TMP 
  170 CONTINUE
C
      DO 300 KOL=1,NP
      DO 270 I=1,N
      RSI=SW(I)/SIG
c-    RSID=(Y(I)-TMP)/DS 
      WRS=PSY(RSI)
      SUM2=DBLE(WRS)
C     IF (IPSI.GE.1.AND.IPSI.LE.3.AND.ABS(RSI).GE.DCHI) WRS=0.0
      IF (DELTA(I).EQ.1.0) THEN
        IF (RSI.LE.TL.OR.RSI.GE.TU) GOTO 270
        WRI=PSP(RSI)
        SUM1=DBLE(WRI)
        JJ=I
        GOTO 250
      ELSE 
        IF (RSI.GE.TU) GOTO 270        
        CALL NTRP0L(SW(I),NU,SZ,II)      
        IF (II.GE.NU-1) THEN
          JJ=IT(NU)
          RSI=SZ(NU)/SIG
          WRI=PSP(RSI)
          SUM1=DBLE(WRI)
          WRS=PSY(RSI)
          SUM2=DBLE(WRS)
          GOTO 250
        ENDIF
        AI=YY(II+1)  
        IF (AI.LT.1.E-5) AI=1.0/FLOAT(N)
        SUMP=0.D0
        SUMR=0.D0
        SUM1=0.D0
        SUM2=0.D0
        DO 200 J=II+1,NU
        TMP=SZ(J)/SIG 
        IF (TMP.GE.TU.AND.SUMP.EQ.0.D0) GOTO 270
        JJ=IT(J)
        WRI=PSP(TMP)
        PIJ=DD(J)
        SUMP=SUMP+DBLE(WRI*PIJ*X(JJ,KOL))/DS
        SUMR=SUMR+DBLE(WRI*PIJ*TMP)/DS
        WRS=PSY(TMP)
        SUM1=SUM1+DBLE(WRS*PIJ*X(JJ,KOL))/DS
        SUM2=SUM2+DBLE(WRS*PIJ*TMP)/DS        
  200   CONTINUE
        SUMP=SUMP/DBLE(AI)
        SUMR=SUMR/DBLE(AI)
        SUM1=SUM1/DBLE(AI)
        SUM2=SUM2/DBLE(AI)
        DO 210 J=1,NP
        RES(J,KOL)=RES(J,KOL)-DBLE(X(I,J))*SUMP
        IF (KOL.EQ.NP) RES(J,NP1)=RES(J,NP1)-DBLE(X(I,J))*SUMR
  210   CONTINUE
        RES(NP1,KOL)=RES(NP1,KOL)-SUM1
        IF (KOL.EQ.NP) RES(NP1,NP1)=RES(NP1,NP1)-SUM2
        GOTO 270
      ENDIF    
  250   DO 260 J=1,NP
        RES(J,KOL)=RES(J,KOL)-DBLE(X(I,J))*SUM1*DBLE(X(JJ,KOL))/DS
        IF (KOL.EQ.NP) 
     *  RES(J,NP1)=RES(J,NP1)-DBLE(X(I,J))*SUM1*DBLE(RSI)/DS
  260   CONTINUE
        RES(NP1,KOL)=RES(NP1,KOL)-SUM2*DBLE(X(JJ,KOL))/DS
        IF (KOL.EQ.NP) 
     *  RES(NP1,NP1)=RES(NP1,NP1)-SUM2*DBLE(RSI)/DS
  270 CONTINUE
  300 CONTINUE
      DO 400 I=1,NP1
      DO 400 J=1,NP1
  400 RES(I,J)=RES(I,J)/DFLOAT(N)
      RETURN
      END    
C
      SUBROUTINE NRMRRES(DBETA,DS,X,DELTA,Y,N,NP,MDX,RES)
C
C.......................................................................
C   PROGRAMMER : A. RANDRIAMIHARISOA
C.......................................................................
C
      DOUBLE PRECISION DBETA(NP+1),X(MDX,NP),Y(N),DS,RSID,AI,
     +                 RES(NP+1),TMP,TTL,TTU,SUM,SUM1 
      DIMENSION DELTA(N)
      EXTERNAL PSY,RHO
      COMMON/PSIPR/IPSI,CPSI,H1,H2,H3,XK,DCHI
      TU=1.E9
      IF (IPSI.EQ.2) TU=H3
      IF (IPSI.EQ.3) TU=1.0
      IF (IPSI.EQ.4) TU=XK
      TL=-TU
      TTU=DBLE(TU)
      TTL=DBLE(TL)
      DO 100 I=1,NP+1
 100  RES(I)=0.D0
      DO 270 I=1,N
      TMP=Y(I)
      DO 120 J=1,NP
      TMP=TMP-DBETA(J)*X(I,J)
  120 CONTINUE
      RSID=TMP/DS
      RSI=SNGL(RSID)      
      IF (DELTA(I).EQ.1.0) THEN
        IF (RSI.GE.TU) GOTO 160
        IF (RSI.LE.TL) GOTO 160
        WRI=PSY(RSI)
        DO 150 J=1,NP
        RES(J)=RES(J)+DBLE(WRI)*X(I,J)
  150   CONTINUE
  160   WRI=RHO(RSI)
        RES(NP+1)=RES(NP+1)+DBLE(WRI)
      ELSE
        CALL GAUSSD(1,RSID,TMP)
        AI=1.D0-TMP
        IF (AI.LT.1.D-6) AI=1.D-6
        CALL REFSNRM(RSID,AI,0,SUM,SUM1)
        DO 250 J=1,NP
        RES(J)=RES(J)+SUM*X(I,J)/AI
  250   CONTINUE
        RES(NP+1)=RES(NP+1)+SUM1/AI
      ENDIF    
  270 CONTINUE
      DO 300 I=1,NP+1
  300 RES(I)=RES(I)/DFLOAT(N)
      RES(NP+1)=RES(NP+1)-0.5D0
      RETURN
      END
C
      SUBROUTINE NRM1RES(DBETA,DS,X,DELTA,Y,N,NP,MDX,RES)
C
C.......................................................................
C   PROGRAMMER : A. RANDRIAMIHARISOA
C.......................................................................
C
      DOUBLE PRECISION DBETA(NP+1),X(MDX,NP),Y(N),DS,RSID,AI,
     +       RES(NP),TMP,TTL,TTU,SUM,SUM1
      DIMENSION DELTA(N)
      EXTERNAL PSY
      COMMON/PSIPR/IPSI,CPSI,H1,H2,H3,XK,DCHI
      TU=1.E9
      IF (IPSI.EQ.2) TU=H3
      IF (IPSI.EQ.3) TU=1.0
      IF (IPSI.EQ.4) TU=XK
      TL=-TU
      TTU=DBLE(TU)
      TTL=DBLE(TL)
      DO 100 I=1,NP 
 100  RES(I)=0.D0
      DO 270 I=1,N
      TMP=Y(I)
      DO 120 J=1,NP
      TMP=TMP-DBETA(J)*X(I,J)
  120 CONTINUE
      RSID=TMP/DS
      RSI=SNGL(RSID)      
      IF (DELTA(I).EQ.1.0) THEN
        IF (RSI.GE.TU) GOTO 270
        IF (RSI.LE.TL) GOTO 270
        WRI=PSY(RSI)
        DO 150 J=1,NP
        RES(J)=RES(J)+DBLE(WRI)*X(I,J)
  150   CONTINUE
c 160   WRI=RHO(RSI)
c       RES(NP+1)=RES(NP+1)+DBLE(WRI)
      ELSE 
        CALL GAUSSD(1,RSID,TMP)
        AI=1.D0-TMP
        IF (AI.LT.1.D-6) AI=1.D-6
        CALL REFSNRM(RSID,AI,1,SUM,SUM1)
        DO 250 J=1,NP
        RES(J)=RES(J)+SUM*X(I,J)/AI
  250   CONTINUE
c       RES(NP+1)=RES(NP+1)+SUM1/AI
      ENDIF    
  270 CONTINUE
      DO 300 I=1,NP+1
  300 RES(I)=RES(I)/DFLOAT(N)
c     RES(NP+1)=RES(NP+1)-0.5D0
      RETURN
      END
C
      SUBROUTINE NRM2RES(DS,DBETA,X,DELTA,Y,N,NP,MDX,RES)
C
C.......................................................................
C   PROGRAMMER : A. RANDRIAMIHARISOA
C.......................................................................
C
      DOUBLE PRECISION DBETA(NP),X(MDX,NP),Y(N),DS,RSID,AI,
     +       RES,TMP,TTL,TTU,SUM,SUM1
      DIMENSION DELTA(N)
      EXTERNAL RHO
      COMMON/PSIPR/IPSI,CPSI,H1,H2,H3,XK,DCHI
      TU=1.E9
      IF (IPSI.EQ.2) TU=H3
      IF (IPSI.EQ.3) TU=1.0
      IF (IPSI.EQ.4) TU=XK
      TL=-TU
      TTU=DBLE(TU)
      TTL=DBLE(TL)
      RES=0.D0
      DO 270 I=1,N
      TMP=Y(I)
      DO 120 J=1,NP
      TMP=TMP-DBETA(J)*X(I,J)
  120 CONTINUE
      RSID=TMP/DS
      RSI=SNGL(RSID)      
      IF (DELTA(I).EQ.1.0) THEN
        WRI=RHO(RSI)
        RES=RES+DBLE(WRI)
      ELSE 
        CALL GAUSSD(1,RSID,TMP)
        AI=1.D0-TMP
        IF (AI.LT.1.D-6) AI=1.D-6
        CALL REFSNRM(RSID,AI,2,SUM,SUM1)
        RES=RES+SUM1/AI
      ENDIF    
  270 CONTINUE
      RES=RES/DFLOAT(N-NP)-0.5D0
      RETURN
      END
C
C==========================================================================
C
      DOUBLE PRECISION FUNCTION FUNORM(DX,WGT,N,EXU,EXV)
      DIMENSION WGT(N)
      DOUBLE PRECISION ANS,EXU,DX
      EXTERNAL EXU,EXV
      COMMON/PSIPR/IPSI,C,H1,H2,H3,XK,D
C
C  Initializations
C
      TMP=WGT(1)
      ANS=EXU(DX)
      FUNORM=0.D0
      IF (ANS.EQ.0.D0) RETURN
      TMP=SNGL(DX)
      TMP=EXV(TMP)
      FUNORM=DBLE(TMP)*ANS
      RETURN
      END
C 
      SUBROUTINE REFIJAC(DPAR,X,Y,DELTA,N,NP,MDX,LINT,RES)     
C
C.......................................................................
C   PROGRAMMER : A. RANDRIAMIHARISOA
C.......................................................................
C
      DIMENSION DELTA(N)
      DOUBLE PRECISION DPAR(NP+1),X(MDX,NP),Y(N),RES(NP+1,NP+1),DS
      COMMON/PSIPR/IPSI,CPSI,H1,H2,H3,XK,DCHI
      DS=DPAR(NP+1)
      IF (DS.LE.1.E-4) DS=1.E-4
      IF (LINT.EQ.0.OR.LINT.EQ.3) THEN 
        CALL NRMRJAC(DPAR,DS,X,DELTA,Y,N,NP,MDX,RES)
      ELSE
        CALL GMBRJAC(DPAR,DS,X,DELTA,Y,LINT,N,NP,MDX,RES)            
      ENDIF
      RETURN
      END
C
      SUBROUTINE REFSNRM(RSID,AI,IOPT,SUMPSY,SUMRHO)
      DIMENSION WGT(1)
      DOUBLE PRECISION RSID,AI,SUMPSY,SUMRHO,TMP,TMPS,TTL,TTU,LOW,HI,
     +       TILD,ERRSTD,WORK,FUNORM,DGAUSS
      EXTERNAL RHO,PSY,FUNORM,DGAUSS
      COMMON/PSIPR/IPSI,CPSI,H1,H2,H3,XK,DCHI
      COMMON/INTEGN/AINTEG(4),IWORK(80),WORK(320)
      DATA KEY,LIMIT,TILD/1,80,0.001D0/
      TU=1.E9
      IF (IPSI.EQ.2) TU=H3
      IF (IPSI.EQ.3) TU=1.0
      IF (IPSI.EQ.4) TU=XK
      TL=-TU
      TTU=DBLE(TU)
      TTL=DBLE(TL)
      RHOU=RHO(TU)
      RHOL=RHO(TL)
      WGT(1)=0.0
C
      LOW=RSID
      IF (RSID.LT.TTL) LOW=TTL
      SUMPSY=0.D0
      IF (IOPT.EQ.2) GOTO 200
      HI=10.D0
      IF (TU.LT.10.0) HI=TTU
      IF (TTL.LE.RSID.AND.RSID.LE.TTU.AND.LOW.LT.HI) THEN
        CALL INTGRD(FUNORM,WGT,1,DGAUSS,PSY,LOW,HI,TILD,TILD,KEY,
     +  LIMIT,SUMPSY,ERRSTD,NEVAL,IERR,WORK,IWORK)
        IF (IERR.NE.0) CALL MESSGE(400+IERR,'FUNORM',0)
      ENDIF
  200 SUMRHO=0.D0
      IF (IOPT.EQ.1) RETURN
      IF (TTU.LE.RSID) THEN
        SUMRHO=AI*DBLE(RHOU)
        RETURN
      ENDIF
      TMPS=0.D0
      LOW=RSID
      IF (LOW.LT.TTL) THEN
        CALL GAUSSD(1,TTL,TMP)
        TMPS=(TMP-1.D0+AI)*DBLE(RHOL)
        LOW=TTL
      ENDIF
      HI=TTU
      IF (HI.GT.10.D0) HI=10.D0
      IF (LOW.LT.TTU.AND.LOW.LT.HI) THEN
        CALL INTGRD(FUNORM,WGT,1,DGAUSS,RHO,LOW,HI,TILD,TILD,KEY,
     +  LIMIT,TMP,ERRSTD,NEVAL,IERR,WORK,IWORK)
        IF (IERR.NE.0) CALL MESSGE(401+IERR,'FUNORM',0)
        TMPS=TMPS+TMP
        LOW=TTU
      ENDIF
      HI=10.D0
      IF (TTU.LT.HI) THEN
       CALL GAUSSD(1,TTU,TMP)
       TMPS=TMPS+(1.D0-TMP)*DBLE(RHOU)
      ENDIF
      SUMRHO=TMPS
      RETURN
      END
C
      SUBROUTINE NRMRJAC(DBETA,DS,X,DELTA,Y,N,NP,MDX,RES)                   
C
C.......................................................................
C   PROGRAMMER : A. RANDRIAMIHARISOA
C.......................................................................
C
      DIMENSION DELTA(N)
      DOUBLE PRECISION DS,DBETA(NP+1),X(MDX,NP),Y(N),SUM1,SUM2,
     +       RSID,AI,FACT,TMP,TMP0,TMPAI,TTL,TTU,RES(NP+1,NP+1),
     +       DGAUSS,FUNORM 
      EXTERNAL PSP,PSY,RHO,FUNORM,DGAUSS
      COMMON/PSIPR/IPSI,CPSI,H1,H2,H3,XK,DCHI
      TU=1.E9
      IF (IPSI.EQ.2) TU=H3
      IF (IPSI.EQ.3) TU=1.0
      IF (IPSI.EQ.4) TU=XK
      TL=-TU
      TTU=DBLE(TU)
      TTL=DBLE(TL)
      NP1=NP+1
      DO 100 I=1,NP1
      DO 100 J=1,NP1
 100  RES(I,J)=0.D0
      DO 300 KOL=1,NP
      DO 270 I=1,N
      TMP=0.D0
      DO 120 J=1,NP
      TMP=TMP+DBETA(J)*X(I,J)
  120 CONTINUE
      RSID=(Y(I)-TMP)/DS 
      RSI=SNGL(RSID)
      WRS=PSY(RSI)
C     IF (IPSI.GE.1.AND.IPSI.LE.3.AND.ABS(RSI).GE.DCHI) WRS=0.0
      IF (DELTA(I).EQ.1.0) THEN
        WRI=PSP(RSI)
        DO 140 J=1,NP
        RES(J,KOL)=RES(J,KOL)-X(I,J)*DBLE(WRI)*X(I,KOL)/DS
        IF (KOL.EQ.NP) 
     *  RES(J,NP1)=RES(J,NP1)-X(I,J)*DBLE(WRI)*RSID/DS
  140   CONTINUE
        RES(NP1,KOL)=RES(NP1,KOL)-DBLE(WRS)*X(I,KOL)/DS
        IF (KOL.EQ.NP) 
     *  RES(NP1,NP1)=RES(NP1,NP1)-DBLE(WRS)*RSID/DS
      ELSE 
        CALL GAUSSD(1,RSID,TMP)
        TMP0=DGAUSS(RSID)
        AI=1.D0-TMP
        TMPAI=TMP0/AI
        IF (AI.LT.1.D-5) THEN
          AI=1.D-5
          TMPAI=RSID
        ENDIF
        CALL REFSNRM(RSID,AI,0,SUM1,SUM2)
        FACT=TMPAI*(SUM1/AI+DBLE(WRS))
        DO 170 J=1,NP
        RES(J,KOL)=RES(J,KOL)+FACT*X(I,J)*X(I,KOL)/DS
        IF (KOL.EQ.NP) RES(J,NP1)=RES(J,NP1)+FACT*X(I,J)*RSID/DS
  170   CONTINUE
        WRI=RHO(RSI)
        FACT=TMPAI*(SUM2/AI+DBLE(WRI))
        RES(NP1,KOL)=RES(NP1,KOL)+FACT*X(I,KOL)/DS
        IF (KOL.EQ.NP) RES(NP1,NP1)=RES(NP1,NP1)+FACT*RSID/DS
      ENDIF    
  270 CONTINUE
  300 CONTINUE
      DO 400 I=1,NP1
      DO 400 J=1,NP1
  400 RES(I,J)=RES(I,J)/DFLOAT(N)
      RETURN
      END    
C      
      SUBROUTINE REFSGMB(ITYP,RSID,AI,IOPT,SUMPSY,SUMRHO)
      DIMENSION WGT(2)
      DOUBLE PRECISION RSID,AI,SUMPSY,SUMRHO,TMP,TMPS,TTL,TTU,LOW,MMU,
     +       HI,TILD,ERRSTD,WORK,XKINT,FGUMBL,FUGMBL,DGUMBL 
      EXTERNAL RHO,PSY,FGUMBL,FUGMBL,DGUMBL
      REAL MU0
      COMMON/PSIPR/IPSI,CPSI,H1,H2,H3,XK,DCHI
      COMMON/INTEGN/AINTEG(4),IWORK(80),WORK(320)
      DATA KEY,LIMIT,TILD/1,80,0.001D0/  !,3.5D0,2.8D1/
      DATA XKI,XKINT/0.0,0.D0/
      TU=1.E9
      MU0=0.1351788
      IF (ITYP.EQ.2) MU0=-MU0
      IF (IPSI.EQ.2) TU=H3
      IF (IPSI.EQ.3) TU=1.0
      IF (IPSI.EQ.4) TU=XK
C     IF (IPSI.EQ.2.OR.IPSI.EQ.3) DCHI=TU
      TL=-TU
      TTU=DBLE(TU)
      TTL=DBLE(TL)
      MMU=DBLE(MU0)
      RHOU=RHO(TU)
      RHOL=RHO(TL)
      WGT(1)=FLOAT(ITYP)
      WGT(2)=MU0
      IF (XKI.NE.XK) THEN
        XKI=XK
        HI=DBLE(XK)
        LOW=-HI
        CALL INTGRD(FUGMBL,WGT,2,DGUMBL,PSY,LOW,HI,TILD,TILD,KEY,
     +  LIMIT,XKINT,ERRSTD,NEVAL,IERR,WORK,IWORK)
      ENDIF
C
      LOW=RSID-MMU
      SUMPSY=XKINT
      IF (LOW.LE.TTL) GOTO 100
      IF (IOPT.EQ.2) GOTO 100
c     HI=DGMAX(ITYP)
C     IF (TTU.LT.HI) HI=TTU
      HI=TTU
      IF (TTL.LE.LOW.AND.LOW.LE.TTU) THEN
        CALL INTGRD(FUGMBL,WGT,2,DGUMBL,PSY,LOW,HI,TILD,TILD,KEY,
     +  LIMIT,SUMPSY,ERRSTD,NEVAL,IERR,WORK,IWORK)
        IF (IERR.NE.0) CALL MESSGE(400+IERR,'FUGMBL',0)
      ENDIF
      IF (TTU.LE.LOW) SUMPSY=0.D0
  100 SUMRHO=AI
      IF (TTU.LE.LOW) RETURN
      IF (IOPT.EQ.1) RETURN
      TMPS=0.D0
      IF (LOW.LT.TTL) THEN
        HI=FGUMBL(DBLE(-XK+MU0),ITYP)
        TMP=FGUMBL(RSID,ITYP)
        TMPS=(HI-TMP)*DBLE(RHOL)
        LOW=TTL
      ENDIF    
      HI=TTU
c     IF (HI.GT.DGMAX(ITYP)) HI=DGMAX(ITYP)
      IF (LOW.LT.TTU.AND.LOW.LT.HI) THEN
        CALL INTGRD(FUGMBL,WGT,2,DGUMBL,RHO,LOW,HI,TILD,TILD,KEY,
     +  LIMIT,TMP,ERRSTD,NEVAL,IERR,WORK,IWORK)
        IF (IERR.NE.0) CALL MESSGE(401+IERR,'FUGMBL',0)
        TMPS=TMPS+TMP
        LOW=TTU
      ENDIF
      TMP=FGUMBL(DBLE(XK+MU0),ITYP)
        TMPS=TMPS+(1.D0-TMP)*DBLE(RHOU)
      SUMRHO=TMPS
      RETURN
      END
C
C==========================================================================
C
      DOUBLE PRECISION FUNCTION F0GMBL(DX,WGT,N,EXU,EXV)
      DIMENSION WGT(N)
      DOUBLE PRECISION ANS,EXU,DX
      EXTERNAL EXU,EXV
      COMMON/PSIPR/IPSI,CPSI,H1,H2,H3,XK,DCHI
C
C  Initializations
C
      ITYP=WGT(1)
      ANS=EXU(DX,ITYP)
      F0GMBL=0.D0
      IF (ANS.EQ.0.D0) RETURN
      F0GMBL=DX*ANS
      RETURN
      END
C     
C------------------------------------------------------------------------
C 
      SUBROUTINE SUMXNRM(RES,X,Y,DELTA,BETA,SINI,NL,N,NP,MDX,
     *                   THETA,RS,DD,YY,SZ,SW,IT,NUR0,SUM)        
C
C.......................................................................
C   PROGRAMMER : A. RANDRIAMIHARISOA
C.......................................................................
C
      DIMENSION X(MDX,NP),Y(N),DELTA(N),BETA(1,NP),THETA(N),
     +          RS(N),DD(N),YY(N),SZ(N),SW(N),DUMMY(1,2)
      DOUBLE PRECISION RES(NL),SUM(NL),SUMJ,SINI
      INTEGER IT(N),NUR0(2)
C
      IF (N.LE.0.OR.MDX.LT.N.OR.NP.LE.0) CALL MESSGE(500,'SUMXNRM',1)
      NB=1
      LINT=0
      CALL KMEDMAD(X,Y,IT,DELTA,BETA,N,NP,MDX,NB,NB,LINT,DUMMY,THETA,
     *             RS,DD,YY,SZ,SW)
      NUR0(1)=IT(1)
      NUR0(2)=IT(1)+1
      NU=IT(1)
C
      DO 370 I=1,NL
        SUMJ=0.D0
        RSI=SNGL(RES(I)*SINI)
        CALL NTRP0L(RSI,NU,SZ,II)
        IF (II.GE.NU) THEN
         SUM(I)=RES(I)
         GOTO 370
        ENDIF
        AI=YY(II+1)
        IF (ABS(AI).LT.1.E-5) THEN
           SUM(I)=RES(I)
           GOTO 370
        ENDIF
       DO 350 J=II+1,NU
        TMP=SZ(J)
        PIJ=DD(J)
        SUMJ=SUMJ+DBLE(TMP)*DBLE(PIJ)
  350   CONTINUE
        SUM(I)=SUMJ/(DBLE(AI)*SINI)
  370 CONTINUE
      RETURN
      END
C     
      SUBROUTINE INTXGMB(RS,NL,ITYP,SUM)
      DIMENSION WGT(1)
      DOUBLE PRECISION RS(NL),SUM(NL),LOW,HI,
     +       TILD,ERRSTD,DGMAX(2),WORK,F0GMBL,DGUMBL
      EXTERNAL PSY,F0GMBL,DGUMBL
      COMMON/PSIPR/IPSI,CPSI,H1,H2,H3,XK,DCHI
      COMMON/INTEGN/AINTEG(4),IWORK(80),WORK(320)
      DATA KEY,LIMIT,TILD,DGMAX/1,80,0.001D0,3.5D0,2.8D1/
C
      WGT(1)=FLOAT(ITYP)
      HI=DGMAX(ITYP)
      DO 100 I=1,NL
      LOW=RS(I)
      SUM(I)=0.D0
      IF (LOW.LT.HI) THEN
        CALL INTGRD(F0GMBL,WGT,1,DGUMBL,PSY,LOW,HI,TILD,TILD,KEY,
     +  LIMIT,SUM(I),ERRSTD,NEVAL,IERR,WORK,IWORK)
        IF (IERR.NE.0) CALL MESSGE(400+IERR,'F0GMBL',0)
      ENDIF
  100 CONTINUE
      RETURN
      END
C
C==========================================================================
C
      DOUBLE PRECISION FUNCTION FUGMBL(DX,WGT,N,EXU,EXV)
      DIMENSION WGT(N)
      DOUBLE PRECISION ANS,EXU,DX,XMU
      EXTERNAL EXU,EXV
      COMMON/PSIPR/IPSI,C,H1,H2,H3,XK,D
C
C  Initializations
C
      ITYP=WGT(1)
      XMU=DX+DBLE(WGT(2))
      ANS=EXU(XMU,ITYP)
      FUGMBL=0.D0
      IF (ANS.EQ.0.D0) RETURN
      TMP=SNGL(DX)
      TMP=EXV(TMP)
      FUGMBL=DBLE(TMP)*ANS
      RETURN
      END
C 
      SUBROUTINE GMBRRES(DBETA,DS,X,DELTA,Y,LINT,N,NP,MDX,RES) 
C
C.......................................................................
C   PROGRAMMER : A. RANDRIAMIHARISOA
C   LINT: 1=Gumbel, 2=LogWeibull
C.......................................................................
C
      DOUBLE PRECISION DBETA(NP+1),X(MDX,NP),Y(N),DS,RSID,AI,
     +       SUM,SUM1,RES(NP+1),TMP,TTL,TTU,FGUMBL    
      REAL DELTA(N),MU0
      EXTERNAL PSY,RHO,FGUMBL
      COMMON/PSIPR/IPSI,CPSI,H1,H2,H3,XK,DCHI
      TU=1.E9
      MU0=0.1351788
      IF (LINT.EQ.2) MU0=-MU0
      IF (IPSI.EQ.2) TU=H3
      IF (IPSI.EQ.3) TU=1.0
      IF (IPSI.EQ.4) TU=XK
C     IF (IPSI.EQ.2.OR.IPSI.EQ.3) DCHI=TU
      TL=-TU
      TTU=DBLE(TU)
      TTL=DBLE(TL)
      NP1=NP+1
      DO 100 I=1,NP1
 100  RES(I)=0.D0
      DO 270 I=1,N
      TMP=Y(I)
      DO 120 J=1,NP
      TMP=TMP-DBETA(J)*X(I,J)
  120 CONTINUE
      RSID=TMP/DS
      RSI=SNGL(RSID)      
      IF (DELTA(I).EQ.1.0) THEN
        IF (RSI-MU0.GE.TU) GOTO 160
        IF (RSI-MU0.LE.TL) GOTO 160
        WRI=PSY(RSI-MU0)
        DO 150 J=1,NP
        RES(J)=RES(J)+DBLE(WRI)*X(I,J)
  150   CONTINUE
  160   WRI=RHO(RSI-MU0)
        RES(NP1)=RES(NP1)+DBLE(WRI)
      ELSE 
        TMP=FGUMBL(RSID,LINT)
        AI=1.D0-TMP
        IF (AI.LT.1.D-5) AI=1.D-5
        CALL REFSGMB(LINT,RSID,AI,0,SUM,SUM1)
        DO 250 J=1,NP
        RES(J)=RES(J)+SUM*X(I,J)/AI
  250   CONTINUE
        RES(NP1)=RES(NP1)+SUM1/AI
      ENDIF    
  270 CONTINUE
      DO 300 I=1,NP+1
  300 RES(I)=RES(I)/DFLOAT(N)
      RES(NP1)=RES(NP1)-0.5D0
      RETURN
      END
C
      SUBROUTINE GMB1RES(DBETA,DS,X,DELTA,Y,LINT,N,NP,MDX,RES) 
C
C.......................................................................
C   PROGRAMMER : A. RANDRIAMIHARISOA
C   LINT: 1=Gumbel, 2=LogWeibull
C.......................................................................
C
      DOUBLE PRECISION DBETA(NP),X(MDX,NP),Y(N),DS,RSID,AI,
     +       SUM,SUM1,RES(NP),TMP,TTL,TTU,FGUMBL    
      REAL DELTA(N),MU0
      EXTERNAL PSY,FGUMBL
      COMMON/PSIPR/IPSI,CPSI,H1,H2,H3,XK,DCHI
      TU=1.E9
      MU0=0.1351788
      IF (LINT.EQ.2) MU0=-MU0
      IF (IPSI.EQ.2) TU=H3
      IF (IPSI.EQ.3) TU=1.0
      IF (IPSI.EQ.4) TU=XK
C     IF (IPSI.EQ.2.OR.IPSI.EQ.3) DCHI=TU
      TL=-TU
      TTU=DBLE(TU)
      TTL=DBLE(TL)
      DO 100 I=1,NP
 100  RES(I)=0.D0
      DO 270 I=1,N
      TMP=Y(I)
      DO 120 J=1,NP
      TMP=TMP-DBETA(J)*X(I,J)
  120 CONTINUE
      RSID=TMP/DS
      RSI=SNGL(RSID)      
      IF (DELTA(I).EQ.1.0) THEN
        IF (RSI-MU0.GE.TU) GOTO 270
        IF (RSI-MU0.LE.TL) GOTO 270
        WRI=PSY(RSI-MU0)
        DO 150 J=1,NP
        RES(J)=RES(J)+DBLE(WRI)*X(I,J)
  150   CONTINUE
c 160   WRI=RHO(RSI-MU0)
c       RES(NP1)=RES(NP1)+DBLE(WRI)
      ELSE 
        TMP=FGUMBL(RSID,LINT)
        AI=1.D0-TMP
        IF (AI.LT.1.D-5) AI=1.D-5
        CALL REFSGMB(LINT,RSID,AI,1,SUM,SUM1)
        DO 250 J=1,NP
        RES(J)=RES(J)+SUM*X(I,J)/AI
  250   CONTINUE
c       RES(NP1)=RES(NP1)+SUM1/AI
      ENDIF    
  270 CONTINUE
      DO 300 I=1,NP+1
  300 RES(I)=RES(I)/DFLOAT(N)
c     RES(NP1)=RES(NP1)-0.5D0
      RETURN
      END
C 
      SUBROUTINE GMB2RES(DS,DBETA,X,DELTA,Y,LINT,N,NP,MDX,RES) 
C
C.......................................................................
C   PROGRAMMER : A. RANDRIAMIHARISOA
C   LINT: 1=Gumbel, 2=LogWeibull
C.......................................................................
C
      DOUBLE PRECISION DBETA(NP),X(MDX,NP),Y(N),DS,RSID,AI,
     +       SUM,SUM1,RES,TMP,TTL,TTU,FGUMBL    
      REAL DELTA(N),MU0
      EXTERNAL RHO,FGUMBL
      COMMON/PSIPR/IPSI,CPSI,H1,H2,H3,XK,DCHI
      TU=1.E9
      MU0=0.1351788
      IF (LINT.EQ.2) MU0=-MU0
      IF (IPSI.EQ.2) TU=H3
      IF (IPSI.EQ.3) TU=1.0
      IF (IPSI.EQ.4) TU=XK
C     IF (IPSI.EQ.2.OR.IPSI.EQ.3) DCHI=TU
      TL=-TU
      TTU=DBLE(TU)
      TTL=DBLE(TL)
      NP1=NP+1
      RES=0.D0
      DO 270 I=1,N
      TMP=Y(I)
      DO 120 J=1,NP
      TMP=TMP-DBETA(J)*X(I,J)
  120 CONTINUE
      RSID=TMP/DS
      RSI=SNGL(RSID)      
      IF (DELTA(I).EQ.1.0) THEN
        WRI=RHO(RSI-MU0)
        RES=RES+DBLE(WRI)
      ELSE 
        TMP=FGUMBL(RSID,LINT)
        AI=1.D0-TMP
        IF (AI.LT.1.D-5) AI=1.D-5
        CALL REFSGMB(LINT,RSID,AI,2,SUM,SUM1)
        RES=RES+SUM1/AI
      ENDIF    
  270 CONTINUE
      RES=RES/DFLOAT(N)-0.5D0
      RETURN
      END
C
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      SUBROUTINE GMBRJAC(DBETA,DS,X,DELTA,Y,LINT,N,NP,MDX,RES)        
C
C.......................................................................
C   PROGRAMMER : A. RANDRIAMIHARISOA
C   Jacobian is computed numerically if delta[I]=0 
C.......................................................................
C
      DIMENSION DELTA(N)
      DOUBLE PRECISION DS,DBETA(NP+1),X(MDX,NP),Y(N),SUM1,SUM2,
     +       RSID,AI,FACT,TMP0,TMPAI,TMP,TTL,TTU,RES(NP+1,NP+1),
     +       FGUMBL,DGUMBL,FUGMBL,XEXPD
      REAL MU0
      EXTERNAL PSP,PSY,RHO,FGUMBL,FUGMBL,DGUMBL,XEXPD
      COMMON/PSIPR/IPSI,CPSI,H1,H2,H3,XK,DCHI
      TU=1.E9
      MU0=0.1351788
      IF (LINT.EQ.2) MU0=-MU0
      IF (IPSI.EQ.2) TU=H3
      IF (IPSI.EQ.3) TU=1.0
      IF (IPSI.EQ.4) TU=XK
C     IF (IPSI.EQ.2.OR.IPSI.EQ.3) DCHI=TU
      TL=-TU
      TTU=DBLE(TU)
      TTL=DBLE(TL)
      NP1=NP+1
      DO 100 I=1,NP1
      DO 100 J=1,NP1
  100 RES(I,J)=0.D0
      DO 300 KOL=1,NP
      DO 270 I=1,N
      TMP=0.D0
      DO 120 J=1,NP
      TMP=TMP+DBETA(J)*X(I,J)
  120 CONTINUE
      RSID=(Y(I)-TMP)/DS 
      RSI=SNGL(RSID)
      WRS=PSY(RSI-MU0)
C       IF (IPSI.GE.1.AND.IPSI.LE.3.AND.ABS(RSI).GE.DCHI) WRS=0.0
      IF (DELTA(I).EQ.1.0) THEN
        WRI=PSP(RSI-MU0)
        DO 140 J=1,NP
        RES(J,KOL)=RES(J,KOL)-X(I,J)*DBLE(WRI)*X(I,KOL)/DS
        IF (KOL.EQ.NP) 
     *  RES(J,NP1)=RES(J,NP1)-X(I,J)*DBLE(WRI)*RSID/DS
  140   CONTINUE
        RES(NP1,KOL)=RES(NP1,KOL)-DBLE(WRS)*X(I,KOL)/DS
        IF (KOL.EQ.NP) 
     *  RES(NP1,NP1)=RES(NP1,NP1)-DBLE(WRS)*RSID/DS
      ELSE 
        TMP=FGUMBL(RSID,LINT)
        TMP0=DGUMBL(RSID,LINT)
        AI=1.D0-TMP
        TMPAI=TMP0/AI
        IF (AI.LT.1.D-5) THEN
          AI=1.D-5
          TMPAI=XEXPD(RSID)-1.D0
          IF (LINT.EQ.1) TMPAI=1.D0-XEXPD(-RSID)
        ENDIF
        CALL REFSGMB(LINT,RSID,AI,0,SUM1,SUM2)
        FACT=TMPAI*(SUM1/AI+DBLE(WRS))
        DO 170 J=1,NP
        RES(J,KOL)=RES(J,KOL)+FACT*X(I,J)*X(I,KOL)/DS
        IF (KOL.EQ.NP) 
     *  RES(J,NP1)=RES(J,NP1)+FACT*X(I,J)*RSID/DS
  170   CONTINUE
        WRI=RHO(RSI-MU0)
        FACT=TMPAI*(SUM2/AI+DBLE(WRI))
        RES(NP1,KOL)=RES(NP1,KOL)+FACT*X(I,KOL)/DS
        IF (KOL.EQ.NP) RES(NP1,NP1)=RES(NP1,NP1)+FACT*RSID/DS
      ENDIF    
  270 CONTINUE
  300 CONTINUE
      DO 400 I=1,NP1
      DO 400 J=1,NP1
 400  RES(I,J)=RES(I,J)/DFLOAT(N)
      RETURN
      END
C///////////////// ROUTINES COMPLEMTAIRES POUR SPLUS 7.0 ET R ///////////////
C
C-----------------------------------------------------------------------
C
      FUNCTION RHO(S)
C.......................................................................
C
C   COPYRIGHT 1992 Alfio Marazzi
C
C   AUTHORS : A. MARAZZI / A. RANDRIAMIHARISOA
C.......................................................................
C
C  PURPOSE
C  -------
C  GIVES THE VALUE OF THE INTEGRAL FROM 0 TO S OF PSI(T)
C
      COMMON/PSIPR/IPSI,C,H1,H2,H3,XK,D
      IPS=IABS(IPSI)
      ABST=ABS(S)
      S2=S*S
      IF (IPS.EQ.0) GOTO 100
      IF (IPS.EQ.1) GOTO 200
      IF (IPS.EQ.2) GOTO 300
      IF (IPS.EQ.3) GOTO 400
      IF (IPS.EQ.4) GOTO 500
C      RHO=URHO(S)
C      RETURN
  100 RHO=S2/2.
      RETURN
  200 TMP=S2/2.
      IF (ABST.GT.C) TMP=C*(ABST-C/2.)
      GOTO 600
  300 IF (ABST.GT.H2) GOTO 350
      TMP=S2/2.
      IF (ABST.GT.H1) TMP=H1*(ABST-H1/2.)
      GOTO 600
  350 TMP=0.5*H1*(H3+H2-H1)
      IF (ABST.LT.H3) TMP=TMP-.5*H1*(H3-ABST)**2/(H3-H2)
      GOTO 600
  400 TMP=1./6.
      IF (ABST.GE.1.) GOTO 600
      TMP=(S2*(S2-3)+3)*S2/6.
      GOTO 600
  500 TMP=1.
      IF (ABST.GE.XK) GOTO 600
      S2=S2/(XK**2)
      TMP=(S2*(S2-3.0)+3.0)*S2
  600 RHO=TMP
      RETURN
      END
C
C-----------------------------------------------------------------------
C
      FUNCTION CHI(S)
C.......................................................................
C
C   COPYRIGHT 1992 Alfio Marazzi
C
C   AUTHORS : A. MARAZZI / A. RANDRIAMIHARISOA
C.......................................................................
C
C  PURPOSE
C  -------
C  GIVES THE VALUE OF THE FUNCTION CHI(S)=S*S/2 IF IPSI=0,
C  CHI(S)=HUBER'S CHI FUNCTION IF ABS(IPSI) IS LESS THAN 4,
C  AND CHI(S)=CHIK(S) FOR S-ESTIMATES IF IPS=4.
C
      COMMON/PSIPR/IPSI,C,H1,H2,H3,XK,D
      IF (IPSI.EQ.0) GOTO 100
      IPS=IABS(IPSI)
      IF (IPS.LE.3) GOTO 200
      IF (IPS.EQ.4) GOTO 300
C      CHI=UCHI(S)
C      RETURN
  100 CHI=S*S/2.
      RETURN
  200 ABST=ABS(S)
      PS=AMIN1(D,ABST)
      CHI=PS*PS/2.
      RETURN
  300 TMP=1.
      ABST=ABS(S)
      IF (ABST.GE.XK) GOTO 400
      S2=(S/XK)**2
      TMP=(S2*(S2-3.0)+3.0)*S2
  400 CHI=TMP
      RETURN
      END
C
      FUNCTION PSY(S)
C.......................................................................
C
C   COPYRIGHT 1992 Alfio Marazzi
C
C   AUTHORS : A. MARAZZI / A. RANDRIAMIHARISOA
C.......................................................................
C
C  PURPOSE
C  -------
C  GIVES THE VALUE IN THE POINT T OF THE FUNCTION PSI
C
      COMMON/PSIPR/IPSI,C,H1,H2,H3,XK,D
      IPS=IABS(IPSI)
      ABST=ABS(S)
      IF (IPS.EQ.0) GOTO 100
      IF (IPS.EQ.1) GOTO 200
      IF (IPS.EQ.2) GOTO 300
      IF (IPS.EQ.3) GOTO 400
      IF (IPS.EQ.4) GOTO 500
C
C  PSI(S)=USER PSI FUNCTION
C
C      PSY=UPSI(S)
C      RETURN
C
C  PSI(S)=S
C
  100 PSY=S
      RETURN
C
C  PSI(S,C)=MAX(-C,MIN(C,S))
C
  200 TMP=AMIN1(C,ABST)
      IF (S.LT.0.) TMP=-TMP
      GOTO 600
C
C  PSI(S,H1,H2,H3)=-PSI(-S,H1,H2,H3)
C                 =S FOR 0 .LE. S .LE. H1
C                 =H1 FOR H1 .LE. S .LE. H2
C                 =H1*(H3-S)/(H3-H2) FOR H2 .LE. S .LE. H3
C                 =0 FOR S .GT. H3
C
  300 TMP=0
      IF (ABST.GE.H3) GOTO 600
      IF (ABST.LE.H2) TMP=AMIN1(H1,ABST)
      IF (ABST.GT.H2) TMP=H1*(H3-ABST)/(H3-H2)
      IF (S.LT.0.) TMP=-TMP
      GOTO 600
C
C  PSI(S)=S*[MAX(1-S**2,0)]**2
C
  400 TMP=0.
      IF (ABST.GE.1.) GOTO 600
      TMP=S*(1.-S*S)*(1.-S*S)
      GOTO 600
C
C  PSI(S)=(6/K)*(S/K)*[MAX(1-(S/K)**2,0)]**2
C
  500 TMP=0.
      IF (ABST.GE.XK) GOTO 600
      SK=S/XK
      TMP=(6.*SK/XK)*(1.-SK*SK)*(1.-SK*SK)
  600 PSY=TMP
      RETURN
      END
C
      FUNCTION PSP(S)
C.......................................................................
C
C   COPYRIGHT 1992 Alfio Marazzi
C
C   AUTHORS : A. MARAZZI / A. RANDRIAMIHARISOA
C.......................................................................
C
C  PURPOSE
C  -------
C  GIVES THE VALUE AT THE POINT S OF THE FIRST DERI-
C  VATIVE OF THE FUNCTION PSI .
C
      COMMON/PSIPR/IPSI,C,H1,H2,H3,XK,D
      IPS=IABS(IPSI)
      ABST=ABS(S)
      IF (IPS.EQ.0) GOTO 100
      IF (IPS.EQ.1) GOTO 200
      IF (IPS.EQ.2) GOTO 300
      IF (IPS.EQ.3) GOTO 400
      IF (IPS.EQ.4) GOTO 500
C      PSP=UPSP(S)
C      RETURN
  100 PSP=1.
      RETURN
  200 TMP=0.
      IF (ABST.LE.C) TMP=1.
      GOTO 600
  300 TMP=1.
      IF (ABST.LE.H1) GOTO 600
      TMP=0.
      IF ((ABST.LT.H2).OR.(ABST.GT.H3)) GOTO 600
      TMP=H1/(H2-H3)
      GOTO 600
  400 TMP=0.
      IF (ABST.GE.1.) GOTO 600
      S2=S*S
      TMP=(1.-S2)*(1.-5.*S2)
      GOTO 600
  500 TMP=0.
      IF (ABST.GE.XK) GOTO 600
      S2=(S/XK)**2
      TMP=(6./XK)*(1.-S2)*(1.-5.*S2)/XK
  600 PSP=TMP
      RETURN
      END
C
      SUBROUTINE GAUSSD (KODE,X,P)
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
      IF (KODE.NE.1.AND.KODE.NE.2) CALL MESSGE(500,'GAUSSD',1)
      CALL CERFD(-X*SQR1D2,CD)
      P = .5D0 * CD
      IF (KODE.EQ.2) P=1.D0-P
      RETURN
      END
C
C-----------------------------------------------------------------------
C
      SUBROUTINE CERFD(X,F)
C.......................................................................
C
C   COPYRIGHT 1992 Alfio Marazzi
C
C   AUTHOR : A. RANDRIAMIHARISOA
C                
C.......................................................................
C
      DOUBLE PRECISION   F,X,XEXPD
      DIMENSION          P(5),Q(4),P1(9),Q1(8),P2(6),Q2(5)
      DOUBLE PRECISION   P,Q,P1,Q1,P2,Q2,XMIN,XLARGE,SQRPI,XX,
     *                   RES,XSQ,XNUM,XDEN,XI,XBIG
      INTEGER            ISW,I
      EXTERNAL XEXPD
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
      Y=X
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
   60 RES = RES*XEXPD(-XSQ)
      IF (ISW.EQ.-1) RES = 2.0D0-RES
      GO TO 70
   65 RES = 0.0D0
   70 F = RES
      RETURN
      END
C
      SUBROUTINE MACH(I,X)
C.......................................................................
C
C   COPYRIGHT 1992 Alfio Marazzi
C
C   AUTHOR : A. MARAZZI
C.......................................................................
C
C     SINGLE PRECISION VERSION (DOUBLE PRECISION VERSION FOLLOWS)
C     ***********************************************************
C
C  MACHINE PARAMETERS : TO  ALTER  THIS  SUBROUTINE  FOR  A PARTICULAR 
C  ++++++++++++++++++   ENVIRONMENT, THE DESIRED SET OF DATA STATEMENT
C  SHOULD BE ACTIVATED BY REMOVING THE "C" FROM COLUMN ONE  AND ADDING
C  THE "C" FOR THE TWO LINES AFTER "... VAX FORTRAN (V5) compiler".
C
C  RADIX IS ALWAYS EQUAL TO 2
C  PREC CAN BE FOUND BY CALLING THE ROBETH SUBROUTINE "PRECS"
C  EPMACH IS APPROXIMATELY EQUAL TO THE EXPONENT PART OF PREC
C  EXMIN, XLGMN, YLGMN AND XBIG CAN BE FOUND BY TRIAL AND ERROR
C
C  for VAX AND DEC-ALPHA FORTRAN (V5) compiler
C     DATA RADIX,PREC,EXMIN,XLGMN,YLGMN,XBIG,EPMACH
C    *     /2.,5.960465E-8,-88.722,0.2939E-38,-88.7227,1.7E38,1.0E-7/
C  for IBM-PC F77L compiler
C     DATA RADIX,PREC,EXMIN,XLGMN,YLGMN,XBIG,EPMACH
C    *     /2.,5.43E-20,-87.4,0.118E-37,-87.3327,3.4E38,1.0E-19/
C  for IBM-PC MICROSOFT FORTRAN compiler
C     DATA RADIX,PREC,EXMIN,XLGMN,YLGMN,XBIG,EPMACH
C    *     /2.,5.47522E-18,-103.972,0.701E-45,-103.279,.E.,1.0E-17/
C  for WATCOM F77 compiler (32 bits version)
C     DATA RADIX,PREC,EXMIN,XLGMN,YLGMN,XBIG,EPMACH
C    *     /2.,6.020E-8,-87.336,0.0588E-37,-88.029,34.02E37,1.0E-7/
C  for ULTRIX DEC and ALPHA/OSF1 FORTRAN-77 compiler
      DATA RADIX,PREC,EXMIN,XLGMN,YLGMN,XBIG,EPMACH
     *     /2.,6.02007E-8,-87.336,0.1176E-37,-87.3361,3.401E38,1.0E-7/
C  for SUN FORTRAN Compiler 
c      DATA RADIX,PREC,EXMIN,XLGMN,YLGMN,XBIG,EPMACH
c    *     /2.,6.02007E-8,-103.972,0.1401E-44,-103.279,3.402E38,1.0E-7/
C  for SILICON GRAPHICS MIPS FORTRAN 77 Compiler 
C      DATA RADIX,PREC,EXMIN,XLGMN,YLGMN,XBIG,EPMACH
C    *     /2.,6.02007E-8,-102.88,4.757E-43,-88.0297,3.395E38,1.0E-7/
C  for HP-UX FORTRAN 77 Compiler 
C      DATA RADIX,PREC,EXMIN,XLGMN,YLGMN,XBIG,EPMACH
C    *     /2.,6.02007E-8,-87.33,0.989E-42,-88.0296,3.401E38,1.0E-7/
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
      SUBROUTINE MACHD(I,X)
C.......................................................................
C
C   COPYRIGHT 1992 Alfio Marazzi
C
C   AUTHOR : A. MARAZZI
C.......................................................................
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
C      DATA RADIX,PREC,EXMIN,XLGMN,YLGMN,XBIG,EPMACH
C    *  /2.D0,1.38778D-17,-88.722D0,2.939D-39,-88.7227D0,1.7D38,1.D-17/
C  for IBM-PC F77L compiler
C      DATA RADIX,PREC,EXMIN,XLGMN,YLGMN,XBIG,EPMACH
C    *  /2.D0,5.422D-20,-708.D0,1.D-307,-706.591D0,1.D308,1.D-19/
C  for IBM-PC MICROSOFT FORTRAN compiler
C      DATA RADIX,PREC,EXMIN,XLGMN,YLGMN,XBIG,EPMACH
C    *  /2.D0,5.47522D-18,-745.133D0,0.9D-48,-110.629D0,1.D308,1.D-17/
C  for WATCOM F77 Compiler (32 bits version)
C      DATA RADIX,PREC,EXMIN,XLGMN,YLGMN,XBIG,EPMACH
C    */2.,0.1121D-15,-709.782D0,0.974D-312,-718.433D0,1.797D308,1.0D-17/
C  for ULTRIX DEC and ALPHA/OSF1 FORTRAN-77 compiler
      DATA RADIX,PREC,EXMIN,XLGMN,YLGMN,XBIG,EPMACH
     */2.,1.12133D-16,-707.9D0,2.226D-308,-708.396D0,1.796D308,1.0D-17/
C  for SUN FORTRAN compiler
c     DATA RADIX,PREC,EXMIN,XLGMN,YLGMN,XBIG,EPMACH
c    */2.,1.12133D-16,-745.13D0,0.494D-323,-744.44D0,1.797D308,1.0D-17/
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
      SUBROUTINE MESSGE(NUMBER,ITEXT,ISTOP)
C.......................................................................
C
C   COPYRIGHT 1992 Alfio Marazzi
C
C   AUTHOR : A. RANDRIAMIHARISOA
C.......................................................................
C
      CHARACTER *6 ITEXT, CC*34
      IF (ISTOP.EQ.1) THEN
         CC='Input parameter error(s) in '//ITEXT
C
C Error exit from R
C
       CC='Input parameter error(s) in '//ITEXT
       CALL REXIT(CC)
      ELSE
       CC='Warning message in '//ITEXT
       CALL INTPR(CC,25,NUMBER,1)
      ENDIF
      RETURN
      END
C
C-----------------------------------------------------------------------
C
      SUBROUTINE SRT1(A,N,K1,K2)
C.......................................................................
C
C   COPYRIGHT 1992 Alfio Marazzi
C
C   AUTHOR : A. MARAZZI
C.......................................................................
C
      REAL A(N)
      LOGICAL NPRCHK
C
      NPRCHK=K1.GE.1.AND.K2.GT.K1.AND.K2.LE.N
      IF (.NOT.NPRCHK) CALL MESSGE(500,'SRT1  ',1)
      N1=K2-K1+1
c      I=1
c   10 I=I+I
c      IF (I.LE.N1) GOTO 10
      M=N1
   20 M=M/2
      IF (M.EQ.0) GOTO 90
      K=N1-M
      DO 40 J=1,K
      L=J
   50 IF (L.LT.1) GOTO 40
      LPM=L+M
      LPM1=LPM+K1-1
      L1=L+K1-1
      IF (A(LPM1).GE.A(L1)) GOTO 40
      X=A(LPM1)
      A(LPM1)=A(L1)
      A(L1)=X
      L=L-M
      GOTO 50
   40 CONTINUE
      GOTO 20
   90 CONTINUE
      END
C
C-----------------------------------------------------------------------
C
      SUBROUTINE SRT2(A,B,N,K1,K2)
C.......................................................................
C
C   COPYRIGHT 1992 Alfio Marazzi
C
C   AUTHOR : A. MARAZZI
C.......................................................................
C
      REAL A(N),B(N)
      LOGICAL NPRCHK
C
      NPRCHK=N.GT.0.AND.K1.GE.1.AND.K2.GE.K1.AND.K2.LE.N
      IF (.NOT.NPRCHK) CALL MESSGE(500,'SRT2  ',1)
      N1=K2-K1+1
c      I=1
c   10 I=I+I
c      IF (I.LE.N1) GOTO 10
      M=N1
   20 M=M/2
      IF (M.EQ.0) GOTO 90
      K=N1-M
      DO 40 J=1,K
      L=J
   50 IF (L.LT.1) GOTO 40
      LPM=L+M
      LPM1=LPM+K1-1
      L1=L+K1-1
      IF (A(LPM1).GE.A(L1)) GOTO 40
      X=A(LPM1)
      Y=B(LPM1)
      A(LPM1)=A(L1)
      B(LPM1)=B(L1)
      A(L1)=X
      B(L1)=Y
      L=L-M
      GOTO 50
   40 CONTINUE
      GOTO 20
   90 CONTINUE
      END
C
      SUBROUTINE NRM2(X,N,INCX,MDX,XNRM)
C.......................................................................
C
C   COPYRIGHT 1979 SOCIETY FOR INDUSTRIAL AND APPLIED MATHEMATICS.
C   ALL RIGHTS RESERVED.
C
C   AUTHOR :     LINPACK (SUBROUTINE SNRM2)
C                REPRINTED WITH PERMISSION FROM 
C                LINPACK USER'S GUIDE.
C                ADAPTED FOR ROBETH BY A. MARAZZI
C                ASSIGN STATEMENTS CANCELED (07/10/2010)
C.......................................................................
C
      LOGICAL NPRCHK
      REAL X(MDX)
      INTEGER LABEL
      DOUBLE PRECISION SUM,ZERO,ONE,XMAX,DXI
C     DATA ZERO,ONE/0.0D0,1.0D0/,CUTLO,CUTHI/4.441E-16,1.304E19/
      PARAMETER (ZERO=0.0D0,ONE=1.0D0,CUTLO=4.441E-16,CUTHI=1.304E19)
C
C  PARAMETER CHECK
C
      NPRCHK=INCX.GT.0.AND.INCX*(N-1)+1.LE.MDX
      IF (.NOT.NPRCHK) CALL MESSGE(500,'NRM2  ',1)
C
      IF (N.GT.0) GOTO 10
      XNRM=0.
      GOTO 300
C
   10 CONTINUE
C     ASSIGN 30 TO NEXT
      LABEL=1
      SUM=ZERO
      NN=N*INCX
C
C  BEGIN MAIN LOOP
C
      I=1
   20 DXI=DBLE(X(I))
C     GOTO NEXT,(30,50,70,110)
      GOTO (30,50,70,110), LABEL
   30 IF (ABS(X(I)).GT.CUTLO) GOTO 85
C     ASSIGN 50 TO NEXT
      LABEL=2
      XMAX=ZERO
C
C  PHASE1.  SUM IS ZERO
C
   50 IF (X(I).EQ.0.) GOTO 200
      IF (ABS(X(I)).GT.CUTLO) GOTO 85
C
C  PREPARE FOR PHASE 2.
C
C     ASSIGN 70 TO NEXT
      LABEL=3
      GOTO 105
C
C  PREPARE FOR PHASE 4.
C
  100 I=J
      DXI=DBLE(X(I))
C     ASSIGN 110 TO NEXT
      LABEL=4
      SUM=(SUM/DXI)/DXI
  105 XMAX=DABS(DXI)
      GOTO 115
C
C  PHASE 2.  SUM IS SMALL. SCALE TO AVOID DESTRUCTIVE UNDERFLOW.
C
   70 IF (ABS(X(I)).GT.CUTLO) GOTO 75
C
C  COMMON CODE FOR PHASE 2 AND 4.
C  IN PHASE 4 SUM IS LARGE.  SCALE TO AVOID OVERFLOW.
C
  110 IF (DABS(DXI).LE.XMAX) GOTO 115
      SUM=ONE+SUM*(XMAX/DXI)**2
      XMAX=DABS(DXI)
      GOTO 200
C
  115 SUM=SUM+(DXI/XMAX)**2
      GOTO 200
C
C  PREPARE FOR PHASE 3.
C
   75 SUM=(SUM*XMAX)*XMAX
C
C  SET HITEST=CUTHI/N
C
   85 HITEST=CUTHI/FLOAT(N)
C
C  PHASE3.  SUM IS MID-RANGE.  NO SCALING.
C
      DO 95 J=I,NN,INCX
      IF (ABS(X(J)).GE.HITEST) GOTO 100
   95 SUM=SUM+X(J)*DBLE(X(J))
      XNRM=DSQRT(SUM)
      GOTO 300
C
  200 CONTINUE
      I=I+INCX
      IF (I.LE.NN) GOTO 20
C
C  END MAIN LOOP
C
      XNRM=XMAX*DSQRT(SUM)
  300 CONTINUE
      RETURN
      END
C
C-----------------------------------------------------------------------
C
      FUNCTION XEXP(X)
C.......................................................................
C
C   COPYRIGHT 1992 Alfio Marazzi
C
C   AUTHOR : A. RANDRIAMIHARISOA
C.......................................................................
      DATA NCALL,DMIN,DMAX,XBIG/0,0.,0.,0./
      IF (NCALL.EQ.0) THEN
        CALL MACH(3,DMIN)
        CALL MACH(6,XBIG)
        XBIG=XBIG/10.
        DMAX=ALOG(XBIG)
        NCALL=1
      ENDIF
C
C  EXTENDED EXPONENTIAL FUNCTION
C
      IF (X.LE.DMIN) THEN
        XEXP=0.
      ELSEIF (X.GE.DMAX) THEN
        XEXP=XBIG
      ELSE
        XEXP=EXP(X)
      ENDIF
      RETURN
      END
C
      DOUBLE PRECISION FUNCTION XEXPD(X)
C.......................................................................
C
C   COPYRIGHT 1992 Alfio Marazzi
C
C   AUTHOR : A. RANDRIAMIHARISOA
C.......................................................................
      DOUBLE PRECISION X,DMIN,DMAX,XBIG
      DATA NCALL,DMIN,DMAX,XBIG/0,0.D0,0.D0,0.D0/
      IF (NCALL.EQ.0) THEN
        CALL MACHD(3,DMIN)
        CALL MACHD(6,XBIG)
        XBIG=XBIG/10.D0
        DMAX=DLOG(XBIG)
        NCALL=1
      ENDIF
C
C  EXTENDED EXPONENTIAL FUNCTION
C
      IF (X.LE.DMIN) THEN
        XEXPD=0.D0
      ELSEIF (X.GE.DMAX) THEN
        XEXPD=XBIG
      ELSE
        XEXPD=DEXP(X)
      ENDIF
      RETURN
      END
C
C
C***********************************************************************
C********************* INTEGRATION SUBROUTINES *************************
C
      SUBROUTINE INTGRD(F,FARR,N,FEXT,GEXT,A,B,EPSABS,EPSREL,KEY,LIMIT,
     1           RESULT,ABSERR,NEVAL,IER,WORK,IWORK)
C.......................................................................
C
C   R O B E T H  -  R O B S Y S   RELEASE 3.0 (COPYRIGHT) 1985, 1990
C
C   PROGRAMMER : QUADPACK
C                ADAPTED FOR ROBETH BY A. RANDRIAMIHARISOA
C.......................................................................
C
      DOUBLE PRECISION A,ABSERR,B,EPSABS,EPSREL,F,RESULT,FEXT,WORK
      INTEGER IER,KEY,LAST,LIMIT,NEVAL,ALIST,BLIST,ELIST,RLIST
C
      DIMENSION FARR(N),WORK(4*LIMIT),IWORK(LIMIT)
C
      EXTERNAL F,FEXT,GEXT
C
C         LIMIT IS THE MAXIMUM NUMBER OF SUBINTERVALS ALLOWED IN THE
C         SUBDIVISION PROCESS OF QAGE1D. TAKE CARE THAT LIMIT.GE.1.
C
C**** DATA LIMIT/500/
C
      IF ((EPSABS.LT.0..AND.EPSREL.LT.0.).OR.LIMIT.LE.1
     1   .OR.LIMIT.GT.500)  CALL MESSGE(500,'INTGRD',1)
      ALIST=1
      BLIST=ALIST+LIMIT
      RLIST=BLIST+LIMIT
      ELIST=RLIST+LIMIT
      CALL QAGE1D(F,FARR,N,FEXT,GEXT,A,B,EPSABS,EPSREL,KEY,LIMIT,
     1     RESULT,ABSERR,NEVAL,IER,
     2     WORK,WORK(BLIST),WORK(RLIST),WORK(ELIST),IWORK,LAST)
C
      RETURN
      END
C
C-----------------------------------------------------------------------
C
      SUBROUTINE QAGE1D(F,FARR,N,FEXT,GEXT,A,B,EPSABS,EPSREL,KEY,LIMIT,
     *  RESULT,ABSERR,NEVAL,IER,ALIST,BLIST,RLIST,ELIST,IORD,LAST)
C.......................................................................
C
C   R O B E T H  -  R O B S Y S   RELEASE 3.0 (COPYRIGHT) 1985, 1990
C
C   PROGRAMMER : QUADPACK
C                ADAPTED FOR ROBETH BY A. RANDRIAMIHARISOA
C.......................................................................
C
      DOUBLE PRECISION A,ABSERR,ALIST,AREA,AREA1,AREA12,AREA2,A1,A2,B,
     *  BLIST,B1,B2,C,DABS,DEFABS,DEFAB1,DEFAB2,DMAX1,ELIST,EPMACH,
     *  EPSABS,EPSREL,ERRBND,ERRMAX,ERROR1,ERROR2,ERRO12,ERRSUM,F,OFLOW,
     *  RESABS,RESULT,RLIST,UFLOW,FEXT
      INTEGER IER,IORD,IROFF1,IROFF2,K,KEY,KEYF,LAST,LIMIT,MAXERR,NEVAL,
     *  NRMAX
C
      DIMENSION ALIST(LIMIT),BLIST(LIMIT),ELIST(LIMIT),IORD(LIMIT),
     *  RLIST(LIMIT),FARR(N)
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
      CALL MACHD(7,EPMACH)
      CALL MACHD(4,UFLOW)
      CALL MACHD(6,OFLOW)
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
      IF ((EPSABS.LT.0..AND.EPSREL.LT.0.).OR.LIMIT.LE.1
     1   .OR.LIMIT.GT.500)  CALL MESSGE(500,'QAGE1D',1)
      IER = 0
C
C           FIRST APPROXIMATION TO THE INTEGRAL
C           -----------------------------------
C
      KEYF = KEY
      IF(KEY.LE.0) KEYF = 1
      IF(KEY.GE.7) KEYF = 6
      C = KEYF
      NEVAL = 0
      IF (KEYF.EQ.1)
     *  CALL Q1K15D(F,FARR,N,FEXT,GEXT,A,B,RESULT,ABSERR,DEFABS,RESABS)
CC    IF (KEYF.EQ.2)
CC   *  CALL Q1K21D(F,FARR,N,FEXT,GEXT,A,B,RESULT,ABSERR,DEFABS,RESABS)
CC    IF (KEYF.EQ.3)
CC   *  CALL Q1K31D(F,FARR,N,FEXT,GEXT,A,B,RESULT,ABSERR,DEFABS,RESABS)
CC    IF (KEYF.EQ.4)
CC   *  CALL Q1K41D(F,FARR,N,FEXT,GEXT,A,B,RESULT,ABSERR,DEFABS,RESABS)
CC    IF (KEYF.EQ.5)
CC   *  CALL Q1K51D(F,FARR,N,FEXT,GEXT,A,B,RESULT,ABSERR,DEFABS,RESABS)
CC    IF (KEYF.EQ.6)
CC   *  CALL Q1K61D(F,FARR,N,FEXT,GEXT,A,B,RESULT,ABSERR,DEFABS,RESABS)
      LAST = 1
      RLIST(1) = RESULT
      ELIST(1) = ABSERR
      IORD(1) = 1
C
C           TEST ON ACCURACY.
C
      ERRBND = DMAX1(EPSABS,EPSREL*DABS(RESULT))
      IF(ABSERR.LE.5.0D+01*EPMACH*DEFABS.AND.ABSERR.GT.ERRBND) IER = 2
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
        A1 = ALIST(MAXERR)
        B1 = 5.0D-01*(ALIST(MAXERR)+BLIST(MAXERR))
        A2 = B1
        B2 = BLIST(MAXERR)
        IF (KEYF.EQ.1)
     *  CALL Q1K15D(F,FARR,N,FEXT,GEXT,A1,B1,AREA1,ERROR1,RESABS,DEFAB1)
CC      IF (KEYF.EQ.2)
CC   *  CALL Q1K21D(F,FARR,N,FEXT,GEXT,A1,B1,AREA1,ERROR1,RESABS,DEFAB1)
CC      IF (KEYF.EQ.3)
CC   *  CALL Q1K31D(F,FARR,N,FEXT,GEXT,A1,B1,AREA1,ERROR1,RESABS,DEFAB1)
CC      IF (KEYF.EQ.4)
CC   *  CALL Q1K41D(F,FARR,N,FEXT,GEXT,A1,B1,AREA1,ERROR1,RESABS,DEFAB1)
CC      IF (KEYF.EQ.5)
CC   *  CALL Q1K51D(F,FARR,N,FEXT,GEXT,A1,B1,AREA1,ERROR1,RESABS,DEFAB1)
CC      IF (KEYF.EQ.6)
CC   *  CALL Q1K61D(F,FARR,N,FEXT,GEXT,A1,B1,AREA1,ERROR1,RESABS,DEFAB1)
        IF (KEYF.EQ.1)
     *  CALL Q1K15D(F,FARR,N,FEXT,GEXT,A2,B2,AREA2,ERROR2,RESABS,DEFAB2)
CC      IF (KEYF.EQ.2)
CC   *  CALL Q1K21D(F,FARR,N,FEXT,GEXT,A2,B2,AREA2,ERROR2,RESABS,DEFAB2)
CC      IF (KEYF.EQ.3)
CC   *  CALL Q1K31D(F,FARR,N,FEXT,GEXT,A2,B2,AREA2,ERROR2,RESABS,DEFAB2)
CC      IF (KEYF.EQ.4)
CC   *  CALL Q1K41D(F,FARR,N,FEXT,GEXT,A2,B2,AREA2,ERROR2,RESABS,DEFAB2)
CC      IF (KEYF.EQ.5)
CC   *  CALL Q1K51D(F,FARR,N,FEXT,GEXT,A2,B2,AREA2,ERROR2,RESABS,DEFAB2)
CC      IF (KEYF.EQ.6)
CC   *  CALL Q1K61D(F,FARR,N,FEXT,GEXT,A2,B2,AREA2,ERROR2,RESABS,DEFAB2)
C
C           IMPROVE PREVIOUS APPROXIMATIONS TO INTEGRAL
C           AND ERROR AND TEST FOR ACCURACY.
C
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
        IF(DMAX1(DABS(A1),DABS(B2)).LE.(1.0D+00+C*1.0D+03*
     *  EPMACH)*(DABS(A2)+1.0D+04*UFLOW)) IER = 3
C
C           APPEND THE NEWLY-CREATED INTERVALS TO THE LIST.
C
    8   IF(ERROR2.GT.ERROR1) GO TO 10
        ALIST(LAST) = A2
        BLIST(MAXERR) = B1
        BLIST(LAST) = B2
        ELIST(MAXERR) = ERROR1
        ELIST(LAST) = ERROR2
        GO TO 20
   10   ALIST(MAXERR) = A2
        ALIST(LAST) = A1
        BLIST(LAST) = B1
        RLIST(MAXERR) = AREA2
        RLIST(LAST) = AREA1
        ELIST(MAXERR) = ERROR2
        ELIST(LAST) = ERROR1
C
C           CALL SUBROUTINE QSORTD TO MAINTAIN THE DESCENDING ORDERING
C           IN THE LIST OF ERROR ESTIMATES AND SELECT THE SUBINTERVAL
C           WITH THE LARGEST ERROR ESTIMATE (TO BE BISECTED NEXT).
C
   20   CALL QSORTD(LIMIT,LAST,MAXERR,ERRMAX,ELIST,IORD,NRMAX)
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
  999 IF (IER.NE.0) CALL MESSGE(400+IER,'QAGE1 ',0)
      RETURN
      END
C
C-----------------------------------------------------------------------
C
      SUBROUTINE Q1K15D
     *  (F,FARR,N,FEXT,GEXT,A,B,RESULT,ABSERR,RESABS,RESASC)
C.......................................................................
C
C   R O B E T H  -  R O B S Y S   RELEASE 3.0 (COPYRIGHT) 1985, 1990
C
C   PROGRAMMER : QUADPACK
C                ADAPTED FOR ROBETH BY A. RANDRIAMIHARISOA
C.......................................................................
C
      DOUBLE PRECISION A,ABSC,ABSERR,B,CENTR,DABS,DHLGTH,DMAX1,DMIN1,
     *  EPMACH,F,FC,FSUM,FVAL1,FVAL2,FV1,FV2,HLGTH,OFLOW,RESABS,RESASC,
     *  RESG,RESK,RESKH,RESULT,UFLOW,WG,WGK,XGK,FEXT
      INTEGER J,JTW,JTWM1
      EXTERNAL F,FEXT,GEXT
C
      DIMENSION FV1(7),FV2(7),WG(4),WGK(8),XGK(8),FARR(N)
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
      CALL MACHD(7,EPMACH)
      CALL MACHD(4,UFLOW)
      CALL MACHD(6,OFLOW)
C
      CENTR = 5.0D-01*(A+B)
      HLGTH = 5.0D-01*(B-A)
      DHLGTH = DABS(HLGTH)
C
C           COMPUTE THE 15-POINT KRONROD APPROXIMATION TO
C           THE INTEGRAL, AND ESTIMATE THE ABSOLUTE ERROR.
C
      FC = F(CENTR,FARR,N,FEXT,GEXT)
      RESG = FC*WG(4)
      RESK = FC*WGK(8)
      RESABS = DABS(RESK)
      DO 10 J=1,3
        JTW = J*2
        ABSC = HLGTH*XGK(JTW)
        FVAL1 = F(CENTR-ABSC,FARR,N,FEXT,GEXT)
        FVAL2 = F(CENTR+ABSC,FARR,N,FEXT,GEXT)
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
        FVAL1 = F(CENTR-ABSC,FARR,N,FEXT,GEXT)
        FVAL2 = F(CENTR+ABSC,FARR,N,FEXT,GEXT)
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
C
C-----------------------------------------------------------------------
C
      SUBROUTINE QSORTD(LIMIT,LAST,MAXERR,ERMAX,ELIST,IORD,NRMAX)
C.......................................................................
C
C   R O B E T H  -  R O B S Y S   RELEASE 3.0 (COPYRIGHT) 1985, 1990
C
C   PROGRAMMER : QUADPACK
C                ADAPTED FOR ROBETH BY A. RANDRIAMIHARISOA
C.......................................................................
C
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
      SUBROUTINE DFCOMN2(IPSI,C,H1,H2,H3,XK,D,BTA,BT0,IUCV,A2,B2,CHK,
     +                   CKW,BB,BT,CW,EM,CR,VK,NP,ENU,V7,IWWW)
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


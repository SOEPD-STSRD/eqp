      PROGRAM PSLOC
C     VERSIONS 2   PHIVOLCS PROGRAM FOR HYPOCENTRAL LOCATION
C     PROGRAM CALCULATES HYPOCENTRAL PARAMETERS , ADAPTED FROM THE ORIGINAL
C     PROGRAM OF ROBERT Mc COFFREY. MODIFIED TO USE S-PHASE ARRIVAL TIMES BY
C     BY BCB (1990)
C     MODIFIED BY DSS (0696) TO BE USER AND PROGRAMMER FRIENDLY.
      PARAMETER (NSTMAX=150,NWSMAX=200,NPH=2)
      DIMENSION TTMTAB(120,15,2),DEP(15,2),PVSURF(2),DEP1ST(2),DDELT(2),
     1          NUMDTA(2),NUMDEP(2),D(60,4),TRES(60,2),RES(300),B(4,4),
     2          BDIAG(4),Z(4),INDEX(4,3),SS(NSTMAX),TOBS(90,2),WT(90,2),
     3          WTFAC(90,2),PTOA(90),KFM(90,2),DELTA(90),TCALC(90,2),
     4          AZIM(90),ISTA(90),FMS(6),WWSS(NWSMAX),ADLAT(NWSMAX),
     5          ADLON(NWSMAX),SELEV(NWSMAX),SCORR(NWSMAX),HEADER(2),
     6          TTAB(120),TC1(90),TC2(90)
      CHARACTER *4 SS,WWSS,AS,YS
      CHARACTER *12 OUT,TTMTBL,STAT,INPUT
      DOUBLE PRECISION D,Z,RPD,B,RES,WP,WT,WTA,RMS,WS
C    PVSURF =  P WAVE AND S WAVE VELOCITY ON SURFACE
      PII   = 3.1415926
      DPR   = 180./PII
      RPD   = 1D0/DPR
      ITRUE  = 1
      IFALSE = 0
C==== READ CONSTRAINT DATA
C     WRITE(*,'(''CUT OFF FOR ZERO WEIGHTING OF RESIDUAL
C    1 (SECONDS)'')')
      R=20
      IF(R.GT.0.) THEN
         RESCUT = R
           ELSE
         RESCUT = 60.
      END IF
C     WRITE(*,'(''ENTER THETA (LEAST SQUARE DAMPER)'')')
      THETA=0.1
C     WRITE(*,'(''ENTER MINIMUM ADJUSMENT TO HYPOCENTER'')')
      ADJHYP=1
C     WRITE(*,'(''ENTER NUMBER OF ITERATIONS FOR EVENT'')')
      ITERAN = 50
C
C     WRITE(*,'(''DEPTH FIXED ? (Y/N)'')')
      YS='N'
      NPARS = 4
C     Count Total number of Stations
      NUMSTA = 0
      STAT='PHILSTA.DAT'
      OPEN(2,FILE=STAT)
       DO
         READ (2,*, END=2)
         NUMSTA = NUMSTA + 1
       END DO
2     CLOSE (2)
  
C     WRITE(*,'(''ENTER P WAVE VELOCITY ON SURFACE '')')
      PVSURF(1)=5
      PVSURF(2) = PVSURF(1)/1.75
C==== BEGIN READ TRAVEL TIME TABLE
C     WRITE(*,'(''ENTER FILE OF TRAVEL TIME TABLE'')')
      TTMTBL='JBTAB.PS'
      OPEN(1,FILE=TTMTBL)
C      READ IN FILE OF STATIOMS AND TABLES
C      NUMDEP = NO. OF DEPTH   ,  DEP1ST = FIRST DEPTH
C      DDELT = DELTA INCREMENT ,  NUMDTA = NO. OF DELTAS
C      DEP() = DEPTH
C      TTMTAB(I,J,K) , I=DELTA ,J=DEPTH ,K=PHASE
C
C      WRITE(*,'(''READING TABLES'')')
       NPA = 0
100    NPA = NPA +1
       READ(1,'(A2)') HEADER(NPA)
C      WRITE(*,'(A2)') HEADER(NPA)
       READ(1,555)NUMDEP(NPA),DEP1ST(NPA),DDELT(NPA),NUMDTA(NPA)
555    FORMAT(I5,2F10.4,I5)
       DO 35 J=1,NUMDEP(NPA)
          READ(1,'(F8.2)')DEP(J,NPA)
          READ(1,'(9F8.2)')(TTAB(IM),IM=1,NUMDTA(NPA))
          DO 40 LL=1,NUMDTA(NPA)
              TTMTAB(LL,J,NPA) = TTAB(LL)
40        CONTINUE
35     CONTINUE
       IF(NPA.LT.NPH) GO TO 100
C==== END READ TRAVEL TIME TABLE
C==== BEGIN READ STATION LOCATION FROM PHILSTA.DAT
      STAT='PHILSTA.DAT'
      OPEN(2,FILE=STAT)
      READ(2,10)(WWSS(I),ADLAT(I),ADLON(I),SELEV(I),SCORR(I),I=1,NUMSTA)
10      FORMAT(A4,2F10.5,2F6.3)
! C     CHANGE SELEV(I) FROM METER TO KILOMETER
! 1     DO 11 I=1,NUMSTA
! 11      SELEV(I) = SELEV(I)/1000.
      CLOSE(2)
C==== END READ STATION LOCATION FROM PHILSTA.DAT
C     WRITE(*,'('' ENTER FILE OF ARRIVAL DATA '')')
      INPUT='PLOT.DAT'
      OPEN(3,FILE=INPUT)
C     WRITE(*,'(''ENTER FILE OF OUTPUT'')')
      OUT='PLOT.OUT'
      OPEN(4,FILE=OUT)
      IEQ = 0
15    IEQ = IEQ + 1
C==== BEGIN TO READ TRIAL MODEL AND ARRIVAL DATA
       READ(3,112,END=850)IY,MON,IDAY,IH,MINT,SECD,OLAT,OLON,DEPP
112    FORMAT(5I2,F5.1,2F7.3,F5.1)
       WRITE(4,'(3I2)')IY,MON,IDAY
       NSTA = 0
       HREF  = FLOAT(IH) * 3600.
       OT = HREF  + FLOAT(MINT) * 60. + SECD
       DO 105 KK = 1,NUMSTA
         READ(3,104)AS,IWP,IMINP,TSECP,IWS,IMINS,TSECS
104        FORMAT(A4,I2,I3,F5.2,I2,I3,F5.2)
         IF(AS.EQ.'XXX') GO TO 106
         SS(KK) = AS
	   DO  111 N = 1, NUMSTA
           IF(AS.EQ.WWSS(N)) GO TO 301
111      CONTINUE

         N = 0
301      ISTA(KK) = N
C
         WP = 1D0 - FLOAT(IWP)/4D0
         WS = 1D0 - FLOAT(IWS)/4D0
         WT(KK,1) = WP
         WT(KK,2) = WS
         TP = FLOAT(IMINP)*60. + TSECP
         IF(TP.EQ.0.) THEN
              TOBS(KK,1) = TP
            ELSE
              TOBS(KK,1) = TP + HREF
         END IF
         TS = FLOAT(IMINS)*60. + TSECS
         IF(TS.EQ.0.) THEN
              TOBS(KK,2) = TS
            ELSE
              TOBS(KK,2) = TS + HREF
         END IF
105    CONTINUE
       NSTA = NUMSTA
       GO TO 103
106    NSTA = KK -1
C====  END OF READING DATA
      IF(NSTA.EQ.0) GO TO 15
C     WRITE HEADER
C 103 WRITE(4,'(I3)') IEQ
103   WRITE(4,'(''IT    LAT     LON      DEPTH HR MIN SEC    RMS '')')
C     WRITE(*,*)IEQ
      ITSTOP = IFALSE
C===== BEGIN LOOP TROUGH ITERATION
      DO 820 ITER=1,ITERAN
C        ......................FLAG LAST ITERATION
         IF(ITER.EQ.ITERAN.OR.ITSTOP.NE.0) THEN
              ITLAST = ITRUE
            ELSE
              ITLAST = IFALSE
         END IF
         SUMWT = 0.
         RMS   = 0D0
C---------------------------------------------------------
C        PRINT CURRENT SOLUTION
C
C        CONVERT OT TO HR ,MIN,SEC AND PRINT SOLUTION
C---------------------------------------------------------
         IAH = INT(OT/3600.)
         AH  = FLOAT(IAH)
         TTIME = OT -AH*3600.
         IAM = INT(TTIME/60.)
         AM  = FLOAT(IAM)
         ASEC = TTIME - AM*60.
C        FIX DEPTH FOR THE FIRST ITERATION
         IF(ITER.EQ.1) THEN
              NPAR = 3
            ELSE
              NPAR = NPARS
         END IF
         Z(4) = 0.
         NOBS = 0
C=====   BEGIN LOOP THROUGH STATIONS
         DO 680 KSTA=1,NSTA
C           IST = WHERE THE STATION IS IN THE LIST
            IST = ISTA(KSTA)
            IF(IST.EQ.0) GO TO 680
C               STATION NOT IN THE LIST SO SKIP
C               GET DELTA AND AZIMUTH
            CALL DELAZ(ADLAT(IST),ADLON(IST),OLAT,OLON,DELTB,DAZ)
            DELTA(KSTA) = DELTB
            AZIM(KSTA) = DAZ
            PTOA(KSTA) = 0.
C =====     BEGIN LOOP THROUGH PHASES
            DO 650 IPH = 1,NPH
C
C           GET THEORITICAL ARRIVAL TIME AND DERIVATIVE DT/DH , DT/DZ
C
               IDEL = INT((DELTB-DEP1ST(IPH))/DDELT(IPH)) +1
               IF(IDEL.GE.NUMDTA(IPH).OR.IDEL.LT.1) THEN
                  TTIME = 0.
                  GO TO 650
               END IF
               INDT = INT(DELTB/DDELT(IPH))
               DDEL = DELTB/DDELT(IPH) - FLOAT(INDT)
               IDEP = NUMDEP(IPH) - 1
               DO 26 I =1,NUMDEP(IPH)-1
                  IF(DEPP.GE.DEP(I,IPH).AND.DEPP.LT.DEP(I+1,IPH)) IDEP=I
26             CONTINUE
               DDEP = (DEPP-DEP(IDEP,IPH))/(DEP(IDEP+1,IPH)-
     1           DEP(IDEP,IPH))
               T0 = TTMTAB(IDEL,IDEP,IPH)
               T2 = TTMTAB(IDEL+1,IDEP,IPH) - T0
               T3 = TTMTAB(IDEL,IDEP+1,IPH) - T0
               T4 = TTMTAB(IDEL+1,IDEP+1,IPH) - T0
               TTIME =T0 + T2*DDEL + T3*DDEP + ( T4-T3-T2 )*DDEL*DDEP
               DTDH = T2/DDELT(IPH) / 111.19
               DTDZ = T3/(DEP(IDEP+1,IPH) -DEP(IDEP,IPH))
C
C      SKIP IF OUT OF DELTA RANGE
C
C      GET TAKE OFF ANGLE FOR P WAVE
C
               VELP = 1./SQRT(DTDH*DTDH + DTDZ*DTDZ)
               IF(IPH.EQ.1) THEN
                  PTOA(KSTA) = FNASIN(VELP*DTDH)*DPR
                  IF(DTDZ.GT.0.) PTOA(KSTA) = 180.-PTOA(KSTA)
               END IF
C ....................SKIP IF NO ARRIVAL TIME
               IF(TOBS(KSTA,IPH).EQ.0.) THEN
                  WTFAC(KSTA,IPH) = 0.
                  TCALC(KSTA,IPH) = 0.
                  TRES(KSTA,IPH)  = 0.
                  GO TO 650
               END IF
C
C ................CORRECT FOR STATION DELAY SCORR AND STATION ELEVATION SELEV
C
               TTIME = TTIME + SCORR(IST) + SELEV(IST) * SQRT(1./(PVSURF
     1         (IPH)* PVSURF(IPH))-DTDH*DTDH)
               TCALC(KSTA,IPH) = TTIME
C
C ...............TRAVEL TIME RESIDUAL
C
               TRES(KSTA,IPH) = TOBS(KSTA,IPH)- TTIME- OT
               TRESA = TRES(KSTA,IPH)
               ABRES = ABS(TRESA)
C
C.................DETERMINE WEIGHTING FOR THIS ARRIVAL
C.................WEIGHTING FACTOR  TO BE MULTIPLIED BY ASSIGNED
C.................WEIGHT USING COSINE BELL TAPER
C
               WTFAC(KSTA,IPH) = 1.
               IF(ITER.GT.2) THEN
                    RCUT = RESCUT
                  ELSE
                    RCUT = 5.*RESCUT
               END IF
               IF(ABRES.GT.RCUT) WTFAC(KSTA,IPH) = 0.5*(1+COS(PII*(
     1           ABRES-RCUT)/RCUT))
               IF(ABRES.GT.2.*RCUT) WTFAC(KSTA,IPH) =0.
               WTA = WT(KSTA,IPH) * WTFAC(KSTA,IPH)
               IF(WTA.EQ.0D0) GO TO 650
C................SUM OF WEIGHTS AND WEIGHTED RMS
               SUMWT = SUMWT + WTA
               RMS = RMS + TRESA*TRESA*WTA*WTA
               IF(ITLAST.NE.0) GO TO 650
C................NO DERIVATIVES FOR LAST ITERATION
C
C...................GET DERIVATIVES
C                   NUMBER OF OBSERVATION
               NOBS = NOBS + 1
C            WEIGHTED RESIDUAL
               RES(NOBS) = TRES(KSTA,IPH)*WTA
C
C     CALCULATE DERIVATIVES 1=ORIGIN TIME, 2 =LATITUDE ,3 = LONGITUDE
C              4 = DEPTH AND WEIGHT THEM
C
             D(NOBS,1) = 1D0 *WTA
             D(NOBS,2) = -COS(DAZ*RPD)*DTDH*WTA
             D(NOBS,3) = -SIN(DAZ*RPD)*DTDH*WTA
             D(NOBS,4) = DTDZ*WTA
650       CONTINUE
C =====   END LOOP THROUGH PHASES
680    CONTINUE
C===== END LOOP THROUGH STATIONS
       IF(SUMWT.EQ.0.) THEN
            WRITE(4,'(''SUM OF WEIGHT = 0 '')')
C           WRITE(*,'(''SUM OF WEIGHT = 0 '')')
          ELSE
            DOF = SUMWT - NPAR
            IF(DOF.LE.0.) DOF = 1.
            RMS = SQRT(RMS/DOF)
C           WRITE(*,'(F7.2)') RMS
            IF(ITLAST.NE.0) THEN
               WRITE(4,500) ITER,OLAT,OLON,DEPP,IAH,IAM,ASEC,RMS
500            FORMAT(I3,2F9.3,F7.2,2I3,F6.2,F7.2)
C
C        WRITE ERROR ELLEPTIC
C
               WRITE(4,'('' '')')
               WRITE(4,'(''ERR   DEG     DEG      KM           SEC '')')
               ABADA = SQRT(BDIAG(2))*RMS/111.2
               ABADE = SQRT(BDIAG(3))*RMS/111.2
               ABADI = SQRT(BDIAG(4))*RMS
               ABADO = SQRT(BDIAG(1))*RMS
               WRITE(4,690)ABADA,ABADE,ABADI,ABADO
690            FORMAT(6X,F6.3,3X,F6.3,2X,F5.2,4X,F8.2)
C              WRITE(4,'(2F9.3,F7.2,6X,F6.2)')ABADA,ABADE,ABADI,ABADO
               GO TO 825
            END IF
C
C       SOLVE FOR NEW ADJUSMENT
C       BUILD MATRIX B = D (TRANSPOSE) D
C
            IF(ITER.EQ.1) THEN
               DO 115 IKL=1,4
                  B(IKL,4) = 0.
                  B(4,IKL) = 0.
115            CONTINUE
            END IF
            CALL BUILDB(NPAR,NOBS,D,RES,Z,B,THETA)
C
C       SOLVE FOR ADJUSMENTS Z() TO MODEL BY GAUSS-JORDAN INVERSION
C
            CALL GJINV(NPAR,B,Z,DETERM,IPMI,ISING)
            DO 101 I=1,NPAR
101             BDIAG(I) = B(I,I)
            IF(ISING.NE.0) THEN
               ITSTOP = ITRUE
               WRITE(4,'('' '')')
               GO TO 820
            END IF
C
C        PRINT ADJUSMENTS
C
C           WRITE(*,'(4F7.2)')Z(1),Z(2),Z(3),Z(4)
C
C           WRITE(4,500) ITER,OLAT,OLON,DEPP,IAH,IAM,ASEC,RMS
C 500         FORMAT(I3,2F9.3,F7.2,2I3,F6.2,F7.2)
C
C        CHECK SIZE OF ADJUSMENT AND SET FLAG IT BELOW
C                        CRITICAL SIZE
C
            ADJ = SQRT(Z(2)*Z(2) + Z(3)*Z(3) + Z(4)*Z(4))
            IF(ITER.GT.2) THEN
               IF(ADJ.LT.ADJHYP) ITSTOP = ITRUE
               GO TO 499
            END IF
C
C       ADJUSMENT TO HYPOCENTER  DETERMINATION
C
499         OT = OT + Z(1)
            OLON = OLON + Z(3)/111.19
            OLAT = OLAT + Z(2)/111.19
            IF(NPAR.EQ.4) DEPP = DEPP +Z(4)
            IF(DEPP.LT.0.) DEPP = 1.
            END IF
820    CONTINUE
C
C                   PRINT RESULTS
C
C           PRINT RESULTS AND SET UP STATE
C
825     WRITE(4,'('' '')')
        WRITE(4,'(''STA  TOA  AZIM  DELTA RESp  RESs  OBSp  OBSs  CALCp
     1CALCs'')')
        DO 66 I=1,NSTA
           IF(ISTA(I).EQ.0) THEN
                WRITE(4,'(''STATION NOT IN LIST'')')
              ELSE
           IF(TOBS(I,1).EQ.0.) THEN
                TC1(I) = 0.0
              ELSE
                TC1(I) = TOBS(I,1) - OT
           ENDIF
           IF(TOBS(I,2).EQ.0.) THEN
                TC2(I) = 0.0
              ELSE
                TC2(I) = TOBS(I,2) - OT
           END IF
           WT(I,1) = WT(I,1) * WTFAC(I,1)
           WT(I,2) = WT(I,2) * WTFAC(I,2)
           WRITE(4,77)SS(I),PTOA(I),AZIM(I),DELTA(I),TRES(I,1),
     1        TRES(I,2),TC1(I),TC2(I),TCALC(I,1),TCALC(I,2)
C          WRITE(4,77)SS(I),PTOA(I),AZIM(I),DELTA(I),TC1(I),TCALC(I,1),
C    1        TRES(I,1),WT(I,1),TC2(I),TCALC(I,2),TRES(I,2),WT(I,2)
C
77         FORMAT(A3,2F6.1,F6.2,8F6.2)
           END IF
66      CONTINUE
C
840    GO TO 15
850    WRITE(*,'('' '')')
       CLOSE(4)
       STOP
       END
C
C =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-
C         SUBROUTINES
C =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-
C
      FUNCTION FNASIN(X)
C---------------------------------------------------------------
      IF(ABS(X).LT.1.) THEN
         FNASIN = ATAN(X/SQRT(1-X*X))
          GO TO 1111
           ELSE
          FNASIN = 1.570796
       END IF
1111   RETURN
      END
C -------------------------------------------------------------
C
      SUBROUTINE BUILDB(NPAR,NINV,D,R,Z,B,THETA)
C--------------------------------------------------------------
C     ROUTINE TO BUILD B = D TRANSPOSE D AND Z = D TRANSPOSE R
C     ONLY LOWER TRIANGLE OF B STORED
C
      DIMENSION D(60,4),B(4,4),Z(4),R(60)
      DOUBLE PRECISION D,Z,B,R
C
      DO 1 K=1,NPAR
         Z(K) = 0D0
         DO 2 L=1,K
            B(K,L) = 0D0
            DO 3 I =1,NINV
3              B(K,L) = B(K,L) + D(I,L)*D(I,K)
2        CONTINUE
         DO 4 I=1,NINV
4           Z(K) = Z(K) + R(I)*D(I,K)
1     CONTINUE
      DO 5 I=1,NPAR
5        B(I,I) = B(I,I) + THETA*THETA
      RETURN
      END
C--------------------------------------------------------------
      SUBROUTINE GJINV(NCOL,B,Z,DETERM,IPMI,ISING)
C --------------------------------------------------------------
C     ROUTINE FOR GAUSS-JORDAN MATRIX INVERSION
C     SOLVES BX = Z ,OUTPUT IS B INV IN B, AND X IN Z
C     IPMI PRINT MATRIX INVERSE IF TRUE
C     ISING BECOMES TRUE (-1) FOR SINGULAR MATRIX
C
      DIMENSION B(4,4),Z(4),INDEX(4,3)
      DOUBLE PRECISION Z,B
      ISING = 0
C
C     MAKE B THE FULL MATRIX
C
      DO 6 K=1,NCOL
         DO 6 L=1,K
            BB=B(K,L)
         B(L,K) = BB
6     CONTINUE
      DO 7 I=1,NCOL
7        INDEX(I,3) = 0
C     DETERM = 1
C     FIND THE LARGEST ELEMENT
      DO 8 I=1,NCOL
         B1=0.
         DO 9 J=1,NCOL
            IF(INDEX(J,3).EQ.1) GO TO 9
            DO 11 IT=1,NCOL
               IF(INDEX(IT,3).GT.1) GO TO 20
                 IF(INDEX(IT,3).EQ.1)GO TO 11
                   IF(B1.GE.ABS(B(J,IT))) GO TO 11
                 I3 = J
               I4 = IT
              B1 = ABS(B(J,IT))
 11          CONTINUE
  9      CONTINUE
         I43 = INDEX(I4,3)
         INDEX(I4,3) = I43 + 1
         INDEX(I,1) = I3
         INDEX(I,2) = I4
C
C     INTERCHANGE ROWS TO PUT PIVOT O DIAGONAL
C
         IF(I3.EQ.I4) GO TO 40
C
C     DETERM = - DETERM
C
         DO 13 L=1,NCOL
            T1 = B(I3,L)
              T2 = B(I4,L)
            B(I3,L) = T2
13       B(I4,L) = T1
         T1 = Z(I3)
         T2 = Z(I4)
         Z(I3) = T2
         Z(I4) = T1
C
C     DIVIDE PIVOT ROW BY PIVOT ELEMENT
C
40     PIVOT = B(I4,I4)
C
C      DETERM = DETERM*PIVOT
C
         B(I4,I4) = 1.
         DO 55 L=1,NCOL
            BI4L = B(I4,L)
55       B(I4,L) = BI4L / PIVOT
         ZI4 = Z(I4)
         Z(I4) = ZI4 / PIVOT
C
C     REDUCE NONPIVOT ROWS
C
         DO 16 L1=1,NCOL
            IF(L1.EQ.I4) GO TO 16
             T = B(L1,I4)
              B(L1,I4) = 0D0
               DO 17 L=1,NCOL
                  BI4L = B(I4,L)
                 BL1L = B(L1,L)
17             B(L1,L) = BL1L - BI4L * T
             ZI4 = Z(I4)
           ZL1 = Z(L1)
         Z(L1) = ZL1 - ZI4 * T
16      CONTINUE
8     CONTINUE
C
C     INTERCHANGE COLUMNS
C
      DO 18 I=1,NCOL
         L = NCOL-I+1
           IF(INDEX(L,1).EQ.INDEX(L,2)) GO TO 18
            I3 = INDEX(L,1)
           I4 = INDEX(L,2)
         DO 19 IT=1,NCOL
           T1 = B(IT,I3)
             T2 = B(IT,I4)
           B(IT,I3) = T2
19       B(IT,I4) = T1
18    CONTINUE
      DO 107 IT=1,NCOL
         IF(INDEX(IT,3).NE.1) GO TO 20
107   CONTINUE
       GO TO 29
20     ISING = -1
       WRITE(*,'('' SINGULAR MATRIX '')')
29      RETURN
      END
C---------------------------------------------------------------
      SUBROUTINE DELAZ(SLAT,SLON,EPLAT,EPLON,DELTA,AZIM)
C-----------------------------------------------------------------
C     ROUTINE TO DETERMINE DISTANT BETWEEN STATION AND EPICENTER
C     AND DETERMINE AZIMUTH OF STATION
C     SLAT & SLON - STATION LAT AND LONG
C     EPLAT & EPLON - EPICENTER LAT AND LONG
C
      F = 1.745329E-02
      TSLAT = SLAT * F
      TSLON = SLON * F
       GLAT = ATAN(0.993277* TAN(TSLAT))
          DCA = COS(GLAT)*COS(TSLON)
          DCB = COS(GLAT)*SIN(TSLON)
          DCC = SIN(GLAT)
          SLAT = TSLAT/F
       GEPLAT = ATAN(0.993277* TAN(EPLAT*F))
              EPDCA = COS(GEPLAT) * COS(EPLON*F)
              EPDCB = COS(GEPLAT) * SIN(EPLON*F)
              EPDCC = SIN(GEPLAT)
C
C      CALCULATE AZIMUTH AND DISTANCE FROM EQ TO STATION
C
      PAR1=(1.-0.5*((DCA-EPDCA)**2+(DCB-EPDCB)**2+(DCC-EPDCC)**2))
      DELTA = 57.29578 * ATAN(SQRT((1./(PAR1*PAR1))-1.))
      IF(PAR1.LT.0.) DELTA = 180.-DELTA
        IF(DELTA.EQ.0.) DELTA = 0.001
          T1= DCA - SIN(EPLON*F)
          T2 = DCB + COS(EPLON*F)
          T3 = SIN(DELTA*F)
      PAR2 = ( 0.5 * ( T1*T1 + T2*T2 + DCC*DCC )-1. ) / T3
          T1 = DCA - SIN(GEPLAT) * COS(EPLON*F)
          T2 = DCB - SIN(GEPLAT) * SIN(EPLON*F)
          T4 = DCC + COS(GEPLAT)
          APAR2 = ABS(PAR2)
       IF(APAR2.LT.1.) GO TO 5
          AZIM = 90.
       GO TO 10
5       IF(PAR2.EQ.0.) PAR2 = 0.001
         IF(PAR2.EQ.1) PAR2 = 0.999
          AZIM = 57.29578 * ATAN(1./SQRT((1./PAR2**2)-1.))
10     PAR3 = ( 0.5 * ( T1*T1 + T2*T2 + T4*T4) -1. ) / T3
       IF(PAR2.GE.0.) GO TO 15
         IF(PAR3.LT.0.) GO TO 20
           AZIM = 360. - AZIM
           GO TO 30
20       AZIM = AZIM + 180.
         GO TO 30
15     IF(PAR3.LT.0.) AZIM = 180.- AZIM
30     RETURN
       END
C =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=

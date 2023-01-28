!******************************************************************************
! Calcula el nivel de significacion de la distribucion F de Fisher con N1 y N2
! grados de libertad (area de la cola derecha calculada a partir de la
! abcisa X).
        REAL FUNCTION FFFISHER(N1,N2,X)
        IMPLICIT NONE
        INTEGER N1,N2
        REAL X
        REAL BETAI
!
        FFFISHER=BETAI(REAL(N2)/2.,REAL(N1)/2.,REAL(N2)/(REAL(N2)+REAL(N1)*X))
        END
!******************************************************************************
! Calcula el nivel de significacion de la distribucion t de Student con N
! grados de libertad (area de la cola derecha calculada a partir de la
! abcisa X).
        REAL FUNCTION FTSTUDENT(N,X)
        IMPLICIT NONE
        INTEGER N
        REAL X
        REAL BETAI
!
        FTSTUDENT=0.5*BETAI(REAL(N)/2.,0.5,REAL(N)/(REAL(N)+X*X))
        END
!******************************************************************************
! Calcula la abcisa a cuya derecha la distribucion t de Student con N grados de
! libertad cubre un area (nivel de significacion) alpha
        REAL FUNCTION FTSTUDENTI(N,ALPHA)
        IMPLICIT NONE
        INTEGER N
        REAL ALPHA
        REAL INVIERTE_FTSTUDENT
!
        FTSTUDENTI=INVIERTE_FTSTUDENT(N,ALPHA)
        END
!******************************************************************************
! Calcula el area de la cola derecha de la curva normal tipificada a partir de
! la abcisa X.
        REAL FUNCTION FNORTIP(X)
        IMPLICIT NONE
        REAL X
        REAL GAMMQ
!
        FNORTIP=0.5*GAMMQ(0.5,X*X/2.)
        END
!******************************************************************************
! Calcula la abcisa a cuya derecha la distribucion normal tipificada cubre un 
! area (nivel de significacion) alpha
        REAL FUNCTION FNORTIPI(ALPHA)
        IMPLICIT NONE
        REAL ALPHA
        REAL INVIERTE_FNORTIP
!
        FNORTIPI=INVIERTE_FNORTIP(ALPHA)
        END
!******************************************************************************
! Calcula el nivel de significacion de la distribucion chi-cuadrado con N
! grados de libertad (area de la cola derecha calculada a partir de la 
! abcisa X).
        REAL FUNCTION FCHISQR(N,X)
        IMPLICIT NONE
        INTEGER N
        REAL X
        REAL GAMMQ
!
        FCHISQR=GAMMQ(REAL(N)/2.,X/2.)
        END
!******************************************************************************
! Calcula la abcisa a cuya derecha la distribucion chi-cuadrado con N
! grados de libertad cubre un area (nivel de significacion) alpha
!
! Nota: me temo que esta funcion no funciona demasiado bien cuando N es muy
! pequen~o o muy grande
!
        REAL FUNCTION FCHISQRI(N,ALPHA)
        IMPLICIT NONE
        INTEGER N
        REAL ALPHA
        REAL INVIERTE_FCHISQR
!
        FCHISQRI=INVIERTE_FCHISQR(N,ALPHA)
        END
!
!******************************************************************************
!******************************************************************************
      FUNCTION GAMMP(A,X)
      IF(X.LT.0..OR.A.LE.0.)THEN
        INCLUDE 'deallocate_arrays.inc'
        STOP 'FATAL ERROR: PAUSE'
      END IF
      IF(X.LT.A+1.)THEN
        CALL GSER(GAMSER,A,X,GLN)
        GAMMP=GAMSER
      ELSE
        CALL GCF(GAMMCF,A,X,GLN)
        GAMMP=1.-GAMMCF
      ENDIF
      RETURN
      END
!******************************************************************************
      FUNCTION GAMMQ(A,X)
      IF(X.LT.0..OR.A.LE.0.)THEN
        INCLUDE 'deallocate_arrays.inc'
        STOP 'FATAL ERROR: PAUSE'
      END IF
      IF(X.LT.A+1.)THEN
        CALL GSER(GAMSER,A,X,GLN)
        GAMMQ=1.-GAMSER
      ELSE
        CALL GCF(GAMMCF,A,X,GLN)
        GAMMQ=GAMMCF
      ENDIF
      RETURN
      END
!******************************************************************************
      SUBROUTINE GSER(GAMSER,A,X,GLN)
      PARAMETER (ITMAX=100,EPS=3.E-7)
      GLN=GAMMLN(A)
      IF(X.LE.0.)THEN
        IF(X.LT.0.)THEN
          INCLUDE 'deallocate_arrays.inc'
          STOP 'FATAL ERROR: PAUSE'
        END IF
        GAMSER=0.
        RETURN
      ENDIF
      AP=A
      SUM=1./A
      DEL=SUM
      DO 11 N=1,ITMAX
        AP=AP+1.
        DEL=DEL*X/AP
        SUM=SUM+DEL
        IF(ABS(DEL).LT.ABS(SUM)*EPS)GO TO 1
11    CONTINUE
      INCLUDE 'deallocate_arrays.inc'
      STOP 'FATAL ERROR: A too large, ITMAX too small'
1     GAMSER=SUM*EXP(-X+A*LOG(X)-GLN)
      RETURN
      END
!******************************************************************************
      SUBROUTINE GCF(GAMMCF,A,X,GLN)
      PARAMETER (ITMAX=100,EPS=3.E-7)
      G=0.0 !evita un WARNING de compilación
      GLN=GAMMLN(A)
      GOLD=0.
      A0=1.
      A1=X
      B0=0.
      B1=1.
      FAC=1.
      DO 11 N=1,ITMAX
        AN=FLOAT(N)
        ANA=AN-A
        A0=(A1+A0*ANA)*FAC
        B0=(B1+B0*ANA)*FAC
        ANF=AN*FAC
        A1=X*A0+ANF*A1
        B1=X*B0+ANF*B1
        IF(A1.NE.0.)THEN
          FAC=1./A1
          G=B1*FAC
          IF(ABS((G-GOLD)/G).LT.EPS)GO TO 1
          GOLD=G
        ENDIF
11    CONTINUE
      INCLUDE 'deallocate_arrays.inc'
      STOP 'FATAL ERROR: A too large, ITMAX too small'
1     GAMMCF=EXP(-X+A*ALOG(X)-GLN)*G
      RETURN
      END
!******************************************************************************
      real FUNCTION GAMMLN(XX)
      REAL*8 COF(6),STP,HALF,ONE,FPF,X,TMP,SER
      DATA COF,STP/76.18009173D0,-86.50532033D0,24.01409822D0,-1.231739516D0,.120858003D-2,-.536382D-5,2.50662827465D0/
      DATA HALF,ONE,FPF/0.5D0,1.0D0,5.5D0/
      X=XX-ONE
      TMP=X+FPF
      TMP=(X+HALF)*LOG(TMP)-TMP
      SER=ONE
      DO 11 J=1,6
        X=X+ONE
        SER=SER+COF(J)/X
11    CONTINUE
      GAMMLN=real(TMP+LOG(STP*SER))
      RETURN
      END
!******************************************************************************
      real FUNCTION BETAI(A,B,X)
      real a,b,x,gammln,betacf
      IF(X.LT.0..OR.X.GT.1.)THEN
        INCLUDE 'deallocate_arrays.inc'
        STOP 'FATAL ERROR: bad argument X in BETAI'
      END IF
      IF(X.EQ.0..OR.X.EQ.1.)THEN
        BT=0.
      ELSE
        BT=GAMMLN(A+B)-GAMMLN(A)-GAMMLN(B)+A*ALOG(X)+B*ALOG(1.-X)
        IF(BT.LT.-69.)THEN !evitamos: Floating point exception 5, underflow
          BT=0.
        ELSE
          BT=EXP(BT)
        END IF
      ENDIF
      IF(X.LT.(A+1.)/(A+B+2.))THEN
        BETAI=BT*BETACF(A,B,X)/A
        RETURN
      ELSE
        BETAI=1.-BT*BETACF(B,A,1.-X)/B
        RETURN
      ENDIF
      END
!******************************************************************************
      real FUNCTION BETACF(A,B,X)
      PARAMETER (ITMAX=100,EPS=3.E-7)
      AM=1.
      BM=1.
      AZ=1.
      QAB=A+B
      QAP=A+1.
      QAM=A-1.
      BZ=1.-QAB*X/QAP
      DO 11 M=1,ITMAX
        EM=M
        TEM=EM+EM
        D=EM*(B-M)*X/((QAM+TEM)*(A+TEM))
        AP=AZ+D*AM
        BP=BZ+D*BM
        D=-(A+EM)*(QAB+EM)*X/((A+TEM)*(QAP+TEM))
        APP=AP+D*AZ
        BPP=BP+D*BZ
        AOLD=AZ
        AM=AP/BPP
        BM=BP/BPP
        AZ=APP/BPP
        BZ=1.
        IF(ABS(AZ-AOLD).LT.EPS*ABS(AZ))GO TO 1
11    CONTINUE
      INCLUDE 'deallocate_arrays.inc'
      STOP 'FATAL ERROR: A or B too big, or ITMAX too small'
1     BETACF=AZ
      RETURN
      END
!******************************************************************************
      real FUNCTION ERFfunction(X)
      IF(X.LT.0.)THEN
        ERFfunction=-GAMMP(.5,X**2)
      ELSE
        ERFfunction=GAMMP(.5,X**2)
      ENDIF
      RETURN
      END
!******************************************************************************
!******************************************************************************
        REAL FUNCTION INVIERTE_FTSTUDENT(N,ALPHA)
        IMPLICIT NONE
        INTEGER N
        REAL ALPHA
!
        REAL FUNK_IFTSTUDENT
        EXTERNAL FUNK_IFTSTUDENT
!
        INTEGER NDEG_OF_FREEDOM
        INTEGER NEVAL
        REAL X0,DX0
        REAL XF,DXF
        REAL SIGNIFICACION
!
        COMMON/BLKIFTSTUDENT1/NDEG_OF_FREEDOM
        COMMON/BLKIFTSTUDENT2/SIGNIFICACION
!------------------------------------------------------------------------------
! Segun el valor de alfa, hacemos una estimacion del punto inicial para
! DOWNHILL, usando el valor de la funcion t de student para infinitos grados 
! de libertad (distribucion normal).
        IF(ALPHA.GT.0.50)THEN
          X0=0.0
        ELSEIF(ALPHA.GT.0.40)THEN
          X0=0.253
        ELSEIF(ALPHA.GT.0.30)THEN
          X0=0.524
        ELSEIF(ALPHA.GT.0.20)THEN
          X0=0.842
        ELSEIF(ALPHA.GT.0.10)THEN
          X0=1.282
        ELSEIF(ALPHA.GT.0.050)THEN
          X0=1.645
        ELSEIF(ALPHA.GT.0.025)THEN
          X0=1.960
        ELSEIF(ALPHA.GT.0.010)THEN
          X0=2.326
        ELSEIF(ALPHA.GT.0.005)THEN
          X0=2.576
        ELSEIF(ALPHA.GT.0.001)THEN
          X0=3.090
        ELSE
          X0=3.291
        END IF
        DX0=.1
!------------------------------------------------------------------------------
        NDEG_OF_FREEDOM=N
        SIGNIFICACION=ALPHA
        CALL DOWNHILL(1,X0,DX0,FUNK_IFTSTUDENT,1.0,0.5,2.0,1.E-7,XF,DXF,NEVAL,500)
        INVIERTE_FTSTUDENT=XF
        RETURN
        END
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
        REAL FUNCTION FUNK_IFTSTUDENT(X)
        IMPLICIT NONE
        REAL X
!
        REAL FTSTUDENT
!
        INTEGER N
        REAL ALPHA
        COMMON/BLKIFTSTUDENT1/N
        COMMON/BLKIFTSTUDENT2/ALPHA
!------------------------------------------------------------------------------        
        FUNK_IFTSTUDENT=ABS(ALPHA-FTSTUDENT(N,X))
        END
!******************************************************************************
!******************************************************************************
        REAL FUNCTION INVIERTE_FNORTIP(ALPHA)
        IMPLICIT NONE
        REAL ALPHA
!
        REAL FUNK_IFNORTIP
        EXTERNAL FUNK_IFNORTIP
!
        INTEGER NEVAL
        REAL X0,DX0
        REAL XF,DXF
        REAL SIGNIFICACION
!
        COMMON/BLKIFTNORTIP1/SIGNIFICACION
!------------------------------------------------------------------------------
! Segun el valor de alfa, hacemos una estimacion del punto inicial para
! DOWNHILL
        IF(ALPHA.GT.0.50)THEN
          X0=0.0
        ELSEIF(ALPHA.GT.0.40)THEN
          X0=0.253
        ELSEIF(ALPHA.GT.0.30)THEN
          X0=0.524
        ELSEIF(ALPHA.GT.0.20)THEN
          X0=0.842
        ELSEIF(ALPHA.GT.0.10)THEN
          X0=1.282
        ELSEIF(ALPHA.GT.0.050)THEN
          X0=1.645
        ELSEIF(ALPHA.GT.0.025)THEN
          X0=1.960
        ELSEIF(ALPHA.GT.0.010)THEN
          X0=2.326
        ELSEIF(ALPHA.GT.0.005)THEN
          X0=2.576
        ELSEIF(ALPHA.GT.0.001)THEN
          X0=3.090
        ELSE
          X0=3.291
        END IF
        DX0=.1
!------------------------------------------------------------------------------
        SIGNIFICACION=ALPHA
        CALL DOWNHILL(1,X0,DX0,FUNK_IFNORTIP,1.0,0.5,2.0,1.E-7,XF,DXF,NEVAL,500)
        INVIERTE_FNORTIP=XF
        RETURN
        END
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
        REAL FUNCTION FUNK_IFNORTIP(X)
        IMPLICIT NONE
        REAL X
!
        REAL FNORTIP
!
        REAL ALPHA
        COMMON/BLKIFTNORTIP1/ALPHA
!------------------------------------------------------------------------------        
        FUNK_IFNORTIP=ABS(ALPHA-FNORTIP(X))
        END
!******************************************************************************
!******************************************************************************
        REAL FUNCTION INVIERTE_FCHISQR(N,ALPHA)
        IMPLICIT NONE
        INTEGER N
        REAL ALPHA
!
        REAL FUNK_IFCHISQR
        EXTERNAL FUNK_IFCHISQR
!
        INTEGER NDEG_OF_FREEDOM
        INTEGER NEVAL
        REAL X0,DX0
        REAL XF,DXF
        REAL SIGNIFICACION
!
        COMMON/BLKIFCHISQR1/NDEG_OF_FREEDOM
        COMMON/BLKIFCHISQR2/SIGNIFICACION
!------------------------------------------------------------------------------
! Valores iniciales par DOWNHILL
        X0=REAL(N)
        DX0=REAL(N)/10.
!------------------------------------------------------------------------------
        NDEG_OF_FREEDOM=N
        SIGNIFICACION=ALPHA
        CALL DOWNHILL(1,X0,DX0,FUNK_IFCHISQR,1.0,0.5,2.0,1.E-7,XF,DXF,NEVAL,500)
        INVIERTE_FCHISQR=XF
        RETURN
        END
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
        REAL FUNCTION FUNK_IFCHISQR(X)
        IMPLICIT NONE
        REAL X
!
        REAL FCHISQR
!
        INTEGER N
        REAL ALPHA
        COMMON/BLKIFCHISQR1/N
        COMMON/BLKIFCHISQR2/ALPHA
!------------------------------------------------------------------------------        
        FUNK_IFCHISQR=ABS(ALPHA-FCHISQR(N,X))
        END
!******************************************************************************


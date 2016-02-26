C Calcula, en la posicion dada como parametro de entrada, la distorsion local
C realizando un ajuste mediante polinomios de segundo grado a la distorsion
C en la region del boundary en la cual se encuentra dicho punto.
C El retorno de la subrutina son dos polinomios, COEFFBL00 y COEFFBA00, de
C grado NDEGBL00 y NDEGBA00 respectivamente.
C NCOL1: si es distinto de -1, se dibuja el polinomios definido por COEFFBL00
C NCOL2: si es distinto de -1, se dibuja el polinomios definido por COEFFBA00
        SUBROUTINE WHEREAMI(XC,YC,COEFFBL00,COEFFBA00,
     +   NCOL1,NCOL2,LDEBUGLOCAL)
        IMPLICIT NONE
        REAL XC,YC
        REAL COEFFBL00(20),COEFFBA00(20)   !coef. polinomio que pasa por XC, YC
        INTEGER NCOL1,NCOL2
        LOGICAL LDEBUGLOCAL
C
        INCLUDE 'largest.inc'
C
        REAL THRESHOLD
        PARAMETER (THRESHOLD=1.E-5)
        REAL THRESHOLDBIS
        PARAMETER (THRESHOLDBIS=1.E-3)
        INTEGER NITERMAX
        PARAMETER (NITERMAX=100)
C
        REAL FPOLY
        REAL ARCLENGTH
        LOGICAL LROMBO
C
        INTEGER I,L,LMIN,K,L1,L2,LL1,LL2
        INTEGER NDEGBL(NXYMAX),NDEGBA(NXYMAX)  !polynomial degrees for boundary
        INTEGER NDEGBL00,NDEGBA00
        INTEGER NLINBL,NLINBA                  !number of arc lines in boundary
        INTEGER NDEGTRAS                        !grado del polinomio trasladado
        INTEGER NITER  !num. de iteraciones con los polinomios de segundo grado
        REAL COEFFBL(20,NXYMAX),COEFFBA(20,NXYMAX)      !pol. coef. in boundary
        REAL COEFFP1(3),COEFFP2(3) !coeficientes de polinomios de segundo grado
        REAL COEFFTRAS(20)                 !coeficiente de polinomio trasladado
        REAL XMINBL(NXYMAX),XMAXBL(NXYMAX)
        REAL YMINBL(NXYMAX),YMAXBL(NXYMAX)
        REAL XMINBA(NXYMAX),XMAXBA(NXYMAX)
        REAL YMINBA(NXYMAX),YMAXBA(NXYMAX)
        REAL X0,Y0
        REAL XN,YN,XM,YM,XQ,YQ,XR,YR
        REAL XN_,YN_,XM_,YM_,XQ_,YQ_,XR_,YR_
        REAL S0,S1,T0,T1                                   !arcos de polinomios
        REAL S0_,S1_,T0_,T1_                               !arcos de polinomios
        REAL N0,N1,M0,M1
        REAL DIST,DISTMIN
        REAL XP(NXYMAX),YP(NXYMAX)
        REAL XF(NXYMAX),YF(NXYMAX)
        REAL XDMIN,YDMIN
        REAL XMIN,XMAX,YMIN,YMAX
        REAL XMIN1,XMAX1,YMIN1,YMAX1
        REAL XMIN2,XMAX2,YMIN2,YMAX2
        REAL VARDUM(20),CHISQR,SR2,COVAR(20,20)
        REAL XINI,XFIN,YINI,YFIN
        LOGICAL LUP,LRIGHT
        LOGICAL LOOP
        LOGICAL IFCOEFF(20)
C
        COMMON/BLKBOUND1/NLINBL,NDEGBL,NDEGBL00
        COMMON/BLKBOUND1B/NLINBA,NDEGBA,NDEGBA00
        COMMON/BLKBOUND2/COEFFBL
        COMMON/BLKBOUND2B/COEFFBA
        COMMON/BLKBOUND4/XMINBL,YMINBL,XMAXBL,YMAXBL
        COMMON/BLKBOUND5/XMINBA,YMINBA,XMAXBA,YMAXBA
C------------------------------------------------------------------------------
C------------------------------------------------------------------------------
        LMIN=0 !evita WARNING de compilacion
        LUP=.FALSE. !evita WARNING de compilacion
        LRIGHT=.FALSE. !evita WARNING de compilacion
C La primera parte del problema es averiguar cual es la region calibrada mas
C proxima.
C------------------------------------------------------------------------------
C------------------------------------------------------------------------------
C Buscamos el punto mas cercano en el grid de lineas. Empezamos buscando
C entre que dos lineas BA se encuentra el punto
        L=0
        DISTMIN=10.*REAL(NXYMAX*NXYMAX)
        DO L=1,NLINBA
c borde inferior
          DIST=(XC-XMINBA(L))*(XC-XMINBA(L))+
     +         (YC-YMINBA(L))*(YC-YMINBA(L))
          IF(DIST.LT.DISTMIN)THEN
            LMIN=L
            DISTMIN=DIST
            LUP=.FALSE.
          END IF
c borde superior
          DIST=(XC-XMAXBA(L))*(XC-XMAXBA(L))+
     +     (YC-YMAXBA(L))*(YC-YMAXBA(L))
          IF(DIST.LT.DISTMIN)THEN
            LMIN=L
            DISTMIN=DIST
            LUP=.TRUE.
          END IF
        END DO
C si no hemos encontrado nada, es que hay algun error de bulto
        IF(LMIN.EQ.0)THEN
          WRITE(*,101) '***FATAL ERROR***'
          WRITE(*,101) '=> LMIN=0 in WHEREAMI.'
          STOP
        END IF
C si estamos en un borde, no hay mucho que decidir; el problema es averiguar
C que ocurre en los casos intermedios
        L1=0
        L2=0
        IF(LMIN.EQ.1)THEN !borde izquierdo
          L1=1
          L2=2
        ELSEIF(LMIN.EQ.NLINBA)THEN !borde derecho
          L1=NLINBA-1
          L2=NLINBA
        ELSE !vamos con los casos intermedios
          L1=LMIN-1  !rombo definido con la linea anterior
          L2=LMIN
          IF(.NOT.LROMBO(XC,YC,
     +     XMINBA(L1),YMINBA(L1),
     +     XMINBA(L2),YMINBA(L2),
     +     XMAXBA(L2),YMAXBA(L2),
     +     XMAXBA(L1),YMAXBA(L1)))THEN
            L1=LMIN !rombo definido con la linea posterior
            L2=LMIN+1
            IF(.NOT.LROMBO(XC,YC,
     +       XMINBA(L1),YMINBA(L1),
     +       XMINBA(L2),YMINBA(L2),
     +       XMAXBA(L2),YMAXBA(L2),
     +       XMAXBA(L1),YMAXBA(L1)))THEN
              L1=0 !no ha habido suerte; todavia no sabemos cual es el rombo
              L2=0
            END IF
          END IF
        END IF
C Si no hemos tenido suerte, decidimos tomar el rombo mas proximo. Para ello
C calculamos el punto del polinomio mas cercano a XC,YC.
        IF(L1*L2.EQ.0)THEN
          IF(LUP)THEN                                           !borde superior
            CALL DISTPP(XC,YC,NDEGBL(NLINBL),COEFFBL(1,NLINBL),
     +       XC,DIST,XDMIN,YDMIN)
            IF(XDMIN.LE.XMAXBA(LMIN))THEN
              L1=LMIN-1
              L2=LMIN
            ELSE
              L1=LMIN
              L2=LMIN+1
            END IF
          ELSE                                                  !borde inferior
            CALL DISTPP(XC,YC,NDEGBL(1),COEFFBL(1,1),
     +       XC,DIST,XDMIN,YDMIN)
            IF(XDMIN.LE.XMINBA(LMIN))THEN
              L1=LMIN-1
              L2=LMIN
            ELSE
              L1=LMIN
              L2=LMIN+1
            END IF
          END IF
        END IF
C bueno, ya tenemos un rombo grande; vamos a dibujarlo
        IF(LDEBUGLOCAL)THEN
          CALL PGSCI(7)
          CALL PGMOVE(XMINBA(L1),YMINBA(L1))
          CALL PGDRAW(XMAXBA(L1),YMAXBA(L1))
          CALL PGDRAW(XMAXBA(L2),YMAXBA(L2))
          CALL PGDRAW(XMINBA(L2),YMINBA(L2))
          CALL PGDRAW(XMINBA(L1),YMINBA(L1))
          CALL PGPOINT(1,XMINBA(LMIN),YMINBA(LMIN),23)
          CALL PGPOINT(1,XMAXBA(LMIN),YMAXBA(LMIN),23)
          CALL PGSCI(1)
        END IF
C------------------------------------------------------------------------------
C Ahora nos falta encontrar en que subrombo se encuentra (es decir, entre que
C dos lineas BL)
        L=0
        DISTMIN=10.*REAL(NXYMAX*NXYMAX)
        DO L=1,NLINBL
c borde izquierdo
          CALL INTERSEC(NDEGBA(L1),COEFFBA(1,L1),NDEGBL(L),COEFFBL(1,L),
     +     XC,XMIN,YMIN)
          DIST=(XC-XMIN)*(XC-XMIN)+(YC-YMIN)*(YC-YMIN)
          IF(DIST.LT.DISTMIN)THEN
            LMIN=L
            DISTMIN=DIST
            LRIGHT=.FALSE.
          END IF
c borde derecho
          CALL INTERSEC(NDEGBA(L2),COEFFBA(1,L2),NDEGBL(L),COEFFBL(1,L),
     +     XC,XMAX,YMAX)
          DIST=(XC-XMAX)*(XC-XMAX)+(YC-YMAX)*(YC-YMAX)
          IF(DIST.LT.DISTMIN)THEN
            LMIN=L
            DISTMIN=DIST
            LRIGHT=.TRUE.
          END IF
        END DO
C si no hemos encontrado nada, es que hay algun error de bulto
        IF(LMIN.EQ.0)THEN
          WRITE(*,101) '***FATAL ERROR***'
          WRITE(*,101) '=> LMIN=0 in WHEREAMI.'
          STOP
        END IF
C si estamos en un borde, no hay mucho que decidir; el problema es averiguar
C que ocurre en los casos intermedios (en este caso, tiene que estar en alguno
C de los dos rombos adyacentes por narices; si los polinomios tienen mucha
C curvatura, es posible que asignemos un punto al rombo equivocado, aunque
C solo va a ocurrir esto muy cerca de los bordes, con lo que la posible
C extrapolacion tiene que ser buena)
        LL1=0
        LL2=0
        IF(LMIN.EQ.1)THEN !borde inferior
          LL1=1
          LL2=2
        ELSEIF(LMIN.EQ.NLINBL)THEN !borde superior
          LL1=NLINBL-1
          LL2=NLINBL
        ELSE !vamos con los casos intermedios
          LL1=LMIN-1  !rombo definido con la linea inferior
          LL2=LMIN
          CALL INTERSEC(NDEGBA(L1),COEFFBA(1,L1),
     +     NDEGBL(LL1),COEFFBL(1,LL1),XC,XMIN1,YMIN1)
          CALL INTERSEC(NDEGBA(L2),COEFFBA(1,L2),
     +     NDEGBL(LL1),COEFFBL(1,LL1),XC,XMIN2,YMIN2)
          CALL INTERSEC(NDEGBA(L2),COEFFBA(1,L2),
     +     NDEGBL(LL2),COEFFBL(1,LL2),XC,XMAX2,YMAX2)
          CALL INTERSEC(NDEGBA(L1),COEFFBA(1,L1),
     +     NDEGBL(LL2),COEFFBL(1,LL2),XC,XMAX1,YMAX1)
          IF(.NOT.LROMBO(XC,YC,
     +     XMIN1,YMIN1,
     +     XMIN2,YMIN2,
     +     XMAX2,YMAX2,
     +     XMAX1,YMAX1))THEN
            LL1=LMIN  !rombo definido con la linea superior
            LL2=LMIN+1
            CALL INTERSEC(NDEGBA(L1),COEFFBA(1,L1),
     +       NDEGBL(LL1),COEFFBL(1,LL1),XC,XMIN1,YMIN1)
            CALL INTERSEC(NDEGBA(L2),COEFFBA(1,L2),
     +       NDEGBL(LL1),COEFFBL(1,LL1),XC,XMIN2,YMIN2)
            CALL INTERSEC(NDEGBA(L2),COEFFBA(1,L2),
     +       NDEGBL(LL2),COEFFBL(1,LL2),XC,XMAX2,YMAX2)
            CALL INTERSEC(NDEGBA(L1),COEFFBA(1,L1),
     +       NDEGBL(LL2),COEFFBL(1,LL2),XC,XMAX1,YMAX1)
            IF(.NOT.LROMBO(XC,YC,
     +       XMIN1,YMIN1,
     +       XMIN2,YMIN2,
     +       XMAX2,YMAX2,
     +       XMAX1,YMAX1))THEN
              LL1=0 !no ha habido suerte; todavia no sabemos cual es el rombo
              LL2=0
            END IF
          END IF
        END IF
C Si no hemos tenido suerte, decidimos tomar el subrombo mas proximo. Para ello
C calculamos el punto del polinomio mas cercano a XC,YC.
        IF(LL1*LL2.EQ.0)THEN
          IF(LRIGHT)THEN                                         !borde derecho
            CALL INTERSEC(NDEGBA(NLINBA),COEFFBA(1,NLINBA),
     +       NDEGBL(LMIN),COEFFBL(1,LMIN),XC,XMAX,YMAX)
            CALL DISTPP(YC,XC,NDEGBA(NLINBA),COEFFBA(1,NLINBA),
     +       YC,DIST,YDMIN,XDMIN)
            IF(YDMIN.LE.YMAX)THEN
              LL1=LMIN-1
              LL2=LMIN
            ELSE
              LL1=LMIN
              LL2=LMIN+1
            END IF
          ELSE                                                 !borde izquierdo
            CALL INTERSEC(NDEGBA(1),COEFFBA(1,1),
     +       NDEGBL(LMIN),COEFFBL(1,LMIN),XC,XMIN,YMIN)
            CALL DISTPP(YC,XC,NDEGBA(1),COEFFBA(1,1),
     +       YC,DIST,YDMIN,XDMIN)
            IF(YDMIN.LE.YMIN)THEN
              LL1=LMIN-1
              LL2=LMIN
            ELSE
              LL1=LMIN
              LL2=LMIN+1
            END IF
          END IF
        END IF
C calculamos intersecciones definitivas
        CALL INTERSEC(NDEGBA(L1),COEFFBA(1,L1),
     +   NDEGBL(LL1),COEFFBL(1,LL1),XC,XMIN1,YMIN1)
        CALL INTERSEC(NDEGBA(L2),COEFFBA(1,L2),
     +   NDEGBL(LL1),COEFFBL(1,LL1),XC,XMIN2,YMIN2)
        CALL INTERSEC(NDEGBA(L2),COEFFBA(1,L2),
     +   NDEGBL(LL2),COEFFBL(1,LL2),XC,XMAX2,YMAX2)
        CALL INTERSEC(NDEGBA(L1),COEFFBA(1,L1),
     +   NDEGBL(LL2),COEFFBL(1,LL2),XC,XMAX1,YMAX1)
C bueno, ya tenemos el rombo pequen~o; vamos a dibujarlo
        IF(LDEBUGLOCAL)THEN
          CALL PGSCI(6)
          CALL PGMOVE(XMIN1,YMIN1)
          CALL PGDRAW(XMAX1,YMAX1)
          CALL PGDRAW(XMAX2,YMAX2)
          CALL PGDRAW(XMIN2,YMIN2)
          CALL PGDRAW(XMIN1,YMIN1)
          CALL PGSCI(1)
        END IF
C------------------------------------------------------------------------------
C------------------------------------------------------------------------------
C Ahora ya podemos calcular la distorsion en el punto seleccionado. Para ello
C vamos a calcular los polinomios con la distorsion correspondiente que pasan
C por dicho punto. Esos polinomios se calculan de forma iterativa, usando
C los polinomios del rombo seleccionado anteriormente.
C------------------------------------------------------------------------------
C------------------------------------------------------------------------------
C polinomio que pasa por el punto XC,YC: tomamos el polinomio que pasa por
C la esquina mas cercana del rombo seleccionado
        X0=XC !evaluamos el polinomio en XC
        IF(LMIN.EQ.LL1)THEN                   !el borde inferior esta mas cerca
          Y0=FPOLY(NDEGBL(LL1),COEFFBL(1,LL1),X0)
        ELSE                                  !el borde superior esta mas cerca
          Y0=FPOLY(NDEGBL(LL2),COEFFBL(1,LL2),X0)
        END IF
        IF(LDEBUGLOCAL) CALL PGPOINT(1,X0,Y0,24)
C coeficientes del polinomio trasladado (son todos iguales menos el primero, en
C que introducimos el offset)
        IF(LMIN.EQ.LL1)THEN
          NDEGTRAS=NDEGBL(LL1)
          DO K=1,NDEGTRAS+1
            COEFFTRAS(K)=COEFFBL(K,LL1)
          END DO
          COEFFTRAS(1)=COEFFTRAS(1)+(YC-Y0)
        ELSE
          NDEGTRAS=NDEGBL(LL2)
          DO K=1,NDEGTRAS+1
            COEFFTRAS(K)=COEFFBL(K,LL2)
          END DO
          COEFFTRAS(1)=COEFFTRAS(1)+(YC-Y0)
        END IF
C calculamos la interseccion de ese polinomio con la linea de la izquierda
        CALL INTERSEC(NDEGBA(L1),COEFFBA(1,L1),NDEGTRAS,COEFFTRAS,
     +   XC,XR,YR)
C calculamos la interseccion de ese polinomio con la linea de la derecha
        CALL INTERSEC(NDEGBA(L2),COEFFBA(1,L2),NDEGTRAS,COEFFTRAS,
     +   XC,XQ,YQ)
        IF(LDEBUGLOCAL)THEN
          CALL PGSCI(3)
          CALL PGPOINT(1,XR,YR,17)
          CALL PGPOINT(1,XQ,YQ,17)
          CALL PGSCI(1)
        END IF
C calculamos relacion de arcos s1/s0
        S0=ARCLENGTH(NDEGTRAS,COEFFTRAS,XR,XQ,100)
        S1=ARCLENGTH(NDEGTRAS,COEFFTRAS,XR,XC,100)
C------------------------------------------------------------------------------
C Ya tenemos el valor inicial para comenzar a iterar el asunto. Aqui empieza
C el meollo.
C
C calculamos el punto N
        CALL X2ARC(NDEGBL(LL2),COEFFBL(1,LL2),XMAX1,XMAX2,S1/S0,XN)
        YN=FPOLY(NDEGBL(LL2),COEFFBL(1,LL2),XN)
C calculamos el punto M
        CALL X2ARC(NDEGBL(LL1),COEFFBL(1,LL1),XMIN1,XMIN2,S1/S0,XM)
        YM=FPOLY(NDEGBL(LL1),COEFFBL(1,LL1),XM)
        IF(LDEBUGLOCAL)THEN
          CALL PGSCI(6)
          CALL PGPOINT(1,XN,YN,17)
          CALL PGPOINT(1,XM,YM,17)
          CALL PGSCI(1)
        END IF
        LOOP=.TRUE.
        NITER=0
        DO WHILE(LOOP)
C calculamos los arcos n0, n1, m0 y m1
          N0=ARCLENGTH(NDEGBL(LL2),COEFFBL(1,LL2),XMAX1,XMAX2,100)
          N1=ARCLENGTH(NDEGBL(LL2),COEFFBL(1,LL2),XMAX1,XN,100)
          M0=ARCLENGTH(NDEGBL(LL1),COEFFBL(1,LL1),XMIN1,XMIN2,100)
          M1=ARCLENGTH(NDEGBL(LL1),COEFFBL(1,LL1),XMIN1,XM,100)
          IF(LDEBUGLOCAL)THEN
            print*,' '
            print*,'NITER=',NITER
            print*,'Hemos impuesto  s1/s0: ',s1/s0
            print*,'Obtenemos n1/n0,m1/m0:',n1/n0,m1/m0
          END IF
C calculamos el polinomio x=g(y) de segundo grado que pasa por M, P y N
          CALL POLY2(YM,XM,YC,XC,YN,XN,COEFFP1)
C calculamos relacion de arcos t1/t0
          T0=ARCLENGTH(2,COEFFP1,YM,YN,100)
          T1=ARCLENGTH(2,COEFFP1,YM,YC,100)
C calculamos el punto R'
          CALL X2ARC(NDEGBA(L1),COEFFBA(1,L1),YMIN1,YMAX1,T1/T0,YR_)
          XR_=FPOLY(NDEGBA(L1),COEFFBA(1,L1),YR_)
C calculamos el punto Q'
          CALL X2ARC(NDEGBA(L2),COEFFBA(1,L2),YMIN2,YMAX2,T1/T0,YQ_)
          XQ_=FPOLY(NDEGBA(L2),COEFFBA(1,L2),YQ_)
C calculamos el polinomio y=f(x) de segundo grado que pasa por R', P y Q'
          CALL POLY2(XR_,YR_,XC,YC,XQ_,YQ_,COEFFP2)
C calculamos la nueva relacion de arcos s1/s0
          S0_=ARCLENGTH(2,COEFFP2,XR_,XQ_,100)
          S1_=ARCLENGTH(2,COEFFP2,XR_,XC,100)
          IF(LDEBUGLOCAL)THEN
            print*,'s0 ,s1 ,s1/s0:',s0,s1,s1/s0
            print*,'s0_,s1_,s1/s0:',s0_,s1_,s1_/s0_
            print*,'error........:',s1/s0-s1_/s0_
          END IF
C calculamos el punto N'
          CALL X2ARC(NDEGBL(LL2),COEFFBL(1,LL2),XMAX1,XMAX2,
     +     S1_/S0_,XN_)
          YN_=FPOLY(NDEGBL(LL2),COEFFBL(1,LL2),XN_)
C calculamos el punto M'
          CALL X2ARC(NDEGBL(LL1),COEFFBL(1,LL1),XMIN1,XMIN2,
     +     S1_/S0_,XM_)
          YM_=FPOLY(NDEGBL(LL1),COEFFBL(1,LL1),XM_)
C calculamos el polinomio x=g(y) de segundo grado que pasa por M', P y N'
          CALL POLY2(YM_,XM_,YC,XC,YN_,XN_,COEFFP1)
C calculamos relacion de arcos t1/t0
          T0_=ARCLENGTH(2,COEFFP1,YM_,YN_,100)
          T1_=ARCLENGTH(2,COEFFP1,YM_,YC,100)
          IF(LDEBUGLOCAL)THEN
            print*,'t0 ,t1 ,t1/t0:',t0,t1,t1/t0
            print*,'t0_,t1_,t1/t0:',t0_,t1_,t1_/t0_
            print*,'error........:',t1/t0-t1_/t0_
          END IF
C decidimos si volvemos a iterar o ya nos vale
          LOOP=((ABS(S1/S0-S1_/S0_).GE.THRESHOLD).OR.
     +     (ABS(T1/T0-T1_/T0_).GE.THRESHOLD))
          IF(LOOP)THEN
            NITER=NITER+1
            IF(NITER.LE.NITERMAX)THEN
              S1=S1_
              S0=S0_
              XN=XN_
              YN=YN_
              XM=XM_
              YM=YM_
            ELSE
              LOOP=.FALSE.
            END IF
          END IF
        END DO
        IF(NITER.GT.NITERMAX)THEN
          WRITE(*,101) '***WARNING***'
          WRITE(*,101) '=> NITERMAX in WHEREAMI'
        END IF
        IF(LDEBUGLOCAL)THEN
          XMIN=XR_
          XMAX=XQ_
          DO I=1,NXYMAX
            XP(I)=REAL(I-1)/REAL(NXYMAX-1)*(XMAX-XMIN)+XMIN
            YP(I)=FPOLY(2,COEFFP2,XP(I))
          END DO
          CALL PGSCI(3)
          CALL PGLINE(NXYMAX,XP,YP)
          CALL PGSCI(1)
          YMIN=YM_
          YMAX=YN_
          DO I=1,NXYMAX
            YP(I)=REAL(I-1)/REAL(NXYMAX-1)*(YMAX-YMIN)+YMIN
            XP(I)=FPOLY(2,COEFFP1,YP(I))
          END DO
          CALL PGSCI(3)
          CALL PGLINE(NXYMAX,XP,YP)
          CALL PGSCI(1)
        END IF
C------------------------------------------------------------------------------
C Ya conocemos la relacion s1/s0 y t1/t0. Ahora vamos a calcular el espectro
C que pasa por ese punto, asi como la linea espectral. Para ello imponemos
C los valores s1/s0 y t1/t0 que acabamos de calcular.
C------------------------------------------------------------------------------
C calculemos el polinomio en la direccion espectral
        DO L=1,NLINBA
          CALL INTERSEC(NDEGBA(L),COEFFBA(1,L),
     +     NDEGBL(LL1),COEFFBL(1,LL1),
     +     0.5*(XMINBA(L)+XMAXBA(L)),XMIN,YMIN)
          CALL INTERSEC(NDEGBA(L),COEFFBA(1,L),
     +     NDEGBL(LL2),COEFFBL(1,LL2),
     +     0.5*(XMINBA(L)+XMAXBA(L)),XMAX,YMAX)
          CALL X2ARC(NDEGBA(L),COEFFBA(1,L),YMIN,YMAX,T1_/T0_,YF(L))
          XF(L)=FPOLY(NDEGBA(L),COEFFBA(1,L),YF(L))
          IF(LDEBUGLOCAL)THEN
            CALL PGPOINT(1,XMIN,YMIN,22)
            CALL PGPOINT(1,XMAX,YMAX,22)
            CALL PGPOINT(1,XF(L),YF(L),22)
          END IF
        END DO
        XF(NLINBA+1)=XC !el polinomio tambien debe pasar por XC,YC
        YF(NLINBA+1)=YC
c
        DO K=1,NDEGBL00+1
          IFCOEFF(K)=.TRUE.
        END DO
        CALL FITPOL(NLINBA+1,XF,YF,NDEGBL00,IFCOEFF,COEFFBL00,
     +   VARDUM,CHISQR,SR2,COVAR)
        IF(LDEBUGLOCAL)THEN
          CALL FINDMML(NLINBA+1,1,NLINBA+1,XF,XMIN,XMAX)
          DO I=1,NXYMAX
            XP(I)=REAL(I-1)/REAL(NXYMAX-1)*(XMAX-XMIN)+XMIN
            YP(I)=FPOLY(NDEGBL00,COEFFBL00,XP(I))
          END DO
          CALL PGSCI(7)
          CALL PGLINE(NXYMAX,XP,YP)
          CALL PGSCI(1)
        END IF
        YINI=FPOLY(NDEGBL00,COEFFBL00,XC)
        LOOP=(ABS(YC-YINI).GE.THRESHOLDBIS)
        NITER=0
        DO WHILE(LOOP)
          COEFFBL00(1)=COEFFBL00(1)-(YINI-YC)
          IFCOEFF(1)=.FALSE. !no ajustamos el primer coeficiente
          CALL FITPOL(NLINBA+1,XF,YF,NDEGBL00,IFCOEFF,COEFFBL00,
     +     VARDUM,CHISQR,SR2,COVAR)
          IF(LDEBUGLOCAL)THEN
            DO I=1,NXYMAX
              XP(I)=REAL(I-1)/REAL(NXYMAX-1)*(XMAX-XMIN)+XMIN
              YP(I)=FPOLY(NDEGBL00,COEFFBL00,XP(I))
            END DO
            CALL PGSCI(2)
            CALL PGLINE(NXYMAX,XP,YP)
            CALL PGSCI(1)
          END IF
          YFIN=FPOLY(NDEGBL00,COEFFBL00,XC)
          LOOP=(ABS(YC-YFIN).GE.THRESHOLDBIS)
          IF(LOOP)THEN
            NITER=NITER+1
            LOOP=(NITER.LT.NITERMAX)
            YINI=YFIN
          END IF
        END DO
        IF(NCOL1.GT.-1)THEN
          XMIN=XMINBL(1)
          XMAX=XMAXBL(1)
          DO I=1,NXYMAX
            XP(I)=REAL(I-1)/REAL(NXYMAX-1)*(XMAX-XMIN)+XMIN
            YP(I)=FPOLY(NDEGBL00,COEFFBL00,XP(I))
          END DO
          CALL PGSCI(NCOL1)
          CALL PGLINE(NXYMAX,XP,YP)
          CALL PGSCI(1)
        END IF
C------------------------------------------------------------------------------
C calculemos el polinomio en la direccion espacial
        DO L=1,NLINBL
          CALL INTERSEC(NDEGBA(L1),COEFFBA(1,L1),
     +     NDEGBL(L),COEFFBL(1,L),
     +     0.5*(XMINBA(L1)+XMAXBA(L1)),XMIN,YMIN)
          CALL INTERSEC(NDEGBA(L2),COEFFBA(1,L2),
     +     NDEGBL(L),COEFFBL(1,L),
     +     0.5*(XMINBA(L2)+XMAXBA(L2)),XMAX,YMAX)
          CALL X2ARC(NDEGBL(L),COEFFBL(1,L),XMIN,XMAX,S1/S0,XF(L))
          YF(L)=FPOLY(NDEGBL(L),COEFFBL(1,L),XF(L))
          IF(LDEBUGLOCAL)THEN
            CALL PGPOINT(1,XMIN,YMIN,22)
            CALL PGPOINT(1,XMAX,YMAX,22)
            CALL PGPOINT(1,XF(L),YF(L),22)
          END IF
        END DO
        XF(NLINBL+1)=XC !el polinomio tambien debe pasar por XC,YC
        YF(NLINBL+1)=YC
c
        DO K=1,NDEGBA00+1
          IFCOEFF(K)=.TRUE.
        END DO
        CALL FITPOL(NLINBL+1,YF,XF,NDEGBA00,IFCOEFF,COEFFBA00,
     +   VARDUM,CHISQR,SR2,COVAR)
        IF(LDEBUGLOCAL)THEN
          CALL FINDMML(NLINBL+1,1,NLINBL+1,YF,YMIN,YMAX)
          DO I=1,NXYMAX
            YP(I)=REAL(I-1)/REAL(NXYMAX-1)*(YMAX-YMIN)+YMIN
            XP(I)=FPOLY(NDEGBA00,COEFFBA00,YP(I))
          END DO
          CALL PGSCI(7)
          CALL PGLINE(NXYMAX,XP,YP)
          CALL PGSCI(1)
        END IF
        XINI=FPOLY(NDEGBA00,COEFFBA00,YC)
        LOOP=(ABS(XC-XINI).GE.THRESHOLDBIS)
        NITER=0
        DO WHILE(LOOP)
          COEFFBA00(1)=COEFFBA00(1)-(XINI-XC)
          IFCOEFF(1)=.FALSE. !no ajustamos el primer coeficiente
          CALL FITPOL(NLINBL+1,YF,XF,NDEGBA00,IFCOEFF,COEFFBA00,
     +     VARDUM,CHISQR,SR2,COVAR)
          IF(LDEBUGLOCAL)THEN
            DO I=1,NXYMAX
              YP(I)=REAL(I-1)/REAL(NXYMAX-1)*(YMAX-YMIN)+YMIN
              XP(I)=FPOLY(NDEGBA00,COEFFBA00,YP(I))
            END DO
            CALL PGSCI(2)
            CALL PGLINE(NXYMAX,XP,YP)
            CALL PGSCI(1)
          END IF
          XFIN=FPOLY(NDEGBA00,COEFFBA00,YC)
          LOOP=(ABS(XC-XFIN).GE.THRESHOLDBIS)
          IF(LOOP)THEN
            NITER=NITER+1
            LOOP=(NITER.LT.NITERMAX)
            XINI=XFIN
          END IF
        END DO
        IF(NCOL2.GT.-1)THEN
          CALL FINDMML(NLINBL+1,1,NLINBL+1,YF,YMIN,YMAX)
          DO I=1,NXYMAX
            YP(I)=REAL(I-1)/REAL(NXYMAX-1)*(YMAX-YMIN)+YMIN
            XP(I)=FPOLY(NDEGBA00,COEFFBA00,YP(I))
          END DO
          CALL PGSCI(NCOL2)
          CALL PGLINE(NXYMAX,XP,YP)
          CALL PGSCI(1)
        END IF
C------------------------------------------------------------------------------
101     FORMAT(A)
        END

! Version 06-February-2001
! Determina la fraccion de luz que atraviesa una rendija con una anchura y
! orientacion determinadas. Tambien realiza simulaciones de Monte Carlo para
! estimar el error.
        SUBROUTINE SLIT_APERTURE(NCBUFF)
        IMPLICIT NONE
        INTEGER NCBUFF
!
        INCLUDE 'dimensions.inc'
!
        INTEGER MXCOR
        PARAMETER (MXCOR=8) !numero maximo de intersecciones de dos rombos
        INTEGER NSIMULMAX
        PARAMETER (NSIMULMAX=1000) !numero maximo de simulaciones para errores
        REAL PI
        PARAMETER (PI=3.141592654)
!
        INTEGER READILIM
        REAL READF
        REAL RANDOMNUMBER
        REAL FMEAN0
!
        INTEGER I,J,K
        INTEGER NX1,NX2,NY1,NY2
        INTEGER NX1_,NX2_,NY1_,NY2_
        INTEGER NAXIS(2,NMAXBUFF)
        INTEGER IXC1,IXC2,IYC1,IYC2
        INTEGER NCOR
        INTEGER NSIMUL
        INTEGER NSEED
        REAL XC,YC
        REAL IMAGEN(NXMAX,NYMAX,NMAXBUFF)
        REAL XC_SLIT,YC_SLIT,ERR_XY_SLIT
        REAL WIDTH_SLIT,LENGTH_SLIT
        REAL PA_SLIT,OFFPA_SLIT,ERR_PA_SLIT
        REAL X1,Y1,X2,Y2,X3,Y3,X4,Y4
        REAL XCOR1(2,4),XCOR2(2,4),XCOM(2,MXCOR)
        REAL COSA,SINA
        REAL SUM_TOTAL,SUM_INSLIT,FAREA
        REAL SUM_INSLIT_(NSIMULMAX)
        REAL SUM_INSLIT_MEAN,SUM_INSLIT_SIGMA
        REAL R1,R2,ERRORX,ERRORY,ERRORPA
        CHARACTER*1 CH
        LOGICAL MASKPIXEL(NXMAXB9,NYMAXB9)
        LOGICAL LOOP
!
        COMMON/BLKIMAGEN1/IMAGEN
        COMMON/BLKNAXIS/NAXIS
        COMMON/BLKXYLIMPLOT/NX1,NX2,NY1,NY2
        COMMON/BLKMASKBOX9/MASKPIXEL !NOTE: aqui estamos usando la misma
                                     !direccion de memoria que MASKBOX9 para
                                     !evitar sobrecargar mucho el programa.
!------------------------------------------------------------------------------
! Proteccion
        IF((NX2-NX1+1).GT.NXMAXB9)THEN
          WRITE(*,100) 'NX2-NX1+1, NXMAXB9: '
          WRITE(*,*) NX2-NX1+1,NXMAXB9
          WRITE(*,101) 'ERROR: NX2-NX1+1.GT.NXMAXB9'
          WRITE(*,100) 'Press <CR> to continue...'
          READ(*,*)
          RETURN
        END IF
        IF((NY2-NY1+1).GT.NYMAXB9)THEN
          WRITE(*,100) 'NY2-NY1+1, NYMAXB9: '
          WRITE(*,*) NY2-NY1+1,NXMAXB9
          WRITE(*,101) 'ERROR: NY2-NY1+1.GT.NYMAXB9'
          WRITE(*,100) 'Press <CR> to continue...'
          READ(*,*)
          RETURN
        END IF
!------------------------------------------------------------------------------
! Definimos los pixels que definen el objeto que queremos medir
        WRITE(*,100) 'Use the mouse to indicate the image regions'
        WRITE(*,101) 'that define the object'
! inicializamos la mascara
        DO I=NY1,NY2
          DO J=NX1,NX2
            MASKPIXEL(J-NX1+1,I-NY1+1)=.FALSE.
          END DO
        END DO
! definimos regiones a utilizar
        LOOP=.TRUE.
        DO WHILE(LOOP)
          CALL RPGBAND(0,0,0.,0.,XC,YC,CH)
          IF(CH.EQ.'X')THEN
            LOOP=.FALSE.
          END IF
          IF(LOOP)THEN
            IXC1=INT(XC+0.5)
            IF(IXC1.LT.NX1) IXC1=NX1
            IF(IXC1.GT.NX2) IXC1=NX2
            IYC1=INT(YC+0.5)
            IF(IYC1.LT.NY1) IYC1=NY1
            IF(IYC1.GT.NY2) IYC1=NY2
            WRITE(*,112) 'Cursor at ',IXC1,IYC1
            CALL PGSCI(5)
            CALL RPGBAND(2,0,REAL(IXC1),REAL(IYC1),XC,YC,CH)
            CALL PGSCI(1)
            IF(CH.EQ.'X')THEN
              LOOP=.FALSE.
            END IF
            IF(LOOP)THEN
              IXC2=INT(XC+0.5)
              IF(IXC2.LT.NX1) IXC2=NX1
              IF(IXC2.GT.NX2) IXC2=NX2
              IYC2=INT(YC+0.5)
              IF(IYC2.LT.NY1) IYC2=NY1
              IF(IYC2.GT.NY2) IYC2=NY2
              WRITE(*,112) 'Cursor at ',IXC2,IYC2
              NX1_=MIN0(IXC1,IXC2)
              NX2_=MAX0(IXC1,IXC2)
              NY1_=MIN0(IYC1,IYC2)
              NY2_=MAX0(IYC1,IYC2)
              DO I=NY1_,NY2_
                DO J=NX1_,NX2_
                  MASKPIXEL(J-NX1+1,I-NY1+1)=.TRUE.
                END DO
              END DO
              CALL SUBLOOK(.TRUE.,NCBUFF,.FALSE.)
              CALL PGBBUF
              CALL PGSCI(4)
              DO I=NY1,NY2
                DO J=NX1,NX2
                  IF(MASKPIXEL(J-NX1+1,I-NY1+1))THEN
                    CALL PGMOVE(REAL(J)-0.5,REAL(I)-0.5)
                    CALL PGDRAW(REAL(J)+0.5,REAL(I)-0.5)
                    CALL PGDRAW(REAL(J)+0.5,REAL(I)+0.5)
                    CALL PGDRAW(REAL(J)-0.5,REAL(I)+0.5)
                    CALL PGDRAW(REAL(J)-0.5,REAL(I)-0.5)
                  END IF
                END DO
              END DO
              CALL PGSCI(1)
              CALL PGEBUF
            END IF
          END IF
        END DO
!------------------------------------------------------------------------------
! Definimos la rendija
        XC_SLIT    =READF('X centroid of the slit (pixels)','@')
        YC_SLIT    =READF('Y centroid of the slit (pixels)','@')
        WIDTH_SLIT =READF('Slit width ............(pixels)','@')
        LENGTH_SLIT=READF('Slit length ...........(pixels)','@')
        PA_SLIT    =READF('PA of the slit .......(degrees)','@')
        OFFPA_SLIT =READF('PA offset ............(degrees)','@')
        SINA=SIN((PA_SLIT+OFFPA_SLIT)*PI/180.0)
        COSA=COS((PA_SLIT+OFFPA_SLIT)*PI/180.0)
! Dibujamos la rendija
        X1=XC_SLIT-LENGTH_SLIT/2.0*COSA+WIDTH_SLIT/2.0*SINA
        Y1=YC_SLIT-LENGTH_SLIT/2.0*SINA-WIDTH_SLIT/2.0*COSA
        X2=XC_SLIT+LENGTH_SLIT/2.0*COSA+WIDTH_SLIT/2.0*SINA
        Y2=YC_SLIT+LENGTH_SLIT/2.0*SINA-WIDTH_SLIT/2.0*COSA
        X3=XC_SLIT+LENGTH_SLIT/2.0*COSA-WIDTH_SLIT/2.0*SINA
        Y3=YC_SLIT+LENGTH_SLIT/2.0*SINA+WIDTH_SLIT/2.0*COSA
        X4=XC_SLIT-LENGTH_SLIT/2.0*COSA-WIDTH_SLIT/2.0*SINA
        Y4=YC_SLIT-LENGTH_SLIT/2.0*SINA+WIDTH_SLIT/2.0*COSA
        CALL PGSCI(7)
        CALL PGMOVE(X1,Y1)
        CALL PGDRAW(X2,Y2)
        CALL PGMOVE(X4,Y4)
        CALL PGDRAW(X3,Y3)
        CALL PGSCI(1)
!------------------------------------------------------------------------------
! Para el calculo de las intersecciones, cambiamos el origen de la escala
! para no tener valores demasiado grandes
        XCOR1(1,1)=X1-XC_SLIT
        XCOR1(2,1)=Y1-YC_SLIT
        XCOR1(1,2)=X2-XC_SLIT
        XCOR1(2,2)=Y2-YC_SLIT
        XCOR1(1,3)=X3-XC_SLIT
        XCOR1(2,3)=Y3-YC_SLIT
        XCOR1(1,4)=X4-XC_SLIT
        XCOR1(2,4)=Y4-YC_SLIT
!------------------------------------------------------------------------------
! Calculamos el area dentro de la rendija y el area total del objeto
! seleccionado.
! NOTA: el procedimiento calcula incluso las fracciones de pixel dentro de la
! rendija, por lo que funciona en un caso general, cualesquiera que sean las
! dimensiones relativas de los pixels y de la rendija considerada.
        SUM_TOTAL=0.0
        SUM_INSLIT=0.0
        DO I=NY1,NY2
          DO J=NX1,NX2
            IF(MASKPIXEL(J-NX1+1,I-NY1+1))THEN !si el pixel puede utilizarse
              SUM_TOTAL=SUM_TOTAL+IMAGEN(J,I,NCBUFF)
              XCOR2(1,1)=REAL(J)-0.5-XC_SLIT
              XCOR2(2,1)=REAL(I)-0.5-YC_SLIT
              XCOR2(1,2)=REAL(J)+0.5-XC_SLIT
              XCOR2(2,2)=REAL(I)-0.5-YC_SLIT
              XCOR2(1,3)=REAL(J)+0.5-XC_SLIT
              XCOR2(2,3)=REAL(I)+0.5-YC_SLIT
              XCOR2(1,4)=REAL(J)-0.5-XC_SLIT
              XCOR2(2,4)=REAL(I)+0.5-YC_SLIT
              CALL PPOLY4(XCOR2(1,1)+XC_SLIT,XCOR2(2,1)+YC_SLIT, &
                          XCOR2(1,2)+XC_SLIT,XCOR2(2,2)+YC_SLIT, &
                          XCOR2(1,3)+XC_SLIT,XCOR2(2,3)+YC_SLIT, &
                          XCOR2(1,4)+XC_SLIT,XCOR2(2,4)+YC_SLIT,8)
              CALL SRCPAT(4,XCOR1,4,XCOR2,NCOR,XCOM)
              IF(NCOR.GT.0)THEN !hay interseccion
                CALL PPOLYN(NCOR,XCOM,FAREA)
                SUM_INSLIT=SUM_INSLIT+FAREA*IMAGEN(J,I,NCBUFF)
                CALL PGSCI(5)
                CALL PGMOVE(XCOM(1,1)+xc_slit,XCOM(2,1)+yc_slit)
                DO K=2,NCOR
                  CALL PGDRAW(XCOM(1,K)+xc_slit,XCOM(2,K)+yc_slit)
                END DO
                CALL PGDRAW(XCOM(1,1)+XC_SLIT,XCOM(2,1)+YC_SLIT)
                CALL PGSCI(1)
              END IF
            END IF
          END DO
        END DO
        WRITE(*,100) '>>> Flux (inside slit, total): '
        WRITE(*,*) SUM_INSLIT,SUM_TOTAL
        WRITE(*,100) '>>> Fraction of light in slit: '
        WRITE(*,*) SUM_INSLIT/SUM_TOTAL
!------------------------------------------------------------------------------
! Simulaciones para estimar errores
        NSIMUL=READILIM('No. of simulations to estimate errors','100',0,NSIMULMAX)
        IF(NSIMUL.GT.0)THEN
          ERR_XY_SLIT=READF('Error in [x,y] slit centroid (pixels)','@')
          ERR_PA_SLIT=READF('Error in PA of the slit ....(degrees)','@')
          NSEED=-1
          DO K=1,NSIMUL
            R1=RANDOMNUMBER(NSEED)
            R2=RANDOMNUMBER(NSEED)
            ERRORX=1.41421356*ERR_XY_SLIT*SQRT(-1.*LOG(1.-R1))*COS(2.*PI*R2)
            R1=RANDOMNUMBER(NSEED)
            R2=RANDOMNUMBER(NSEED)
            ERRORY=1.41421356*ERR_XY_SLIT*SQRT(-1.*LOG(1.-R1))*COS(2.*PI*R2)
            R1=RANDOMNUMBER(NSEED)
            R2=RANDOMNUMBER(NSEED)
            ERRORPA=1.41421356*ERR_PA_SLIT*SQRT(-1.*LOG(1.-R1))*COS(2.*PI*R2)
            SINA=SIN((PA_SLIT+ERRORPA+OFFPA_SLIT)*PI/180.0)
            COSA=COS((PA_SLIT+ERRORPA+OFFPA_SLIT)*PI/180.0)
            X1=XC_SLIT+ERRORX-LENGTH_SLIT/2.0*COSA+WIDTH_SLIT/2.0*SINA
            Y1=YC_SLIT+ERRORY-LENGTH_SLIT/2.0*SINA-WIDTH_SLIT/2.0*COSA
            X2=XC_SLIT+ERRORX+LENGTH_SLIT/2.0*COSA+WIDTH_SLIT/2.0*SINA
            Y2=YC_SLIT+ERRORY+LENGTH_SLIT/2.0*SINA-WIDTH_SLIT/2.0*COSA
            X3=XC_SLIT+ERRORX+LENGTH_SLIT/2.0*COSA-WIDTH_SLIT/2.0*SINA
            Y3=YC_SLIT+ERRORY+LENGTH_SLIT/2.0*SINA+WIDTH_SLIT/2.0*COSA
            X4=XC_SLIT+ERRORX-LENGTH_SLIT/2.0*COSA-WIDTH_SLIT/2.0*SINA
            Y4=YC_SLIT+ERRORY-LENGTH_SLIT/2.0*SINA+WIDTH_SLIT/2.0*COSA
            CALL PGSCI(0)
            CALL PGMOVE(X1,Y1)
            CALL PGDRAW(X2,Y2)
            CALL PGMOVE(X4,Y4)
            CALL PGDRAW(X3,Y3)
            CALL PGSCI(1)
            XCOR1(1,1)=X1-XC_SLIT
            XCOR1(2,1)=Y1-YC_SLIT
            XCOR1(1,2)=X2-XC_SLIT
            XCOR1(2,2)=Y2-YC_SLIT
            XCOR1(1,3)=X3-XC_SLIT
            XCOR1(2,3)=Y3-YC_SLIT
            XCOR1(1,4)=X4-XC_SLIT
            XCOR1(2,4)=Y4-YC_SLIT
            SUM_INSLIT_(K)=0.0
            DO I=NY1,NY2
              DO J=NX1,NX2
                IF(MASKPIXEL(J-NX1+1,I-NY1+1))THEN
                  XCOR2(1,1)=REAL(J)-0.5-XC_SLIT
                  XCOR2(2,1)=REAL(I)-0.5-YC_SLIT
                  XCOR2(1,2)=REAL(J)+0.5-XC_SLIT
                  XCOR2(2,2)=REAL(I)-0.5-YC_SLIT
                  XCOR2(1,3)=REAL(J)+0.5-XC_SLIT
                  XCOR2(2,3)=REAL(I)+0.5-YC_SLIT
                  XCOR2(1,4)=REAL(J)-0.5-XC_SLIT
                  XCOR2(2,4)=REAL(I)+0.5-YC_SLIT
                  CALL SRCPAT(4,XCOR1,4,XCOR2,NCOR,XCOM)
                  IF(NCOR.GT.0)THEN !hay interseccion
                    CALL PPOLYN(NCOR,XCOM,FAREA)
                    SUM_INSLIT_(K)=SUM_INSLIT_(K)+FAREA*IMAGEN(J,I,NCBUFF)
                  END IF
                END IF
              END DO
            END DO
          END DO
          SUM_INSLIT_MEAN=FMEAN0(NSIMUL,SUM_INSLIT_,SUM_INSLIT_SIGMA)
          WRITE(*,100) '>>> Flux inside slit & rms...........: '
          WRITE(*,*) SUM_INSLIT_MEAN,' +/-',SUM_INSLIT_SIGMA
          WRITE(*,100) '>>> Fraction of light in slit & error: '
          WRITE(*,*) SUM_INSLIT_MEAN/SUM_TOTAL,' +/-',SUM_INSLIT_SIGMA/SUM_TOTAL
        END IF
!------------------------------------------------------------------------------
100     FORMAT(A,$)
101     FORMAT(A)
112     FORMAT(A,I5,2X,I5)
        END

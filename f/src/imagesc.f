      PROGRAM IMAGESC
      IMPLICIT NONE
      INTEGER :: STATUS,NAXES(3),XMIN,XMAX,YMIN,YMAX,NARGS,K,ITF,FRAME
      REAL,PARAMETER :: DEFAULTSC=0.0429718,
     &  GL(2)=(/0.0, 1.0/),
     &  GR(2)=(/0.0, 1.0/),
     &  GG(2)=(/0.0, 1.0/),
     &  GB(2)=(/0.0, 1.0/),
     &  RL(9)=(/-0.5, 0.0, 0.17, 0.33, 0.50, 0.67, 0.83, 1.0, 1.7/),
     &  RR(9)=(/ 0.0, 0.0,  0.0,  0.0,  0.6,  1.0,  1.0, 1.0, 1.0/),
     &  RG(9)=(/ 0.0, 0.0,  0.0,  1.0,  1.0,  1.0,  0.6, 0.0, 1.0/),
     &  RB(9)=(/ 0.0, 0.3,  0.8,  1.0,  0.3,  0.0,  0.0, 0.0, 1.0/),
     &  HL(5)=(/0.0, 0.2, 0.4, 0.6, 1.0/),
     &  HR(5)=(/0.0, 0.5, 1.0, 1.0, 1.0/),
     &  HG(5)=(/0.0, 0.0, 0.5, 1.0, 1.0/),
     &  HB(5)=(/0.0, 0.0, 0.0, 0.3, 1.0/)
      REAL :: TR(6),VPX1,VPX2,VPY1,VPY2,D,SC,CMIN,CMAX
      DOUBLE PRECISION,ALLOCATABLE :: DIMG(:,:)
      CHARACTER(LEN=256) :: IMGFILE,ARG,TITLE,XLAB,YLAB,CLAB,COLOR,
     &  DEVICE,MAP
C
      INTERFACE
      SUBROUTINE READIMAGE(FILENAME,FPIXELS,LPIXELS,DIMG)
      INTEGER, INTENT(IN) :: FPIXELS(3),LPIXELS(3)
      DOUBLE PRECISION, INTENT(OUT) :: DIMG(*)
      CHARACTER(LEN=*), INTENT(IN) :: FILENAME
      END SUBROUTINE READIMAGE
      SUBROUTINE IMAGESIZE(FILENAME,NAXES)
      INTEGER, INTENT(OUT) :: NAXES(3)
      CHARACTER(LEN=*) :: FILENAME
      END SUBROUTINE IMAGESIZE
      FUNCTION PGOPEN(DEVICE)
      CHARACTER(LEN=*) :: DEVICE
      INTEGER :: PGOPEN
      END FUNCTION PGOPEN
      SUBROUTINE PGENV(XMIN,XMAX,YMIN,YMAX,JUST,AXIS)
      REAL :: XMIN,XMAX,YMIN,YMAX
      INTEGER :: JUST,AXIS
      END SUBROUTINE PGENV
      SUBROUTINE PGLAB(XLBL,YLBL,TOPLBL)
      CHARACTER(LEN=*) :: XLBL,YLBL,TOPLBL
      END SUBROUTINE PGLAB
      SUBROUTINE PGCLOS
      END SUBROUTINE
      SUBROUTINE PGIMAG(A,NI,NJ,I1,I2,J1,J2,A1,A2,TR)
      INTEGER :: NI,NJ,I1,I2,J1,J2
      REAL :: A(NI,NJ),A1,A2,TR(6)
      END SUBROUTINE PGIMAG
      SUBROUTINE PGSVP(XLEFT,XRIGHT,YBOT,YTOP)
      REAL :: XLEFT,XRIGHT,YBOT,YTOP
      END SUBROUTINE PGSVP
      SUBROUTINE PGPAGE
      END SUBROUTINE PGPAGE
      SUBROUTINE PGCTAB(L,R,G,B,NC,CONTRA,BRIGHT)
      INTEGER :: NC
      REAL :: L(NC),R(NC),G(NC),B(NC),CONTRA,BRIGHT
      END SUBROUTINE PGCTAB
      SUBROUTINE PGWNAD(X1,X2,Y1,Y2)
      REAL :: X1, X2, Y1, Y2
      END SUBROUTINE PGWNAD
      SUBROUTINE PGBOX(XOPT,XTICK,NXSUB,YOPT,YTICK,NYSUB)
      CHARACTER(LEN=*) :: XOPT, YOPT
      REAL :: XTICK, YTICK
      INTEGER :: NXSUB, NYSUB
      END SUBROUTINE PGBOX
      SUBROUTINE PGMTXT(SIDE,DISP,COORD,FJUST,TEXT)
      CHARACTER(LEN=*) :: SIDE, TEXT
      REAL :: DISP, COORD, FJUST
      END SUBROUTINE PGMTXT
      SUBROUTINE PGSCH(SZ)
      REAL :: SZ
      END SUBROUTINE PGSCH
      SUBROUTINE PGVSTD
      END SUBROUTINE PGVSTD
      SUBROUTINE PGQVP(UNITS, X1, X2, Y1, Y2)
      INTEGER :: UNITS
      REAL :: X1, X2, Y1, Y2
      END SUBROUTINE PGQVP
      SUBROUTINE PGWEDG(SIDE, DISP, WIDTH, FG, BG, LABEL)
      CHARACTER(LEN=*) :: SIDE,LABEL
      REAL             :: DISP, WIDTH, FG, BG
      END SUBROUTINE PGWEDG
      SUBROUTINE PGSCF(FONT)
      INTEGER :: FONT
      END SUBROUTINE PGSCF
      SUBROUTINE PGSITF(ITF)
      INTEGER :: ITF
      END SUBROUTINE PGSITF
      END INTERFACE
C
      STATUS=0
      CALL GET_COMMAND_ARGUMENT(1,ARG)
      IF(INDEX(ARG,'-help').GT.0)THEN
        PRINT *,'Usage:'
        PRINT *,'======'
        PRINT *,'imagesc filename [-xmin=x_min] [-xmax=x_max]'
        PRINT *,'  [-ymin=y_min] [-ymax=y_max] [-scale=sc]'
        PRINT *,'  [-max=c_max] [-min=c_min] [-xlabel=x_label]'
        PRINT *,'  [-ylabel=y_label] [-title=title] [-clabel=c_label]'
        PRINT *,'  [-color=color_scheme] [-device=device_name]'
        PRINT *,'  [-map=linear,log,sqrt] [-frame=n]'
        STOP
      END IF
      NARGS=COMMAND_ARGUMENT_COUNT()
      CALL GET_COMMAND_ARGUMENT(1,IMGFILE)
      CALL IMAGESIZE(IMGFILE,NAXES)
      ALLOCATE(DIMG(NAXES(1),NAXES(2)),STAT=STATUS)
      IF(STATUS.NE.0)THEN
        WRITE(*,*)'error: out of memory.'
        STOP
      END IF
      FRAME=1
      DO K=2,NARGS
        CALL GET_COMMAND_ARGUMENT(K,ARG)
        IF(INDEX(ARG,'-frame=').GT.0)THEN
          READ(ARG(INDEX(ARG,'-frame=')+7:),*)FRAME
          FRAME=MIN(FRAME,NAXES(3))
          WRITE(*,'(A,I5)') ' frame: ',FRAME
        END IF
      END DO
      CALL READIMAGE(IMGFILE,(/1,1,FRAME/),
     &  (/NAXES(1),NAXES(2),FRAME/),DIMG)
      XMIN=1
      YMIN=1
      XMAX=NAXES(1)
      YMAX=NAXES(2)
      SC=DEFAULTSC
      CMIN=REAL(MINVAL(DIMG))
      CMAX=REAL(MAXVAL(DIMG))
      XLAB='x/arcsec'
      YLAB='y/arcsec'
      CLAB='Intensity'
      COLOR='gray'
      MAP='linear'
      DEVICE='?'
      TITLE=IMGFILE
      DO K=2,NARGS
        CALL GET_COMMAND_ARGUMENT(K,ARG)
        IF(INDEX(ARG,'-xmin=').GT.0)THEN
          READ(ARG(INDEX(ARG,'-xmin=')+6:),*)XMIN
        ELSE IF(INDEX(ARG,'-frame=').GT.0)THEN
        ELSE IF(INDEX(ARG,'-xmax=').GT.0)THEN
          READ(ARG(INDEX(ARG,'-xmax=')+6:),*)XMAX
        ELSE IF(INDEX(ARG,'-ymin=').GT.0)THEN
          READ(ARG(INDEX(ARG,'-ymin=')+6:),*)YMIN
        ELSE IF(INDEX(ARG,'-ymax=').GT.0)THEN
          READ(ARG(INDEX(ARG,'-ymax=')+6:),*)YMAX
        ELSE IF(INDEX(ARG,'-min=').GT.0)THEN
          READ(ARG(INDEX(ARG,'-min=')+5:),*)CMIN
        ELSE IF(INDEX(ARG,'-max=').GT.0)THEN
          READ(ARG(INDEX(ARG,'-max=')+5:),*)CMAX
        ELSE IF(INDEX(ARG,'-scale=').GT.0)THEN
          READ(ARG(INDEX(ARG,'-scale=')+7:),*)SC
        ELSE IF(INDEX(ARG,'-xlabel=').GT.0)THEN
          XLAB=ARG(INDEX(ARG,'-xlabel=')+8:)
        ELSE IF(INDEX(ARG,'-ylabel=').GT.0)THEN
          YLAB=ARG(INDEX(ARG,'-ylabel=')+8:)
        ELSE IF(INDEX(ARG,'-title=').GT.0)THEN
          TITLE=ARG(INDEX(ARG,'-title=')+7:)
        ELSE IF(INDEX(ARG,'-clabel=').GT.0)THEN
          CLAB=ARG(INDEX(ARG,'-clabel=')+8:)
        ELSE IF(INDEX(ARG,'-color=').GT.0)THEN
          COLOR=ARG(INDEX(ARG,'-color=')+7:)
        ELSE IF(INDEX(ARG,'-device=').GT.0)THEN
          DEVICE=ARG(INDEX(ARG,'-device=')+8:)
        ELSE IF(INDEX(ARG,'-map=').GT.0)THEN
          MAP=ARG(INDEX(ARG,'-map=')+5:)
        ELSE
          PRINT *,'unknown argument '//ARG
          STOP
        END IF
      END DO
      WRITE(*,'(A,I3,A,I3,A,I3,A,I3,A)')
     &  ' subimage: (',XMIN,', ',YMIN,') to (',XMAX,', ',YMAX,')'
      IF(PGOPEN(TRIM(DEVICE)) .LT. 1)THEN
        STOP
      END IF
      TR=SC*(/0.0, 1.0, 0.0, 0.0, 0.0, 1.0/)
      IF(INDEX(MAP,'linear') .GT. 0)THEN
        ITF=0
      ELSE IF(INDEX(MAP,'log') .GT. 0)THEN
        ITF=1
      ELSE IF(INDEX(MAP,'sqrt') .GT. 0)THEN
        ITF=2
      ELSE
        PRINT *,'error: unknown image transfer function.'
        ITF=0
      END IF
      CALL PGPAGE
      CALL PGSVP(0.0, 1.0, 0.0, 1.0)
      CALL PGQVP(1, VPX1, VPX2, VPY1, VPY2)
      D=MIN(VPX2-VPX1, VPY2-VPY1)/40.0
      VPX1 = VPX1 + 6.0*D
      VPX2 = VPX2 - 2.0*D
      VPY1 = VPY1 + 9.0*D
      VPY2 = VPY2 - 2.0*D
      CALL PGVSIZ(VPX1, VPX2, VPY1, VPY2)
      CALL PGWNAD(SC*REAL(XMIN-1),SC*REAL(XMAX+1),
     &  SC*REAL(YMIN-1),SC*REAL(YMAX+1))
      IF(INDEX(COLOR,'gray').GT.0)THEN
        CALL PGCTAB(GL,GR,GG,GB,2,1.0,0.5)
      ELSE IF(INDEX(COLOR,'heat').GT.0)THEN
        CALL PGCTAB(HL,HR,HG,HB,5,1.0,0.5)
      ELSE IF(INDEX(COLOR,'rainbow').GT.0)THEN
        CALL PGCTAB(RL,RR,RG,RB,9,1.0,0.5)
      ELSE
        PRINT *,'error: unknown color table.'
        CALL PGCTAB(GL,GR,GG,GB,2,1.0,0.5)
      END IF
      CALL PGSITF(ITF)
      CALL PGIMAG(REAL(DIMG),NAXES(1),NAXES(2),XMIN,XMAX,YMIN,YMAX
     &  ,CMIN,CMAX,TR)
      CALL PGSCH(1.0)
      CALL PGMTXT('T',1.0,0.5,0.5,TRIM(TITLE))
      CALL PGBOX('BCNTSI',0.0,0,'BCNTSIV',0.0,0)
      CALL PGMTXT('B',3.0,0.5,0.5,TRIM(XLAB))
      CALL PGMTXT('L',3.0,0.5,0.5,TRIM(YLAB))
      CALL PGWEDG('BI',4.0,5.0,CMIN,CMAX,TRIM(CLAB))
      CALL PGCLOS
      DEALLOCATE(DIMG)
      STOP
      END PROGRAM IMAGESC

      PROGRAM CENTER
C  Purpose:
C  ========
C  align the centroid of the given image to its center gradually. then using the
C  pixels locate on the border of the shifted image to determine the background.
C
C  Usage:
C  ======
C  center filename_input [-prefix=output] [-border=n] [-fit=fit_method]
C
C  Arguments:
C  ==========
C  filename_input - input filename.
C  prefix         - optional, prefix of output filename.
C  border         - optional, border width of rectangular signal region.
C  fit            - optional, fitting method. permitted values: -1, 0, ..., n.
C                   -1 indicates maximum value method, while 0 indicates
C                   constant average background, and positive values indicate
C                   polynomials fitting methods.
C
C  Variables:
C  ==========
C  NARGS  - Number of command arguments.
C  NAXES  - Lengths of axes of image.
C  STATUS - Status.
C  K      - Index of command arguments.
C  DIMG   - Input image to be centred and output.
C  INFILE - Input filename.
C  PREFIX - Prefix of output filename.
C  ARG    - Command argument.
      IMPLICIT NONE
      INTEGER :: NARGS,NAXES(3),STATUS,K,HW,HH,X,Y,XC,YC,SHIFT(2),BORDER
     &  ,FITN
      INTEGER, PARAMETER :: DEFAULTBORDER=2,DEFAULTFITN=-1
      DOUBLE PRECISION, ALLOCATABLE :: DIMG(:,:),DFLAG(:,:),DBG(:,:),
     &  DPARAMS(:),DTMP(:,:)
      CHARACTER(LEN=256) :: INFILE,PREFIX,ARG,BASENAME,EXTNAME
      INTERFACE
      SUBROUTINE IMAGESIZE(FILENAME,NAXES)
      INTEGER, INTENT(IN) :: NAXES(3)
      CHARACTER(LEN=*), INTENT(IN) :: FILENAME
      END SUBROUTINE IMAGESIZE
      SUBROUTINE CENTERIMAGE(NX,NY,DIMG,SHIFT,HW,HH)
      INTEGER, INTENT(IN) :: NX,NY
      INTEGER, INTENT(OUT) :: SHIFT(2),HW,HH
      DOUBLE PRECISION, INTENT(INOUT) :: DIMG(NX,NY)
      END SUBROUTINE CENTERIMAGE
      SUBROUTINE RESOLVEPATH(PATH,BASENAME,EXTNAME)
      CHARACTER(LEN=*), INTENT(IN) :: PATH
      CHARACTER(LEN=*), INTENT(OUT) :: BASENAME,EXTNAME
      END SUBROUTINE RESOLVEPATH
      SUBROUTINE READIMAGE(FILENAME,FPIXELS,LPIXELS,DIMG)
      INTEGER, INTENT(IN) :: FPIXELS(3),LPIXELS(3)
      DOUBLE PRECISION, INTENT(OUT) :: DIMG(*)
      CHARACTER(LEN=*), INTENT(IN) :: FILENAME
      END SUBROUTINE READIMAGE
      SUBROUTINE WRITEIMAGE(FILENAME,FPIXELS,LPIXELS,DIMG)
      INTEGER, INTENT(IN) :: FPIXELS(3),LPIXELS(3)
      DOUBLE PRECISION, INTENT(IN) :: DIMG(*)
      CHARACTER(LEN=*), INTENT(IN) :: FILENAME
      END SUBROUTINE WRITEIMAGE
      SUBROUTINE BGFIT2PN(NX,NY,DFLAG,DIMG,N,DBG,DB)
      INTEGER, INTENT(IN) :: NX,NY,N
      DOUBLE PRECISION, INTENT(IN) :: DFLAG(NX,NY),DIMG(NX,NY)
      DOUBLE PRECISION, INTENT(OUT) :: DBG(NX,NY),DB(*)
      END SUBROUTINE BGFIT2PN
      END INTERFACE
      STATUS=0
      NARGS=COMMAND_ARGUMENT_COUNT()
      CALL GET_COMMAND_ARGUMENT(1,ARG)
      IF(INDEX(ARG,'-help').GT.0)THEN
        PRINT *,'Usage:'
        PRINT *,'======'
        PRINT *,'center filename_input [-prefix=output] [-border=n]'//
     &    ' [-fit=fit_method]'
        PRINT *,'Arguments:'
        PRINT *,'=========='
        PRINT *,'filename_input - input filename.'
        PRINT *,'prefix         - optional, prefix of output filename.'
        PRINT *,'border         '//
     &    '- optional, border width of rectangular signal region.'
        PRINT *,'fit            '//
     &    '- optional, fitting method. permitted values: -1, 0, ..., n.'
        PRINT *,'                 '//
     &    '-1 indicates maximum value method, while 0 indicates'
        PRINT *,'                 '//
     &    'constant average background, and positive values indicate'
        PRINT *,'                 '//
     &    'polynomials fitting methods.'
        STOP
      END IF
      BORDER=DEFAULTBORDER
      FITN=DEFAULTFITN
      CALL GET_COMMAND_ARGUMENT(1,INFILE)
      PRINT *,'input: ',TRIM(INFILE)
      CALL IMAGESIZE(INFILE,NAXES)
      WRITE(*,'(A,I4,A,I4)')' size (width x height):',NAXES(1),' x',
     &  NAXES(2)
      IF(NAXES(3).GT.1)THEN
        PRINT *,'input multi-frame FITS file.'
        PRINT *,'only process the first frame.'
      END IF
      CALL RESOLVEPATH(INFILE,BASENAME,EXTNAME)
      PREFIX=TRIM(BASENAME)
      DO K=2,NARGS
        CALL GET_COMMAND_ARGUMENT(K,ARG)
        IF(INDEX(ARG,'-border=').GT.0)THEN
          READ(ARG(INDEX(ARG,'-border=')+8:),*)BORDER
        ELSE IF(INDEX(ARG,'-fit=').GT.0)THEN
          READ(ARG(INDEX(ARG,'-fit=')+6:),*)FITN
        ELSE IF(INDEX(ARG,'-prefix=').GT.0)THEN
          PREFIX=TRIM(ARG(INDEX(ARG,'-prefix=')+8:))
        ELSE
          PRINT *,'unknown argument '//ARG
          STOP
        END IF
      END DO
      ALLOCATE(DIMG(NAXES(1),NAXES(2)),STAT=STATUS)
      IF(STATUS.NE.0)THEN
        PRINT *,'out of memory.'
        RETURN
      END IF
      CALL READIMAGE(INFILE,(/1,1,1/),(/NAXES(1),NAXES(2),1/),DIMG)
      CALL CENTERIMAGE(NAXES(1),NAXES(2),DIMG,SHIFT,HW,HH)
      WRITE(*,'(A,I3.1,A,I3.1,A)')' total shift: (',SHIFT(1),', ',
     &  SHIFT(2),')'
      CALL WRITEIMAGE(TRIM(PREFIX)//'_cen.fits',(/1,1,1/),
     &  (/NAXES(1),NAXES(2),1/),DIMG)
      PRINT *,'output: '//TRIM(PREFIX)//'_cen.fits'
      ALLOCATE(DFLAG(NAXES(1),NAXES(2)),STAT=STATUS)
      IF(STATUS.NE.0)THEN
        PRINT *,'out of memory.'
        RETURN
      END IF
      XC=INT(FLOOR(0.5*REAL(1+NAXES(1))))
      YC=INT(FLOOR(0.5*REAL(1+NAXES(2))))
      DO X=1,NAXES(1)
        DO Y=1,NAXES(2)
          IF(((IABS(X-XC).LT.HW).AND.(IABS(Y-YC).LT.HH)).AND.
     &      ((IABS(X-XC).GE.HW-BORDER).OR.(IABS(Y-YC).GE.HH-BORDER)))
     &      THEN
            DFLAG(X,Y)=1.0D0
          ELSE
            DFLAG(X,Y)=0.0D0
          END IF
        END DO
      END DO
      ALLOCATE(DBG(NAXES(1),NAXES(2)),STAT=STATUS)
      IF(STATUS.NE.0)THEN
        PRINT *,'out of memory.'
        RETURN
      END IF
      ALLOCATE(DPARAMS(NINT(SUM(DFLAG))),STAT=STATUS)
      IF(STATUS.NE.0)THEN
        PRINT *,'out of memory.'
        RETURN
      END IF
      PRINT *,'before fitting.'
      CALL BGFIT2PN(NAXES(1),NAXES(2),DFLAG,DIMG,FITN,DBG,DPARAMS)
      PRINT *,'after fitting.'
      ALLOCATE(DTMP(NAXES(1),NAXES(2)),STAT=STATUS)
      IF(STATUS.NE.0)THEN
        PRINT *,'out of memory.'
        RETURN
      END IF
      DTMP=DIMG*DFLAG
      CALL WRITEIMAGE(TRIM(PREFIX)//'_cen_flag.fits',(/1,1,1/),
     &  (/NAXES(1),NAXES(2),1/),DFLAG)
      PRINT *,'background flag: '//TRIM(PREFIX)//'_cen_flag.fits'
      CALL WRITEIMAGE(TRIM(PREFIX)//'_cen_bg.fits',(/1,1,1/),
     &  (/NAXES(1),NAXES(2),1/),DTMP)
      PRINT *,'defined background: '//TRIM(PREFIX)//'_cen_bg.fits'
      CALL WRITEIMAGE(TRIM(PREFIX)//'_cen_bgfit.fits',(/1,1,1/),
     &  (/NAXES(1),NAXES(2),1/),DBG)
      PRINT *,'fitted background: '//TRIM(PREFIX)//'_cen_bgfit.fits'
      DTMP=(DIMG-DBG)*DFLAG
      CALL WRITEIMAGE(TRIM(PREFIX)//'_cen_bgres.fits',(/1,1,1/),
     &  (/NAXES(1),NAXES(2),1/),DTMP)
      PRINT *,'fitting residual: '//TRIM(PREFIX)//'_cen_bgres.fits'
      DTMP=DIMG-DBG
      DO X=1,NAXES(1)
        DO Y=1,NAXES(2)
          IF(DTMP(X,Y).LT.0.0D0)THEN
            DTMP(X,Y)=0.0D0
          END IF
        END DO
      END DO
      CALL WRITEIMAGE(TRIM(PREFIX)//'_cen_sig.fits',(/1,1,1/),
     &  (/NAXES(1),NAXES(2),1/),DTMP)
      PRINT *,'output signal: '//TRIM(PREFIX)//'_cen_sig.fits'
      DEALLOCATE(DFLAG)
      DEALLOCATE(DBG)
      DEALLOCATE(DTMP)
      DEALLOCATE(DPARAMS)
      DEALLOCATE(DIMG)
      STOP
      END PROGRAM CENTER

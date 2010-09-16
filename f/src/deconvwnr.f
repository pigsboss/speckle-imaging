      PROGRAM MAIN
C  Wiener deconvolution executable.
C
C  Usage:
C  ======
C  deconvwnr filename_g -psf=filename_h -snr=SNR -o=filename_f
C
C  Purpose:
C  ========
C  g=conv(f,h). estimates f given g and h.
C
C  Arguments:
C  ==========
C  filename_g - input file. contains matrix g.
C  filename_h - input file. contains matrix h.
C  SNR        - signal-to-noise ratio.
C  filename_f - output file. contains the estimate of matrix f.
C
C  Declaration:
C  ============
      IMPLICIT NONE
      INTEGER :: NAXES(3),NARGS,STATUS,K
      DOUBLE PRECISION :: DSNR
      DOUBLE PRECISION, ALLOCATABLE :: DF(:,:),DG(:,:),DH(:,:)
      CHARACTER(LEN=256) :: ARG,FILEG,FILEH,FILEF,BASENAME,EXTNAME
      INTERFACE
      SUBROUTINE IMAGESIZE(FILENAME,NAXES)
      INTEGER, INTENT(OUT) :: NAXES(3)
      CHARACTER(LEN=*) :: FILENAME
      END SUBROUTINE IMAGESIZE
      SUBROUTINE RESOLVEPATH(PATH,BASENAME,EXTNAME)
      CHARACTER(LEN=*), INTENT(IN) :: PATH
      CHARACTER(LEN=*), INTENT(OUT) :: BASENAME,EXTNAME
      END SUBROUTINE RESOLVEPATH
      SUBROUTINE DECONVWNR(NX,NY,DG,DF,DH,DSNR)
      INTEGER, INTENT(IN) :: NX,NY
      DOUBLE PRECISION, INTENT(IN) :: DG(NX,NY),DH(NX,NY),DSNR
      DOUBLE PRECISION, INTENT(OUT) :: DF(NX,NY)
      END SUBROUTINE DECONVWNR
      END INTERFACE
      STATUS=0
      NARGS=COMMAND_ARGUMENT_COUNT()
      CALL GET_COMMAND_ARGUMENT(1,FILEG)
      CALL GET_COMMAND_ARGUMENT(1,ARG)
      IF(INDEX(ARG,'-help').GT.0)THEN
        PRINT *,'Usage:'
        PRINT *,'======'
        PRINT *,'deconvwnr filename_g -psf=filename_h'//
     &    ' -snr=SNR -o=filename_f'
        STOP
      END IF
      CALL RESOLVEPATH(FILEG,BASENAME,EXTNAME)
      FILEF=TRIM(BASENAME)//'_wnr.fits'
      CALL IMAGESIZE(FILEG,NAXES)
      DO K=2,NARGS
        CALL GET_COMMAND_ARGUMENT(K,ARG)
        IF(INDEX(ARG,'-snr=').GT.0)THEN
          READ(ARG(INDEX(ARG,'-snr=')+5:),*)DSNR
        ELSE IF(INDEX(ARG,'-psf=').GT.0)THEN
          FILEH=ARG(INDEX(ARG,'-psf=')+5:)
        ELSE IF(INDEX(ARG,'-o=').GT.0)THEN
          FILEF=ARG(INDEX(ARG,'-o=')+3:)
        ELSE
          PRINT *,'unknown argument '//ARG
          STOP
        END IF
      END DO
      ALLOCATE(DG(NAXES(1),NAXES(2)),STAT=STATUS)
      IF(STATUS.NE.0)THEN
        PRINT *,'out of memory.'
        RETURN
      END IF
      ALLOCATE(DH(NAXES(1),NAXES(2)),STAT=STATUS)
      IF(STATUS.NE.0)THEN
        PRINT *,'out of memory.'
        RETURN
      END IF
      CALL READIMAGE(FILEG,(/1,1,1/),(/NAXES(1),NAXES(2),1/),DG)
      CALL READIMAGE(FILEH,(/1,1,1/),(/NAXES(1),NAXES(2),1/),DH)
      ALLOCATE(DF(NAXES(1),NAXES(2)),STAT=STATUS)
      IF(STATUS.NE.0)THEN
        PRINT *,'out of memory.'
        RETURN
      END IF
      CALL DECONVWNR(NAXES(1),NAXES(2),DG,DF,DH,DSNR)
      WRITE(*,'(A,ES7.1)')' SNR: ',DSNR
      CALL WRITEIMAGE(FILEF,(/1,1,1/),(/NAXES(1),NAXES(2),1/),DF)
      PRINT *,'output: '//TRIM(FILEF)
      DEALLOCATE(DG)
      DEALLOCATE(DH)
      DEALLOCATE(DF)
      STOP
      END PROGRAM MAIN

      PROGRAM SNR
C  estimate the signal-to-noise ratio of given image.
C
C  Usage:
C  ======
C  snr filename_img filename_flag [-fit=fit_method] [-prefix=output]
C
      IMPLICIT NONE
      INTEGER :: NARGS,K,FITN
      INTEGER, PARAMETER :: DEFAULTFITN=4
      DOUBLE PRECISION :: DSNR
      CHARACTER(LEN=256) :: IMGFILE,FLAGFILE,PREFIX,ARG,BASENAME,EXTNAME
      INTERFACE
      SUBROUTINE ESTSNR(IMGFILE,FLAGFILE,FITN,DSNR,PREFIX)
      INTEGER, INTENT(IN) :: FITN
      DOUBLE PRECISION, INTENT(OUT) :: DSNR
      CHARACTER(LEN=*), INTENT(IN) :: IMGFILE,FLAGFILE,PREFIX
      END SUBROUTINE ESTSNR
      SUBROUTINE RESOLVEPATH(PATH,BASENAME,EXTNAME)
      CHARACTER(LEN=*), INTENT(IN) :: PATH
      CHARACTER(LEN=*), INTENT(OUT) :: BASENAME,EXTNAME
      END SUBROUTINE RESOLVEPATH
      END INTERFACE
      NARGS=COMMAND_ARGUMENT_COUNT()
      CALL GET_COMMAND_ARGUMENT(1,ARG)
      IF(INDEX(ARG,'-help').GT.0)THEN
        PRINT *,'Usage:'
        PRINT *,'======'
        PRINT *,'snr filename_img filename_flag'//
     &    ' [-fit=fit_method] [-prefix=output]'
        STOP
      END IF
      CALL GET_COMMAND_ARGUMENT(1,IMGFILE)
      CALL GET_COMMAND_ARGUMENT(2,FLAGFILE)
      FITN=DEFAULTFITN
      CALL RESOLVEPATH(IMGFILE,BASENAME,EXTNAME)
      PREFIX=TRIM(BASENAME)
      DO K=3,NARGS
        CALL GET_COMMAND_ARGUMENT(K,ARG)
        IF(INDEX(ARG,'-fit=').GT.0)THEN
          READ(ARG(INDEX(ARG,'-fit=')+6:),*)FITN
          WRITE(*,'(A,I2)')' fitting order: ',FITN
        ELSE IF(INDEX(ARG,'-prefix=').GT.0)THEN
          PREFIX=TRIM(ARG(INDEX(ARG,'-prefix=')+8:))
        ELSE
          PRINT *,'unknown argument '//ARG
          STOP
        END IF
      END DO
      CALL ESTSNR(IMGFILE,FLAGFILE,FITN,DSNR,PREFIX)
      WRITE(*,'(A,ES9.2)')' estimated SNR: ',DSNR
      STOP
      END PROGRAM SNR

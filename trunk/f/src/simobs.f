      PROGRAM SIMOBS
C  Generate simulated observation based on real reference star observation and
C  simulated objects.
C
C  Usage:
C  ======
C  ./simobs ref_file [-range=m,n] -obj=obj_file [-bg=cts_bg] [-noise]
C    [-prefix=...]
C
      IMPLICIT NONE
      INCLUDE 'fftw3.f'
      INTEGER :: STATUS,NARGS,K,RNG(2),NAXES(3),PLANF,PLANB,INFO,X,Y
      LOGICAL :: SIMNOISE
      DOUBLE PRECISION :: DBGL
      DOUBLE PRECISION, ALLOCATABLE :: DREF(:,:),DOBJ(:,:),DOBS(:,:)
      DOUBLE COMPLEX, ALLOCATABLE :: ZIN(:,:),ZOUT(:,:),ZSPE(:,:)
      CHARACTER(LEN=256) :: REFFILE,OBJFILE,PREFIX,ARG,BASENAME,EXTNAME
      INTERFACE
      SUBROUTINE IMAGESIZE(FILENAME,NAXES)
      INTEGER, INTENT(OUT) :: NAXES(3)
      CHARACTER(LEN=*) :: FILENAME
      END SUBROUTINE IMAGESIZE
      SUBROUTINE RESOLVEPATH(PATH,BASENAME,EXTNAME)
      CHARACTER(LEN=*), INTENT(IN) :: PATH
      CHARACTER(LEN=*), INTENT(OUT) :: BASENAME,EXTNAME
      END SUBROUTINE RESOLVEPATH
      SUBROUTINE READIMAGE(FILENAME,FPIXELS,LPIXELS,DIMG)
      INTEGER, INTENT(IN) :: FPIXELS(3),LPIXELS(3)
      DOUBLE PRECISION, INTENT(OUT) :: DIMG(*)
      CHARACTER(LEN=*), INTENT(IN) :: FILENAME
      END SUBROUTINE READIMAGE
      SUBROUTINE DELETEFILE(FILENAME,STATUS)
      INTEGER, INTENT(INOUT) :: STATUS
      CHARACTER*(*), INTENT(IN) :: FILENAME
      END SUBROUTINE DELETEFILE
      SUBROUTINE APPENDIMAGE(FILENAME,FPIXELS,LPIXELS,DIMG)
      INTEGER, INTENT(IN) :: FPIXELS(3),LPIXELS(3)
      DOUBLE PRECISION, INTENT(IN) :: DIMG(*)
      CHARACTER*(*), INTENT(IN) :: FILENAME
      END SUBROUTINE APPENDIMAGE
      SUBROUTINE ZIFFTSHIFT(NX,NY,ZSP)
      INTEGER, INTENT(IN) :: NX,NY
      DOUBLE COMPLEX, INTENT(INOUT) :: ZSP(NX,NY)
      END SUBROUTINE ZIFFTSHIFT
      SUBROUTINE DFFTSHIFT(NX,NY,DIMG)
      INTEGER, INTENT(IN) :: NX,NY
      DOUBLE PRECISION, INTENT(INOUT) :: DIMG(NX,NY)
      END SUBROUTINE DFFTSHIFT
      SUBROUTINE ZFFTSHIFT(NX,NY,ZSP)
      INTEGER, INTENT(IN) :: NX,NY
      DOUBLE COMPLEX, INTENT(INOUT) :: ZSP(NX,NY)
      END SUBROUTINE ZFFTSHIFT
      SUBROUTINE INIT_RANDOM_SEED()
      END SUBROUTINE INIT_RANDOM_SEED
      FUNCTION NORMRND(NX,NY,DMU,DSIGMA)
      INTEGER,INTENT(IN):: NX,NY
      DOUBLE PRECISION,INTENT(IN) :: DMU(NX,NY),DSIGMA(NX,NY)
      DOUBLE PRECISION :: NORMRND(NX,NY)
      END FUNCTION NORMRND
      END INTERFACE
C  Statements:
C  ===========
      STATUS=0
C    Resolve the command line options:
C    =================================
      CALL GET_COMMAND_ARGUMENT(1,ARG)
      IF(INDEX(ARG,'-help').GT.0)THEN
        PRINT *,'Usage:'
        PRINT *,'======'
        PRINT *,'simobs ref_file [-range=m,n] -obj=obj_file'
        PRINT *,'  [-prefix=prefix_of_output] [-bg=cts_bg] [-noise]'
        STOP
      END IF
      NARGS=COMMAND_ARGUMENT_COUNT()
      CALL GET_COMMAND_ARGUMENT(1,REFFILE)
      CALL IMAGESIZE(REFFILE,NAXES)
      CALL RESOLVEPATH(REFFILE,BASENAME,EXTNAME)
      SIMNOISE = .FALSE.
      DBGL = 0.0D0
      OBJFILE=''
      RNG(1)=1
      RNG(2)=NAXES(3)
      PREFIX=TRIM(BASENAME)
      DO K=2,NARGS
        CALL GET_COMMAND_ARGUMENT(K,ARG)
        IF(INDEX(ARG,'-range=').GT.0)THEN
          READ(ARG(INDEX(ARG,'-range=')+7:),*)RNG(1),RNG(2)
        ELSE IF(INDEX(ARG,'-prefix=').GT.0)THEN
          PREFIX=ARG(INDEX(ARG,'-prefix=')+8:)
        ELSE IF(INDEX(ARG,'-obj=').GT.0)THEN
          OBJFILE=ARG(INDEX(ARG,'-obj=')+5:)
        ELSE IF(INDEX(ARG,'-noise').EQ.1)THEN
          SIMNOISE=.TRUE.
        ELSE IF(INDEX(ARG,'-bg=').EQ.1)THEN
          READ(ARG(INDEX(ARG,'-bg=')+4:),*)DBGL
        ELSE
          PRINT *,'Unknown argument '//TRIM(ARG)
          STOP
        END IF
      END DO
      IF(LEN_TRIM(OBJFILE) .LT. 1)THEN
        PRINT *,'the object must be specified.'
        STOP
      END IF
      WRITE(*,'(A,I5,A,I5)')' real reference: '//TRIM(REFFILE)//
     &  ', from frame ',RNG(1),' to frame ',RNG(2)
      WRITE(*,'(A)')' simulated object: '//TRIM(OBJFILE)
      WRITE(*,'(A,I3,A,I3)')' image size (width x height): ',
     &  NAXES(1),' x ',NAXES(2)
C    Allocation:
C    ===========
      ALLOCATE(DREF(NAXES(1),NAXES(2)),STAT=STATUS)
      IF(STATUS.NE.0)THEN
        PRINT *,'error: out of memory.'
        STOP
      END IF
      ALLOCATE(DOBJ(NAXES(1),NAXES(2)),STAT=STATUS)
      IF(STATUS.NE.0)THEN
        PRINT *,'error: out of memory.'
        STOP
      END IF
      ALLOCATE(DOBS(NAXES(1),NAXES(2)),STAT=STATUS)
      IF(STATUS.NE.0)THEN
        PRINT *,'error: out of memory.'
        STOP
      END IF
      ALLOCATE(ZIN(NAXES(1),NAXES(2)),STAT=STATUS)
      IF(STATUS.NE.0)THEN
        PRINT *,'error: out of memory.'
        STOP
      END IF
      ALLOCATE(ZOUT(NAXES(1),NAXES(2)),STAT=STATUS)
      IF(STATUS.NE.0)THEN
        PRINT *,'error: out of memory.'
        STOP
      END IF
      ALLOCATE(ZSPE(NAXES(1),NAXES(2)),STAT=STATUS)
      IF(STATUS.NE.0)THEN
        PRINT *,'error: out of memory.'
        STOP
      END IF
C    Generate the simulated observations:
C    ====================================
      CALL DELETEFILE(TRIM(PREFIX)//'_sim_ref.fits',STATUS)
      CALL DELETEFILE(TRIM(PREFIX)//'_sim_obs.fits',STATUS)
      CALL INIT_RANDOM_SEED()
      CALL DFFTW_INIT_THREADS(INFO)
      IF (INFO .EQ. 0)THEN
        PRINT *,'error: DFFTW_INIT_THREADS failed.'
        STOP
      END IF
      CALL DFFTW_PLAN_WITH_NTHREADS(2)
      CALL DFFTW_IMPORT_SYSTEM_WISDOM(INFO)
      IF (INFO .EQ. 0)THEN
        PRINT *,'warning: DFFTW_IMPORT_SYSTEM_WISDOM failed.'
      END IF
      CALL DFFTW_PLAN_DFT_2D(PLANF,NAXES(1),NAXES(2),ZIN,ZOUT,-1,
     &  FFTW_MEASURE+FFTW_DESTROY_INPUT)
      CALL DFFTW_PLAN_DFT_2D(PLANB,NAXES(1),NAXES(2),ZIN,ZOUT,1,
     &  FFTW_MEASURE+FFTW_DESTROY_INPUT)
      CALL READIMAGE(OBJFILE,(/1,1,1/),(/NAXES(1),NAXES(2),1/),DOBJ)
      ZIN=DCMPLX(DOBJ)
      CALL ZIFFTSHIFT(NAXES(1),NAXES(2),ZIN)
      CALL DFFTW_EXECUTE_DFT(PLANF,ZIN,ZOUT)
      ZSPE=ZOUT/DSQRT(DBLE(NAXES(1)*NAXES(2)))
      DO K=RNG(1),RNG(2)
        CALL READIMAGE(REFFILE,(/1,1,K/),(/NAXES(1),NAXES(2),K/),DREF)
        IF(MOD(K,2) .EQ. 0) THEN
          CALL APPENDIMAGE(TRIM(PREFIX)//'_sim_ref.fits',(/1,1,K/2/),
     &      (/NAXES(1),NAXES(2),K/2/),DREF)
          WRITE(*,'(A,I5,A,I5)')' append frame ',K,
     &      ' to simulated reference as frame ',K/2
        ELSE
          ZIN=DCMPLX(DREF/SUM(DREF))
          CALL ZIFFTSHIFT(NAXES(1),NAXES(2),ZIN)
          CALL DFFTW_EXECUTE_DFT(PLANF,ZIN,ZOUT)
          ZIN=ZSPE*ZOUT
          CALL DFFTW_EXECUTE_DFT(PLANB,ZIN,ZOUT)
          CALL ZFFTSHIFT(NAXES(1),NAXES(2),ZOUT)
          DOBS=ZABS(ZOUT)/DSQRT(DBLE(NAXES(1)*NAXES(2)))
          IF (SIMNOISE) THEN
            DOBS=DABS(NORMRND(NAXES(1),NAXES(2),DOBS,DSQRT(DOBS)))
          END IF
          CALL APPENDIMAGE(TRIM(PREFIX)//'_sim_obs.fits',
     &      (/1,1,(K+1)/2/),(/NAXES(1),NAXES(2),(K+1)/2/),DOBS)
          WRITE(*,'(A,I5,A,I5)')' append frame ',K,
     &      ' to simulated observation as frame ',(K+1)/2
        END IF
      END DO
      DEALLOCATE(DREF)
      DEALLOCATE(DOBJ)
      DEALLOCATE(DOBS)
      DEALLOCATE(ZIN)
      DEALLOCATE(ZOUT)
      DEALLOCATE(ZSPE)
      CALL DFFTW_DESTROY_PLAN(PLANF)
      CALL DFFTW_DESTROY_PLAN(PLANB)
      STOP
      END PROGRAM SIMOBS

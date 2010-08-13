      PROGRAM MAIN
      IMPLICIT NONE
      INCLUDE 'fftw3.f'
      INTEGER :: K,KOBS,KREF,FPIXELS(3),LPIXELS(3),M,N,NFRAMES,INFO
     &  ,STATUS,I,J
      INTEGER*8 :: PLANF,PLANB
      DOUBLE PRECISION :: DTMP,X,Y,XC,YC,DR
      DOUBLE PRECISION, ALLOCATABLE :: DOBJ(:,:),DREF(:,:),DOBS(:,:),
     &  DPUP(:,:),DPSF(:,:)
      DOUBLE COMPLEX, ALLOCATABLE :: ZIN(:,:),ZOUT(:,:),ZOBJ(:,:),
     &  ZOTF(:,:)
      CHARACTER(LEN=256) :: ARG,FILEREF,PREFIX
      STATUS=0
      CALL GETARG(1,FILEREF)
      CALL GETARG(2,ARG)
      READ(ARG,*) FPIXELS(1),FPIXELS(2),FPIXELS(3)
      CALL GETARG(3,ARG)
      READ(ARG,*) LPIXELS(1),LPIXELS(2),LPIXELS(3)
      CALL GETARG(4,ARG)
      READ(ARG,*) DR
      CALL GETARG(5,PREFIX)
      CALL DELETEFILE(TRIM(PREFIX)//'_OBJ.FITS',STATUS)
      CALL DELETEFILE(TRIM(PREFIX)//'_OBS.FITS',STATUS)
      CALL DELETEFILE(TRIM(PREFIX)//'_REF.FITS',STATUS)
      M=LPIXELS(1)-FPIXELS(1)+1
      N=LPIXELS(2)-FPIXELS(2)+1
      ALLOCATE(DOBJ(M,N))
      ALLOCATE(DREF(M,N))
      ALLOCATE(DOBS(M,N))
      ALLOCATE(DPUP(M,N))
      ALLOCATE(DPSF(M,N))
      NFRAMES=LPIXELS(3)-FPIXELS(3)+1
      ALLOCATE(ZIN(M,N))
      ALLOCATE(ZOUT(M,N))
      ALLOCATE(ZOBJ(M,N))
      ALLOCATE(ZOTF(M,N))
      CALL DFFTW_INIT_THREADS(INFO)
      IF (INFO .EQ. 0)THEN
        PRINT *,'DFFTW_INIT_THREADS failed.'
        RETURN
      END IF
      CALL DFFTW_PLAN_WITH_NTHREADS(2)
      CALL DFFTW_IMPORT_SYSTEM_WISDOM(INFO)
      IF (INFO .EQ. 0)THEN
        PRINT *,'DFFTW_IMPORT_SYSTEM_WISDOM failed.'
      END IF
      PRINT *,'Start planning.'
      CALL DFFTW_PLAN_DFT_2D(PLANF,M,N,ZIN,ZOUT,-1,
     &  FFTW_MEASURE+FFTW_DESTROY_INPUT)
      CALL DFFTW_PLAN_DFT_2D(PLANB,M,N,ZIN,ZOUT,1,
     &  FFTW_MEASURE+FFTW_DESTROY_INPUT)
      PRINT *,'Finished planning.'
C  Setup the source:
      DOBJ(1,1)=100.D0
      DOBJ(10,10)=50.D0
      XC=DBLE(N)/2
      YC=DBLE(M)/2
      DO I=1,M
        DO J=1,N
          X=DBLE(J)
          Y=DBLE(I)
          IF ((X-XC)*(X-XC)+(Y-YC)*(Y-YC).GT.DR*DR)THEN
            DPUP(I,J)=DBLE(0)
          ELSE
            DPUP(I,J)=DBLE(1)
          END IF
        END DO
      END DO
      CALL DIFFTSHIFT(M,N,DPUP)
      ZIN=CMPLX(DPUP)
      CALL DFFTW_EXECUTE_DFT(PLANF,ZIN,ZOUT)
      DPSF=DREAL(ZOUT*CONJG(ZOUT))
      CALL DFFTSHIFT(M,N,DPSF)
      CALL WRITEIMAGE(TRIM(PREFIX)//'_PSF.FITS',
     &  (/1,1,1/),(/M,N,1/),DPSF)
      CALL DIFFTSHIFT(M,N,DPSF)
      ZIN=CMPLX(DPSF)
      CALL DFFTW_EXECUTE_DFT(PLANF,ZIN,ZOUT)
      ZOTF=ZOUT
      CALL ZFFTSHIFT(M,N,ZOTF)
      CALL WRITEIMAGE(TRIM(PREFIX)//'_MTF.FITS',
     &  (/1,1,1/),(/M,N,1/),ZABS(ZOTF))
      CALL ZIFFTSHIFT(M,N,ZOTF)
      ZIN=CMPLX(DOBJ)
      CALL DFFTW_EXECUTE_DFT(PLANF,ZIN,ZOUT)
      ZOBJ=ZOUT
      ZIN=ZOBJ*ZOTF
      CALL DFFTW_EXECUTE_DFT(PLANB,ZIN,ZOUT)
      CALL ZFFTSHIFT(M,N,ZOUT)
      ZOUT=ZOUT/SUM(ZABS(ZOUT))*DBLE(M*N)
      CALL WRITEIMAGE(TRIM(PREFIX)//'_OBJ.FITS',
     &  (/1,1,1/),(/M,N,1/),DREAL(ZOUT))
      KOBS=0
      KREF=0
      DO K=1,NFRAMES
        CALL READIMAGE(FILEREF,(/FPIXELS(1),FPIXELS(2),K/),
     &    (/LPIXELS(1),LPIXELS(2),K/),DREF)
        DREF=DREF/SUM(DREF)*DBLE(M*N)
        DREF=DREF-DBLE(1.1)
        DO I=1,M
          DO J=1,N
            IF (DREF(I,J) .LT. DBLE(0))THEN
              DREF(I,J)=DBLE(0)
            END IF
          END DO
        END DO
        IF (MOD(K,2) .EQ. 0)THEN
          KOBS=KOBS+1
          ZIN=CMPLX(DREF)
          CALL ZIFFTSHIFT(M,N,ZIN)
          CALL DFFTW_EXECUTE_DFT(PLANF,ZIN,ZOUT)
          ZIN=ZOUT*ZOBJ*ZOTF
          CALL DFFTW_EXECUTE_DFT(PLANB,ZIN,ZOUT)
          DOBS=DREAL(ZOUT)/DBLE(M*N)
          CALL DFFTSHIFT(M,N,DOBS)
          CALL APPENDIMAGE(TRIM(PREFIX)//'_OBS.FITS',
     &      (/1,1,KOBS/),(/M,N,KOBS/),DOBS)
        ELSE
          KREF=KREF+1
          CALL APPENDIMAGE(TRIM(PREFIX)//'_REF.FITS',
     &      (/1,1,KREF/),(/M,N,KREF/),DREF)
        END IF
      END DO
      DEALLOCATE(ZIN)
      DEALLOCATE(ZOUT)
      DEALLOCATE(ZOBJ)
      DEALLOCATE(DOBS)
      DEALLOCATE(DREF)
      DEALLOCATE(DOBJ)
      DEALLOCATE(DPUP)
      DEALLOCATE(DPSF)
      DEALLOCATE(ZOTF)
      CALL DFFTW_DESTROY_PLAN(PLANF)
      CALL DFFTW_DESTROY_PLAN(PLANB)
      STOP
      END PROGRAM MAIN
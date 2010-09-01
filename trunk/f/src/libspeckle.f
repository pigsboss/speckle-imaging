      SUBROUTINE SHIFTADD(INFILE,NRNG,RNG,PREFIX)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NRNG,RNG(2,NRNG)
      INTEGER, PARAMETER :: BUFFERSIZE=64
      INTEGER :: STATUS,K,L,L1,L2,NR,LBUFFER,NBUFS,NB,NAXES(3),
     &  XM,YM,XC,YC,NFRAMES
      DOUBLE PRECISION, ALLOCATABLE :: DBUF(:,:,:),DIMG(:,:)
      CHARACTER(LEN=*), INTENT(IN) :: INFILE,PREFIX
      INTERFACE
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
      SUBROUTINE IMAGESIZE(FILENAME,NAXES)
      INTEGER, INTENT(OUT) :: NAXES(3)
      CHARACTER(LEN=*) :: FILENAME
      END SUBROUTINE IMAGESIZE
      END INTERFACE
      STATUS=0
      WRITE(*,'(A,I3,A)')' buffer size: ',BUFFERSIZE,'MB'
      CALL IMAGESIZE(INFILE,NAXES)
      XC=INT(FLOOR(0.5*REAL(1+NAXES(1))))
      YC=INT(FLOOR(0.5*REAL(1+NAXES(2))))
      WRITE(*,'(A,I3,A,I3)')' image size (width x height): ',
     &  NAXES(1),' x ',NAXES(2)
      LBUFFER=NINT(DBLE(BUFFERSIZE*1024*1024)/DBLE(NAXES(1)*NAXES(2)*8))
      WRITE(*,'(A,I4,A)')' buffer length: ',LBUFFER,' frames'
      DO NR=1,NRNG
        WRITE(*,'(A,I3,A,I5,A,I5)')' range ',NR,': from ',
     &    RNG(1,NR),' to ',RNG(2,NR)
      END DO
      ALLOCATE(DBUF(NAXES(1),NAXES(2),LBUFFER),STAT=STATUS)
      IF(STATUS.NE.0)THEN
        PRINT *,'out of memory.'
        RETURN
      END IF
      ALLOCATE(DIMG(NAXES(1),NAXES(2)),STAT=STATUS)
      IF(STATUS.NE.0)THEN
        PRINT *,'out of memory.'
        RETURN
      END IF
      DIMG=0.0D0
      NFRAMES=0
      DO NR=1,NRNG
        NBUFS=CEILING(DBLE(RNG(2,NR)+1-RNG(1,NR))/DBLE(LBUFFER))
        DO NB=1,NBUFS
          L1=RNG(1,NR)+(NB-1)*LBUFFER
          L2=MIN(RNG(1,NR)+NB*LBUFFER-1,RNG(2,NR))
          CALL READIMAGE(INFILE,
     &     (/1,1,L1/),(/NAXES(1),NAXES(2),L2/),DBUF)
          DO L=1,L2+1-L1
            NFRAMES=NFRAMES+1
            XM=MAXLOC(MAXVAL(DBUF(:,:,L),2),1)
            YM=MAXLOC(MAXVAL(DBUF(:,:,L),1),2)
            WRITE(*,'(A,I5,A,I3,A,I3,A)')' maximum location ',
     &        NFRAMES,': (',XM,', ',YM,')'
            DIMG=DIMG+EOSHIFT(EOSHIFT(DBUF(:,:,L),
     &        YM-YC,0.0D0,2),XM-XC,0.0D0,1)
          END DO
        END DO
      END DO
      DIMG=DIMG/DBLE(NFRAMES)
      DEALLOCATE(DBUF)
      CALL WRITEIMAGE(TRIM(PREFIX)//'_ssa.fits',(/1,1,1/),
     &  (/NAXES(1),NAXES(2),1/),DIMG)
      PRINT *,'output: '//TRIM(PREFIX)//'_ssa.fits'
      DEALLOCATE(DIMG)
      RETURN
      END SUBROUTINE SHIFTADD
C ******************************************************************************
      SUBROUTINE CENTERIMAGE(NX,NY,DIMG,SHIFT,HW,HH)
C  Variables:
C  ==========
C  NX,NY         - Size of dimension along x axis and y axis.
C  K,KMAX        - Index of iteration, maximum of number of iterations.
C  XC,YC         - Centre of image.
C  DIMG          - Double precision valued image.
C  DXC,DYC       - Centroid of image.
C  HW,HH,DFLAG - Half width, half height, flag.
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NX,NY
      INTEGER, INTENT(OUT) :: SHIFT(2),HW,HH
      INTEGER :: K,XC,YC,X,Y
      INTEGER, PARAMETER :: KMAX=20
      DOUBLE PRECISION, INTENT(INOUT) :: DIMG(NX,NY)
      DOUBLE PRECISION :: DXC,DYC,DFLAG(NX,NY)
      INTERFACE
      SUBROUTINE GETCENTROID(NX,NY,DIMG,DX,DY)
      INTEGER, INTENT(IN) :: NX,NY
      DOUBLE PRECISION, INTENT(IN) :: DIMG(NX,NY)
      DOUBLE PRECISION, INTENT(OUT) :: DX,DY
      END SUBROUTINE GETCENTROID
      END INTERFACE
      XC=INT(FLOOR(0.5*REAL(1+NX)))
      YC=INT(FLOOR(0.5*REAL(1+NY)))
      HW=MIN(NX-XC,XC-1)
      HH=MIN(NY-YC,YC-1)
      SHIFT=0
      DO K=1,KMAX
        CALL GETCENTROID(NX,NY,DIMG,DXC,DYC)
        WRITE(*,'(A,I2,A,F5.1,A,F5.1,A)')
     &    ' loop',K,', centroid: (',DXC,', ',DYC,')'
        IF((NINT(DXC).EQ.XC).AND.(NINT(DYC).EQ.YC))THEN
          EXIT
        ELSE
          HW=MIN(XC+HW-NINT(DXC),NINT(DXC)-XC+HW)
          HH=MIN(YC+HH-NINT(DYC),NINT(DYC)-YC+HH)
          DO X=1,NX
            DO Y=1,NY
              IF((IABS(X-NINT(DXC)).GE.HW).OR.
     &          (IABS(Y-NINT(DYC)).GE.HH))THEN
                DFLAG(X,Y)=0.0D0
              ELSE
                DFLAG(X,Y)=1.0D0
              END IF
            END DO
          END DO
          DIMG=DIMG*DFLAG
          SHIFT=SHIFT+(/NINT(DXC)-XC,NINT(DYC)-YC/)
          DIMG=EOSHIFT(EOSHIFT(DIMG,NINT(DXC)-XC,0.0D0,1),
     &     NINT(DYC)-YC,0.0D0,2)
        END IF
      END DO
      RETURN
      END SUBROUTINE CENTERIMAGE
C ******************************************************************************
      SUBROUTINE GETCENTROID(NX,NY,DIMG,DX,DY)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NX,NY
      INTEGER :: X,Y
      DOUBLE PRECISION, INTENT(IN) :: DIMG(NX,NY)
      DOUBLE PRECISION, INTENT(OUT) :: DX,DY
      DOUBLE PRECISION :: DXI(NX,NY),DYI(NX,NY)
      DO X=1,NX
        DO Y=1,NY
          DXI(X,Y)=DBLE(X)
          DYI(X,Y)=DBLE(Y)
        END DO
      END DO
      DX=SUM(DXI*DIMG)/SUM(DIMG)
      DY=SUM(DYI*DIMG)/SUM(DIMG)
      RETURN
      END SUBROUTINE
C ******************************************************************************
      SUBROUTINE IMAGESIZE(IMGFILE,NAXES)
      IMPLICIT NONE
      INTEGER :: STATUS,UNIT,RWMODE,BLOCKSIZE,NAXIS
      INTEGER, INTENT(OUT) :: NAXES(3)
      CHARACTER*(*), INTENT(IN) :: IMGFILE
      STATUS=0
      CALL FTGIOU(UNIT,STATUS)
      RWMODE=0
      CALL FTOPEN(UNIT,IMGFILE,RWMODE,BLOCKSIZE,STATUS)
      CALL FTGIDM(UNIT,NAXIS,STATUS)
      IF(NAXIS.NE.3)THEN
        PRINT *,'dimensional fault: the image must have 3 axes.'
        RETURN
      END IF
      CALL FTGISZ(UNIT,3,NAXES,STATUS)
      CALL FTCLOS(UNIT, STATUS)
      CALL FTFIOU(UNIT, STATUS)
      IF(STATUS .GT. 0)CALL PRINTERROR(STATUS)
      RETURN
      END SUBROUTINE IMAGESIZE
C ******************************************************************************
      SUBROUTINE SUBIMAGE(IMGFILE,FPIXELS,LPIXELS,SUBFILE,BUFFERSIZE)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: FPIXELS(3),LPIXELS(3)
      CHARACTER*(*), INTENT(IN) :: IMGFILE,SUBFILE
      INTEGER, OPTIONAL, INTENT(INOUT) :: BUFFERSIZE
      INTEGER :: STATUS,K,M,N,NPIXELS,NFRAMES,LBUFFER,NBUFFER
      INTEGER :: L,L1,L2
      DOUBLE PRECISION, ALLOCATABLE :: DBUF(:,:,:)
      INTERFACE
      SUBROUTINE DELETEFILE(FILENAME,STATUS)
      INTEGER, INTENT(INOUT) :: STATUS
      CHARACTER*(*), INTENT(IN) :: FILENAME
      END SUBROUTINE DELETEFILE
      SUBROUTINE READIMAGE(FILENAME,FPIXELS,LPIXELS,DIMG)
      INTEGER, INTENT(IN) :: FPIXELS(3),LPIXELS(3)
      DOUBLE PRECISION, INTENT(OUT) :: DIMG(*)
      CHARACTER*(*), INTENT(IN) :: FILENAME
      END SUBROUTINE READIMAGE
      SUBROUTINE APPENDIMAGE(FILENAME,FPIXELS,LPIXELS,DIMG)
      INTEGER, INTENT(IN) :: FPIXELS(3),LPIXELS(3)
      DOUBLE PRECISION, INTENT(IN) :: DIMG(*)
      CHARACTER*(*), INTENT(IN) :: FILENAME
      END SUBROUTINE APPENDIMAGE
      SUBROUTINE PRINTERROR(STATUS)
      INTEGER, INTENT(IN) :: STATUS
      END SUBROUTINE PRINTERROR
      END INTERFACE
      STATUS=0
      CALL DELETEFILE(SUBFILE,STATUS)
      M=LPIXELS(1)-FPIXELS(1)+1
      N=LPIXELS(2)-FPIXELS(2)+1
      NPIXELS=M*N
      NFRAMES=LPIXELS(3)-FPIXELS(3)+1
      IF (PRESENT(BUFFERSIZE))THEN
        LBUFFER=INT(FLOOR(DBLE(BUFFERSIZE*1024*1024/8)/DBLE(NPIXELS)))
      ELSE
        LBUFFER=INT(FLOOR(DBLE(32*1024*1024/8)/DBLE(NPIXELS)))
      END IF
      NBUFFER=INT(CEILING(DBLE(NFRAMES)/DBLE(LBUFFER)))
      ALLOCATE(DBUF(M,N,LBUFFER),STAT=STATUS)
      IF(STATUS.NE.0)THEN
        PRINT *,'out of memory.'
        RETURN
      END IF
      DO K=1,NBUFFER
        L1=(K-1)*LBUFFER+1
        L2=MIN(K*LBUFFER,NFRAMES)
C       PRINT *,'READIMAGE from',L1,'to',L2
        CALL READIMAGE(IMGFILE,(/FPIXELS(1),FPIXELS(2),L1/),
     &    (/LPIXELS(1),LPIXELS(2),L2/),DBUF)
        CALL APPENDIMAGE(SUBFILE,(/1,1,L1/),(/M,N,L2/),DBUF)
      END DO
      DEALLOCATE(DBUF)
      IF (STATUS .GT. 0)CALL PRINTERROR(STATUS)
      RETURN
      END SUBROUTINE SUBIMAGE
C ******************************************************************************
      SUBROUTINE FLATFIELDCORRECT(TARGFILE,FLATFILE,DARKFILE,
     &  FPIXELS,LPIXELS,CRCTFILE)
C  Flat-field correction subroutine.
C
C  Usage:
C  =====
C  Input direct observation results, flat field file, as well as dark current
C  file. Return flat-field corrected observation results as output.
C
C                 I_original - I_dark
C  I_corrected = ---------------------
C                   I_flat - I_dark
C
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: FPIXELS(3),LPIXELS(3)
      INTEGER :: STATUS,UNIT,RWMODE,BLOCKSIZE,GROUP
      INTEGER :: K,M,N,NPIXELS,NFRAMES,NAXES(3),NFOUND
      DOUBLE PRECISION, ALLOCATABLE :: DFLAT(:),DDARK(:),DTARG(:),
     &  DCRCT(:)
      CHARACTER*(*), INTENT(IN) :: TARGFILE,FLATFILE,DARKFILE,CRCTFILE
      CHARACTER(LEN=256) :: KEYVAL,COMMENT
      INTERFACE
      SUBROUTINE AVERAGE(AVGFILE,FPIXELS,LPIXELS,DAVG,BUFFERSIZE)
      INTEGER, INTENT(IN) :: FPIXELS(3),LPIXELS(3)
      INTEGER, INTENT(INOUT), OPTIONAL :: BUFFERSIZE
      DOUBLE PRECISION :: DAVG(LPIXELS(1)+1-FPIXELS(1),
     &  LPIXELS(2)+1-FPIXELS(2))
      CHARACTER*(*), INTENT(IN) :: AVGFILE
      END SUBROUTINE AVERAGE
      END INTERFACE
      STATUS=0
      CALL FTGIOU(UNIT,STATUS)
      RWMODE=0
      CALL FTOPEN(UNIT,TARGFILE,RWMODE,BLOCKSIZE,STATUS)
      CALL FTGKNJ(UNIT,'NAXIS',1,3,NAXES,NFOUND,STATUS)
      IF (NFOUND .NE. 3)THEN
        PRINT *,'FLATFIELDCORRECT failed to read the NAXIS keywords'//
     &    ' from '//TRIM(TARGFILE)
        RETURN
      ENDIF
      IF ((NAXES(1).LT.LPIXELS(1)) .OR.
     &  (NAXES(2).LT.LPIXELS(2)))THEN
        PRINT *,'Specified subimage size exceeds boundary.'
        RETURN
      END IF
      M=NAXES(1)
      N=NAXES(2)
      CALL FTCLOS(UNIT, STATUS)
      CALL FTOPEN(UNIT,FLATFILE,RWMODE,BLOCKSIZE,STATUS)
      CALL FTGKNJ(UNIT,'NAXIS',1,3,NAXES,NFOUND,STATUS)
      IF (NFOUND .NE. 3)THEN
        PRINT *,'FLATFIELDCORRECT failed to read the NAXIS keywords'//
     &    ' from '//TRIM(FLATFILE)
        RETURN
      ENDIF
      IF ((NAXES(1).LT.LPIXELS(1)) .OR.
     &  (NAXES(2).LT.LPIXELS(2)))THEN
        PRINT *,'Specified subimage size exceeds boundary.'
        RETURN
      END IF
      IF ((NAXES(1).NE.M) .OR. (NAXES(2).NE.N))THEN
        PRINT *,TRIM(TARGFILE)//' and '//TRIM(FLATFILE)//
     &    ' differ in size.'
        RETURN
      END IF
      CALL FTCLOS(UNIT, STATUS)
      ALLOCATE(DFLAT(NAXES(1)*NAXES(2)),STAT=STATUS)
      IF(STATUS.NE.0)THEN
        PRINT *,'out of memory.'
        RETURN
      END IF
      CALL AVERAGE(FLATFILE,(/FPIXELS(1),FPIXELS(2),1/),
     &  (/LPIXELS(1),LPIXELS(2),NAXES(3)/),DFLAT)
      CALL FTOPEN(UNIT,DARKFILE,RWMODE,BLOCKSIZE,STATUS)
      CALL FTGKNJ(UNIT,'NAXIS',1,3,NAXES,NFOUND,STATUS)
      IF (NFOUND .NE. 3)THEN
        PRINT *,'FLATFIELDCORRECT failed to read the NAXIS keywords'//
     &    ' from '//TRIM(DARKFILE)
        RETURN
      ENDIF
      IF ((NAXES(1).LT.LPIXELS(1)) .OR.
     &  (NAXES(2).LT.LPIXELS(2)))THEN
        PRINT *,'Specified subimage size exceeds boundary.'
        RETURN
      END IF
      IF ((NAXES(1).NE.M) .OR. (NAXES(2).NE.N))THEN
        PRINT *,TRIM(TARGFILE)//' and '//TRIM(DARKFILE)//
     &    ' differ in size.'
        RETURN
      END IF
      CALL FTCLOS(UNIT, STATUS)
      ALLOCATE(DDARK(NAXES(1)*NAXES(2)),STAT=STATUS)
      IF(STATUS.NE.0)THEN
        PRINT *,'out of memory.'
        RETURN
      END IF
      CALL AVERAGE(DARKFILE,(/FPIXELS(1),FPIXELS(2),1/),
     &  (/LPIXELS(1),LPIXELS(2),NAXES(3)/),DDARK)
      IF (MAXVAL(DDARK) .GE. MINVAL(DFLAT))THEN
        PRINT *,'Maximum of dark current larger than or equal'
     &    //' to minimum of flat field.'
        RETURN
      END IF
      M=LPIXELS(1)-FPIXELS(1)+1
      N=LPIXELS(2)-FPIXELS(2)+1
      NPIXELS=M*N
      NFRAMES=LPIXELS(3)-FPIXELS(3)+1
      ALLOCATE(DTARG(NPIXELS))
      ALLOCATE(DCRCT(NPIXELS))
      CALL DELETEFILE(CRCTFILE,STATUS)
      DO K=FPIXELS(3),LPIXELS(3)
        CALL READIMAGE(TARGFILE,(/FPIXELS(1),FPIXELS(2),K/),
     &    (/LPIXELS(1),LPIXELS(2),K/),DTARG)
        DCRCT=(DTARG-DDARK)/(DFLAT-DDARK)
        CALL APPENDIMAGE(CRCTFILE,(/1,1,K/),(/M,N,K/),DCRCT)
      END DO
      DEALLOCATE(DTARG)
      DEALLOCATE(DFLAT)
      DEALLOCATE(DDARK)
      DEALLOCATE(DCRCT)
      CALL FTFIOU(UNIT, STATUS)
      IF (STATUS .GT. 0)CALL PRINTERROR(STATUS)
      RETURN
      END SUBROUTINE FLATFIELDCORRECT
C ******************************************************************************
      SUBROUTINE PRINTERROR(STATUS)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: STATUS
      CHARACTER :: ERRTEXT*30,ERRMESSAGE*80
      IF (STATUS .LE. 0)RETURN
      CALL FTGERR(STATUS,ERRTEXT)
      PRINT *,'FITSIO ERROR STATUS =',STATUS,': ',ERRTEXT
      CALL FTGMSG(ERRMESSAGE)
      DO WHILE (ERRMESSAGE .NE. ' ')
        PRINT *,ERRMESSAGE
        CALL FTGMSG(ERRMESSAGE)
      ENDDO
      RETURN
      END SUBROUTINE PRINTERROR
C ******************************************************************************
      SUBROUTINE DELETEFILE(FILENAME,STATUS)
      IMPLICIT NONE
      INTEGER, INTENT(INOUT) :: STATUS
      INTEGER :: UNIT,BLOCKSIZE
      CHARACTER*(*), INTENT(IN) :: FILENAME
      INTERFACE
      SUBROUTINE PRINTERROR(STATUS)
      INTEGER, INTENT(IN) :: STATUS
      END SUBROUTINE PRINTERROR
      END INTERFACE
      IF (STATUS .GT. 0)RETURN
      CALL FTGIOU(UNIT,STATUS)
      CALL FTOPEN(UNIT,FILENAME,1,BLOCKSIZE,STATUS)
      IF (STATUS .EQ. 0)THEN
        CALL FTDELT(UNIT,STATUS)
      ELSE IF (STATUS .EQ. 103)THEN
        STATUS=0
        CALL FTCMSG
      ELSE
        STATUS=0
        CALL FTCMSG
        CALL FTDELT(UNIT,STATUS)
      END IF
      CALL FTFIOU(UNIT, STATUS)
      IF (STATUS .GT. 0)CALL PRINTERROR(STATUS)
      RETURN
      END SUBROUTINE DELETEFILE
C ******************************************************************************
      SUBROUTINE WRITEIMAGE(FILENAME,FPIXELS,LPIXELS,ARRAY)
      IMPLICIT NONE
      INTEGER :: STATUS,UNIT,BLOCKSIZE,BITPIX,NAXIS,GROUP
      INTEGER, INTENT(IN) :: FPIXELS(3),LPIXELS(3)
      INTEGER :: NAXES(3)
      DOUBLE PRECISION, INTENT(IN) :: ARRAY(*)
      CHARACTER*(*), INTENT(IN) :: FILENAME
      INTERFACE
      SUBROUTINE DELETEFILE(FILENAME,STATUS)
      INTEGER, INTENT(INOUT) :: STATUS
      CHARACTER*(*), INTENT(IN) :: FILENAME
      END SUBROUTINE DELETEFILE
      SUBROUTINE PRINTERROR(STATUS)
      INTEGER, INTENT(IN) :: STATUS
      END SUBROUTINE PRINTERROR
      END INTERFACE
      STATUS=0
      CALL FTGIOU(UNIT,STATUS)
      BLOCKSIZE=1
      BITPIX=-64
      NAXIS=3
      NAXES(1)=LPIXELS(1)-FPIXELS(1)+1
      NAXES(2)=LPIXELS(2)-FPIXELS(2)+1
      NAXES(3)=LPIXELS(3)-FPIXELS(3)+1
      CALL DELETEFILE(FILENAME,STATUS)
      CALL FTINIT(UNIT,FILENAME,BLOCKSIZE,STATUS)
      CALL FTPHPS(UNIT,BITPIX,NAXIS,NAXES,STATUS)
      GROUP=1
      CALL FTPSSD(UNIT,GROUP,NAXIS,NAXES,FPIXELS,LPIXELS,ARRAY,STATUS)
      CALL FTCLOS(UNIT,STATUS)
      CALL FTFIOU(UNIT,STATUS)
      IF (STATUS .GT. 0)CALL PRINTERROR(STATUS)
      RETURN
      END SUBROUTINE WRITEIMAGE
C ******************************************************************************
      SUBROUTINE READIMAGE(FILENAME,FPIXELS,LPIXELS,ARRAY)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: FPIXELS(3),LPIXELS(3)
      INTEGER :: STATUS,UNIT,READWRITE,BLOCKSIZE,NAXIS,NFOUND
      INTEGER :: GROUP,INCS(3),NAXES(3)
      DOUBLE PRECISION, INTENT(OUT) :: ARRAY(*)
      DOUBLE PRECISION :: NULLVAL
      LOGICAL :: ANYF
      CHARACTER*(*), INTENT(IN) :: FILENAME
      INTERFACE
      SUBROUTINE PRINTERROR(STATUS)
      INTEGER, INTENT(IN) :: STATUS
      END SUBROUTINE PRINTERROR
      END INTERFACE
      STATUS=0
      CALL FTGIOU(UNIT,STATUS)
      READWRITE=0
      CALL FTOPEN(UNIT,FILENAME,READWRITE,BLOCKSIZE,STATUS)
      CALL FTGKNJ(UNIT,'NAXIS',1,3,NAXES,NFOUND,STATUS)
      IF (NFOUND .NE. 3)THEN
        PRINT *,'READIMAGE failed to read the NAXIS keywords.'
        RETURN
      ENDIF
      GROUP=1
      NAXIS=3
      INCS(1)=1
      INCS(2)=1
      INCS(3)=1
      NULLVAL=-999
      CALL FTGSVD(UNIT,GROUP,NAXIS,NAXES,FPIXELS,LPIXELS,INCS,
     &  NULLVAL,ARRAY,ANYF,STATUS)
      IF (ANYF)THEN
        PRINT *,'One or more pixels are undefined in the FITS image.'
        RETURN
      ENDIF
      CALL FTCLOS(UNIT, STATUS)
      CALL FTFIOU(UNIT, STATUS)
      IF (STATUS .GT. 0)CALL PRINTERROR(STATUS)
      RETURN
      END SUBROUTINE READIMAGE
C ******************************************************************************
      SUBROUTINE AVERAGE(FILENAME,FPIXELS,LPIXELS,ARRAY,BUFFERSIZE)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: FPIXELS(3),LPIXELS(3)
      INTEGER, INTENT(INOUT), OPTIONAL :: BUFFERSIZE
      INTEGER :: L,LF,LT,LR,K,NPIXELS,M,N,NFRAMES,LBUFFER,STATUS
      CHARACTER*(*), INTENT(IN) :: FILENAME
      DOUBLE PRECISION, INTENT(OUT) :: ARRAY(LPIXELS(1)-FPIXELS(1)+1,
     &  LPIXELS(2)-FPIXELS(2)+1)
      DOUBLE PRECISION, ALLOCATABLE :: BUFFER(:,:,:)
      INTERFACE
      SUBROUTINE READIMAGE(FILENAME,FPIXELS,LPIXELS,DIMG)
      INTEGER, INTENT(IN) :: FPIXELS(3),LPIXELS(3)
      DOUBLE PRECISION, INTENT(OUT) :: DIMG(*)
      CHARACTER*(*), INTENT(IN) :: FILENAME
      END SUBROUTINE READIMAGE
      END INTERFACE
      STATUS=0
      M=LPIXELS(1)-FPIXELS(1)+1
      N=LPIXELS(2)-FPIXELS(2)+1
      NFRAMES=LPIXELS(3)-FPIXELS(3)+1
      NPIXELS=M*N
      IF (PRESENT(BUFFERSIZE))THEN
        LBUFFER=INT(FLOOR(DBLE(BUFFERSIZE)*1024*128/DBLE(NPIXELS)))
      ELSE
        LBUFFER=INT(FLOOR(DBLE(4*1024*1024)/DBLE(NPIXELS)))
      END IF
      ALLOCATE(BUFFER(M,N,LBUFFER),STAT=STATUS)
      IF(STATUS.NE.0)THEN
        PRINT *,'out of memory.'
        RETURN
      END IF
      ARRAY=0
      LF=1
      DO L=1,INT(CEILING(DBLE(NFRAMES)/DBLE(LBUFFER)))
        LT=MIN(LPIXELS(3),LF+LBUFFER-1)
        CALL READIMAGE(FILENAME,(/FPIXELS(1),FPIXELS(2),LF/),
     &    (/LPIXELS(1),LPIXELS(2),LT/),BUFFER(1:M,1:N,1:LBUFFER))
        LR=LT-LF+1
        LF=LT+1
        DO K=1,LR
          ARRAY=ARRAY+BUFFER(1:M,1:N,K)
        END DO
      END DO
      ARRAY=ARRAY/DBLE(NFRAMES)
      DEALLOCATE(BUFFER)
      RETURN
      END SUBROUTINE AVERAGE
C ******************************************************************************
      SUBROUTINE RESOLVEPATH(PATH,BASENAME,EXTNAME)
      IMPLICIT NONE
      CHARACTER*(*), INTENT(IN) :: PATH
      CHARACTER*(*), INTENT(OUT) :: BASENAME,EXTNAME
      CHARACTER :: FILENAME*80
      INTEGER :: K,L
      K=SCAN(PATH,'/',.TRUE.)
      L=LEN_TRIM(PATH)
      READ(PATH(K+1:L),*) FILENAME
      K=SCAN(FILENAME,'.')
      L=LEN_TRIM(FILENAME)
      READ(FILENAME(1:K-1),*) BASENAME
      K=SCAN(FILENAME,'.',.TRUE.)
      READ(FILENAME(K+1:L),*) EXTNAME
      RETURN
      END SUBROUTINE RESOLVEPATH
C ******************************************************************************
      SUBROUTINE ESTSNR(IMGFILE,FLAGFILE,FITN,DSNR,PREFIX)
C  Estimate the signal-to-noise ratio of given image.
C
      IMPLICIT NONE
      INCLUDE 'fftw3.f'
      INTEGER, INTENT(IN) :: FITN
      INTEGER :: STATUS,PLAN,NAXES(3),K,X,Y,NSAMPLES
      INTEGER, PARAMETER :: NBIN=20
      DOUBLE PRECISION, INTENT(OUT) :: DSNR
      DOUBLE PRECISION :: DMAX,DMIN,DMU,DSIGMA
      DOUBLE PRECISION, ALLOCATABLE :: DIMG(:,:),DFLAG(:,:),DBG(:,:),
     &  DB(:),DHIST(:)
      DOUBLE COMPLEX, ALLOCATABLE :: ZIN(:,:),ZOUT(:,:)
      CHARACTER(LEN=*), INTENT(IN) :: IMGFILE,FLAGFILE,PREFIX
      INTERFACE
      SUBROUTINE BGFIT2PN(NX,NY,DFLAG,DIMG,N,DBG,DB)
      INTEGER, INTENT(IN) :: NX,NY,N
      DOUBLE PRECISION, INTENT(IN) :: DFLAG(NX,NY),DIMG(NX,NY)
      DOUBLE PRECISION, INTENT(OUT) :: DBG(NX,NY),DB(*)
      END SUBROUTINE BGFIT2PN
      SUBROUTINE IMAGESIZE(FILENAME,NAXES)
      INTEGER, INTENT(IN) :: NAXES(3)
      CHARACTER(LEN=*), INTENT(IN) :: FILENAME
      END SUBROUTINE IMAGESIZE
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
      END INTERFACE
      CALL IMAGESIZE(IMGFILE,NAXES)
      WRITE(*,'(A,I3,A,I3)')' input image size: ',NAXES(1),' x ',
     &  NAXES(2)
      ALLOCATE(DIMG(NAXES(1),NAXES(2)),STAT=STATUS)
      IF(STATUS.NE.0)THEN
        PRINT *,'out of memory.'
        RETURN
      END IF
      CALL READIMAGE(IMGFILE,(/1,1,1/),(/NAXES(1),NAXES(2),1/),DIMG)
      ALLOCATE(DFLAG(NAXES(1),NAXES(2)),STAT=STATUS)
      IF(STATUS.NE.0)THEN
        PRINT *,'out of memory.'
        RETURN
      END IF
      CALL READIMAGE(FLAGFILE,(/1,1,1/),(/NAXES(1),NAXES(2),1/),DFLAG)
      NSAMPLES=NINT(SUM(DFLAG))
      WRITE(*,'(A,I5)')' number of flagged samples:',NSAMPLES
      ALLOCATE(DBG(NAXES(1),NAXES(2)),STAT=STATUS)
      IF(STATUS.NE.0)THEN
        PRINT *,'out of memory.'
        RETURN
      END IF
      ALLOCATE(DB(NAXES(1)*NAXES(2)),STAT=STATUS)
      IF(STATUS.NE.0)THEN
        PRINT *,'out of memory.'
        RETURN
      END IF
      CALL BGFIT2PN(NAXES(1),NAXES(2),DFLAG,DIMG,FITN,DBG,DB)
      CALL WRITEIMAGE(TRIM(PREFIX)//'_snr_bg.fits',(/1,1,1/),
     &  (/NAXES(1),NAXES(2),1/),DBG)
      PRINT *,'fitted background: '//TRIM(PREFIX)//'_snr_bg.fits'
      ALLOCATE(DHIST(NBIN),STAT=STATUS)
      IF(STATUS.NE.0)THEN
        PRINT *,'out of memory.'
        RETURN
      END IF
      DBG=(DBG-DIMG)*DFLAG
      CALL WRITEIMAGE(TRIM(PREFIX)//'_snr_res.fits',(/1,1,1/),
     &  (/NAXES(1),NAXES(2),1/),DBG)
      PRINT *,'background residual: '//TRIM(PREFIX)//'_snr_res.fits'
      DMAX=MAXVAL(DBG)
      DMIN=MINVAL(DBG)
      DHIST=0.0D0
      DO X=1,NAXES(1)
        DO Y=1,NAXES(2)
          DHIST(NINT((DBG(X,Y)-DMIN)/(DMAX-DMIN)*DBLE(NBIN-1))+1)=1.0D0+
     &      DHIST(NINT((DBG(X,Y)-DMIN)/(DMAX-DMIN)*DBLE(NBIN-1))+1)
        END DO
      END DO
      DHIST=DHIST/DBLE(NSAMPLES)
      DO K=1,NBIN
        WRITE(*,'(A,I2,A,ES10.3)')' x=',K,', y=',DHIST(K)
      END DO
      DMU=SUM(DBG)/DBLE(NSAMPLES)
      DBG=DBG-DMU
      DSIGMA=DSQRT(SUM(DBG*DBG)/DBLE(NSAMPLES))
      DBG=DBG/DSIGMA
      WRITE(*,'(A,ES9.2,A,ES9.2,A,ES9.2,A,ES9.2)')
     &  ' min:',DMIN,'max:',DMAX,' N(',DMU,',',DSIGMA,')'
      DEALLOCATE(DIMG)
      DEALLOCATE(DFLAG)
      DEALLOCATE(DBG)
      DEALLOCATE(DB)
      DEALLOCATE(DHIST)
      RETURN
      END SUBROUTINE ESTSNR
C ******************************************************************************
      SUBROUTINE BGFIT2PN(NX,NY,DFLAG,DIMG,N,DBG,DB)
C  2-dimensional background fitting subroutine using N-th polynomials.
C
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NX,NY,N
      INTEGER :: NSAMPLES,LWORK,INFO,LDA,LDB,NPARAMS,X,Y,L,PX,PY,K,
     &  STATUS
      DOUBLE PRECISION, INTENT(IN) :: DIMG(NX,NY),DFLAG(NX,NY)
      DOUBLE PRECISION, INTENT(OUT) :: DBG(NX,NY),DB(*)
      DOUBLE PRECISION :: DXC,DYC
      DOUBLE PRECISION, ALLOCATABLE :: WORK(:),A(:,:),DXI(:,:),DYI(:,:),
     &  DSAMPLEX(:),DSAMPLEY(:)
      STATUS=0
      DXC=0.5D0*(DBLE(NX+1))
      DYC=0.5D0*(DBLE(NY+1))
      NPARAMS=(N+1)*(N+2)/2
      NSAMPLES=NINT(SUM(DFLAG))
      WRITE(*,'(I6,A)'),NSAMPLES,' pixels have been sampled to'//
     &  ' determine the fitted background.'
      ALLOCATE(DSAMPLEX(NSAMPLES),STAT=STATUS)
      IF(STATUS.NE.0)THEN
        PRINT *,'out of memory.'
        RETURN
      END IF
      ALLOCATE(DSAMPLEY(NSAMPLES),STAT=STATUS)
      IF(STATUS.NE.0)THEN
        PRINT *,'out of memory.'
        RETURN
      END IF
      ALLOCATE(DXI(NX,NY),STAT=STATUS)
      IF(STATUS.NE.0)THEN
        PRINT *,'out of memory.'
        RETURN
      END IF
      ALLOCATE(DYI(NX,NY),STAT=STATUS)
      IF(STATUS.NE.0)THEN
        PRINT *,'out of memory.'
        RETURN
      END IF
      NSAMPLES=0
      DO X=1,NX
        DO Y=1,NY
          DXI(X,Y)=DBLE(X)
          DYI(X,Y)=DBLE(Y)
          IF(DFLAG(X,Y) .EQ. 1.0D0)THEN
            NSAMPLES=NSAMPLES+1
            DB(NSAMPLES)=DIMG(X,Y)
            DSAMPLEX(NSAMPLES)=DBLE(X)
            DSAMPLEY(NSAMPLES)=DBLE(Y)
          END IF
        END DO
      END DO
      IF(N.EQ.-1)THEN
        PRINT *,'use maximum value.'
        DB(1)=MAXVAL(DB(1:NSAMPLES))
        DBG=DB(1)
        RETURN
      END IF
      IF(N.EQ.0)THEN
        PRINT *,'use constant average value.'
      ELSE
        WRITE(*,'(A,I2,A)')'use ',N,'-th polynomials.'
      END IF
      ALLOCATE(A(NSAMPLES,NPARAMS),STAT=STATUS)
      IF(STATUS.NE.0)THEN
        PRINT *,'out of memory.'
        RETURN
      END IF
      L=1
      DO K=0,N
        DO PX=0,K
          PY=K-PX
          A(1:NSAMPLES,L)=DEXP(DLOG(DSAMPLEX(1:NSAMPLES))*DBLE(PX))*
     &      DEXP(DLOG(DSAMPLEY(1:NSAMPLES))*DBLE(PY))
          L=L+1
        END DO
      END DO
      DEALLOCATE(DSAMPLEX)
      DEALLOCATE(DSAMPLEY)
      LDA=NSAMPLES
      LDB=MAX(NSAMPLES,NPARAMS)
      LWORK=-1
      ALLOCATE(WORK(1),STAT=STATUS)
      IF(STATUS.NE.0)THEN
        PRINT *,'out of memory.'
        RETURN
      END IF
      CALL DGELS('N',NSAMPLES,NPARAMS,1,A,LDA,DB,LDB,WORK,
     &  LWORK,INFO)
      LWORK=INT(WORK(1))
      DEALLOCATE(WORK)
      ALLOCATE(WORK(LWORK),STAT=STATUS)
      IF(STATUS.NE.0)THEN
        PRINT *,'out of memory.'
        RETURN
      END IF
      CALL DGELS('N',NSAMPLES,NPARAMS,1,A,LDA,DB,LDB,WORK,
     &  LWORK,INFO)
      DEALLOCATE(WORK)
      DEALLOCATE(A)
      IF (INFO.GT.0)THEN
        PRINT *,'The least-squares solution could not be computed ',
     &    'because A does not have full rank.'
        RETURN
      ENDIF
      IF (INFO.LT.0)THEN
        PRINT *,'The argument has illegal value: ',ABS(INFO)
        RETURN
      ENDIF
      DBG=0.0D0
      L=1
      DO K=0,N
        DO PX=0,K
          PY=K-PX
          DBG=DBG+DB(L)*DEXP(DLOG(DXI)*DBLE(PX))*
     &      DEXP(DLOG(DYI)*DBLE(PY))
          L=L+1
        END DO
      END DO
      DEALLOCATE(DXI)
      DEALLOCATE(DYI)
      RETURN
      END SUBROUTINE BGFIT2PN
C ******************************************************************************
      SUBROUTINE BGFIT2P0(M,N,D,IMG,BG,B)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: M,N
      INTEGER :: NSAMPLES,LWORK,INFO,LDA,LDB,NPARAMS
      INTEGER :: I,J,K
      DOUBLE PRECISION, INTENT(IN) :: D,IMG(M,N)
      DOUBLE PRECISION, INTENT(OUT) :: BG(M,N),B(*)
      DOUBLE PRECISION :: X,Y,XC,YC
      NPARAMS=1
      NSAMPLES=0
      B(1)=0
      XC=0.5*(1+DBLE(N))
      YC=0.5*(1+DBLE(M))
      DO J=1,N
        DO I=1,M
          X=DBLE(J)
          Y=DBLE(I)
          IF (SQRT((X-XC)*(X-XC)+(Y-YC)*(Y-YC)).GE.D) THEN
            NSAMPLES=NSAMPLES+1
            B(1)=B(1)+IMG(I,J)
          END IF
        END DO
      END DO
      B(1)=B(1)/DBLE(NSAMPLES)
      BG=B(1)
      RETURN
      END SUBROUTINE BGFIT2P0
C ******************************************************************************
      SUBROUTINE BGFIT2P2(M,N,D,IMG,BG,B)
C  2-dimensional background fitting subroutine using 2nd polynomials.
C
C  Purpose
C  =======
C  z = a_0 + a_1*x + a_2*y + a_3*x^2 + a_4*x*y + a_5*y^2
C  x is column number.
C  y is row number.
C  z is value of IMG at (x,y).
C  Try to determine parameters a_i, i=0,1,2,3,4,5.
C
C  Arguments
C  =========
C  M is the number of rows of matrix IMG.
C  N is the number of columns of matrix IMG.
C  IMG is the matrix the background of which is to be fitted.
C  BG,B are the output of this subroutine.
C  Each row of A is vector (1, x, y, x^2, x*y, y^2) for a specific (x,y).
C  B is (z_1, z_2, z_3, ..., z_n)'. The subscript of z denotes different
C  location.
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: M,N
      INTEGER :: NSAMPLES,LWORK,INFO,LDA,LDB,NPARAMS
      INTEGER :: I,J,K
      DOUBLE PRECISION, INTENT(IN) :: D,IMG(M,N)
      DOUBLE PRECISION, INTENT(OUT) :: BG(M,N),B(*)
      DOUBLE PRECISION :: X,Y,XC,YC
      DOUBLE PRECISION, ALLOCATABLE :: WORK(:),A(:,:)
      NPARAMS=6
      NSAMPLES=0
      XC=0.5*(1+DBLE(N))
      YC=0.5*(1+DBLE(M))
      DO J=1,N
        DO I=1,M
          X=DBLE(J)
          Y=DBLE(I)
          IF (SQRT((X-XC)*(X-XC)+(Y-YC)*(Y-YC)).GE.D) THEN
            NSAMPLES=NSAMPLES+1
          ENDIF
        ENDDO
      ENDDO
      ALLOCATE(A(NSAMPLES,6))
      K=1
      DO J=1,N
        DO I=1,M
          X=DBLE(J)
          Y=DBLE(I)
          IF (SQRT((X-XC)*(X-XC)+(Y-YC)*(Y-YC)).GE.D) THEN
            A(K,:)=(/DBLE(1),X,Y,X*X,X*Y,Y*Y/)
            B(K)=IMG(INT(Y),INT(X))
            K=K+1
          ENDIF
        ENDDO
      ENDDO
      LDA=NSAMPLES
      LDB=MAX(NSAMPLES,6)
      LWORK=-1
      ALLOCATE(WORK(1))
      CALL DGELS('N',NSAMPLES,6,1,A,LDA,B,LDB,WORK,
     &  LWORK,INFO)
      LWORK=INT(WORK(1))
      DEALLOCATE(WORK)
      ALLOCATE(WORK(LWORK))
      CALL DGELS('N',NSAMPLES,6,1,A,LDA,B,LDB,WORK,
     &  LWORK,INFO)
      IF (INFO.GT.0)THEN
        PRINT *,'The least-squares solution could not be computed ',
     &    'because A does not have full rank.'
        RETURN
      ENDIF
      IF (INFO.LT.0)THEN
        PRINT *,'The argument has illegal value: ',ABS(INFO)
        RETURN
      ENDIF
      DO I=1,M
        DO J=1,N
          BG(I,J)=B(1)+B(2)*DBLE(J)+B(3)*DBLE(I)+B(4)*DBLE(J)*DBLE(J)
     &      +B(5)*DBLE(J)*DBLE(I)+B(6)*DBLE(I)*DBLE(I)
        ENDDO
      ENDDO
      DEALLOCATE(WORK)
      DEALLOCATE(A)
      RETURN
      END SUBROUTINE BGFIT2P2
C ******************************************************************************
      SUBROUTINE BGFIT2P4(M,N,D,IMG,BG,B)
C  2-dimensional background fitting subroutine using 4th polynomials.
C
C  Purpose
C  =======
C  z = a_0 + 
C    a_1*x + a_2*y + 
C    a_3*x^2 + a_4*x*y + a_5*y^2 + 
C    a_6*x^3 + a_7*x^2*y + a_8*x*y^2 + a_9*y^3 +
C    a_10*x^4 + a_11*x^3*y + a_12*x^2*y^2 + a_13*x*y^3 + a_14*y^4
C  x is column number.
C  y is row number.
C  z is value of IMG at (x,y).
C  Try to determine parameters a_i, i=0,1,2,3,...,14.
C
C  Arguments
C  =========
C  M is the number of rows of matrix IMG.
C  N is the number of columns of matrix IMG.
C  D is distance. Elements whose distance from the centre of IMG is larger than
C    D are sampled. BG is fitted to the sampled elements.
C  IMG is the matrix the background of which is to be fitted.
C  BG is the output of this subroutine.
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: M,N
      INTEGER :: NSAMPLES,LWORK,INFO,LDA,LDB,NPARAMS
      INTEGER :: I,J,K
      DOUBLE PRECISION, INTENT(IN) :: D,IMG(M,N)
      DOUBLE PRECISION, INTENT(OUT) :: BG(M,N),B(*)
      DOUBLE PRECISION :: X,Y,XC,YC
      DOUBLE PRECISION, ALLOCATABLE :: WORK(:),A(:,:)
      NPARAMS=15
      NSAMPLES=0
      XC=0.5*(1+DBLE(N))
      YC=0.5*(1+DBLE(M))
      DO J=1,N
        DO I=1,M
          X=DBLE(J)
          Y=DBLE(I)
          IF (SQRT((X-XC)*(X-XC)+(Y-YC)*(Y-YC)).GE.D) THEN
            NSAMPLES=NSAMPLES+1
          ENDIF
        ENDDO
      ENDDO
      ALLOCATE(A(NSAMPLES,NPARAMS))
      K=1
      DO J=1,N
        DO I=1,M
          X=DBLE(J)
          Y=DBLE(I)
          IF (DSQRT((X-XC)*(X-XC)+(Y-YC)*(Y-YC)).GE.D) THEN
            A(K,:)=(/DBLE(1),X,Y,X*X,X*Y,Y*Y,
     &        X*X*X,X*X*Y,X*Y*Y,Y*Y*Y,
     &        X*X*X*X,X*X*X*Y,X*X*Y*Y,X*Y*Y*Y,Y*Y*Y*Y/)
            B(K)=IMG(INT(Y),INT(X))
            K=K+1
          ENDIF
        ENDDO
      ENDDO
      LDA=NSAMPLES
      LDB=MAX(NSAMPLES,NPARAMS)
      LWORK=-1
      ALLOCATE(WORK(1))
      CALL DGELS('N',NSAMPLES,NPARAMS,1,A,LDA,B,LDB,WORK,
     &  LWORK,INFO)
      LWORK=INT(WORK(1))
      DEALLOCATE(WORK)
      ALLOCATE(WORK(LWORK))
      CALL DGELS('N',NSAMPLES,NPARAMS,1,A,LDA,B,LDB,WORK,
     &  LWORK,INFO)
      IF (INFO.GT.0)THEN
        PRINT *,'The least-squares solution could not be computed ',
     &    'because A does not have full rank.'
        RETURN
      ENDIF
      IF (INFO.LT.0)THEN
        PRINT *,'The argument has illegal value: ',ABS(INFO)
        RETURN
      ENDIF
      DO I=1,M
        DO J=1,N
          X=DBLE(J)
          Y=DBLE(I)
          BG(I,J)=B(1)+B(2)*X+B(3)*Y+
     &      B(4)*X*X+B(5)*X*Y+B(6)*Y*Y+
     &      B(7)*X*X*X+B(8)*X*X*Y+B(9)*X*Y*Y+B(10)*Y*Y*Y+
     &      B(11)*X*X*X*X+B(12)*X*X*X*Y+B(13)*X*X*Y*Y+B(14)*X*Y*Y*Y+
     &      B(15)*Y*Y*Y*Y
        ENDDO
      ENDDO
      DEALLOCATE(WORK)
      DEALLOCATE(A)
      RETURN
      END SUBROUTINE BGFIT2P4
C ******************************************************************************
      SUBROUTINE GETSNR(M,N,D,SNR,IMG,FTMETHOD)
      IMPLICIT NONE
      INCLUDE 'fftw3.f'
      INTEGER*8 :: PLAN1D,PLAN2D
      INTEGER :: M,N,I,J,K,NPIXS,NSPLS,INFO
      DOUBLE PRECISION :: SNR,D,X,Y,XC,YC,PS,PN
      DOUBLE PRECISION :: IMG(M,N),CLN(M,N),WORK(M,N)
      DOUBLE PRECISION, ALLOCATABLE :: BUFFER(:)
      DOUBLE COMPLEX :: ZIN2D(M,N),ZOUT2D(M,N)
      DOUBLE COMPLEX, ALLOCATABLE :: ZIN1D(:),ZOUT1D(:)
      CHARACTER*(*) :: FTMETHOD
      CALL DFFTW_INIT_THREADS(INFO)
      IF (INFO .EQ. 0)THEN
        PRINT *,'DFFTW_INIT_THREADS failed.'
        RETURN
      END IF
      CALL DFFTW_PLAN_WITH_NTHREADS(2)
      NPIXS=M*N
      NSPLS=0
      XC=0.5*(1+DBLE(N))
      YC=0.5*(1+DBLE(M))
      DO I=1,M
        DO J=1,N
          X=DBLE(J)
          Y=DBLE(I)
          IF (SQRT((X-XC)*(X-XC)+(Y-YC)*(Y-YC)).GE.D) THEN
            NSPLS=NSPLS+1
          END IF
        END DO
      END DO
      ALLOCATE(BUFFER(NSPLS))
      ALLOCATE(ZIN1D(NSPLS))
      ALLOCATE(ZOUT1D(NSPLS))
      CALL DFFTW_PLAN_DFT_1D(PLAN1D,NSPLS,ZIN1D,ZOUT1D,-1,
     &  FFTW_ESTIMATE+FFTW_DESTROY_INPUT)
      CALL DFFTW_PLAN_DFT_2D(PLAN2D,M,N,ZIN2D,ZOUT2D,-1,
     &  FFTW_ESTIMATE+FFTW_DESTROY_INPUT)
      IF (FTMETHOD .EQ. 'P2') THEN
        CALL BGFIT2P2(M,N,D,IMG,WORK,BUFFER)
      ELSE IF (FTMETHOD .EQ. 'P4') THEN
        CALL BGFIT2P4(M,N,D,IMG,WORK,BUFFER)
      ELSE IF (FTMETHOD .EQ. 'P0') THEN
        CALL BGFIT2P0(M,N,D,IMG,WORK,BUFFER(1))
      ELSE
        PRINT *,' Unknown fitting method.'
        RETURN
      END IF
      CLN=IMG-WORK
      K=0
      DO I=1,M
        DO J=1,N
          X=DBLE(J)
          Y=DBLE(I)
          IF (SQRT((X-XC)*(X-XC)+(Y-YC)*(Y-YC)).GE.D) THEN
            K=K+1
            BUFFER(K)=CLN(I,J)
          END IF
        END DO
      END DO
      ZIN2D(1:M,1:N)=CMPLX(CLN(1:M,1:N))
      CALL DFFTW_EXECUTE_DFT(PLAN2D,ZIN2D,ZOUT2D)
      PS=SUM(ZOUT2D*CONJG(ZOUT2D))/DBLE(NPIXS)/DSQRT(DBLE(NPIXS))
      ZIN1D(1:NSPLS)=CMPLX(BUFFER(1:NSPLS))
      CALL DFFTW_EXECUTE_DFT(PLAN1D,ZIN1D,ZOUT1D)
      PN=SUM(ZOUT1D*CONJG(ZOUT1D))/DBLE(NSPLS)/DSQRT(DBLE(NSPLS))
      SNR=PS/PN
      CALL DFFTW_DESTROY_PLAN(PLAN1D)
      CALL DFFTW_DESTROY_PLAN(PLAN2D)
C     CALL DFFTW_CLEANUP_THREADS()
      DEALLOCATE(ZIN1D)
      DEALLOCATE(BUFFER)
      DEALLOCATE(ZOUT1D)
      RETURN
      END SUBROUTINE GETSNR
C ******************************************************************************
      SUBROUTINE DECONVWNR(M,N,DG,DF,DH,DSNR)
C  Wiener deconvolution subroutine.
C
C  Purpose:
C  ========
C  G = conv(F, H). This routine returns F.
C
C  Arguments:
C  ==========
C  M - Number of rows of G.
C  N - Number of columns of G.
C  DG - Input.
C  DF - Output. size(F) = size(G).
C  DH - Input. size(PSF) = size (G)
C  DSNR - Signal-to-noise ratio.
C 
C  Declarations:
C  =============
      IMPLICIT NONE
      INCLUDE 'fftw3.f'
      INTEGER*8 :: PLAN
      INTEGER :: M,N,INFO
      DOUBLE PRECISION :: DG(M,N),DF(M,N),DH(M,N)
      DOUBLE PRECISION :: DSNR
      DOUBLE COMPLEX :: ZG(M,N),ZH(M,N)
      DOUBLE COMPLEX :: ZIN(M,N),ZOUT(M,N),ZDECONV(M,N)
      CALL DFFTW_INIT_THREADS(INFO)
      IF (INFO .EQ. 0)THEN
        PRINT *,'DFFTW_INIT_THREADS failed.'
      END IF
      CALL DFFTW_PLAN_WITH_NTHREADS(2)
      CALL DFFTW_PLAN_DFT_2D(PLAN,M,N,ZIN,ZOUT,-1,
     &  FFTW_ESTIMATE+FFTW_DESTROY_INPUT)
      ZIN=CMPLX(DG)
      CALL DFFTW_EXECUTE_DFT(PLAN,ZIN,ZOUT)
      ZG=ZOUT
      ZIN=CMPLX(DH)
      CALL ZFFTSHIFT(M,N,ZIN)
      CALL DFFTW_EXECUTE_DFT(PLAN,ZIN,ZOUT)
      ZH=ZOUT
      ZDECONV=CONJG(ZH)/(ZH*CONJG(ZH)+CMPLX(DBLE(1)/DSNR))
      CALL DFFTW_PLAN_DFT_2D(PLAN,M,N,ZIN,ZOUT,1,
     &  FFTW_ESTIMATE+FFTW_DESTROY_INPUT)
      ZIN=ZG*ZDECONV
      CALL DFFTW_EXECUTE_DFT(PLAN,ZIN,ZOUT)
      DF=DBLE(ZOUT)
      CALL DFFTW_DESTROY_PLAN(PLAN)
C     CALL DFFTW_CLEANUP_THREADS()
      RETURN
      END SUBROUTINE DECONVWNR
C ******************************************************************************
      SUBROUTINE ZFFTSHIFT(M,N,ZX)
C  Shift zero-frequency component to the centre of spectrum.
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: M,N
      DOUBLE COMPLEX, INTENT(INOUT) :: ZX(M,N)
      ZX=CSHIFT(CSHIFT(ZX,INT(FLOOR(0.5*DBLE(M))),1),
     &  INT(FLOOR(0.5*DBLE(N))),2)
      RETURN
      END SUBROUTINE ZFFTSHIFT
C ******************************************************************************
      SUBROUTINE DFFTSHIFT(M,N,DX)
C  Shift zero-frequency component to the centre of spectrum.
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: M,N
      DOUBLE PRECISION, INTENT(INOUT) :: DX(M,N)
      DX=CSHIFT(CSHIFT(DX,INT(FLOOR(0.5*DBLE(M))),1),
     &  INT(FLOOR(0.5*DBLE(N))),2)
      RETURN
      END SUBROUTINE DFFTSHIFT
C ******************************************************************************
      SUBROUTINE ZIFFTSHIFT(M,N,ZX)
C  Shift zero-frequency component to (1,1) position of spectrum.
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: M,N
      DOUBLE COMPLEX, INTENT(INOUT) :: ZX(M,N)
      ZX=CSHIFT(CSHIFT(ZX,INT(FLOOR(-0.5*DBLE(M))),1),
     &  INT(FLOOR(-0.5*DBLE(N))),2)
      RETURN
      END SUBROUTINE ZIFFTSHIFT
C ******************************************************************************
      SUBROUTINE DIFFTSHIFT(M,N,DX)
C  Shift zero-frequency component to (1,1) position of spectrum.
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: M,N
      DOUBLE PRECISION, INTENT(INOUT) :: DX(M,N)
      DX=CSHIFT(CSHIFT(DX,INT(FLOOR(-0.5*DBLE(M))),1),
     &  INT(FLOOR(-0.5*DBLE(N))),2)
      RETURN
      END SUBROUTINE DIFFTSHIFT
C ******************************************************************************
      SUBROUTINE DCORR2D(M,N,DA,DB,DC)
C  correlation coefficients of 2-dimensional array.
C
C  Purpose:
C  ========
C  c = corr(a, b). calculates c given a and b using fast fourier transform.
C  fft(c) = fft(a) * conj(fft(b)) (periodic extension is implied).
C
      IMPLICIT NONE
      INCLUDE 'fftw3.f'
      INTEGER*8 :: PLAN
      INTEGER :: M,N
      DOUBLE PRECISION, DIMENSION(M,N) :: DA,DB,DC
      DOUBLE COMPLEX, DIMENSION(M,N) :: ZA,ZB,ZIN,ZOUT
      CALL DFFTW_INIT_THREADS()
      CALL DFFTW_PLAN_WITH_NTHREADS(PLAN,2)
      CALL DFFTW_PLAN_DFT_2D(PLAN,M,N,ZIN,ZOUT,-1,
     &  FFTW_ESTIMATE+FFTW_DESTROY_INPUT)
      ZIN=CMPLX(DA)
      CALL DFFTW_EXECUTE_DFT(PLAN,ZIN,ZOUT)
      ZA=ZOUT
      ZIN=CMPLX(DB)
      CALL DFFTW_EXECUTE_DFT(PLAN,ZIN,ZOUT)
      ZB=ZOUT
      CALL DFFTW_PLAN_DFT_2D(PLAN,M,N,ZIN,ZOUT,1,
     &  FFTW_ESTIMATE+FFTW_DESTROY_INPUT)
      ZIN=ZA*CONJG(ZB)
      CALL DFFTW_EXECUTE_DFT(PLAN,ZIN,ZOUT)
      DC=DBLE(ZOUT)
      RETURN
      END SUBROUTINE DCORR2D

C ******************************************************************************
      SUBROUTINE SPECKLEINTERFEROMETRY(TFILE,TRNG,RFILE,RRNG,PX,PY,
     &  PREFIX,BUFFERSIZE)
C  Declarations:
C  =============
      IMPLICIT NONE
      INCLUDE 'fftw3.f'
      INTEGER*8 :: PLAN
      INTEGER, INTENT(IN) :: PX,PY,TRNG(2),RRNG(2)
      INTEGER, OPTIONAL, INTENT(INOUT) :: BUFFERSIZE
      INTEGER :: X,Y,NAXES(3)
      DOUBLE PRECISION, DIMENSION(PX,PY) :: DTPSD,DRPSD,DPSD,DACF
      DOUBLE COMPLEX :: ZIN(PX,PY),ZOUT(PX,PY)
      CHARACTER*(*), INTENT(IN) :: TFILE,RFILE,PREFIX
      INTERFACE
      SUBROUTINE IMAGESIZE(IMGFILE,NAXES)
      CHARACTER*(*), INTENT(IN) :: IMGFILE
      INTEGER, INTENT(OUT) :: NAXES(3)
      END SUBROUTINE IMAGESIZE
      SUBROUTINE GETPSD(FILENAME,FPIXELS,LPIXELS,DX,DY,DPSD,BUFFERSIZE)
      INTEGER, INTENT(IN) :: FPIXELS(3),LPIXELS(3),DX,DY
      INTEGER, OPTIONAL, INTENT(INOUT) :: BUFFERSIZE
      DOUBLE PRECISION, INTENT(OUT) :: DPSD(DX,DY)
      CHARACTER*(*), INTENT(IN) :: FILENAME
      END SUBROUTINE GETPSD
      SUBROUTINE WRITEIMAGE(FILENAME,FPIXELS,LPIXELS,ARRAY)
      INTEGER, INTENT(IN) :: FPIXELS(3),LPIXELS(3)
      DOUBLE PRECISION, INTENT(IN) :: ARRAY(*)
      CHARACTER*(*), INTENT(IN) :: FILENAME
      END SUBROUTINE WRITEIMAGE
      SUBROUTINE DFFTSHIFT(W,H,DX)
      INTEGER, INTENT(IN) :: W,H
      DOUBLE PRECISION, INTENT(INOUT) :: DX(W,H)
      END SUBROUTINE DFFTSHIFT
      SUBROUTINE DIFFTSHIFT(W,H,DX)
      INTEGER, INTENT(IN) :: W,H
      DOUBLE PRECISION, INTENT(INOUT) :: DX(W,H)
      END SUBROUTINE DIFFTSHIFT
      END INTERFACE
C  Statements:
C  ===========
      CALL IMAGESIZE(TFILE,NAXES)
      IF(TRNG(2).LT.NAXES(3))THEN
        PRINT *,'request is out of range.'
        RETURN
      END IF
      IF(PRESENT(BUFFERSIZE))THEN
        CALL GETPSD(TFILE,(/1,1,TRNG(1)/),(/NAXES(1),NAXES(2),TRNG(2)/),
     &    PX,PY,DTPSD,BUFFERSIZE)
      ELSE
        CALL GETPSD(TFILE,(/1,1,TRNG(1)/),(/NAXES(1),NAXES(2),TRNG(2)/),
     &    PX,PY,DTPSD)
      END IF
      CALL DFFTSHIFT(NAXES(1),NAXES(2),DTPSD)
      CALL WRITEIMAGE(TRIM(PREFIX)//'_tgabs.fits',
     &  (/1,1,1/),(/PX,PY,1/),DSQRT(DTPSD))
      PRINT *,'spectral amplitude of target: '//
     &  TRIM(PREFIX)//'_tgabs.fits'
      CALL IMAGESIZE(RFILE,NAXES)
      IF(RRNG(2).LT.NAXES(3))THEN
        PRINT *,'request is out of range.'
        RETURN
      END IF
      IF(PRESENT(BUFFERSIZE))THEN
        CALL GETPSD(RFILE,(/1,1,RRNG(1)/),(/NAXES(1),NAXES(2),RRNG(2)/),
     &    PX,PY,DRPSD,BUFFERSIZE)
      ELSE
        CALL GETPSD(RFILE,(/1,1,RRNG(1)/),(/NAXES(1),NAXES(2),RRNG(2)/),
     &    PX,PY,DRPSD)
      END IF
      CALL DFFTSHIFT(NAXES(1),NAXES(2),DRPSD)
      CALL WRITEIMAGE(TRIM(PREFIX)//'_rfabs.fits',
     &  (/1,1,1/),(/PX,PY,1/),DSQRT(DRPSD))
      PRINT *,'spectral amplitude of reference: '//
     &  TRIM(PREFIX)//'_rfabs.fits'
      DO X=1,PX
        DO Y=1,PY
          IF (DABS(DRPSD(X,Y)) .GT. 0)THEN
            DPSD(X,Y)=DTPSD(X,Y)/DRPSD(X,Y)
          ELSE
            DPSD(X,Y)=DBLE(0)
          END IF
        END DO
      END DO
      CALL WRITEIMAGE(TRIM(PREFIX)//'_dmpsd.fits',
     &  (/1,1,1/),(/PX,PY,1/),DPSD)
      PRINT *,'demodulated power spectral density: '//
     &  TRIM(PREFIX)//'_dmpsd.fits'
      CALL DFFTW_PLAN_DFT_2D(PLAN,PX,PY,ZIN,ZOUT,1,
     &  FFTW_ESTIMATE+FFTW_DESTROY_INPUT)
      CALL DIFFTSHIFT(NAXES(1),NAXES(2),DPSD)
      ZIN=CMPLX(DPSD)
      CALL DFFTW_EXECUTE_DFT(PLAN,ZIN,ZOUT)
      DACF=DREAL(ZOUT/DBLE(NAXES(1))/DBLE(NAXES(2)))
      CALL DFFTSHIFT(PX,PY,DACF)
      CALL WRITEIMAGE(TRIM(PREFIX)//'_acf.fits',
     &  (/1,1,1/),(/PX,PY,1/),DACF)
      PRINT *,'reconstructed auto-correlation function: '//
     &  TRIM(PREFIX)//'_acf.fits'
      CALL DFFTW_DESTROY_PLAN(PLAN)
      RETURN
      END SUBROUTINE SPECKLEINTERFEROMETRY
C ******************************************************************************
      SUBROUTINE GETPSD(FILENAME,FPIXELS,LPIXELS,DX,DY,DPSD,BUFFERSIZE)
C  Purpose:
C  ========
C  Get the mean power spectral density of all frames in the given FITS file.
C
C  Declarations:
C  =============
      IMPLICIT NONE
      INCLUDE 'fftw3.f'
      INTEGER, INTENT(IN):: FPIXELS(3),LPIXELS(3),DX,DY
      INTEGER, OPTIONAL, INTENT(INOUT) :: BUFFERSIZE
      INTEGER :: W,H,INFO,K,NPIXELS,NFRAMES,LBUF,NBUF,L,L1,L2
      INTEGER*8 :: PLAN
      DOUBLE PRECISION, INTENT(OUT) :: DPSD(DX,DY)
      DOUBLE PRECISION, ALLOCATABLE :: DBUF(:,:,:)
      DOUBLE COMPLEX, DIMENSION(DX,DY) :: ZIN,ZOUT
      CHARACTER*(*), INTENT(IN) :: FILENAME
      INTERFACE
      SUBROUTINE READIMAGE(FILENAME,FPIXELS,LPIXELS,DIMG)
      INTEGER, INTENT(IN) :: FPIXELS(3),LPIXELS(3)
      DOUBLE PRECISION, INTENT(OUT) :: DIMG(*)
      CHARACTER*(*), INTENT(IN) :: FILENAME
      END SUBROUTINE READIMAGE
      END INTERFACE
C  Statements:
C  ===========
      W=LPIXELS(1)-FPIXELS(1)+1
      H=LPIXELS(2)-FPIXELS(2)+1
      NPIXELS=W*H
      NFRAMES=LPIXELS(3)-FPIXELS(3)+1
      IF(PRESENT(BUFFERSIZE))THEN
        LBUF=INT(FLOOR(DBLE(BUFFERSIZE*1024*1024)/DBLE(NPIXELS*8)))
      ELSE
        LBUF=INT(FLOOR(DBLE(1024*1024)/DBLE(NPIXELS)))
      END IF
      NBUF=INT(CEILING(DBLE(NFRAMES)/DBLE(LBUF)))
      ALLOCATE(DBUF(W,H,LBUF))
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
C     PRINT *,'Start planning.'
      CALL DFFTW_PLAN_DFT_2D(PLAN,DX,DY,ZIN,ZOUT,-1,
     &  FFTW_MEASURE+FFTW_DESTROY_INPUT)
C     PRINT *,'Finished planning.'
      DPSD=DBLE(0)
      DO K=1,NBUF
        L1=(K-1)*LBUF+1
        L2=MIN(NFRAMES,K*LBUF)
        CALL READIMAGE(FILENAME,(/FPIXELS(1),FPIXELS(2),L1/),
     &    (/LPIXELS(1),LPIXELS(2),L2/),DBUF)
        DO L=1,L2+1-L1
          ZIN=CMPLX(DBLE(0))
          ZIN(1:W,1:H)=DBUF(:,:,L)*DBLE(NPIXELS)/SUM(DBUF(:,:,L))
          CALL DFFTW_EXECUTE_DFT(PLAN,ZIN,ZOUT)
          DPSD=DPSD+DREAL(ZOUT*CONJG(ZOUT))
        END DO
      END DO
      DPSD=DPSD/DBLE(NFRAMES)
      DEALLOCATE(DBUF)
      CALL DFFTW_DESTROY_PLAN(PLAN)
      RETURN
      END SUBROUTINE GETPSD
C ******************************************************************************
      SUBROUTINE DECONVCLEAN(M,N,DG,DF,DH,DBETA,MNUMIT)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: M,N
      INTEGER :: X,Y,XC,YC,K,L,MNUMIT
      DOUBLE PRECISION, INTENT(IN) :: DG(M,N),DH(M,N),DBETA
      DOUBLE PRECISION, INTENT(OUT) :: DF(M,N)
      DOUBLE PRECISION :: DRES(M,N),DPSF(M,N),DSNR,DR
      XC=INT(CEILING(0.5*DBLE(N+1)))
      YC=INT(CEILING(0.5*DBLE(M+1)))
      DR=DBLE(0.5)*DBLE(MIN(M,N))
      X=MAXLOC(MAXVAL(DH,1),2)
      Y=MAXLOC(MAXVAL(DH,2),1)
      DPSF=EOSHIFT(EOSHIFT(DH,Y-YC,DBLE(0),1),X-XC,DBLE(0),2)
      DPSF=DPSF/SUM(DPSF)
      DRES=DG
      L=1
      DO K=1,MNUMIT
        L=L+1
        X=MAXLOC(MAXVAL(DRES,1),2)
        Y=MAXLOC(MAXVAL(DRES,2),1)
        DF(Y,X)=DF(Y,X)+DRES(Y,X)*DBETA
        DRES=DRES-DRES(Y,X)*DBETA*EOSHIFT(EOSHIFT(DPSF,YC-Y,DBLE(0),1),
     &    XC-X,DBLE(0),2)
        IF (SUM(DRES) .LE. 0)THEN
          EXIT
        END IF
        IF (L .GT. 100)THEN
          L=1
          PRINT *,'Loop: ',K
          CALL GETSNR(M,N,DR,DSNR,DRES,'P4')
          WRITE (*,'(A,ES10.3)') ' Residual SNR:',DSNR
          WRITE (*,'(A,ES10.3,A,ES10.3,A,ES10.3)')
     &     ' Residual sum:',SUM(DRES),
     &      ', CLEAN sum:',SUM(DF),', Total sum:',SUM(DF)+SUM(DRES)
          WRITE (*,'(A,I3,A,I3,A,ES10.3,ES10.3)')
     &      ' Maxima: (',X,',',Y,'),',MAXVAL(DRES),DRES(Y,X)
        END IF
      END DO
      CALL WRITEIMAGE('CLEAN_RESIDUAL.FITS',(/1,1,1/),(/M,N,1/),DRES)
      RETURN
      END SUBROUTINE DECONVCLEAN
C ******************************************************************************
      SUBROUTINE APPENDIMAGE(FILENAME,FPIXELS,LPIXELS,ARRAY)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: FPIXELS(3),LPIXELS(3)
      INTEGER :: STATUS,UNIT,BLOCKSIZE,BITPIX,NAXIS,GROUP,RWMODE,EXISTS
      INTEGER :: NAXES(3),NFOUND
      DOUBLE PRECISION, INTENT(IN) :: ARRAY(*)
      CHARACTER*(*), INTENT(IN) :: FILENAME
      INTERFACE
      SUBROUTINE PRINTERROR(STATUS)
      INTEGER, INTENT(IN) :: STATUS
      END SUBROUTINE PRINTERROR
      END INTERFACE
      STATUS=0
      CALL FTGIOU(UNIT,STATUS)
      BLOCKSIZE=1
      BITPIX=-64
      NAXIS=3
      RWMODE=1
      CALL FTEXIST(FILENAME,EXISTS,STATUS)
      IF (EXISTS .EQ. 1)THEN
C       PRINT *,TRIM(FILENAME),' exists.'
        CALL FTIOPN(UNIT,FILENAME,RWMODE,STATUS)
        CALL FTGKNJ(UNIT,'NAXIS',1,3,NAXES,NFOUND,STATUS)
        IF (NFOUND .NE. 3)THEN
          PRINT *,'APPENDIMAGE failed to read the NAXIS keywords.'
          RETURN
        ENDIF
        NAXES(3)=NAXES(3)+LPIXELS(3)-FPIXELS(3)+1
        CALL FTUKYJ(UNIT,'NAXIS3',NAXES(3),'length of data axis 3',
     &    STATUS)
      ELSE
C       PRINT *,TRIM(FILENAME),' does not exist.'
        CALL FTINIT(UNIT,FILENAME,BLOCKSIZE,STATUS)
        NAXES=LPIXELS-FPIXELS+1
        CALL FTPHPS(UNIT,BITPIX,NAXIS,NAXES,STATUS)
      END IF
      GROUP=1
      CALL FTPSSD(UNIT,GROUP,NAXIS,NAXES,FPIXELS,LPIXELS,ARRAY,STATUS)
      CALL FTCLOS(UNIT,STATUS)
      CALL FTFIOU(UNIT,STATUS)
      IF (STATUS .GT. 0)CALL PRINTERROR(STATUS)
      RETURN
      END SUBROUTINE APPENDIMAGE

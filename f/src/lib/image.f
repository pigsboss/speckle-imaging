      SUBROUTINE READBMP(FILENAME,BMOFFSET,NX,NY,PSIZE,DIMG)
      IMPLICIT NONE
      CHARACTER(LEN=*),INTENT(IN) :: FILENAME
      INTEGER,INTENT(IN) :: BMOFFSET,NX,NY,PSIZE
      DOUBLE PRECISION,INTENT(OUT) :: DIMG(NX,NY)
      INTEGER :: STATUS,ROWSIZE,X,Y
      INTEGER,PARAMETER :: UNIT=8
      INTEGER*1,ALLOCATABLE :: BM(:)
      DOUBLE PRECISION,ALLOCATABLE :: DBM(:)
      OPEN(UNIT=UNIT,FILE=TRIM(FILENAME),ACCESS='STREAM',
     & FORM='UNFORMATTED',ACTION='READ',POSITION='REWIND',
     & IOSTAT=STATUS)
      IF(STATUS.NE.0)THEN
        PRINT *,'error: open file failed.'
        RETURN
      END IF
      ROWSIZE=INT(CEILING(FLOAT(NX*PSIZE)/4.0))*4
      ALLOCATE(BM(ROWSIZE),STAT=STATUS)
      IF(STATUS.NE.0)THEN
        PRINT *,'error: out of memory.'
        RETURN
      END IF
      ALLOCATE(DBM(ROWSIZE),STAT=STATUS)
      IF(STATUS.NE.0)THEN
        PRINT *,'error: out of memory.'
        RETURN
      END IF
      DO Y=1,NY
        READ(UNIT,POS=(BMOFFSET+(Y-1)*ROWSIZE)) BM
        DO X=1,ROWSIZE
          IF(BM(X) .LT. 0) THEN
            DBM(X)=DBLE(BM(X))+256.0D0
          ELSE
            DBM(X)=DBLE(BM(X))
          END IF
        END DO
        DO X=1,NX
          DIMG(X,Y)=SUM(DBLE(DBM(1+(X-1)*PSIZE:PSIZE*X)))
C         DIMG(X,Y)=DBLE(DBM(1+(X-1)*PSIZE))
        END DO
      END DO
      DEALLOCATE(BM)
      DEALLOCATE(DBM)
      RETURN
      END SUBROUTINE READBMP
C ******************************************************************************
      SUBROUTINE BMPINFO(FILENAME,BMOFFSET,NX,NY,PBITS,COMPRESS)
      IMPLICIT NONE
      CHARACTER(LEN=*),INTENT(IN) :: FILENAME
      INTEGER*4,INTENT(OUT) :: BMOFFSET,NX,NY,PBITS,COMPRESS
      INTEGER :: STATUS
      INTEGER, PARAMETER :: UNIT=8
      INTEGER*4 :: BMSIZE,HEADSIZE
      INTEGER*2 :: SPBITS
      CHARACTER(LEN=2) :: MAGNUM
      STATUS=0
      PBITS=0
      COMPRESS=0
      OPEN(UNIT=UNIT,FILE=TRIM(FILENAME),ACCESS='STREAM',
     & FORM='UNFORMATTED',ACTION='READ',POSITION='REWIND',
     & IOSTAT=STATUS)
      IF(STATUS.NE.0)THEN
        PRINT *,'open file failed.'
        RETURN
      END IF
      READ(UNIT,POS=1) MAGNUM
      READ(UNIT,POS=3) BMSIZE
      READ(UNIT,POS=11) BMOFFSET
      READ(UNIT,POS=15) HEADSIZE
      READ(UNIT,POS=19) NX
      READ(UNIT,POS=23) NY
      WRITE(*,*) TRIM(FILENAME)
      WRITE(*,*)          'magic number   :          '//MAGNUM
      WRITE(*,'(A,I8,A)')' total size     : ',BMSIZE,   ' bytes'
      WRITE(*,'(A,I8,A)')' raw data offset: ',BMOFFSET, ' bytes'
      WRITE(*,'(A,I8,A)')' header size    : ',HEADSIZE, ' bytes'
      WRITE(*,'(A,I8,A)')' bitmap width   : ',NX,       ' pixels'
      WRITE(*,'(A,I8,A)')' bitmap height  : ',NY,       ' pixels'
      IF(HEADSIZE .GE. 40) THEN
        READ(UNIT,POS=29) SPBITS
        PBITS=SPBITS
        READ(UNIT,POS=31) COMPRESS
        WRITE(*,'(A,I8)')' bits per pixel : ',PBITS
        WRITE(*,'(A,I8)')' compress type  : ',COMPRESS
      END IF
      CLOSE(UNIT)      
      RETURN
      END SUBROUTINE BMPINFO
C ******************************************************************************
      SUBROUTINE SPECTRUMTOIMAGE(NX,NY,ZSP,DIMG)
      IMPLICIT NONE
      INCLUDE 'fftw3.f'
      INTEGER, INTENT(IN) :: NX,NY
      INTEGER :: PLAN
      DOUBLE PRECISION, INTENT(OUT) :: DIMG(NX,NY)
      DOUBLE COMPLEX, INTENT(IN) :: ZSP(NX,NY)
      DOUBLE COMPLEX :: ZIN(NX,NY),ZOUT(NX,NY)
      INTERFACE
      SUBROUTINE DFFTSHIFT(NX,NY,DIMG)
      INTEGER, INTENT(IN) :: NX,NY
      DOUBLE PRECISION, INTENT(INOUT) :: DIMG(NX,NY)
      END SUBROUTINE DFFTSHIFT
      END INTERFACE
C
      CALL DFFTW_PLAN_DFT_2D(PLAN,NX,NY,ZIN,ZOUT,1,
     &  FFTW_ESTIMATE+FFTW_DESTROY_INPUT)
      ZIN=ZSP
      CALL DFFTW_EXECUTE_DFT(PLAN,ZIN,ZOUT)
      DIMG=DREAL(ZOUT)/DSQRT(DBLE(NX*NY))
      CALL DFFTSHIFT(NX,NY,DIMG)
      CALL DFFTW_DESTROY_PLAN(PLAN)
      RETURN
      END SUBROUTINE SPECTRUMTOIMAGE
C ******************************************************************************
      SUBROUTINE IMAGETOSPECTRUM(NX,NY,DIMG,ZSP)
      IMPLICIT NONE
      INCLUDE 'fftw3.f'
      INTEGER, INTENT(IN) :: NX,NY
      INTEGER :: PLAN
      DOUBLE PRECISION, INTENT(IN) :: DIMG(NX,NY)
      DOUBLE COMPLEX, INTENT(OUT) :: ZSP(NX,NY)
      DOUBLE COMPLEX :: ZIN(NX,NY)
      INTERFACE
      SUBROUTINE ZIFFTSHIFT(NX,NY,ZSP)
      INTEGER, INTENT(IN) :: NX,NY
      DOUBLE COMPLEX, INTENT(INOUT) :: ZSP(NX,NY)
      END SUBROUTINE ZIFFTSHIFT
      END INTERFACE
C
      CALL DFFTW_PLAN_DFT_2D(PLAN,NX,NY,ZIN,ZSP,-1,
     &  FFTW_ESTIMATE+FFTW_DESTROY_INPUT)
      
      ZIN=DCMPLX(DIMG)
      CALL ZIFFTSHIFT(NX,NY,ZIN)
      CALL DFFTW_EXECUTE_DFT(PLAN,ZIN,ZSP)
      ZSP=ZSP/DSQRT(DBLE(NX*NY))
      CALL DFFTW_DESTROY_PLAN(PLAN)
      RETURN
      END SUBROUTINE IMAGETOSPECTRUM
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
        L1=(K-1)*LBUFFER+FPIXELS(3)
        L2=MIN(K*LBUFFER,LPIXELS(3))
        PRINT *,'READIMAGE from',L1,'to',L2
        CALL READIMAGE(IMGFILE,(/FPIXELS(1),FPIXELS(2),L1/),
     &    (/LPIXELS(1),LPIXELS(2),L2/),DBUF)
        CALL APPENDIMAGE(SUBFILE,(/1,1,L1+1-FPIXELS(3)/),
     &    (/M,N,L2+1-FPIXELS(3)/),DBUF)
      END DO
      DEALLOCATE(DBUF)
      IF (STATUS .GT. 0)CALL PRINTERROR(STATUS)
      RETURN
      END SUBROUTINE SUBIMAGE
C ******************************************************************************
      SUBROUTINE SPLITIMAGE(IMGFILE,FPIXELS,LPIXELS,PREFIX,BUFFERSIZE)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: FPIXELS(3),LPIXELS(3)
      CHARACTER*(*), INTENT(IN) :: IMGFILE,PREFIX
      CHARACTER*(256) :: ZEROS,NUMSTR,SUBFILE
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
      SUBROUTINE WRITEIMAGE(FILENAME,FPIXELS,LPIXELS,DIMG)
      INTEGER, INTENT(IN) :: FPIXELS(3),LPIXELS(3)
      DOUBLE PRECISION, INTENT(IN) :: DIMG(*)
      CHARACTER*(*), INTENT(IN) :: FILENAME
      END SUBROUTINE WRITEIMAGE
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
      ZEROS='0000000000000'
      DO K=1,NBUFFER
        L1=(K-1)*LBUFFER+FPIXELS(3)
        L2=MIN(K*LBUFFER,LPIXELS(3))
        PRINT *,'READIMAGE from',L1,'to',L2
        CALL READIMAGE(IMGFILE,(/FPIXELS(1),FPIXELS(2),L1/),
     &    (/LPIXELS(1),LPIXELS(2),L2/),DBUF)
        DO L=L1,L2
          WRITE(NUMSTR,*) L
          PRINT *,'WRITEIMAGE ',L
          SUBFILE=TRIM(PREFIX)//
     &      TRIM(ZEROS(1:6-LEN_TRIM(ADJUSTL(NUMSTR))))//
     &      TRIM(ADJUSTL(NUMSTR))//'.fits'
          CALL WRITEIMAGE(SUBFILE,(/1,1,1/),
     &      (/M,N,1/),DBUF(:,:,L+1-FPIXELS(3)-(K-1)*LBUFFER))
        END DO
      END DO
      DEALLOCATE(DBUF)
      IF (STATUS .GT. 0)CALL PRINTERROR(STATUS)
      RETURN
      END SUBROUTINE SPLITIMAGE
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
C ******************************************************************************
      INCLUDE './lib/writeimage.f'
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
C ******************************************************************************

      PROGRAM MAIN
C  Declarations:
C  =============
      INTEGER :: STATUS,UNIT,FPIXELS(3),LPIXELS(3),NAXES(3),NPIXELS
      DOUBLE PRECISION, ALLOCATABLE :: AVG(:),IMG(:,:),BG(:,:),CLN(:)
      DOUBLE PRECISION :: IMGRAD
      CHARACTER :: INFILE*80,OUTFILE*80,AVGFILE*80,ARG*255
      CHARACTER :: BASENAME*80,EXTNAME*80,BGFILE*80,CLNFILE*80
C  Statements:
C  ===========
      STATUS=0
C    Resolve the command line options:
C    =================================
      CALL GETARG(1,INFILE)
      WRITE(*,*) 'Input: ',INFILE
      CALL GETARG(2,ARG)
      READ(ARG,*) FPIXELS(1)
      CALL GETARG(3,ARG)
      READ(ARG,*) FPIXELS(2)
      CALL GETARG(4,ARG)
      READ(ARG,*) FPIXELS(3)
      CALL GETARG(5,ARG)
      READ(ARG,*) LPIXELS(1)
      CALL GETARG(6,ARG)
      READ(ARG,*) LPIXELS(2)
      CALL GETARG(7,ARG)
      READ(ARG,*) LPIXELS(3)
      CALL GETARG(8,ARG)
      READ(ARG,*) IMGRAD
      CALL GETARG(9,OUTFILE)
      CALL RESOLVEPATH(INFILE,BASENAME,EXTNAME)
C    Determine the size of the problem:
C    ==================================
      NAXES(1)=LPIXELS(1)-FPIXELS(1)+1
      NAXES(2)=LPIXELS(2)-FPIXELS(2)+1
      NAXES(3)=LPIXELS(3)-FPIXELS(3)+1
      NPIXELS=NAXES(1)*NAXES(2)
C    Allocate the space:
C    ===================
      ALLOCATE(AVG(NPIXELS))
      ALLOCATE(IMG(NAXES(1),NAXES(2)))
      ALLOCATE(BG(NAXES(1),NAXES(2)))
      ALLOCATE(CLN(NPIXELS))
C    Calculate the average of all the frames:
C    ========================================
      PRINT *,'Start averaging all frames...'
      CALL AVERAGE(INFILE,FPIXELS,LPIXELS,AVG)
      PRINT *,'Done.'
      FPIXELS=(/1,1,1/)
      LPIXELS=(/NAXES(1),NAXES(2),1/)
      AVGFILE=TRIM(BASENAME)//'_'//TRIM(OUTFILE)//'_AVG.FITS'
      CALL WRITEIMAGE(AVGFILE,(/1,1,1/),(/NAXES(1),NAXES(2),1/),AVG)
      PRINT *,'Write the average frame to: ',AVGFILE
C    Fitting the background of the image:
C    ====================================
      IMG=RESHAPE(AVG,(/NAXES(1),NAXES(2)/))
C      Using 2nd polynomials:
C      ======================
      PRINT *,'Start fitting the background of the average frame',
     &  ' using 2nd polynomials.'
      CALL BGFIT2P2(NAXES(1),NAXES(2),IMGRAD,IMG,BG)
      BGFILE=TRIM(BASENAME)//'_'//TRIM(OUTFILE)//'_BG_P2.FITS'
      CALL WRITEIMAGE(BGFILE,FPIXELS,LPIXELS,BG)
      PRINT *,'Write the background to: ',BGFILE
      CALL DCOPY(NPIXELS,AVG,1,CLN,1)
      PRINT *,'Subtract the background from the average frame.'
      CALL DAXPY(NPIXELS,DBLE(-1),BG,1,CLN,1)
      CLNFILE=TRIM(BASENAME)//'_'//TRIM(OUTFILE)//'_CLN_P2.FITS'
      CALL WRITEIMAGE(CLNFILE,FPIXELS,LPIXELS,CLN)
      PRINT *,'Write the clean frame to: ',CLNFILE
C      Using 4th polynomials:
C      ======================
      PRINT *,'Start fitting the background of the average frame',
     &  ' using 4th polynomials.'
      CALL BGFIT2P4(NAXES(1),NAXES(2),IMGRAD,IMG,BG)
      BGFILE=TRIM(BASENAME)//'_'//TRIM(OUTFILE)//'_BG_P4.FITS'
      CALL WRITEIMAGE(BGFILE,FPIXELS,LPIXELS,BG)
      PRINT *,'Write the background to: ',BGFILE
      CALL DCOPY(NPIXELS,AVG,1,CLN,1)
      PRINT *,'Subtract the background from the average frame.'
      CALL DAXPY(NPIXELS,DBLE(-1),BG,1,CLN,1)
      CLNFILE=TRIM(BASENAME)//'_'//TRIM(OUTFILE)//'_CLN_P4.FITS'
      CALL WRITEIMAGE(CLNFILE,FPIXELS,LPIXELS,CLN)
      PRINT *,'Write the clean frame to: ',CLNFILE
C    Deallocate space:
C    =================
      DEALLOCATE(IMG)
      DEALLOCATE(CLN)
      DEALLOCATE(BG)
      DEALLOCATE(AVG)
      STOP
      END
C **********************************************************************
      SUBROUTINE SSAREC(FILENAME,FPIXELS,LPIXELS,ARRAY)
      INTEGER :: FPIXELS(*),LPIXELS(*),K,NPIXELS,F(3),L(3)
      CHARACTER*(*) :: FILENAME
      DOUBLE PRECISION :: ARRAY(*),DA
      DOUBLE PRECISION, ALLOCATABLE :: BUFFER(:)
      INTEGER :: NAXES(3)
      NAXES(1)=LPIXELS(1)-FPIXELS(1)+1
      NAXES(2)=LPIXELS(2)-FPIXELS(2)+1
      NAXES(3)=LPIXELS(3)-FPIXELS(3)+1
      NPIXELS=NAXES(1)*NAXES(2)
      F(1)=FPIXELS(1)
      F(2)=FPIXELS(2)
      L(1)=LPIXELS(1)
      L(2)=LPIXELS(2)
      ALLOCATE(BUFFER(NPIXELS))
      DO K=1,NPIXELS
        ARRAY(K)=0
      ENDDO
      DA=1.0
      DO K=FPIXELS(3),LPIXELS(3)
        F(3)=K
        L(3)=K
        CALL READIMAGE(FILENAME,F,L,BUFFER)
        CALL DAXPY(NPIXELS,DA,BUFFER,1,ARRAY,1)
      ENDDO
      DA=1.0/DBLE(NAXES(3))
      CALL DSCAL(NPIXELS,DA,ARRAY,1)
      DEALLOCATE(BUFFER)
      RETURN
      ENDSUBROUTINE SSAREC
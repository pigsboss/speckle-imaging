      PROGRAM MAIN
C  Declarations:
C  =============
      INTEGER :: STATUS,UNIT,FPIXELS(3),LPIXELS(3),NAXES(3)
      INTEGER :: K,NPIXELS
      DOUBLE PRECISION, ALLOCATABLE :: AVG(:),IMG(:,:),BG(:,:),CLN(:)
      DOUBLE PRECISION, ALLOCATABLE :: WORK(:),A(:,:),B(:)
      DOUBLE PRECISION :: IMGRAD
      CHARACTER :: INFILE*80,SUFFIX*80,AVGFILE*80,ARG*255
      CHARACTER :: BASENAME*80,EXTNAME*80,BGFILE*80,CLNFILE*80
      CHARACTER :: OUTFILE*80
C  Statements:
C  ===========
      STATUS=0
C    Resolve the command line options:
C    =================================
      CALL GETARG(1,INFILE)
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
      CALL GETARG(9,SUFFIX)
      CALL RESOLVEPATH(INFILE,BASENAME,EXTNAME)
C    Determine the size of the problem:
C    ==================================
      NAXES(1)=LPIXELS(1)-FPIXELS(1)+1
      NAXES(2)=LPIXELS(2)-FPIXELS(2)+1
      NAXES(3)=LPIXELS(3)-FPIXELS(3)+1
      NPIXELS=NAXES(1)*NAXES(2)
      PRINT *,'*****************************************************'
      PRINT *,'*                    Definition                     *'
      PRINT *,'*****************************************************'
      PRINT *,'Input: ',INFILE
      WRITE(*,'(A,I6)') ' First frame: ',FPIXELS(3)
      WRITE(*,'(A,I6)') ' Last frame: ',LPIXELS(3)
      WRITE(*,'(A,I5,A,I5,A)') ' First pixel: (',FPIXELS(1),',',
     &  FPIXELS(2),')'
      WRITE(*,'(A,I5,A,I5,A)') ' Last pixel: (',LPIXELS(1),',',
     &  LPIXELS(2),')'
      WRITE(*,'(A,I5,A,I5)') ' Frame size: ',NAXES(1),' x ',NAXES(2)
      WRITE(*,'(A,I6)') ' Number of frames: ',NAXES(3)
C    Allocate the space:
C    ===================
      ALLOCATE(AVG(NPIXELS))
      ALLOCATE(WORK(NPIXELS))
      ALLOCATE(IMG(NAXES(1),NAXES(2)))
      ALLOCATE(BG(NAXES(1),NAXES(2)))
      ALLOCATE(CLN(NPIXELS))
      ALLOCATE(B(NPIXELS))
C    Calculate the average of all the frames:
C    ========================================
      PRINT *,'*'
      PRINT *,'*'
      PRINT *,'*'
      PRINT *,'*****************************************************'
      PRINT *,'*              Average of the frames:               *'
      PRINT *,'*****************************************************'
      CALL AVERAGE(INFILE,FPIXELS,LPIXELS,AVG,WORK)
      AVGFILE=TRIM(BASENAME)//'_'//TRIM(SUFFIX)//'_AVG.FITS'
      CALL WRITEIMAGE(AVGFILE,(/1,1,1/),(/NAXES(1),NAXES(2),1/),AVG)
      PRINT *,'Average: ',AVGFILE
C    Fitting the background of the image:
C    ====================================
      PRINT *,'*'
      PRINT *,'*'
      PRINT *,'*'
      PRINT *,'*****************************************************'
      PRINT *,'*      Calculate the background of the image        *'
      PRINT *,'*****************************************************'
      IMG=RESHAPE(AVG,(/NAXES(1),NAXES(2)/))
C      Using 2nd polynomials:
C      ======================
      PRINT *,'*'
      PRINT *,'2nd Polynomials:'
      PRINT *,'--------------------------------------------'
      CALL BGFIT2P2(NAXES(1),NAXES(2),IMGRAD,IMG,BG,B)
      PRINT *,'Fitting result:'
      PRINT *,'--------------------------------------------'
      DO K=1,6
        WRITE(*,'(A,ES10.3,A,I1)') ' a_i = ',B(K),', i = ',K
      ENDDO
      BGFILE=TRIM(BASENAME)//'_'//TRIM(SUFFIX)//'_BG_P2.FITS'
      CALL WRITEIMAGE(BGFILE,(/1,1,1/),(/NAXES(1),NAXES(2),1/),BG)
      PRINT *,'Background: ',BGFILE
      CALL DCOPY(NPIXELS,AVG,1,CLN,1)
      CALL DAXPY(NPIXELS,DBLE(-1),BG,1,CLN,1)
      CLNFILE=TRIM(BASENAME)//'_'//TRIM(SUFFIX)//'_CLN_P2.FITS'
      CALL DSCAL(NPIXELS,DBLE(NPIXELS)/SUM(CLN),CLN,1)
      CALL WRITEIMAGE(CLNFILE,(/1,1,1/),(/NAXES(1),NAXES(2),1/),CLN)
      PRINT *,'Clean image: ',CLNFILE
C      Using 4th polynomials:
C      ======================
      PRINT *,'*'
      PRINT *,'*'
      PRINT *,'4th Polynomials:'
      PRINT *,'--------------------------------------------'
      CALL BGFIT2P4(NAXES(1),NAXES(2),IMGRAD,IMG,BG)
      BGFILE=TRIM(BASENAME)//'_'//TRIM(SUFFIX)//'_BG_P4.FITS'
      CALL WRITEIMAGE(BGFILE,(/1,1,1/),(/NAXES(1),NAXES(2),1/),BG)
      PRINT *,'Background: ',BGFILE
      CALL DCOPY(NPIXELS,AVG,1,CLN,1)
      CALL DAXPY(NPIXELS,DBLE(-1),BG,1,CLN,1)
      CLNFILE=TRIM(BASENAME)//'_'//TRIM(SUFFIX)//'_CLN_P4.FITS'
      CALL DSCAL(NPIXELS,DBLE(NPIXELS)/SUM(CLN),CLN,1)
      CALL WRITEIMAGE(CLNFILE,(/1,1,1/),(/NAXES(1),NAXES(2),1/),CLN)
      PRINT *,'Clean image: ',CLNFILE
C    Simple shift-and-add:
C    =====================
      PRINT *,'*'
      PRINT *,'*'
      PRINT *,'*'
      PRINT *,'*****************************************************'
      PRINT *,'*              Simple shift-and-add                 *'
      PRINT *,'*****************************************************'
      CALL SSAREC(INFILE,FPIXELS,LPIXELS,IMG,BG,WORK)
      OUTFILE=TRIM(BASENAME)//'_'//TRIM(SUFFIX)//'_SSA.FITS'
      CALL WRITEIMAGE(OUTFILE,(/1,1,1/),(/NAXES(1),NAXES(2),1/),IMG)
      PRINT *,'SSA: ',OUTFILE
C    Deallocate space:
C    =================
      DEALLOCATE(IMG)
      DEALLOCATE(WORK)
      DEALLOCATE(CLN)
      DEALLOCATE(BG)
      DEALLOCATE(AVG)
      DEALLOCATE(B)
      STOP
      END
C *
C *
C *
C **********************************************************************
C *
C *
C *
      SUBROUTINE SSAREC(INFILE,FPIXELS,LPIXELS,IMG,BG,WORK)
      INTEGER :: FPIXELS(*),LPIXELS(*),NAXES(2)
      INTEGER :: K,NPIXELS,NFRAMES,X,Y,XC,YC
      DOUBLE PRECISION :: IMG(*),BG(*),DA,DSUM
      DOUBLE PRECISION :: WORK(LPIXELS(1)-FPIXELS(1)+1,*)
      CHARACTER*(*) :: INFILE
      NAXES=(/LPIXELS(1)-FPIXELS(1)+1,LPIXELS(2)-FPIXELS(2)+1/)
      NFRAMES=LPIXELS(3)-FPIXELS(3)+1
      NPIXELS=NAXES(1)*NAXES(2)
      XC=INT(FLOOR(0.5*DBLE(NAXES(2)+1)))
      YC=INT(FLOOR(0.5*DBLE(NAXES(1)+1)))
      CALL DSCAL(NPIXELS,DBLE(0),IMG,1)
      DO K=FPIXELS(3),LPIXELS(3)
          CALL READIMAGE(INFILE,(/FPIXELS(1),FPIXELS(2),K/),
     &    (/LPIXELS(1),LPIXELS(2),K/),WORK)
        DSUM=SUM(WORK(1:NAXES(1),1:NAXES(2)))
        CALL DSCAL(NPIXELS,DBLE(NPIXELS)/DSUM,WORK,1)
        CALL DAXPY(NPIXELS,DBLE(-1),BG,1,WORK,1)
        X=MAXLOC(MAXVAL(WORK(1:NAXES(1),1:NAXES(2)),1),2)
        Y=MAXLOC(MAXVAL(WORK(1:NAXES(1),1:NAXES(2)),2),1)
        DSUM=SUM(WORK(1:NAXES(1),1:NAXES(2)))
        CALL DSCAL(NPIXELS,DBLE(NPIXELS)/DSUM,WORK,1)
        WORK(1:NAXES(1),1:NAXES(2))=EOSHIFT(
     &    WORK(1:NAXES(1),1:NAXES(2)),X-XC,DBLE(0),2)
        WORK(1:NAXES(1),1:NAXES(2))=EOSHIFT(
     &    WORK(1:NAXES(1),1:NAXES(2)),Y-YC,DBLE(0),1)
        CALL DAXPY(NPIXELS,DBLE(1),WORK,1,IMG,1)
      END DO
      CALL DSCAL(NPIXELS,DBLE(1)/DBLE(NFRAMES),IMG,1)
      RETURN
      END SUBROUTINE SSAREC

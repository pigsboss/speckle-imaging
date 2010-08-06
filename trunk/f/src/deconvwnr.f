      PROGRAM MAIN
C  Wiener deconvolution executable.
C
C  Usage:
C  ======
C  ./deconvwnr filename_g filename_h filename_f M N SNR [IMGRAD]
C
C  Purpose:
C  ========
C  g=conv(f,h). estimates f given g and h.
C
C  Arguments:
C  ==========
C  filename_g - input file. contains matrix g.
C  filename_h - input file. contains matrix h.
C  filename_f - output file. contains the estimate of matrix f.
C  M - number of rows of g, h, and f.
C  N - number of columns of g, h, and f.
C  SNR --- P2 - uses 2nd polynomials to fit the background then estimates the
C       |       SNR.
C       |
C       |- P4 - the same as P2 but uses 4th polymonials to fit the background.
C       |
C       |- Otherwise - take this argument as pre-estimated SNR.
C  IMGRAD - required only if SNR is not given but should be estimate in this
C           program.
C
C  Declaration:
C  ============
      INTEGER :: M,N
      DOUBLE PRECISION :: DSNR,DIMGRAD
      DOUBLE PRECISION, ALLOCATABLE :: DF(:,:),DG(:,:),DH(:,:)
      CHARACTER :: ARG*80,FILEG*80,FILEH*80,FILEF*80
C
C  Main routine:
C  =============
      CALL GETARG(1,FILEG)
      CALL GETARG(2,FILEH)
      CALL GETARG(3,FILEF)
      CALL GETARG(4,ARG)
      READ(ARG,*) M
      CALL GETARG(5,ARG)
      READ(ARG,*) N
      ALLOCATE(DG(M,N))
      ALLOCATE(DH(M,N))
      CALL READIMAGE(FILEG,(/1,1,1/),(/M,N,1/),DG)
      CALL READIMAGE(FILEH,(/1,1,1/),(/M,N,1/),DH)
      CALL GETARG(6,ARG)
      IF (TRIM(ARG).EQ.'P2')THEN
        CALL GETARG(7,ARG)
        READ(ARG,*) DIMGRAD
        CALL GETSNR(M,N,DIMGRAD,DSNR,DG,2)
      ELSE IF (TRIM(ARG).EQ.'P4')THEN
        CALL GETARG(7,ARG)
        READ(ARG,*) DIMGRAD
        CALL GETSNR(M,N,DIMGRAD,DSNR,DG,4)
      ELSE
        READ(ARG,*) DSNR
      END IF
      WRITE(*,'(A,ES10.3)') 'SNR = ',DSNR
      ALLOCATE(DF(M,N))
      CALL DECONVWNR(M,N,DG,DF,DH,DSNR)
      CALL WRITEIMAGE(FILEF,(/1,1,1/),(/M,N,1/),DF)
      WRITE(*,'(2A)') 'Output: ',TRIM(FILEF)
C
C  Deallocations:
C  ==============
      DEALLOCATE(DG)
      DEALLOCATE(DH)
      DEALLOCATE(DF)
      STOP
      END PROGRAM MAIN

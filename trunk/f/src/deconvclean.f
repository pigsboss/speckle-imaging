      PROGRAM MAIN
C  CLEAN deconvolution executable.
C
C  Usage:
C  ======
C  ./deconvclean filename_g filename_h filename_f M N beta max_numit
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
C  M          - number of rows of g, h, and f.
C  N          - number of columns of g, h, and f.
C  beta       - loop gain.
C  max_numit  - maximum number of iterations.
C
C  Declaration:
C  ============
      INTEGER :: M,N,MNUMIT
      DOUBLE PRECISION :: DBETA
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
      CALL GETARG(6,ARG)
      READ(ARG,*) DBETA
      CALL GETARG(7,ARG)
      READ(ARG,*) MNUMIT
      ALLOCATE(DG(M,N))
      ALLOCATE(DH(M,N))
      CALL READIMAGE(FILEG,(/1,1,1/),(/M,N,1/),DG)
      CALL READIMAGE(FILEH,(/1,1,1/),(/M,N,1/),DH)
      ALLOCATE(DF(M,N))
      CALL DECONVCLEAN(M,N,DG,DF,DH,DBETA,MNUMIT)
      CALL WRITEIMAGE(FILEF,(/1,1,1/),(/M,N,1/),DF)
      WRITE(*,'(2A)') ' Output: ',TRIM(FILEF)
C
C  Deallocations:
C  ==============
      DEALLOCATE(DG)
      DEALLOCATE(DH)
      DEALLOCATE(DF)
      STOP
      END PROGRAM MAIN

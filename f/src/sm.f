      PROGRAM MAIN
C  Speckle masking main program.
C
C  Usage:
C  ======
C  sm file_obs fpixels_obs lpixels_obs radius_obs
C    file_ref fpixels_ref lpixels_ref radius_ref
C    fit_method p_size max_width
C
C  Declaration:
C  ============
      IMPLICIT NONE
      INTEGER :: FPOBS(3),LPOBS(3),FPREF(3),LPREF(3),
     &  P,MW,BUFFERSIZE,M,N,NBISP
      DOUBLE PRECISION :: DROBS,DRREF
      DOUBLE PRECISION, ALLOCATABLE :: DAVG(:,:)
      DOUBLE COMPLEX, ALLOCATABLE :: ZBISPOBS(:),ZBISPREF(:)
      CHARACTER(LEN=256) :: FILEOBS,ARG,FILEREF
      CHARACTER(LEN=10) :: FTMETHOD
C
C  Statements:
C  ===========
C
      CALL GETARG(1,FILEOBS)
      CALL GETARG(2,ARG)
      READ(ARG,*) FPOBS(1),FPOBS(2),FPOBS(3)
      CALL GETARG(3,ARG)
      READ(ARG,*) LPOBS(1),LPOBS(2),LPOBS(3)
      CALL GETARG(4,ARG)
      READ(ARG,*) DROBS
      CALL GETARG(5,REFFILE)
      CALL GETARG(6,ARG)
      READ(ARG,*) FPOBS(1),FPOBS(2),FPOBS(3)
      CALL GETARG(7,ARG)
      READ(ARG,*) LPOBS(1),LPOBS(2),LPOBS(3)
      CALL GETARG(8,ARG)
      READ(ARG,*) DROBS
      CALL GETARG(9,FTMETHOD)
      CALL GETARG(10,ARG)
      READ(ARG,*) P
      CALL GETARG(11,ARG)
      READ(ARG,*) MW
C  The padding size must be even:
      P=INT(CEILING(DBLE(P)/DBLE(2))*DBLE(2))
      NBISP=(MW+1)*(2*P-MW)*P*(P+2)/8
      M=LPOBS(1)-FPOBS(1)+1
      N=LPOBS(2)-FPOBS(2)+1
      ALLOCATE(DAVG(M,N))
      CALL AVERAGE(FILEOBS,FPOBS,LPOBS,DAVG)
      ALLOCATE(ZBISPOBS(NBISP))
      CALL BISPECTRUM(OBSFILE,FPOBS,LPOBS,DAVG,DROBS,FTMETHOD,
     &  P,MW,ZBISPREF)
      M=LPREF(1)-FPREF(1)+1
      N=LPREF(2)-FPREF(2)+1
      DEALLOCATE(DAVG)
      ALLOCATE(DAVG(M,N))
      CALL AVERAGE(FILEREF,FPREF,LPREF,DAVG)
      ALLOCATE(ZBISPREF(NBISPREF))
      CALL BISPECTRUM(REFFILE,FPREF,LPREF,DAVG,DRREF,FTMETHOD,
     &  P,MW,ZBISPREF)
      DEALLOCATE(ZBISPOBS)
      DEALLOCATE(ZBISPREF)
      STOP
      END PROGRAM MAIN
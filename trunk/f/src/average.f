      PROGRAM MAIN
      INTEGER :: STATUS,UNIT,FPIXELS(3),LPIXELS(3),M,N,NPIXELS,I
      DOUBLE PRECISION, ALLOCATABLE :: ARRAY(:)
      CHARACTER(LEN=256) :: FILENAME,ARG,AVGFILE,DATESTR
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
      CALL GETARG(1,FILENAME)
      WRITE(*,'(2A)') ' Input: ',TRIM(FILENAME)
      CALL GETARG(2,ARG)
      READ(ARG,*) FPIXELS(1),FPIXELS(2),FPIXELS(3)
      CALL GETARG(3,ARG)
      READ(ARG,*) LPIXELS(1),LPIXELS(2),LPIXELS(3)
      CALL GETARG(4,AVGFILE)
      M=LPIXELS(1)-FPIXELS(1)+1
      N=LPIXELS(2)-FPIXELS(2)+1
      NPIXELS=M*N
      ALLOCATE(ARRAY(NPIXELS))
      I=TIME()
      CALL CTIME(I,DATESTR)
      PRINT *,TRIM(DATESTR),' start averaging.'
      CALL AVERAGE(TRIM(FILENAME),FPIXELS,LPIXELS,ARRAY)
      I=TIME()
      CALL CTIME(I,DATESTR)
      PRINT *,TRIM(DATESTR),' finished averaging.'
      CALL WRITEIMAGE(TRIM(AVGFILE),(/1,1,1/),(/M,N,1/),ARRAY)
      WRITE(*,'(2A)') ' Output: ',TRIM(AVGFILE)
      DEALLOCATE(ARRAY)
      STOP
      END PROGRAM

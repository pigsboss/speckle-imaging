      PROGRAM BISPOSTEST
      IMPLICIT NONE
      INTEGER :: STATUS,C,K,X1,Y1,X2,Y2,NX,NY,Y2MAX,LB,NARGS
      INTEGER, PARAMETER :: UNIT=8,DEFAULTNX=16
      CHARACTER(LEN=256) :: ARG,STRC,STRK
C
      INTERFACE
      FUNCTION BISPOS(X1,Y1,X2,Y2,NX,NY)
      INTEGER, INTENT(IN) :: X1,Y1,X2,Y2,NX,NY
      INTEGER :: BISPOS
      END FUNCTION BISPOS
      END INTERFACE
C
      STATUS=0
      OPEN(UNIT=UNIT,FILE='bispostest.txt',STATUS='REPLACE',
     &  ACTION='WRITE',IOSTAT=STATUS)
      IF(STATUS.NE.0)THEN
        PRINT *,'open file failed.'
        RETURN
      END IF
      NX=DEFAULTNX
      NY=NX
      Y2MAX=NY
      NARGS=COMMAND_ARGUMENT_COUNT()
      IF(NARGS .GE. 1)THEN
        CALL GET_COMMAND_ARGUMENT(1,ARG)
        READ(ARG,*)NX
        NY=NX
        Y2MAX=NY
      END IF
      IF(NARGS .GE. 2)THEN
        CALL GET_COMMAND_ARGUMENT(2,ARG)
        READ(ARG,*)NY
        Y2MAX=NY
      END IF
      IF(NARGS .GE. 3)THEN
        CALL GET_COMMAND_ARGUMENT(3,ARG)
        READ(ARG,*)Y2MAX
      END IF
      WRITE(*,'(A,I3,A,I3)')' size: ',NX,' x ',NY
      WRITE(UNIT,'(A,I3,A,I3)')' size: ',NX,' x ',NY
      WRITE(*,'(A,I3)')' upper limit of y_2: ',Y2MAX
      WRITE(UNIT,'(A,I3)')' levels: ',Y2MAX
      LB=(Y2MAX+1)*(2*NY-Y2MAX)*NX*(NX+2)/8
      WRITE(*,*)'length of bispectrum array: ',LB
      WRITE(UNIT,*)'length of bispectrum array: ',LB
      K=1
      DO Y2=0,Y2MAX
        WRITE(*,'(A,I3,A,I3)')' level ',Y2,' of ',Y2MAX
        DO X2=0,NX-1
          DO Y1=0,NY-1-Y2
            DO X1=0,MIN(X2,NX-1-X2)
              C=BISPOS(X1,Y1,X2,Y2,NX,NY)
              WRITE(STRC,*)C
              WRITE(STRK,*)K
              STRC=ADJUSTL(STRC)
              STRK=ADJUSTL(STRK)
              IF(C .NE. K)THEN
                WRITE(*,'(A,I3,A,I3,A,I3,A,I3,A)')
     &            ' X (',X1,',',Y1,',',X2,',',Y2,'): '//
     &            TRIM(STRC)//'('//TRIM(STRK)//')'
                WRITE(UNIT,'(A,I3,A,I3,A,I3,A,I3,A)')
     &            ' X (',X1,',',Y1,',',X2,',',Y2,'): '//
     &            TRIM(STRC)//'('//TRIM(STRK)//')'
C             ELSE
C               WRITE(*,'(A,I3,A,I3,A,I3,A,I3,A)')
C    &            ' O (',X1,',',Y1,',',X2,',',Y2,'): '//
C    &            TRIM(STRC)//'('//TRIM(STRK)//')'
C               WRITE(UNIT,'(A,I3,A,I3,A,I3,A,I3,A)')
C    &            ' O (',X1,',',Y1,',',X2,',',Y2,'): '//
C    &            TRIM(STRC)//'('//TRIM(STRK)//')'
              END IF
              K=K+1
            END DO
          END DO
        END DO
      END DO
      CLOSE(UNIT)
      STOP
      END PROGRAM BISPOSTEST

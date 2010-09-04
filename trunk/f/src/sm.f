      PROGRAM SM

      ALLOCATE(ZBISP(LBISP),STAT=STATUS)
      IF(STATUS.NE.0)THEN
        WRITE(*,*)'error: out of memory.'
        WRITE(UNIT,*)'error: out of memory.'
        RETURN
      END IF
      WRITE(*,'(A,F7.1)')' size of bispectral array: ',
     &  DBLE(LBISP*16)/DBLE(1024*1024),'MB'
      WRITE(UNIT,'(A,F7.1)')' size of bispectral array: ',
     &  DBLE(LBISP*16)/DBLE(1024*1024),'MB'
      
      END PROGRAM SM
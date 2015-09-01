!=======================================================================

    SUBROUTINE PLBLK(Z,NR,R0,LABR,NLIST)

! print a lower-triangular packed matrix.
! this version prints eight columns at a time.
! parameter ncol and format 10 should be modified to print
! a different number of columns.

! input:
! Z(*) = matrix to be printed.
! NR   = row and column dimension.
! R0   = row number offset.
! LABR = character row and column label.
! NLIST= output unit nubmer.

    IMPLICIT INTEGER(A-Z)

    PARAMETER (NCOL=6)
    10 FORMAT(/8X,6(6X,A4,I4,1X))
    1 FORMAT(1X,A4,I4,6F15.8)

    real*8 :: ZERO
    PARAMETER(ZERO=0D0)

    CHARACTER*(*) LABR
    real*8 :: Z(*)

!    ASSIGN 1 TO FMTZ

    JLAST=0
    DO 400 JSTRT=1,NR,NCOL
        JLAST=MIN(NR,JLAST+NCOL)
    
        JLAB1=JSTRT+R0
        JLAB2=JLAST+R0
        WRITE(NLIST,10)(LABR,J,J=JLAB1,JLAB2)
    
        IJ0=(JSTRT*(JSTRT-1))/2
        DO 300 I=JSTRT,NR
            ILAB=I+R0
            J2=MIN(I,JLAST)
        
        ! print the row if a nonzero element is found.
        
            DO 100 J=JSTRT,J2
                IF(Z(IJ0+J) /= ZERO)THEN
                    WRITE(NLIST,1)LABR,ILAB,(Z(IJ0+JT),JT=JSTRT,J2)
                    GO TO 101
                ENDIF
            100 end do
            101 IJ0=IJ0+I
        300 end do
    400 end do

    RETURN
    end SUBROUTINE PLBLK

! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


!=======================================================================

    SUBROUTINE PRBLK(Z,NRD,NR,NC,R0,C0,LABR,LABC)

! print a sub-block of a rectangular matrix.
! this version prints eight columns at a time.
! parameter ncol and format 10 should be modified to print
! a different number of columns.

! input:
! Z(*) = matrix to be printed.
! NRD  = row dimension.
! NR   = number of rows to print.
! NC   = column dimension.
! R0   = row number offset.
! C0   = column number offset.
! LABR = character row label.
! LABC = character column label.

    IMPLICIT INTEGER(A-Z)

    PARAMETER (NCOL=6)
    10 FORMAT(8X,6(3X,A4,I3))
    1 FORMAT(1X,A4,I4,6F10.5)

    real*8 :: ZERO
    PARAMETER (ZERO=0D0)

    CHARACTER*(*) LABR,LABC
    real*8 :: Z(NRD,NC)

!    ASSIGN 1 TO FMTZ

    JLAST=0
    DO 400 JSTRT=1,NC,NCOL
        JLAST=MIN(NC,JLAST+NCOL)
    
        JLAB1=JSTRT+C0
        JLAB2=JLAST+C0
        write(*,*)' '
        WRITE(*,10)(LABC,J,J=JLAB1,JLAB2)
    
        DO 300 I=1,NR
            ILAB=I+R0
        
        ! print the row if a nonzero element is found.
        
            DO 100 J=JSTRT,JLAST
            !C               IF(Z(I,J).NE.ZERO)THEN
                WRITE(*,1)LABR,ILAB,(Z(I,JT),JT=JSTRT,JLAST)
                GO TO 300
            !C               ENDIF
            100 end do
        300 end do
    400 end do

    RETURN
    end SUBROUTINE PRBLK


!=======================================================================








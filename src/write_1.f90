! ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    subroutine write_1(iunit,filen,norb,ecore,r1int,tol,nsite)

! This routine writes all the nonzero one-electron integrals to the
! specified file.
! ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!    implicit real*8(a-h,o-z)
!    implicit integer*4(i-n)

!    dimension r1int((norb*(norb+1))/2)
IMPLICIT NONE

    integer,intent(in)::norb,nsite,iunit
    real(kind=8),intent(in)::ecore,tol,r1int((norb*(norb+1))/2)
    integer::n12,nzero,ij
    real(kind=8)::x
    character*(*) filen

    open(unit=iunit,file=filen,form='unformatted',status='unknown', &
    err=101)

    n12=(norb*(norb+1))/2

! First a dry run to find how many nonzero elements are there

    nzero=0
    do ij=1,n12
    
        if(dabs(r1int(ij)) >= tol)nzero=nzero+1
    
    end do

    write(iunit)norb,nsite,nzero

    do ij=1,n12
    
        x=r1int(ij)
        if(dabs(x) >= tol)write(iunit)ij,x
    
    end do

    write(iunit)ecore

    close(unit=iunit,status='keep')

    goto 102

    101 print*,'Error in WRITE_1','Error opening one-electron integral file'

    102 continue

    return
    end subroutine write_1

! ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

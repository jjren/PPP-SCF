! ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    subroutine write_2(iunit,filen,norb,r2int,tol)


! This routine writes all the nonzero two-electron integrals to the
! specified file.
! ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
IMPLICIT NONE

    integer,intent(in)::iunit,norb
    real(kind=8),intent(in)::tol,r2int((norb*(norb+1)*(norb*norb+norb+2))/8)
    integer::nij,ijkl,ij,kl
    real(kind=8)::x
    character*(*) filen

    open(unit=iunit,file=filen,form='unformatted',status='unknown', &
    err=101)

    nij=(norb*(norb+1))/2

    ijkl=0
    do ij=1,nij
        do kl=1,ij
        
            ijkl=ijkl+1
            x=r2int(ijkl)
            if(dabs(x) >= tol)write(iunit)ij,kl,x
        
        end do
    end do

    close(unit=iunit,status='keep')

    goto 102

    101 print*,'Error in WRITE_2','Error opening two-electron integral file'

    102 continue

    return
    end subroutine write_2

! ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

! ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    subroutine dipout(iunit,filen,norb,r1int,idir,ndir,tol)

! This routine writes all the dipole integrals to the
! specified file.
! ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
IMPLICIT NONE

    integer,intent(in)::ndir,norb,iunit,idir(ndir)
    real(kind=8),intent(in)::r1int((norb*(norb+1))/2,ndir),tol
    character*(*) filen
    integer::n12,jdir,iaxis,ij,nzero
    real(kind=8)::x

    open(unit=iunit,file=filen,form='unformatted',status='unknown', &
    err=101)

    n12=(norb*(norb+1))/2

    write(iunit)ndir,norb

    do jdir=1,ndir
    
        iaxis=idir(jdir)
    
    ! First a dry run to find how many nonzero elements are there
    
 
        nzero=0
    
        do ij=1,n12
        
  
            if(dabs(r1int(ij,jdir)) >= tol)nzero=nzero+1
        
        end do
    
        write(iunit)iaxis,nzero
    
        do ij=1,n12
        
            x=r1int(ij,jdir)
            if(dabs(x) >= tol)write(iunit)ij,x
 
        
        end do
    
    end do

    close(unit=iunit,status='keep')

    goto 102

    101 print*,'Error in DIPOUT','Error opening dipole integral file'

    102 continue

    return
    end subroutine dipout

! ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

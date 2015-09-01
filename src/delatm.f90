! ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    subroutine delatm(rr_site,idelatm,swap,ndelatm,nsitex)

! This subroutine deletes atoms from the list
! cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!    implicit real*8(a-h,o-z)
!    implicit integer*4(i-n)

!    dimension rr_site(3,nsitex),swap(3,nsitex-ndelatm)

!    dimension idelatm(ndelatm)

IMPLICIT NONE
    integer,intent(in)::idelatm(*),ndelatm,nsitex
    real(kind=8),intent(out)::rr_site(3,*)
    real(kind=8),intent(inout)::swap(3,nsitex-ndelatm)
    integer::isite,iflag,ix,iatm,jsite,nsite,i

    nsite=nsitex-ndelatm

    isite=0
    do jsite=1,nsitex
        iflag=0
    
        do ix=1,ndelatm
            iatm=idelatm(ix)
            if(iatm == jsite)iflag=iflag+1
        end do
    
        if(iflag == 0)then
            isite=isite+1
            swap(1,isite)=rr_site(1,jsite)
            swap(2,isite)=rr_site(2,jsite)
            swap(3,isite)=rr_site(3,jsite)
        end if
    
    end do

    if(isite /= nsite)then
        print*,'isite,nsite',isite,nsite
        print*,'Error:DELATM','mismatch in the site count'
    end if

    do isite=1,nsite
        do i=1,3
            rr_site(i,isite)=swap(i,isite)
        end do
    end do

    return
    end subroutine delatm

! cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

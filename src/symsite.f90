! ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! This routine symmetrizes the coordinates of the atoms to place the origin
! of the molecule at the center of mass of the molecule.

    subroutine symsite(rr_site,nsite)

! ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
implicit none


    integer,intent(in)::nsite
    real(kind=8),intent(inout)::rr_site(3,nsite)
    integer::isite
    real(kind=8)::x,y,z

! Now symmetrize the system, i.e., place the origin at the center of
! mass

    x=0.d0
    y=0.d0
    z=0.d0
    do isite=1,nsite
        x=x+rr_site(1,isite)
        y=y+rr_site(2,isite)
        z=z+rr_site(3,isite)
    end do

    x=x/dble(nsite)
    y=y/dble(nsite)
    z=z/dble(nsite)
! print*,'x,y,z',x,y,z

    do isite=1,nsite
        rr_site(1,isite)=rr_site(1,isite)-x
        rr_site(2,isite)=rr_site(2,isite)-y
        rr_site(3,isite)=rr_site(3,isite)-z
    end do

    return
    end subroutine symsite

! ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

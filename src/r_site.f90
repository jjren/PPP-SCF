! ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! This routine generates the coordinates of all the pi-electron sites in
! the molecule from the coordinates of origins of the unit cells and
! the locations of pi-electron sites within a given cell.

    subroutine r_site(rr_site,rr_atom,tvec,nsite,natom,ncell)

! ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!    implicit real*8(a-h,o-z)
!    implicit integer*4(i-n)

!    dimension rr_site(3,nsite)
!    dimension rr_atom(3,natom)
!    dimension tvec(3,1)
IMPLICIT NONE
   integer,intent(in)::nsite,natom,ncell
   real(kind=8),intent(in)::rr_atom(3,natom),tvec(3,1)
   real(kind=8),intent(out)::rr_site(3,nsite)
   integer::isite,iatom,icell
   real(kind=8)::x,y,z
! Note that the sites are numbered in the order of atomic positions
! and the cell positions provided by the user.

    isite=0
    x=0.d0
    y=0.d0
    z=0.d0
        
    do icell=1,ncell
    
    
        do iatom=1,natom
        
            isite=isite+1
            if(isite > nsite)then
                print*,'isite,nsite',isite,nsite
                print*,'Error:R_site', &
                'No. of generaed sites exceeds the limit'
                stop
            end if
        
            rr_site(1,isite)=x+rr_atom(1,iatom)
            rr_site(2,isite)=y+rr_atom(2,iatom)
            rr_site(3,isite)=z+rr_atom(3,iatom)
        
        end do
    
        x=x+tvec(1,1)
        y=y+tvec(2,1)
        z=z+tvec(3,1)
    
    end do

!    do isite=1,nsite
!      print*,isite,rr_site(1,isite),rr_site(2,isite),rr_site(3,isite)
!    end do

    if(isite /= nsite)then
        print*,'isite,nsite',isite,nsite
        print*,'Error:R_site', &
        'No. of generaed sites not equal to nsite'
        stop
    end if


    return
    end subroutine r_site

! ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

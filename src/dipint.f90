! ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    subroutine dipint(z,r_site,orb,idir,msite,morb,mdir)

! This routine transforms the dipole operator from the site representation
! into the representation of the orbitals contained in the array ORB.
! The transformed dipole matrix is kept in the array Z upon output.
! cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
IMPLICIT NONE



   integer,intent(in)::morb,mdir,idir(mdir),msite
   real(kind=8),intent(in)::r_site(3,msite),orb(msite,morb)
   real(kind=8),intent(out)::z((morb*morb+morb)/2,mdir)
   integer::ij,i,j,iaxis,isite,jdir
   
    ij=0
    do i=1,morb
        do j=1,i
        
            ij=ij+1
        
            do iaxis=1,mdir
                z(ij,iaxis)=0.d0
            end do
        
            do isite=1,msite
                do jdir=1,mdir
                    iaxis=idir(jdir)
                    z(ij,jdir)=z(ij,jdir)+orb(isite,i)* &
                    orb(isite,j)*r_site(iaxis,isite)
                end do
            end do
        
        end do
    end do

    return
    end subroutine dipint

! ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

! cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! Density matrix construction for a closed-shell system.

    subroutine density(den,orb,norb,nsite)

! cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

IMPLICIT NONE

!    implicit real*8(a-h,o-z)
!    implicit integer*4(i-n)

!    dimension den((nsite*(nsite+1))/2)
!    dimension orb(nsite,norb)

     integer::ij,i,j,ii,iorb
     integer,intent(in)::norb,nsite
     real(kind=8),intent(out)::den((nsite*(nsite+1))/2)
     real(kind=8),intent(in)::orb(nsite,*)
     real(kind=8)::fac
     
    ij=0

    do i=1,nsite
    
        do j=1,i
        
            ij=ij+1
            den(ij)=0.d0
        
            fac=4.d0
            if(i == j)fac=2.d0
            do iorb=1,norb
            
                den(ij)=den(ij)+fac*orb(i,iorb)*orb(j,iorb)
            
            end do
        
        end do
    
    end do


    return
    end subroutine density

! cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


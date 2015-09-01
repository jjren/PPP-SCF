! ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! Computes the Fourier transform of 1-d nearest-neighbor-interacting
! tight-binding Hamiltonian.


    subroutine foutra_tb(fkval,tij,tij_rl,irij,fk,matom)

! ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

IMPLICIT NONE
    
    integer,intent(in)::matom,irij((matom*(matom+1))/2)
    real(kind=8),intent(in)::tij((matom*(matom+1))/2),fkval,&
               tij_rl((matom*(matom+1))/2)   
    complex*16,intent(out)::fk((matom*(matom+1))/2)
    integer::nij,ij,i,j
    real(kind=8)::pi,ang,xcos,xsin,x,y
! pi

    pi=4.d0*datan(1.d0)

    nij=(matom*(matom+1))/2
    ij=0
    do i=1,matom
        do j=1,i
        
            ij=ij+1
            ang=pi*fkval*dble(irij(ij))
            xcos=dcos(ang)
            xsin=dsin(ang)
            x=tij(ij)+xcos*tij_rl(ij)
            y=xsin*tij_rl(ij)
            fk(ij)=cmplx(x,y)
         
        end do
    end do

    return
    end subroutine foutra_tb

! cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

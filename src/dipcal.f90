! cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    subroutine dipcal(r_site,orb,dipmat,idip,nbas,norb,ndip)

! This routine computes the dipole moment matrix corresponding to
! the single-particle orbitals for the PPP Hamiltonian
! cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


IMPLICIT NONE
   
   integer,intent(in)::idip(ndip),nbas,norb,ndip
   real(kind=8),intent(in)::r_site(3,nbas),orb(nbas,norb)
   real(kind=8),intent(out)::dipmat((norb*(norb+1))/2,ndip)
   integer::ij,i,j,jdip,k,ii
   real(kind=8)::thr,xxx


!    call prblk(orb,nbas,nbas,norb,0,0,'bas','orb',6)

    thr=1.d-10
    ij=0
    do i=1,norb
        do j=1,i
            ij=ij+1
            do jdip=1,ndip
            
                k=idip(jdip)
                dipmat(ij,jdip)=0.d0
                xxx=0.d0
            
                do ii=1,nbas
                
                    xxx=xxx+orb(ii,i)*orb(ii,j)* &
                    r_site(k,ii)
                
                end do
            
                if(abs(xxx) > thr)dipmat(ij,jdip)=xxx
            
            end do
        end do
    end do



    return
    end subroutine dipcal

! cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

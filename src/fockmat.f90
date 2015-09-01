! cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! This subroutine constructs the Fock matrix and the SCF energy from the
! density matrix and the one- and two-electron integrals.

    subroutine fockmat(t,h,v,den,fock,nsite,escf)

! cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

IMPLICIT NONE

!    implicit real*8(a-h,o-z)
!    implicit integer*4(i-n)

!    dimension t((nsite*(nsite+1))/2)
!    dimension h((nsite*(nsite+1))/2)
!    dimension v((nsite*(nsite+1))/2)
!    dimension g((nsite*(nsite+1))/2)
!    dimension fock(nsite,nsite)
!    dimension den((nsite*(nsite+1))/2)

    integer::ij,i,j,ii,jj
    integer,intent(in)::nsite
    real(kind=8),intent(in)::t((nsite*(nsite+1))/2),h((nsite*(nsite+1))/2),&
                  v((nsite*(nsite+1))/2),den((nsite*(nsite+1))/2)
    real(kind=8),intent(out)::fock((nsite*(nsite+1))/2),escf
!    real(kind=8)escf,xen1
 
      
      fock=0.d0


! Construct the electron repulsion (ER) matrix (=2J- K, where
! J=Coulomb matrix, K=Exchange Matrix)

    do i=1,nsite
    
        ii=(i*(i-1))/2+i
    
    ! i.ne.j part
    
        do j=1,i-1
        
            jj=(j*(j-1))/2+j
            ij=(i*(i-1))/2+j
        
        ! Coulomb contribution
        
            fock(ii)=fock(ii)+v(ij)*den(jj)
            fock(jj)=fock(jj)+v(ij)*den(ii)
        ! goto 3
        
        ! Exchange contribution
        
            fock(ij)=fock(ij)-0.25d0*v(ij)*den(ij)
        ! 51        continue
        
        end do
    
    ! i.eq.j part
    
    ! Coulomb contribution
    
        fock(ii)=fock(ii)+v(ii)*den(ii)
    
    ! goto 4
    
    ! Exchange contribution
    
        fock(ii)=fock(ii)-0.5d0*v(ii)*den(ii)
    
    end do
!
! Add up all contributions to get the Fock matrix
! 
    
      fock=fock+h

!get the SCF energy    
!    escf=dot_product(den,h)+0.5d0*dot_product(den,fock)

  escf=dot_product(den,h)+dot_product(den,fock)
  escf= 0.5d0*escf   
  
!
! Add up all contributions to get the Fock matrix
! 
    
!      fock=fock+h

    end subroutine fockmat

! cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc



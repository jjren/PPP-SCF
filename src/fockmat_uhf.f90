! This subroutine constructs the Fock matrix for open shell systems and the SCF energy 
!from the density matrix and the one- and two-electron integrals.
!
      subroutine fockmat_uhf(t,h,v,den,den_alpha,den_beta,fock_alpha,fock_beta,nsite,escf)
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
IMPLICIT NONE
INTEGER::i,j,ij,ii,jj
INTEGER,INTENT(IN)::nsite
REAL(KIND=8),INTENT(IN)::t((nsite*(nsite+1))/2),h((nsite*(nsite+1))/2),v((nsite*(nsite+1))/2),&
                         den_alpha((nsite*(nsite+1))/2),den_beta((nsite*(nsite+1))/2),den((nsite*(nsite+1))/2)
REAL(KIND=8),INTENT(OUT)::fock_alpha((nsite*(nsite+1))/2),fock_beta((nsite*(nsite+1))/2),escf

fock_alpha=0.d0
fock_beta=0.d0
do i=1,nsite
      ii=(i*(i-1))/2+i
!
!i.ne.j part
!
   do j=1,i-1
        jj=(j*(j-1))/2+j
        ij=(i*(i-1))/2+j
!
!Coulomb contribution
!
       fock_alpha(ii)=fock_alpha(ii)+v(ij)*den(jj)
       fock_alpha(jj)=fock_alpha(jj)+v(ij)*den(ii)
       fock_beta(ii)=fock_beta(ii)+v(ij)*den(jj)
       fock_beta(jj)=fock_beta(jj)+v(ij)*den(ii)
!
!Exchange contribution
!       
       fock_alpha(ij)=fock_alpha(ij)-0.50d0*v(ij)*den_alpha(ij)
       fock_beta(ij)=fock_beta(ij)-0.50d0*v(ij)*den_beta(ij)
    end do
!
!i.eq.j part
!
! Coulomb contribution
! 
        fock_alpha(ii)=fock_alpha(ii)+v(ii)*den(ii)
        fock_beta(ii)=fock_beta(ii)+v(ii)*den(ii)	
!
! Exchange contribution
!
        fock_alpha(ii)=fock_alpha(ii)-v(ii)*den_alpha(ii)
        fock_beta(ii)=fock_beta(ii)-v(ii)*den_beta(ii)
end do
!    
! Add up all contributions to get the Fock matrix
!
             fock_alpha=fock_alpha+h
	     fock_beta=fock_beta+h
!
!and also get the SCF energy
!
          escf=dot_product(den_alpha,fock_alpha)+dot_product(den_beta,fock_beta)+dot_product(den,h)
	  escf=0.5d0*escf
! 
return
end subroutine fockmat_uhf 

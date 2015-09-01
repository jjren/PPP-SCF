! Density matrix construction for a open-shell system.
!
  subroutine density_uhf(den,den_alpha,den_beta,orb_alpha,orb_beta,nalpha,nbeta,nsite)
!
!nalpha = no. of up spins
!nbeta = no. of down spins
!
!##########################################################################################################
!
IMPLICIT NONE
INTEGER::iorb,i,j,ij,ii
INTEGER,INTENT(IN)::nalpha,nbeta,nsite
REAL(KIND=8),INTENT(IN)::orb_alpha(nsite,*),orb_beta(nsite,*)
REAL(KIND=8),INTENT(OUT)::den((nsite*(nsite+1))/2),den_alpha((nsite*(nsite+1))/2),den_beta((nsite*(nsite+1))/2)

ij=0

do i=1,nsite
 do j=1,i
   ij=ij+1
   den_alpha(ij)=0.d0
   den_beta(ij)=0.d0
   do iorb=1,nalpha
      den_alpha(ij)=den_alpha(ij)+2.d0*orb_alpha(i,iorb)*orb_alpha(j,iorb)
   end do
   do iorb=1,nbeta
      den_beta(ij)=den_beta(ij)+2.d0*orb_beta(i,iorb)*orb_beta(j,iorb)
   end do

 end do
 ii=(i*(i-1))/2+i
 den_alpha(ii)=0.5d0*den_alpha(ii)
 den_beta(ii)=0.5d0*den_beta(ii)
end do 
den=den_alpha+den_beta
!
!
return
end subroutine density_uhf 


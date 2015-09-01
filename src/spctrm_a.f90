! ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! This routine computes the absorption spectrum according to the Wigner
! formula if supplied with the energy levels, dipole moments and the
! line width parameter (gamma),lower frequency limit (omega1), upper
! frequency limit (omega2) and the frequency increment (domga).


    subroutine spctrm_a(omega_i,omega_f,nlvl_i,nlvl_f,dip, &
    omega1,omega2,domega,gamma,iwrite,scale,ndim)

! ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
implicit none

   integer,intent(in)::nlvl_i,nlvl_f,ndim,iwrite
   real(kind=8),intent(in)::omega_i(nlvl_i),omega_f(nlvl_f),&
               dip(nlvl_i*nlvl_f,ndim),omega1,omega2,&
	       domega,gamma,scale
   integer::nomega,iomega,ij,i,j,idim
   real(kind=8)::omega,alpha,alphax,alphay,alphaz,dsqij,omegaij,&
                dnum,denom


    nomega=nint((omega2-omega1)/domega)+1

    omega=omega1
    do iomega=1,nomega
        alpha=0.0
        alphax=0.0
        alphay=0.0
        alphaz=0.0
        ij=0
        do i=1,nlvl_i
            do j=1,nlvl_f
                ij=ij+1
		dsqij=0.d0
		do idim=1,ndim
                   dsqij=dsqij+dip(ij,idim)**2
		end do
                omegaij=(omega_f(j)-omega_i(i))/scale
                dnum=2.d0*omegaij*dsqij*gamma
                denom=(omega-omegaij)**2+gamma**2
                alpha=alpha+(dnum/denom)
		if(ndim>1)then
                alphax=alphax+(2.d0*omegaij*(dip(ij,1)**2)*gamma)/denom
                alphay=alphay+(2.d0*omegaij*(dip(ij,2)**2)*gamma)/denom
		if(ndim==3)alphaz=alphaz+(2.d0*omegaij*(dip(ij,3)**2)*gamma)/denom
                end if
            end do
        end do
    
        write(iwrite,101)omega,alpha
        write(33,101)omega,alphax
        write(34,101)omega,alphay
	if(ndim==3)write(35,101)omega,alphaz
        omega=omega+domega
    
    end do

    101 format(' ',f6.2,2x,f9.4)
    
    close(iwrite)
    close(33)
    close(34)
    close(35)
    return
    end subroutine spctrm_a

! cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

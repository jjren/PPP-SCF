! ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! This routine computes the absorption spectrum from a singly-excited state
! (as compared to the ground state). Lorentzian formula is used to compute
! the spectrum.


    subroutine spctrm_1ex(eval,nbas,i_hole,i_part,nocc,dip,ndim, &
    omega1,omega2,domega,gamma,iwrite,scale)

! ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
implicit none

   integer,intent(in)::iwrite,nbas,i_hole,i_part,nocc,&
                   ndim 
   real(kind=8),intent(in)::eval(nbas),dip((nbas*(nbas+1))/2,ndim),&
              omega1,omega2,domega,gamma,scale
   integer::nomega,iomega,ij,i,j,ivirt,idim,iocc
   real(kind=8)::omega,alpha,alphax,alphay,alphaz,dsqij,omegaij,&
                dnum,denom,fac


    nomega=nint((omega2-omega1)/domega)+1


    write(iwrite,*)'#'


    omega=omega1
    do  iomega=1,nomega
        alpha=0.d0
    
    ! Contribution from singly-excited states with the holes in the same state as the
    ! reference state
    
        do  ivirt=i_part+1,nbas
            ij=(ivirt*(ivirt-1))/2+i_part
            dsqij=0.d0
		do idim=1,ndim
                   dsqij=dsqij+dip(ij,idim)**2
		end do	    
            omegaij=(eval(ivirt)-eval(i_part))/scale
            dnum=omegaij*dsqij*gamma
            denom=(omega-omegaij)**2+gamma**2
            alpha=alpha+(dnum/denom)
        end do
    
    ! Contribution from singly-excited states with the particle in the same state as the
    ! reference state
    
        do  iocc=1,i_hole-1
            ij=(i_hole*(i_hole-1))/2+iocc
              dsqij=0.d0
		do idim=1,ndim
                   dsqij=dsqij+dip(ij,idim)**2
		end do	    
            omegaij=(eval(i_hole)-eval(iocc))/scale
            dnum=omegaij*dsqij*gamma
            denom=(omega-omegaij)**2+gamma**2
            alpha=alpha+(dnum/denom)
        end do
    
    ! Contribution from the doubly excited configuration and singlet (and triplet)
    ! intermediate coupling
    
        do ivirt=nocc+1,nbas
            do  iocc=1,nocc
                fac=2.d0
                if(ivirt /= i_part .AND. iocc /= i_hole)fac=2.d0+6.d0
                ij=(ivirt*(ivirt-1))/2+iocc
                dsqij=0.d0
		do idim=1,ndim
                   dsqij=dsqij+dip(ij,idim)**2
		end do	    		
                omegaij=(eval(ivirt)-eval(iocc))/scale
                dnum=omegaij*dsqij*gamma
                denom=(omega-omegaij)**2+gamma**2
                alpha=alpha+(dnum/denom)
            end do
        end do

        write(iwrite,101)omega,alpha
        omega=omega+domega
    
    end do
    write(iwrite,*)' '

    101 format(' ',f6.2,2x,f9.4)

    close(iwrite)
    return
    end subroutine spctrm_1ex

! cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

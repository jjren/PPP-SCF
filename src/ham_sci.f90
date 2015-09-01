! ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! This routine constructs the single-CI Hamiltonian matrix for the
! singlet subspace. A closed-shell reference state is assumed.

    subroutine ham_sci(orb,v,ham,eval,morb,mocc_orb,mvirt_orb,msite,ndimh)

! ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


IMPLICIT NONE
   integer,intent(in)::msite,mocc_orb,mvirt_orb,morb,ndimh
   real(kind=8),intent(in)::eval(*),v((msite*(msite+1))/2),&
                            orb(msite,*)
   real(kind=8),intent(out)::ham((ndimh*(ndimh+1))/2)
   integer::ij,iocc,ivirt,isite,iisite,jsite,jjsite,&
            i,k,ik,jl,j,ok,l,ijsite,ntorb,ii,iii,iik,jjl,ijkl
   real(kind=8)::dir,exc,lmin,kvirt,lvirt,xham
 
! Dimension checks

    ntorb=mocc_orb+mvirt_orb
    ham=0.d0
 
    ij=0
    do iocc=1,mocc_orb
        do ivirt=mocc_orb+1,morb
        
            ij=ij+1
        
            dir=0.d0
            exc=0.d0
        
            do isite=1,msite
            
                iisite=(isite*(isite-1))/2
            
                do jsite=1,isite
                
                    ijsite=iisite+jsite
                
                    dir=dir+orb(isite,iocc)*orb(isite,ivirt)* &
                    orb(jsite,iocc)*orb(jsite,ivirt)*v(ijsite)
                
                    exc=exc+orb(isite,iocc)*orb(isite,iocc)* &
                    orb(jsite,ivirt)*orb(jsite,ivirt)*v(ijsite)

                end do
            
                do jsite=isite+1,msite
                
                    ijsite=(jsite*(jsite-1))/2+isite
                
                    dir=dir+orb(isite,iocc)*orb(isite,ivirt)* &
                    orb(jsite,iocc)*orb(jsite,ivirt)*v(ijsite)
                
                    exc=exc+orb(isite,iocc)*orb(isite,iocc)* &
                    orb(jsite,ivirt)*orb(jsite,ivirt)*v(ijsite)
                
                end do
            
            end do
	    
	    ii=ij+1
            iii=(ii*(ii-1))/2+ii
            ham(iii)=eval(ivirt)-eval(iocc)+2.d0*dir-exc
        
        end do
    end do

! Off-diagonal elements (Eq. 2.14)

    ik=0
    do i=1,mocc_orb
        do k=1,mvirt_orb
        
            ik=ik+1
            kvirt=k+mocc_orb
        
            do j=i,mocc_orb
            
                lmin=k+1
                if(j /= i)lmin=1
            
                do l=lmin,mvirt_orb
                
                    jl=(j-1)*mvirt_orb+l
                    lvirt=l+mocc_orb
                
                    dir=0.d0
                    exc=0.d0
                
                    do isite=1,msite
                        iisite=(isite*(isite-1))/2
                        do jsite=1,isite
                        
                            ijsite=iisite+jsite
                            dir=dir+orb(isite,i)*orb(isite,kvirt)* &
                            orb(jsite,j)*orb(jsite,lvirt)*v(ijsite)
                        
                            exc=exc+orb(isite,i)*orb(isite,j)* &
                            orb(jsite,kvirt)*orb(jsite,lvirt)*v(ijsite)
                        
                        end do
                    
                        do jsite=isite+1,msite
                        
                            ijsite=(jsite*(jsite-1))/2+isite
                            dir=dir+orb(isite,i)*orb(isite,kvirt)* &
                            orb(jsite,j)*orb(jsite,lvirt)*v(ijsite)
                        
                            exc=exc+orb(isite,i)*orb(isite,j)* &
                            orb(jsite,kvirt)*orb(jsite,lvirt)*v(ijsite)
                        
                        end do
                    
                    end do
                
                    xham=2.d0*dir-exc
                    iik=max(ik+1,jl+1)
		    jjl=min(ik+1,jl+1)
		    ijkl=(iik*(iik-1))/2+jjl
                    ham(ijkl)=xham
                    
                
                end do
            end do
        
        end do
    end do

    return
    end subroutine ham_sci

! ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

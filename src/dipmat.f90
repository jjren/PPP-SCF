! ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! This routine computes the dipole matrix with respect to the single-CI
! CSFs.

    subroutine dipmat(dmu_csf,dmu_scf,mocc_orb,mvirt_orb,morb, &
    mdimh,ndim)

! ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
implicit none


    integer,intent(in)::mocc_orb,mvirt_orb,morb,mdimh,ndim
    real(kind=8),intent(in)::dmu_scf((morb*(morb+1))/2,ndim)
    real(kind=8),intent(out)::dmu_csf((mdimh*(mdimh+1))/2,ndim)
    integer::iaxis,ij,iocc,ivirt,ik,ij1,i,k,kvirt,j,ji,lmin,&
          l,jl,ijkl,lvirt,ll,kk,kl
    real(kind=8)::xx2
    
    dmu_csf=0.d0

    xx2=dsqrt(2.d0)

! Only off-diagonal elements are nonzero

! First the matrix elements w.r.t. the SCF ground state

    do iaxis=1,ndim
        ij=1
        do iocc=1,mocc_orb
            do ivirt=mocc_orb+1,mocc_orb+mvirt_orb
            
                ik=(ivirt*(ivirt-1))/2+iocc
                ij=ij+1
                ij1=(ij*(ij-1))/2+1
                dmu_csf(ij1,iaxis)=xx2*dmu_scf(ik,iaxis)
            end do
        end do
    end do

! Now the other off-diagonal elements

    ik=1
    do i=1,mocc_orb
        do k=1,mvirt_orb
        
            ik=ik+1
            kvirt=k+mocc_orb
        
            do j=i,mocc_orb
            
                ji=(j*(j-1))/2+i
                lmin=k+1
                if(j /= i)lmin=1
            
                do l=lmin,mvirt_orb
                
                    jl=(j-1)*mvirt_orb+l+1
                    ijkl=(jl*(jl-1))/2+ik
                    lvirt=mocc_orb+l
                
                    if(i == j .AND. k /= l)then
                        kk=max(kvirt,lvirt)
                        ll=min(kvirt,lvirt)
                        kl=(kk*(kk-1))/2+ll
                    
                        do iaxis=1,ndim
                            dmu_csf(ijkl,iaxis)=dmu_scf(kl,iaxis)
                        end do
                    
                    elseif(i /= j .AND. k == l)then
                    
                        do iaxis=1,ndim
                            dmu_csf(ijkl,iaxis)=-dmu_scf(ji,iaxis)
                        end do
                    
                    end if
                
                end do
            end do
        
        end do
    end do

    return
    end subroutine dipmat

! ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

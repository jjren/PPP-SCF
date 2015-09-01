! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! This routine four-index transforms the 2-electron
! integrals in the site representation into those in the
! representation of orbitals stored in the array orb.

    subroutine fourind(v,r2int,buf,orb,msite,morb)

! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

IMPLICIT NONE
   
    integer,intent(in)::msite,morb
    real(kind=8),intent(in)::v((msite*(msite+1))/2),orb(msite,morb)
    real(kind=8),intent(out)::r2int((morb*(morb+1)*(morb*morb+morb+2))/8),&
                         buf((morb*(morb+1))/2)
    integer::nmorb,nijklo,ijkl,k,i,klorb,korb,lorb,ijorb,ijklor,&
                iorb,jorb,lomx
    real(kind=8)::efac,q

    nmorb=(morb*(morb+1))/2
    nijklo=(nmorb*(nmorb+1))/2

 
        r2int=0.d0



    ijkl=0
    do i=1,msite
      
            buf=0.d0
    
        do k=1,i
        
            ijkl=ijkl+1
        
            efac=1.d0
            if(i == k)efac=2.d0
        
            q=v(ijkl)
        
        ! transform the index k
        
            klorb=0
            do korb=1,morb
                do lorb=1,korb
                    klorb=klorb+1
                    buf(klorb)=buf(klorb)+(orb(k,korb)*orb(k,lorb)*q) &
                    /efac
                end do
            end do
        
        end do
    
    ! finally transform the index i and manipulate the transformed integrals
    
        ijorb=0
        ijklor=0
        do iorb=1,morb
            do jorb=1,iorb
                ijorb=ijorb+1
                klorb=0
                do korb=1,iorb
                    lomx=korb
                    if(iorb == korb)lomx=jorb
                    do lorb=1,lomx
                    
                        klorb=klorb+1
                        ijklor=ijklor+1
                    
                        r2int(ijklor)=r2int(ijklor) &
                        +orb(i,iorb)*orb(i,jorb)*buf(klorb) &
                        +orb(i,korb)*orb(i,lorb)*buf(ijorb)

                    end do
                end do
            end do
        end do
    
    end do


    if(ijklor /= nijklo)then
        write(*,*)'sr fourind: ijklor.ne.nijklo'
        write(*,*)'something wrong in the orb count'
        write(*,*)'ijklor,nijklo',ijklor,nijklo
        stop
    end if

    RETURN

! END! fourind

    end subroutine fourind
! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC







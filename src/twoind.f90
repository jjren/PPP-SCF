! ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! This routine performs the two-index transformation on the one-electron
! matrix elements.

    subroutine twoind(r1int,h,orb,buf,msite,morb)

! ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
IMPLICIT NONE

   
   integer,intent(in)::msite,morb
   real(kind=8),intent(in)::h((msite*(msite+1))/2),orb(msite,morb)
   real(kind=8),intent(out)::buf(morb),r1int((morb*(morb+1))/2)
   integer::nmorb,ijorb,iorb,i,j,ij,jorb
   real(kind=8)::efac

! Now perform the two-index transformation

    nmorb=(morb*(morb+1))/2

    
        r1int=0.d0
    

    ij=0
    do i=1,msite
    
        do iorb=1,morb
            buf(iorb)=0.d0
        end do
    
        do j=1,i
        
            ij=ij+1
            efac=1.d0
            if(i == j)efac=0.5d0
        
            do iorb=1,morb
                buf(iorb)=buf(iorb)+efac*h(ij)*orb(j,iorb)
            end do
        
        end do
    
        ijorb=0
        do iorb=1,morb
            do jorb=1,iorb
                ijorb=ijorb+1
                r1int(ijorb)=r1int(ijorb)+buf(iorb)*orb(i,jorb) &
                +buf(jorb)*orb(i,iorb)
            end do
        end do
    
    end do

    if(ijorb /= nmorb)then
        write(*,*)'ijorb,nmorb',ijorb,nmorb
        print*,'Error in Twoind','ijorb not equal to nmorb'
    end if

    return
    end subroutine twoind

! ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

! ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! This subroutine computes the Coulomb matrix elements Vij.

    subroutine vij_cal(vij,rr_site,msite,natmic,natmjc,itypec,&
             jtypec,iham,ihampar,ichngcl,r0new,&
                unew,dielnew,r0,u,vv,diel,ifcorc)

! ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
IMPLICIT NONE

!    implicit real*8(a-h,o-z)
!    implicit integer*4(i-n)
    
    integer::ok,ij,i,ii,j,jj,iii,jjj
    integer,intent(in)::msite,natmic,natmjc,itypec(*),jtypec(*),&
                   iham,ihampar,ichngcl,ifcorc
    real(kind=8),intent(out)::vij((msite*(msite+1))/2)
    real(kind=8),intent(in)::rr_site(3,msite),diel,r0new,&
                unew,dielnew,r0,u,vv

    real(kind=8)::rinv2,x1,y1,z1,x2,y2,z2,rr,rinv
    ! jjren define
    real(kind=8),allocatable :: hubbardU(:)
    real(kind=8) :: averageU

!
 vij=0.d0
! 

! Hubbard model

    if(iham == 3)then
    
        do i=1,msite
            ij=(i*(i-1))/2+i
            vij(ij)=u
        end do
    
    end if

! PPP Hamiltonian

    if(iham == 1)then
    
        ij=0
    
    ! Ohno parameterization
    
        if(ihampar == 1)then
        
		! jjren different site have different hubbardU
		open(unit=1001,file="hubbard.inp",status="old")
		allocate(hubbardU(msite))
		do i=1,msite,1
			read(1001,*) hubbardU(i)
			ii=(i*(i-1))/2+i
			vij(ii)=hubbardU(i)
		end do
		close(1001)
!
            do i=1,msite
            
                x1=rr_site(1,i)
                y1=rr_site(2,i)
                z1=rr_site(3,i)
            
                do j=1,i-1
                
                    ij=(i*(i-1))/2+j
                
                    x2=rr_site(1,j)
                    y2=rr_site(2,j)
                    z2=rr_site(3,j)
                
                    rr=(abs(x1-x2))**2 &
                    +(abs(y1-y2))**2 &
                    +(abs(z1-z2))**2
			! jjren
			averageU=(hubbardU(i)+hubbardU(j))/2.0D0
			vij(ij)=averageU/sqrt(1+averageU*averageU*rr/14.397D0/14.397D0)
                
                end do
            end do
		deallocate(hubbardU)
!        
!            if(ichngcl > 0)then
!                rinv2=1.d0/(r0new*r0new)
!            
!                do i=1,natmic
!                    ii=itypec(i)
!                    x1=rr_site(1,ii)
!                    y1=rr_site(2,ii)
!                    z1=rr_site(3,ii)
!                    do j=1,natmjc
!                        jj=jtypec(j)
!                        iii=max(ii,jj)
!                        jjj=min(ii,jj)
!                        ij=(iii*(iii-1))/2+jjj
!                        x2=rr_site(1,jj)
!                        y2=rr_site(2,jj)
!                        z2=rr_site(3,jj)
!                    
!                        rr=(abs(x1-x2))**2 &
!                        +(abs(y1-y2))**2 &
!                        +(abs(z1-z2))**2
!                    
!                        vij(ij)=unew/(dielnew*sqrt(1.d0+rr*rinv2))
!
!                    end do
!                end do
!            end if

  
        ! Mataga-Nishimoto parameterization
        
        elseif(ihampar == 2)then
        
            rinv=1.d0/r0
        
            do i=1,msite
            
                x1=rr_site(1,i)
                y1=rr_site(2,i)
                z1=rr_site(3,i)
            
                do j=1,i
                
                    ij=ij+1
                
                    x2=rr_site(1,j)
                    y2=rr_site(2,j)
                    z2=rr_site(3,j)
                
                    rr=dsqrt((x1-x2)*(x1-x2) &
                    +(y1-y2)*(y1-y2) &
                    +(z1-z2)*(z1-z2))*rinv
                
                    vij(ij)=u/(1.d0+rr)
                
                end do
            
            end do
        
        ! Exponential parameterization
        
        elseif(ihampar == 3)then
        
            rinv=1.d0/r0
        
            do i=1,msite
            
                x1=rr_site(1,i)
                y1=rr_site(2,i)
                z1=rr_site(3,i)
            
                do j=1,i
                
                    ij=ij+1
                
                    x2=rr_site(1,j)
                    y2=rr_site(2,j)
                    z2=rr_site(3,j)
                
                    rr=dsqrt((x1-x2)*(x1-x2) &
                    +(y1-y2)*(y1-y2) &
                    +(z1-z2)*(z1-z2))*rinv
                
                    vij(ij)=u*exp(-rr)
                
                end do
            
            end do
        
        end if
    
    end if


! The nearest-neighborextended Hubbard (U-V) model

    if(iham == 2)then
    
        do  i=1,msite
            ij=(i*(i-1))/2+i
            vij(ij)=u
            ij=ij-1
            if(i > 1)vij(ij)=vv
        end do
        ij=0
        do i=1,msite
            do j=1,i
                ij=ij+1
            end do
        end do
    
    end if

    return
    end subroutine vij_cal

! ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

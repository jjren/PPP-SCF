! ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    subroutine read_ei(t,nsite,nunit,natom)

! This subroutine reads the site energies, if any.
! cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
IMPLICIT NONE

!    implicit real*8(a-h,o-z)
!    implicit integer*4(i-n)

    integer::i,nitem,j,jj,iunit,ii,kk
    integer,intent(in)::nsite,nunit,natom
    real(kind=8)::xx
    real(kind=8),intent(out)::t((nsite*(nsite+1))/2)

    character(len=80)::card,line 


    read(*,*,err=1001)line
    read(*,*)card
    card=adjustl(card)
    if(card(1:4) == 'SITE')then
    
    
        read(*,*,err=1002)nitem
    
        write(*,*)'Number of Nonzero Site Energies/cell =',nitem
    
        do i=1,nitem
        
            read(*,*,err=1003)j,xx
            if(j <= nsite)then
                ii=(j*(j-1))/2+j
                t(ii)=xx
            
            ! If more than one unit cell is there, generate rest of the hoppings
            ! here
            
                do iunit=2,nunit
                
                    jj=j+(iunit-1)*natom
                    if(jj <= nsite)then
		      kk=(jj*(jj-1))/2+jj
		      t(kk)=xx
                    end if
                end do
            
            else
            
                write(*,*)'For item, site energy',i,xx
                write(*,*)'j',j
                print*,'Error:Read_ei','j.gt.nsite'
            
            end if

        end do
    
   
    end if

    goto 1004

    1001 print*,'Read_ei','Error reading the label'
    1002 print*,'Read_ei','Error reading no. of site energies to read'
    1003 write(*,*)'For i=',i
    print*,'Read_ei','Error read the site energy'

    1004 return
    end subroutine read_ei

! cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

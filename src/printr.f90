! cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! This routine truncates the position coordinates to places to right
! of the decimal. In addition, it also prints out the coordinates of
! all the sites to the specified logical unit

    subroutine printr(rr_site,nsite)

! cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

implicit none

!    implicit real*8(a-h,o-z)
!    implicit integer*4(i-n)

!    dimension rr_site(3,nsite)

    integer,intent(in)::nsite
    real(kind=8),intent(inout)::rr_site(3,nsite)
    integer::i
    real(kind=8)::x,y,z,xfac,yfac,zfac

    print*,' '
    write(*,101)
    101 format('Site #',5x,'X-Coord',5x,'Y-Coord',4x,'Z-coord')

    write(*,*)' '
    do i=1,nsite
    
        x=rr_site(1,i)
        y=rr_site(2,i)
        z=rr_site(3,i)
        xfac=1.d+7
        yfac=1.d+7
        zfac=1.d+7
        if(dabs(x) > 200.d0)xfac=1.d0
        if(dabs(y) > 200.d0)yfac=1.d0
        if(dabs(z) > 200.d0)zfac=1.d0
        x=nint(rr_site(1,i)*xfac)
        y=nint(rr_site(2,i)*yfac)
        z=nint(rr_site(3,i)*zfac)
        rr_site(1,i)=x/xfac
        rr_site(2,i)=y/yfac
        rr_site(3,i)=z/zfac

        write(*,102)i,rr_site(1,i),rr_site(2,i),rr_site(3,i)
    
    end do

    102 format(1x,i4.1,2x,f11.5,1x,f11.5,1x,f11.5)

    write(*,*)' '
    
!      do i=1,nsite
!        print*,i,rr_site(1,i),rr_site(2,i),rr_site(3,i)
!      end do

!creating xcrysden file for lattice coordinates

    open(10,file='atomic_coord.xsf',status='unknown',action='write')
      write(10,*)'ATOM'
      do i=1,nsite
       write(10,103)6,rr_site(1,i),rr_site(2,i),rr_site(3,i)
      end do
      
         103 format(1x,i4.1,2x,f11.5,1x,f11.5,1x,f11.5)

    close(10)
    
    return
    end subroutine printr

! cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

! cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    subroutine stline(r_atom,blngth,rorg,nrotat,iaxis,angle)

! This routine generates a bond (2 carbon atoms) with the first of
! the bond at the specified location
! cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!    implicit real*8(a-h,o-z)
!    implicit integer*4(i-n)

!    dimension r_atom(3,2),rorg(3),iaxis(3),angle(3)

! Beginning of the bond is taken to be the origin and initially the
! bond is assumed to be along the x-axis.
! Rotations, if desired, are performed keeping
! the origin fixed. Finally, the center is translated to the location
! specified in array rorg.

IMPLICIT NONE
    integer,intent(in)::iaxis(3),nrotat
    real(kind=8),intent(in)::rorg(3)
    integer::i,irotat,jaxis,j,iatom
    real(kind=8)::pi,r_atom(3,2),blngth,&
                  ang,tx,ty,tz,x,y,z,angle(3)

! Value of pi

    pi=4.d0*datan(1.d0)

    do i=1,3
        angle(i)=(angle(i)*pi)/180.d0
    end do

! Start assigning coordinates

! Atom 1

    r_atom(1,1)=0.d0
    r_atom(2,1)=0.d0
    r_atom(3,1)=0.d0

! Atom 2

    r_atom(1,2)=blngth
    r_atom(2,2)=0.d0
    r_atom(3,2)=0.d0


! Rotate the bond, if needed
! For the rotation, +ve angles are considered to be
! anti-clockwise while the -ve angles are treated clockwise

    if(nrotat == 0)goto 4

    do irotat=1,nrotat
    
        jaxis=iaxis(irotat)
        ang=angle(irotat)
    ! write(iwrite,*)'jaxis,ang',jaxis,ang
    
        do iatom=1,2
        
            tx=r_atom(1,iatom)
            ty=r_atom(2,iatom)
            tz=r_atom(3,iatom)
        

            if(jaxis == 1)then
                x=tx
                y=ty*cos(ang)-tz*sin(ang)
                z=ty*sin(ang)+tz*cos(ang)
            elseif(jaxis == 2)then
                x=tz*sin(ang)+tx*cos(ang)
                y=ty
                z=tz*cos(ang)-tx*sin(ang)
            elseif(jaxis == 3)then
                x=tx*cos(ang)-ty*sin(ang)
                y=tx*sin(ang)+ty*cos(ang)
                z=tz
            end if
        
            r_atom(1,iatom)=x
            r_atom(2,iatom)=y
            r_atom(3,iatom)=z
        
        end do
    
    end do

    4 continue

! Finally the translation operation

    do iatom=1,2
    
        do j=1,3
            r_atom(j,iatom)=r_atom(j,iatom)+rorg(j)
        end do
    
    end do

    return
    end subroutine stline

! ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

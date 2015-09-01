! cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    subroutine benperp(r_atom,blngth,rorg,iplane,nrotat,iaxis, &
    angle)

! This routine generates the coordianates of six carbon atoms
! forming the benzene backbone. The orientation is perpendicular to the
! conventional orientation.
! cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!    implicit real*8(a-h,o-z)
!    implicit integer*4(i-n)

!    dimension r_atom(3,6),rorg(3),iaxis(3),angle(3)

! Note that the origin of the benzene skeleton is considered to
! center of hexagon. Rotations, if desired, are performed keeping
! the origin fixed. Finally, the center is translated to the location
! specified in array rorg.

IMPLICIT NONE
 
   integer,intent(in)::iaxis(3),iplane,nrotat
   real(kind=8),intent(in)::rorg(3),blngth
   integer::i,i1,i2,i3,irotat,jaxis,iatom,j
   real(kind=8)::pi,thirty,tx,ty,tz,x,y,z,r_atom(3,6),ang,angle(3)
 


! Value of pi

    pi=4.d0*datan(1.d0)

! All the angles are converted to radians

    thirty=pi/6.d0

    do  i=1,3
        angle(i)=(angle(i)*pi)/180.d0
    end do

! See in which plane is the molecule located
! iplane =  1 -> x-y plane (atoms 1 and 6 displaced along x-axis)
! =  2 -> y-z plane (atoms 1 and 6 displaced along y-axis)
! =  3 -> z-x plane (atoms 1 and 6 displaced along z-axis)

    if(iplane == 1)then
    
        i1=1
        i2=2
        i3=3
    
    elseif(iplane == 2)then
    
        i1=2
        i2=3
        i3=1
    
    elseif(iplane == 3)then
    
        i1=3
        i2=1
        i3=2
    
    else
    
        print*,'iplane=',iplane
        print*,'Benzene','Wrong value of iplane'
    
    end if

! Start assigning coordinates

! Atom 1

    r_atom(i1,1)=-blngth*dcos(thirty)
    r_atom(i2,1)=blngth*dsin(thirty)
    r_atom(i3,1)=0.d0

! Atom 2

    r_atom(i1,2)=-blngth*dcos(thirty)
    r_atom(i2,2)=-blngth*dsin(thirty)
    r_atom(i3,2)=0.d0

! Atom 3

    r_atom(i1,3)=0.d0
    r_atom(i2,3)=blngth
    r_atom(i3,3)=0.d0


! Atom 4

    r_atom(i1,4)=0.d0
    r_atom(i2,4)=-blngth
    r_atom(i3,4)=0.d0

! Atom 5

    r_atom(i1,5)=blngth*dcos(thirty)
    r_atom(i2,5)=blngth*dsin(thirty)
    r_atom(i3,5)=0.d0

! Atom 6

    r_atom(i1,6)=blngth*dcos(thirty)
    r_atom(i2,6)=-blngth*dsin(thirty)
    r_atom(i3,6)=0.d0

! Rotate the molecule, if needed
! For the rotation, +ve angles are considered to be
! anti-clockwise while the -ve angles are treated clockwise

    if(nrotat == 0)goto 4

    do  irotat=1,nrotat
    
        jaxis=iaxis(irotat)
        ang=angle(irotat)
    
        do  iatom=1,6
        
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

    do  iatom=1,6
    
        do  j=1,3
            r_atom(j,iatom)=r_atom(j,iatom)+rorg(j)
        end do
    
    end do

    return
    end subroutine benperp

! ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

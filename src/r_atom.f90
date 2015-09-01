! ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! This subroutine reads the coordinates of the atoms of a unit cell
! with respect to its (user specified) origin. It can also generate
! coordinates for the phenyl group (benzene) and the double bond.

    subroutine r_atom(rr_atom,angle,rorg,iaxis,matom)

! ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!    implicit real*8(a-h,o-z)
!    implicit integer*4(i-n)

!    character*(*) line,card

!    dimension rr_atom(3,matom)
!    dimension rorg(3),angle(3),iaxis(3)
IMPLICIT NONE

    integer,intent(in)::matom
    integer,intent(out)::iaxis(3)
    integer::iatom,nrotat,iplane,irotat,k,i,icoord=0
    real(kind=8),intent(inout)::rr_atom(3,matom)
    real(kind=8),intent(out)::rorg(3),angle(3)
    real(kind=8)::blngth,single,double
    character(80)::line,card
  
     
    iatom=0

    read(*,*)line

    1 continue

    read*,card
    card=adjustl(card)

    if(card(1:4) == 'BENP')then
    
        if(card(6:7) == 'XY')then
            iplane=1
        elseif(card(6:7) == 'YZ')then
            iplane=2
        elseif(card(6:7) == 'ZX')then
            iplane=3
        else
            print*,'Specified benzene plane',card(5:6)
            print*,'is wrong'
            stop
        end if
    
        read(*,*,err=1002)(rorg(i),i=1,3)
        read(*,*,err=1003)blngth,nrotat
    
        if(nrotat > 0)then
        
            if(nrotat > 3)then
	        print*,'R_atom','For benzene nrotat cant be more than 3'  
		stop
            endif       
            do  irotat=1,nrotat
                read(*,*,err=1004)iaxis(irotat),angle(irotat)
             end do
        
        end if
    
        call benperp(rr_atom(1,iatom+1),blngth,rorg,iplane,nrotat,iaxis, &
        angle)
        iatom=iatom+6
        goto 1
    
    elseif(card(1:3) == 'BEN')then
    
        if(card(5:6) == 'XY')then
            iplane=1
        elseif(card(5:6) == 'YZ')then
            iplane=2
        elseif(card(5:6) == 'ZX')then
            iplane=3
        else
            print*,'Specified benzene plane',card(5:6)
            print*,'is wrong'
            stop
        end if
    
        read(*,*,err=1002)(rorg(i),i=1,3)
        read(*,*,err=1003)blngth,nrotat
    
        if(nrotat > 0)then
        
            if(nrotat > 3)then
	      print*,'R_atom','For benzene nrotat cant be more than 3'
              stop
	    endif
	    
            do  irotat=1,nrotat
                read(*,*,err=1004)iaxis(irotat),angle(irotat)
            end do
        
        end if
    
        call benzen(rr_atom(1,iatom+1),blngth,rorg,iplane,nrotat,iaxis, &
        angle)
        iatom=iatom+6
        goto 1
    
    elseif(card(1:4) == 'BOND')then
    
        read(*,*,err=1005)(rorg(i),i=1,3)
        read(*,*,err=1006)blngth,nrotat
    
        if(nrotat > 0)then
        
            if(nrotat > 3)then
	      print*,'R_atom','For a bond nrotat cant be more than 3'
              stop
	    endif
	    
            do irotat=1,nrotat
                read(*,*,err=1007)iaxis(irotat),angle(irotat)
            end do
        
        end if
    
        call bond(rr_atom(1,iatom+1),blngth,rorg,nrotat,iaxis,angle)
    
        iatom=iatom+2
        goto 1
    
    elseif(card(1:4) == 'LINE')then
    
        read(*,*,err=1005)(rorg(i),i=1,3)
        read(*,*,err=1006)blngth,nrotat
    
        if(nrotat > 0)then
        
            if(nrotat > 3)then
	        print*,'R_atom','For a bond nrotat cant be more than 3'
                stop
            endif
	    
            do  irotat=1,nrotat
                read(*,*,err=1007)iaxis(irotat),angle(irotat)
            end do
        
        end if
    
        call stline(rr_atom(1,iatom+1),blngth,rorg,nrotat,iaxis,angle)
    
        iatom=iatom+2
        goto 1
    
    elseif(card(1:4) == 'ATOM')then
    
        read(*,*,err=1008)(rr_atom(i,iatom+1),i=1,3)
        iatom=iatom+1
        goto 1
    
    elseif(card(1:3) == 'C60')then
        read*,single,double
        read*,(rorg(k),k=1,3)
        call c60_gen(single,double,rr_atom(1,iatom+1),rorg)
        iatom=iatom+60
        goto 1
	
    elseif(card(1:5) == 'COORD')then
        icoord=1
    	do iatom=1,matom
            read(*,*,err=1008)(rr_atom(i,iatom),i=1,3)
        end do
!	print*,'hello'
!	print*,rr_atom
	goto 1
	
    elseif(card(1:4) == 'ENDA')then
      if(icoord==0)then
        if(iatom /= matom)then
            print*,'iatom,matom',iatom,matom
            print*,'R_atom','Wrong atom count with symbolic atom input'
        end if
      end if
        return
    
    endif


    1001 print*,'R_atom','Error reading symbolic atom card'
    stop
    1002 print*,'R_atom', &
    'Error reading the origin coordinates of the benzene unit'
    stop
    1003 print*,'R_atom', &
    'Error reading bond length/rotation info of the benzene unit'
    stop
    1004 print*,'R_atom', &
    'Error reading rotation angles for the benzene unit'
    stop
    1005 print*,'R_atom', &
    'Error reading the coordinates of the bond center'
    stop
    1006 print*,'R_atom', &
    'Error reading rotation angles for the bond'
    stop
    1007 print*,'R_atom', &
    'Error reading rotation angles for the bond'
    stop
    1008 print*,'R_atom','Error reading atomic coordinates'
    stop

    return
    end subroutine r_atom

! ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

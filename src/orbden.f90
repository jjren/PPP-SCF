! cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! This routine computes the charge density on specified sites and
! also on the remaining sites

    subroutine orbden(orb,iosite,nosite,norb,nsite,natom,nunit,iflag)

! cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
implicit none



    integer,intent(in)::iflag,natom,nunit,nsite,nosite,norb,iosite(nosite)
    real(kind=8),intent(in)::orb(nsite,*)
    integer::incr,nxsite,iunit,ix,ii,iorb,isite,ok
    real(kind=8)::xden

    integer,allocatable::ibuff(:)
    
        allocate(ibuff(nsite),stat=ok)
	if(ok/=0)then
          print*,'NOT ENOUGH MEMORY IS AVAILABLE FOR ALLOCATING ibuff'
	  stop
        end if 



    if(iflag > 0)then
        if(nunit > 1)then
            incr=0
            nxsite=0
            do iunit=1,nunit
                do ix=1,nosite
                    ii=iosite(ix)+incr
                    if(ii <= nsite)then
                        nxsite=nxsite+1
                        ibuff(nxsite)=ii
                    end if
                end do
            
                incr=incr+natom
            
            end do
        
        else
        
            print*,'Erron in Orbden','iflag.gt.0.and.nunit<=1 not allowed'
        
        end if
    
    else
        nxsite=nosite
        do ix=1,nosite
            ibuff(ix)=iosite(ix)
        end do
    
    end if

    write(*,*)'Sites identified in orbital density analysis:'
    write(*,*)(ibuff(ix),ix=1,nxsite)

    write(*,*)' '
    write(*,*) &
    'Orbital #    Charge Den (Chosen)   Charge Den (Remaining)'
    write(*,*)' '

    do iorb=1,norb
        xden=0.d0
        do ix=1,nxsite
            isite=ibuff(ix)
            xden=xden+orb(isite,iorb)**2
        end do
        write(*,101)iorb,xden,1.d0-xden
    end do

    101 format(1x,i4.1,14x,f8.4,14x,f8.4)
    deallocate(ibuff)
    return
    end subroutine orbden

! ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

! Revision: 1.1 $
! cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    subroutine writorb(iorbu,orbfile,orb,norb,ncom)

! This subroutine writes the starting orbitals to a file whose name
! is stored in the variable orbfile.
! cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    implicit real*8(a-h,o-z)

    dimension orb(ncom,norb)
    character*(*) orbfile

    open(unit=iorbu,file=orbfile,access='sequential', &
    form='formatted',status='unknown',err=111)
    goto 112

    111 print*,'sr orbout: error opening output orb file ',orbfile
    stop

    112 continue

    rewind(iorbu)

    istage=1
    write(iorbu,*,err=102)norb,ncom


    istage=2

    do 1 iorb=1,norb
        write(iorbu,101,err=103)iorb
        write(iorbu,*,err=104)(orb(i,iorb),i=1,ncom)
    1 end do

    101 format(1x,'Orbital #:',1x,i5.1)

    goto 105

    102 print*,'error writing to the orbital file,istage',istage
    stop

    103 print*,'error writing to the header of orbital #: ',iorb
    stop

    104 print*,'error writing orbital #: ',iorb
    stop

    105 continue


    close(iorbu)

    return
    end subroutine writorb
! ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

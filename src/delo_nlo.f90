! cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    subroutine delo_nlo(orb,eval,nbas,norb,nocc,ndelo1,ndelo2,ndelv1, &
    ndelv2)

! This subroutine deletes occupied orbitals ranging from ndelo1 to
! ndelo2 and virtual orbitals ndelv1,ndelv2 from the orbital array orb.
! cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


IMPLICIT NONE

    integer,intent(in)::nbas,ndelo1,ndelo2,ndelv1,ndelv2
    integer,intent(inout)::norb,nocc
    real(kind=8),intent(inout)::orb(nbas,norb),eval(norb)
    integer::norb_o,ndelo,ndelv,nocc_o,nvrt_o,nvrt,&
          istrt,iorb,ibas



    norb_o=norb
    ndelo=ndelo2-ndelo1+1
    ndelv=ndelv2-ndelv1+1
    norb=norb_o-ndelo-ndelv

    nocc_o=nocc
    nvrt_o=norb_o-nocc_o
    nvrt=nvrt_o-ndelv
    nocc=nocc_o-ndelo
    norb=nvrt+nocc
    if(nocc <= 0)then
        write(*,*)'sr delo_nlo: new norb.le.0'
        write(*,*)'norb=',norb
        stop
    end if

! First move the occupied orbitals

    istrt=ndelo1
    do iorb=ndelo2+1,nocc_o
        eval(istrt)=eval(iorb)
        do ibas=1,nbas
            orb(ibas,istrt)=orb(ibas,iorb)
        end do
        istrt=istrt+1
    end do

! Now move the virtual orbitals

    istrt=nocc+1
    do iorb=nocc_o+1,ndelv1-1
        eval(istrt)=eval(iorb)
        do ibas=1,nbas
            orb(ibas,istrt)=orb(ibas,iorb)
        end do
        istrt=istrt+1
    end do

    do iorb=ndelv2+1,norb_o
        eval(istrt)=eval(iorb)
        do ibas=1,nbas
            orb(ibas,istrt)=orb(ibas,iorb)
        end do
        istrt=istrt+1
    end do

    istrt=istrt-1
    if(istrt /= norb)then
        write(*,*)'sr delo_nlo: istrt.ne.norb'
        write(*,*)'istrt,norb',istrt,norb
        stop
    end if

    return
    end subroutine delo_nlo

! cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


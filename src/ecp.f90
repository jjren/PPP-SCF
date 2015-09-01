! cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    subroutine ecp(ecore,ep,orb,v,h,norb,nsite)

! This program creates the effective one-electron potential
! corresponding to the frozen orbitals, and puts them in the array
! ep. The energy contribution of the frozen orbitals is put in
! the variable ecore.
! cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

IMPLICIT NONE

    integer,intent(in)::nsite,norb
    real(kind=8),intent(in)::v((nsite*(nsite+1))/2),h((nsite*(nsite+1))/2)
    real(kind=8),intent(out)::ep((nsite*(nsite+1))/2),&
               orb(nsite,norb),ecore	
    real(kind=8),allocatable::den(:)
    integer::nsite0,ij,i,ii,j,jj,ok
 !
    allocate(den((nsite*(nsite+1))/2),stat=ok)
      if(ok/=0)then
       print*,'ECP: Allocation for den fails'
       stop
      end if
  !



    
        den=0.d0
        ep=0.d0
    


! The following routine reads the orbitals one-by-one from the file
! and accumulates their contribution to density matrix on the fly.
! Note that whole orbital set is never in the memory.

    call density(den,orb,norb,nsite)

! Contruct the ecp

    do i=1,nsite
    
        ii=(i*(i-1))/2+i
    
    ! i.ne.j part
    
        do j=1,i-1
        
            jj=(j*(j-1))/2+j
            ij=(i*(i-1))/2+j
        
        ! Coulomb contribution
        
            ep(ii)=ep(ii)+v(ij)*den(jj)
            ep(jj)=ep(jj)+v(ij)*den(ii)
        
        ! Exchange contribution
        
            ep(ij)=ep(ij)-0.25d0*v(ij)*den(ij)
        
        end do
    
    ! i.eq.j part
    
    ! Coulomb contribution
    
        ep(ii)=ep(ii)+v(ii)*den(ii)
    
    ! Exchange contribution
    
        ep(ii)=ep(ii)-0.5d0*v(ii)*den(ii)
    
    end do

! Now get the core energy.

    ecore=0.d0
    ij=0
    do i=1,nsite
        do j=1,i
        
            ij=ij+1
            ecore=ecore+den(ij)*h(ij)+0.5d0*den(ij)*ep(ij)
        
        end do
    
    end do

    return
    end subroutine ecp

! ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

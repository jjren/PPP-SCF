! cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! This is the master routine meant for generating the one-electron
! (tij) and two-electron (vij) matrix elements.

    subroutine matel(iham,ibands,enuc,tij,hij,ihij,vij,&
    rr_site,nunit,natom,nsite,icorfnc,xcor,natmi,natmj,natmic,natmjc, &
    itype,jtype,iefield,efield,ihampar,ichngcl,itypec,&
             jtypec,r0new,unew,dielnew,r0,u,vv,diel,ifcorc)

! cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
IMPLICIT NONE

!    implicit real*8(a-h,o-z)
!    implicit integer*4(i-n)
     
      integer::nij,ij,ihopgen,i,ii,j,jj,jjj,iii
      integer,intent(in)::ibands,iham,natom,nunit,icorfnc,iefield,nsite,&
                 natmic,natmjc,natmi,natmj,&
		 itype(*),jtype(*),itypec(*),jtypec(*),&
		 ihampar,ichngcl,ifcorc
      integer,intent(out)::ihij((nsite*(nsite+1))/2)
		
      real(kind=8),intent(in)::rr_site(3,nsite),diel,r0new,&
                unew,dielnew,r0,u,vv,efield(3),xcor
      real(kind=8),intent(out)::tij((nsite*(nsite+1))/2),hij((nsite*(nsite+1))/2),&
          vij((nsite*(nsite+1))/2)
      real(kind=8),intent(inout)::enuc
	! jjren defined every site core positive charge
	real(kind=8),allocatable :: nuclQ(:)
          
 
    nij=(nsite*(nsite+1))/2

! Zero out the arrays

    do ij=1,nij
        tij(ij)=0.d0
        hij(ij)=0.d0
        vij(ij)=0.d0
        ihij(ij)=0
    end do

! Read the Tij matrix elements

    if(ibands > 0 .AND. iham == 4)then
    
        call tij_readb(tij,hij,ihij,natom)
        call read_ei(tij,nsite,nunit,natom)

        return
    
    else
 
        ihopgen=0

        call tij_read(tij,nsite,nunit,natom, &
        ihopgen)
        if(ihopgen>0)call tij_gen(tij,nsite,rr_site)
        call read_ei(tij,nsite,nunit,natom)

    end if
! c
! Compute the Vij matrix elements
! Skip if Hueckel Hamiltonian

    if(iham == 4)goto 101

      call vij_cal(vij,rr_site,nsite,natmic,natmjc,itypec,&
             jtypec,iham,ihampar,ichngcl,r0new,&
                unew,dielnew,r0,u,vv,diel,ifcorc)

!     jjren suited for different core charge condition
	open(unit=1002,file="nuclQ.inp",status="old")
	allocate(nuclQ(nsite))
	do i=1,nsite,1
		read(1002,*) nuclQ(i)
	end do
	close(1002)
! Construct the full one-electron matrix Hij.

! First the on-site energies

    do i=1,nsite
    
        ii=(i*(i-1))/2+i
    
        do j=1,i-1
        
            ij=(i*(i-1))/2+j
            hij(ii)=hij(ii)-vij(ij)*nuclQ(j)
        
        end do
    
        do j=i+1,nsite
        
            ij=(j*(j-1))/2+i
            hij(ii)=hij(ii)-vij(ij)*nuclQ(j)
        
        end do
    
    end do

    101 continue


! Add the hopping and the onsite energies, if any

        hij=hij+tij
  
! Add the contribution due to the external electric field, if needed
! nucl energy

    if(iefield > 0)then
        do i=1,nsite
            ii=(i*(i-1))/2+i
            do j=1,3
                hij(ii)=hij(ii)+efield(j)*rr_site(j,i)
            end do
        end do
    end if

! call plblk(hij(1),nsite,0,'Site',6)

! Calculate the nucleus-nucleus interaction

    enuc=0.d0
    do i=1,nsite
        do j=1,i-1
            ij=(i*(i-1))/2+j
            enuc=enuc+vij(ij)*nuclQ(i)*nuclQ(j)
        end do
    end do

! Add the contribution of the e-field to enuc


    if(iefield > 0)then
        do i=1,nsite
            do j=1,3
                enuc=enuc-efield(j)*rr_site(j,i)*nuclQ(i)
            end do
        end do
    end if

! Add the contribution corresponding to the correlation functions,
! if needed 
! Should be used only for RHF

    if(icorfnc > 0)then
      print*,'CORRELATION FUNCTION CALCULATIONS WILL BE DONE'
      print*,' '
      print*,'\sum_i,j <n_i*n_j> can be computed using the data energies'
      write(*,208)xcor
      print*,'i-type atoms are: ',(itype(i),i=1,natmi)
      print*,'j-type atoms are: ',(jtype(i),i=1,natmj)
      print*,' '
   
   208 format(1x,'Small amount by which vij is incremented:',1x,f8.4)
      
        do i=1,natmi
            ii=itype(i)
            do j=1,natmj
                jj=jtype(j)
                iii=max(ii,jj)
                jjj=min(ii,jj)
                ij=(iii*(iii-1))/2+jjj
                vij(ij)=vij(ij)+xcor
            end do
        end do
    
        do i=1,natmi
            ii=itype(i)
            iii=(ii*(ii-1))/2+ii
            hij(iii)=hij(iii)-xcor*natmj
        end do
    
        do j=1,natmj
            ii=jtype(j)
            iii=(ii*(ii-1))/2+ii
            hij(iii)=hij(iii)-xcor*natmi
        end do
    
        enuc=enuc+dble(natmi)*dble(natmj)*xcor
    
    end if
!    print*,'enuc=',enuc
    deallocate(nuclQ)
    return
    end subroutine matel

! cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

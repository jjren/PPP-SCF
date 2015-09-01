! ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    subroutine ci_drv(eval,orb,r_site,h,v,msite,norbt,&
    norbf,mcell,nfreez,ifreez,idipint,ndelete,idelete,enuc,idir,ndir,&
    ifrez)

! This is the driver routine which controls the preparations for the
! subsequent CI calculations. This routine calls routines to:

! (1) Read the occupied orbitals
! (2) Read the virtual orbitals
! (3) Generate the orbitals of the whole system
! by translation operation if the
! orbitals read were of one unit cell
! and the system on hand has multiple
! unit cells
! (4) Generate the one- and two-electron integrals in the
! basis of these orbitals (exciton basis) by transforming
! the corresponding integrals in the site representaions.
! They are written to different files so that they
! can be read and used by a separate CI program.

! cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
IMPLICIT NONE

    integer,intent(in)::msite,norbt,norbf,mcell,ifrez,idipint,&
                  nfreez,ifreez(*),ndir,idir(*),ndelete,idelete(*)
    real(kind=8),intent(in)::v((msite*(msite+1))/2), h((msite*(msite+1))/2),&
             r_site(3,msite),orb(msite,norbt),eval(norbt),enuc  
    integer::n12,iorbu,istrt,nsite0,norbt0,nbuf,ok,nij,nn12,norb_a,nsite,&
          nm12,nm18,morb,iorb,j,n1234,korb,icheck,jorb,nx12
   
   real(kind=8),allocatable::ep(:),r1int(:),r2int(:),forb(:,:),orb_a(:,:),buf(:),&
                    eval_a(:),dipoles(:,:)
    real(kind=8)::ecore

    
    nsite=msite
    morb=norbt 

      if(nfreez>0)then
        allocate(ep((nsite*(nsite+1))/2),stat=ok)
        if(ok /= 0)then
           WRITE(*,*)'CI_DRV:Allocation fails for ep'
           stop
        end if
       end if
 
	
    if(nfreez>0)then
      allocate(forb(nsite,nfreez),stat=ok)
      !
      if(ok/=0)then
        print*,'CI_drv: Allocation for forb fails'
        stop
      end if
        !
      do iorb=1,nfreez
       do j=1,nsite
         forb(j,iorb)=orb(j,ifreez(iorb))
       end do
      end do
  !
    end if
!	
   if(nfreez>0)then
      norb_a=morb-nfreez
   else
      norb_a=morb
   end if
   
   if(ndelete>0)norb_a=norb_a-ndelete
  !
   write(*,3)norb_a
  !
 
    nij=(nsite*(nsite+1))/2
    nn12=(norb_a*(norb_a+1))/2
    n1234=(nn12*(nn12+1))/2
    nbuf=max(2*nij+nsite,nn12,norb_a*nsite)+500
    
    allocate(buf(nbuf),stat=ok)
      if(ok/=0)then
        print*,'NOT ENOUGH MEMORY IS AVAILABLE FOR ALLOCATING buf'
        stop
      end if
 !
  nm12=(norb_a*(norb_a+1))/2
  nm18=(nm12*(nm12+1))/2
  !   
   
 3 format(1x,'Total number of orbitals in the CI calculations:',1x,i3.1)
   if(nfreez>0)then 
      print*,' '
      write(*,1)nfreez
 1    format(1x,'Total Number of Frozen Core Orbitals=',i3.1)
   end if
   if(ndelete>0)then 
      print*,' '
      write(*,2)ndelete
 2    format(1x,'Total Number of Deleted Virtual Orbitals=',i3.1)
   end if	
  !
    allocate(orb_a(nsite,norb_a),eval_a(norb_a),stat=ok)
      if(ok/=0)then
       print*,'CI_DRV: Allocation for orb_a and eval_a fails'
       stop
      end if
  !	   
      allocate(r1int(nm12),r2int(nm18),stat=ok)
        if(ok /= 0)then
           WRITE(*,*)'CI_DRV:Allocation fails for r1int and r2int'
           stop
        end if	
	
      
    n12=(msite*(msite+1))/2
    iorbu=20
    ecore=enuc

    

! If some orbitals need to be frozen, generate the effective potential
! corresponding to them
! w
  !
  if(nfreez>0.or.ndelete>0)then
  
      korb=0
      do iorb=1,norbt
        icheck=0
        !
        do jorb=1,nfreez
           if(iorb==ifreez(jorb))icheck=icheck+1
        end do
        !
        do jorb=1,ndelete
           if(iorb==idelete(jorb))icheck=icheck+1
        end do
        !
        if(icheck==0)then
           korb=korb+1
           do j=1,nsite
              orb_a(j,korb)=orb(j,iorb)
           end do
	   eval_a(korb)=eval(iorb)
        end if
      end do
   else
      orb_a=orb
      eval_a=eval
  end if


      if(nfreez>0)then
      

        call ecp(ecore,ep,forb,v,h,norbf,nsite)
		
        ecore=ecore+enuc
        print*,'Core + Nuc-Nuc Energy=',ecore
!    end if

! Transform the one-electron integrals


!    if(i > 0)then
    
        ep=ep+h

        call twoind(r1int,ep,orb_a,buf,nsite,norb_a)
    
    else
    
        call twoind(r1int,h,orb_a,buf,nsite,norb_a)
    
    end if

! Write out the one-electron integrals

    call write_1(iorbu,'ONEINT001.DAT',norb_a,ecore, &
    r1int,1.d-8,nsite)

! Transform the two-electron integrals


    call fourind(v,r2int,buf,orb_a,nsite,norb_a)



! Spit out the two-electron integrals

    call write_2(iorbu,'TWOINT001.DAT',norb_a,r2int, &
    1.d-8)

! Finally, if the dipole matrix elements are needed, compute them here.

    if(idipint > 0)then
    
      !	   
       allocate(dipoles(nm12,ndir),stat=ok)
        if(ok /= 0)then
           WRITE(*,*)'CI_DRV:Allocation fails for dipoles'
           stop
        end if
	
        nx12=(norb_a*norb_a+norb_a)/2

        call dipint(dipoles,r_site,orb_a,idir,nsite,norb_a,ndir)
        call dipout(iorbu,'DIPINT001.DAT',norb_a,dipoles, &
        idir,ndir,1.d-8)
    
    end if

    return
    end subroutine ci_drv

! ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

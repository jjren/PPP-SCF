! ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! This routine solves the SCF eigenvalue problem by diagonalizing the HF
! matrix using subroutines like DSPEV and DSPEVX from the LAPACK/BLAS library.

    subroutine scf_rhf(t,h,v,msite,idamp,xdamp,maxitr,norb,orb,&
    iprorb,conv,eval,enuc,escf,iham)

! ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
IMPLICIT NONE

   
    integer::ok,iorbu,isite,i,j,itr,iflag1,iflag2,mevlfnd,&
       info
    integer,intent(in)::msite,idamp,maxitr,norb,iprorb,iham
    real(kind=8),intent(in)::t((msite*(msite+1))/2),enuc,&
               h((msite*(msite+1))/2),v((msite*(msite+1))/2),xdamp,conv
    real(kind=8),intent(out)::orb(msite,msite),eval(msite),escf
    real(kind=8),allocatable::fock(:),fock_o(:),&
                 den(:),eval_o(:)      
   
    real(kind=8)::dif1,dif2,xdiff,escf1,abstol,VL,VU,dlamch
             
    real(kind=8),allocatable::work(:)
    integer,allocatable::ifail(:),iwork(:)  
    
    

    
      abstol=2.d0*dlamch('S')
      allocate(work(8*msite),stat=ok)
       if(ok/=0)then
        print*,"Allocation for work fails"
        stop
       end if
       work=0.d0
  !
       allocate(iwork(5*msite),stat=ok)
         if(ok/=0)then
          print*,"Allocation of iwork fails"
         stop
       end if
       iwork=0
  !
  !
       allocate(ifail(msite),stat=ok)
        if(ok/=0)then
          print*,"Allocation of ifail ceases"
        stop
       end if
       ifail=0

    allocate (fock((msite*(msite+1))/2),stat=ok)
      if(ok/=0)then
        print*,'NOT ENOUGH MEMORY IS AVAILABLE FOR ALLOCATING fock'
        stop
      end if
    allocate (den((msite*(msite+1))/2),eval_o(norb),stat=ok)
      if(ok/=0)then
        print*,'NOT ENOUGH MEMORY IS AVAILABLE FOR ALLOCATING den,g, or eval'
        stop
      end if

      if(idamp>0)then
         allocate (fock_o((msite*(msite+1))/2),stat=ok)
         if(ok/=0)then
           print*,'NOT ENOUGH MEMORY IS AVAILABLE FOR ALLOCATING fock_o'
           stop
         end if
      end if

!      call plblk(v,msite,0,'site',6)

    write(*,*)' '
    write(*,255)
    255 format(1x,'Itr(i)',5x,'E(scf)',8x,'Ediff1(i-1)',7x,'Ediff2(i)', &
    7x,'Evldiff(i)')
    write(*,*)' '


! Starting occupied vectors are zero

           orb=0.d0
           fock=0.d0
	   eval=0.d0
	   eval_o=0.d0

! dif1=E(i)-E(i-1), dif2=E(i-1)-E(i-2), where i=present iteration
! E is the total energy

    dif1=1000.d0
    dif2=1000.d0
    xdiff=1000.d0
    escf=20.d0
    itr=0

    3 continue
    itr=itr+1
    if(itr > (maxitr+1))then
        write(*,*)'Present iteration=',itr
        write(*,*)'Maximum allowed iterations=',maxitr
        print*,'SCF','No convergence achieved in Max iterations'
	stop
    end if

! Construct the density matrix. Presently it assumes that the
! orbitals are doubly occupied. For open-shell systems modifications
! will be needed

    call density(den,orb,norb,msite)

! Construct the Fock matrix and compute the energy.

    dif2=dif1
    escf1=escf

! First iteration with Hueckel Hamiltonian

    if(itr == 1)then
      fock=fock+t
  
      else
    
    ! Construct the Fock matrix and compute the energy.
    
        call fockmat(t,h,v,den,fock,msite,escf)
   
    end if

!
    dif1=escf-escf1

! Quit if the convergence has been achieved

    iflag1=0
    iflag2=0
    if((dabs(dif1) < conv) .AND. (dabs(dif2) < conv))iflag1=1

! Damp the Fock matrix, if needed. 

    if(idamp == 1 .AND. itr > 1)then
        fock=xdamp*fock+(1.d0-xdamp)*fock_o
	fock_o=fock
      end if
   
!   call plblk(fock,msite,0,'site',6)

! Diagonalize the Fock matrix

    call dspevx('V','I','U',msite,fock,VL,VU,1,norb,abstol,mevlfnd,&
             eval,orb,msite,work,iwork,ifail,info)
	     
    
      if(info/=0)then
        print*,"SCF_RHF;matrix diagonalization failed by DSPEVX"
        print*,"info=",info
        stop
      end if


! Quit after first iteration if the Hueckel Hamiltonian was chosen
!
      if(iham==4)goto 101

! eigenvalue convergence check

    if(itr >= 2)then
    
        xdiff=-2.d0

        do i=1,norb
                
                xdiff=max(xdiff,dabs(eval(i)-eval_o(i)))

        end do
    
        if(xdiff <= (conv*1000.d0))iflag2=1
    
    end if

    write(*,11)itr,escf+enuc,dif1,dif2,xdiff
    11 format(1x,i6.1,3x,e14.7,2x,e14.7,2x,e15.7,2x,e14.7)

    if(iflag1 == 1 .AND. iflag2 == 1)goto 101
    eval_o=eval

!
    goto 3

    
! Convergence achieved


    101 continue

    write(*,*)' '
    write(*,102)itr
    102 format(1x,'Convergence achieved after ',i6.1,' iterations')
    write(*,*)' '
    write(*,104)escf+enuc
    104 format(1x,'Final SCF energy: ',2x,f18.12)
    write(*,*)' '
    write(*,204)enuc
    204 format(1x,'Nuc-Nuc energy: ',2x,f18.12)
    write(*,*)' '

   ! Construct the Fock matrix again and diagonalize it to obtain all eigenvalues and eigenvectors (occupied and virtual)  .
    
        call fockmat(t,h,v,den,fock,msite,escf)
	deallocate(eval_o,iwork,ifail)

        call DSPEV('V','U',msite,fock,eval,orb,msite,work,info)
    !
       if(info/=0)then
        print*,"SCF_RHF;matrix diagonalization failed by DSPEV"
        print*,"info=",info
        stop
       
       end if
    !
!writing the eigen values and orbitals to a binary file
!
    iorbu=11
    open(unit=iorbu,file='eval_orb.dat',form='unformatted',status='unknown')

       write(iorbu)msite,msite
       write(iorbu)(eval(i),i=1,msite)
       write(iorbu)((orb(j,i),j=1,msite),i=1,msite)
       
    close(unit=iorbu,status='keep')
! Output the eigenvalues

    write(*,*)'Eigen values of the Fock Matrix'
    write(*,*)' '
    call prblk(eval,msite,msite,1,0,0,'#   ','E')
    write(*,*)' '
  

  
    iorbu=20
    call writorb(iorbu,'ORB001.DAT',orb,msite,msite)
    if(iprorb == 0)goto 67
    write(*,*)'  '
    write(*,*)'Occupied Orbitals'
    write(*,*)'  '
    call prblk(orb,msite,msite,norb,0,0,'Site','ORB')
    write(*,*)'  '
    write(*,*)'Virtual orbitals'
    write(*,*)'  '
    call prblk(orb(1,norb+1),msite,msite,norb,0,norb,'Site','ORB')
    67 continue

    deallocate(fock,den,work)
    if(idamp>0)deallocate(fock_o)
    return
    end subroutine scf_rhf

! ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

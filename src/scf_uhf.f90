! ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! This routine solves the SCF eigenvalue problem by diagonalizing the HF
! matrix using a Householder diagonalizer.

    subroutine scf_uhf(t,h,v,msite,idamp,xdamp,maxitr,iprorb,conv,nalpha,&
                   nbeta,orb_alpha,orb_beta,eval_alpha,eval_beta,enuc)

! ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
IMPLICIT NONE

   
    integer::ok,iorbu,isite,i,j,itr,iflag1,iflag2,mevlfnd,&
       info
    integer,intent(in)::msite,idamp,maxitr,iprorb,nalpha,nbeta
    real(kind=8),intent(in)::t((msite*(msite+1))/2),enuc,&
               h((msite*(msite+1))/2),v((msite*(msite+1))/2),xdamp,conv
    real(kind=8),intent(out)::orb_alpha(msite,msite),orb_beta(msite,msite),eval_alpha(msite),&
                           eval_beta(msite)
    real(kind=8),allocatable::fock_alpha(:),&
                         fock_beta(:),fock_alpha_o(:),fock_beta_o(:),&
                         den(:),den_alpha(:),den_beta(:),&
			 eval_alpha_o(:),eval_beta_o(:)       
   
    real(kind=8)::dif1,dif2,xdiff,escf,escf1,abstol,VL,VU,dlamch
             
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


    allocate (den((msite*(msite+1))/2),den_alpha((msite*(msite+1))/2),&
           den_beta((msite*(msite+1))/2),&
	   eval_alpha_o(nalpha),eval_beta_o(nbeta),stat=ok)
      if(ok/=0)then
        print*,'NOT ENOUGH MEMORY IS AVAILABLE FOR ALLOCATING den or eval'
        stop
      end if
      
      allocate (fock_alpha((msite*(msite+1))/2),fock_beta((msite*(msite+1))/2),stat=ok)
         if(ok/=0)then
           print*,'NOT ENOUGH MEMORY IS AVAILABLE FOR ALLOCATING fock_alpha and fock_beta'
           stop
         end if
   


      if(idamp>0)then
         allocate (fock_alpha_o((msite*(msite+1))/2),fock_beta_o((msite*(msite+1))/2),stat=ok)
         if(ok/=0)then
           print*,'NOT ENOUGH MEMORY IS AVAILABLE FOR ALLOCATING fock_o'
           stop
         end if
      end if

    write(*,*)' '
    write(*,255)
    255 format(1x,'Itr(i)',5x,'E(scf)',8x,'Ediff1(i-1)',7x,'Ediff2(i)', &
    7x,'Evldiff(i)')
    write(*,*)' '


! Starting occupied vectors are zero

           orb_alpha=0.d0
	   orb_beta=0.d0
           fock_alpha=0.d0
	   fock_beta=0.d0
	   eval_alpha=0.d0
	   eval_beta=0.d0
	   eval_alpha_o=0.d0
	   eval_beta_o=0.d0

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
 
   
    call density_uhf(den,den_alpha,den_beta,orb_alpha,orb_beta,nalpha,nbeta,msite)
 
! Construct the Fock matrix and compute the energy.

    dif2=dif1
    escf1=escf

! First iteration with Hueckel Hamiltonian

    if(itr == 1)then
      fock_alpha=fock_alpha+t
      fock_beta=fock_beta+t
  
      else
    
    ! Construct the Fock matrix and compute the energy.
 
        call fockmat_uhf(t,h,v,den,den_alpha,den_beta,fock_alpha,fock_beta,msite,escf)
 

    end if

!
    dif1=escf-escf1

! Quit if the convergence has been achieved

    iflag1=0
    iflag2=0
    if((dabs(dif1) < conv) .AND. (dabs(dif2) < conv))iflag1=1

! Damp the Fock matrix, if needed. 

      if(idamp == 1 .AND. itr > 1)then
        fock_alpha=xdamp*fock_alpha+(1.d0-xdamp)*fock_alpha_o
	fock_alpha_o=fock_alpha
        fock_beta=xdamp*fock_beta+(1.d0-xdamp)*fock_beta_o
	fock_beta_o=fock_beta	
      end if
      

! Diagonalize the Fock matrix

    

!    call plblk(fock_alpha,msite,0,'site',6)
    
    call dspevx('V','I','U',msite,fock_alpha,VL,VU,1,nalpha,abstol,mevlfnd,&
             eval_alpha,orb_alpha,msite,work,iwork,ifail,info)

      if(info/=0)then
        print*,"SCF_UHF;matrix diagonalization failed by DSPEVX for alpha sector"
        print*,"info=",info
        stop
      end if	     
 	     
    call dspevx('V','I','U',msite,fock_beta,VL,VU,1,nbeta,abstol,mevlfnd,&
             eval_beta,orb_beta,msite,work,iwork,ifail,info)	     
    
      if(info/=0)then
        print*,"SCF_UHF;matrix diagonalization failed by DSPEVX for beta sector"
        print*,"info=",info
        stop
      end if

! Orbital convergence check

    if(itr >= 2)then
    

	xdiff=-2.d0
	do i=1,nalpha
		xdiff=max(xdiff,dabs(eval_alpha(i)-eval_alpha_o(i)))
	end do
	do i=1,nbeta
		xdiff=max(xdiff,dabs(eval_beta(i)-eval_beta_o(i)))
	end do
    
        if(xdiff <= (conv*1000.d0))iflag2=1
    
    end if
 
    write(*,11)itr,escf+enuc,dif1,dif2,xdiff
    11 format(1x,i6.1,3x,e14.7,2x,e14.7,2x,e15.7,2x,e14.7)

    if(iflag1 == 1 .AND. iflag2 == 1)goto 101
    
    eval_alpha_o=eval_alpha
    eval_beta_o=eval_beta
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
    
        call fockmat_uhf(t,h,v,den,den_alpha,den_beta,fock_alpha,fock_beta,msite,escf)
	
	deallocate(eval_alpha_o,eval_beta_o,iwork,ifail)
	

        call DSPEV('V','U',msite,fock_alpha,eval_alpha,orb_alpha,msite,work,info)
   !
       if(info/=0)then
        print*,"SCF_RHF;matrix diagonalization failed by DSPEV for alpha sector"
        print*,"info=",info
        stop
       end if
       
      !
    !The eigen values and eigen vectors
  !
  iorbu=11
  open(unit=iorbu,file='eval_orb.dat',form='unformatted',status='unknown')
  !     
  !
      write(iorbu)msite,msite,msite
      write(iorbu)(eval_alpha(i),i=1,msite)
      write(iorbu)((orb_alpha(j,i),j=1,msite),i=1,msite)
!      
    !	
	call DSPEV('V','U',msite,fock_beta,eval_beta,orb_beta,msite,work,info)
    !
       if(info/=0)then
        print*,"SCF_RHF;matrix diagonalization failed by DSPEV for beta sector"
        print*,"info=",info
        stop
       end if
    !
       write(iorbu)(eval_beta(i),i=1,msite)
       write(iorbu)((orb_beta(j,i),j=1,msite),i=1,msite)
 !
   close(unit=iorbu,status='keep')       
! Output the eigenvalues

    write(*,*)'Eigen values from alpha spins of the Fock Matrix'
    write(*,*)' '
    call prblk(eval_alpha,msite,msite,1,0,0,'#   ','E')
    write(*,*)' '
 
    write(*,*)'Eigen values from beta spins of the Fock Matrix'
    write(*,*)' '
    call prblk(eval_beta,msite,msite,1,0,0,'#   ','E')
    write(*,*)' '
     
    iorbu=20
    call writorb(iorbu,'ORBALPHA.DAT',orb_alpha,msite,msite)
    call writorb(iorbu,'ORBBETA.DAT',orb_beta,msite,msite)
        
    if(iprorb == 0)goto 67
    write(*,*)'  '
    write(*,*)'Occupied Orbitals for alpha spins'
    write(*,*)'  '
    call prblk(orb_alpha,msite,msite,nalpha,0,0,'Site','ORB_A')
    write(*,*)'  '
    write(*,*)'Virtual orbitals for alpha spins'
    write(*,*)'  '
    call prblk(orb_alpha(1,nalpha+1),msite,msite,nalpha,0,nalpha,'Site','ORB_A')
    write(*,*)'  '
    write(*,*)'Occupied Orbitals for beta spins'
    write(*,*)'  '
    call prblk(orb_beta,msite,msite,nbeta,0,0,'Site','ORB_B')
    write(*,*)'  '
    write(*,*)'Virtual orbitals for beta spins'
    write(*,*)'  '
    call prblk(orb_beta(1,nbeta+1),msite,msite,nbeta,0,nbeta,'Site','ORB_B')
    67 continue

    deallocate(fock_alpha,fock_beta,den,&
            den_alpha,den_beta,work)
    if(idamp>0)deallocate(fock_alpha_o,fock_beta_o)        
    return
    end subroutine scf_uhf

! ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

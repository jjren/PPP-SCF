! ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! This routine computes the optical absorption both for the SCF and
! the SCI states.

    subroutine optics(r_site,orb,e_ci,eval,ci_vec, &
    msite,mocc_orb,mvirt_orb,morb,mdimh,mdimd,isci,nabsorb,&
    nexcite,iabsorb,ndim,idim,omega1,omega2,domega,gamma,scale,&
    ihole,ipart)

! ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

IMPLICIT NONE

    integer,intent(in)::mocc_orb,mvirt_orb,isci,nabsorb,msite,&
                   nexcite(*),mdimh,iabsorb(*),morb,&
		   mdimd,ihole(*),ipart(*),ndim,idim(ndim)
    real(kind=8),intent(in)::r_site(3,msite),orb(msite,*),&
                   eval(msite),ci_vec(mdimh,mdimh),e_ci(mdimh),&
		   omega1,omega2,domega,gamma,scale
    integer::ntorb,ispecu,istate,ij,isite,iaxis,jscf,icsf,jcsf,&
             ijcsf,jabsorb,icount,ncount,jstate,ok,ndimh,i,j,&
	     nsite
    real(kind=8)::xcut,dx,dy,dxij,dyij,dzij,dz
    real(kind=8),allocatable::dmu_scf(:,:),dmu_csf(:,:),buff(:,:)
    character(80) :: spcfil

    ndimh=mdimh
    allocate(dmu_scf(mdimd,ndim),stat=ok)
      if(ok/=0)then
        print*,'NOT ENOUGH MEMORY IS AVAILABLE FOR ALLOCATING dmu_scf'
        stop
      end if     	  

   if(isci>0)then
    allocate(dmu_csf((mdimh*(mdimh+1))/2,ndim),buff(mdimh,ndim),stat=ok)
      if(ok/=0)then
        print*,'NOT ENOUGH MEMORY IS AVAILABLE FOR ALLOCATING dmu_csf,buff'
        stop
      end if  
    end if
    
    ntorb=mocc_orb+mvirt_orb
    
    xcut=0.1d0


! Absorption spectrum will be written in a separate file

    ispecu=91
    spcfil='spec001.dat'
    open(unit=ispecu,file=spcfil,form='formatted', &
    status='unknown',err=31)
    goto 32

    31 print*,'Error in Optics','Error opening Spectrum File'

    32 continue

    write(*,*)' '
    write(*,42)xcut
    42 format(1x,'Single particle dipoles with cutoff: ',f5.2)
    write(*,*)' '

 
    if(isci == 0)then
        write(ispecu,*)'#'
        write(ispecu,41)
        41 format(1x,'#Single particle optical absorption spectrum')
        write(ispecu,*)'# '
        write(ispecu,*)'# Omega   Probability'
        write(ispecu,*)'# '
        do istate=1,nabsorb
            if(nexcite(istate) == 0)then
                ij=0
                do i=1,mocc_orb
                   do j=mocc_orb+1,ntorb

                        ij=ij+1
                        dx=0.d0
                        dy=0.d0
			
                        do isite=1,msite
                            dx=dx+orb(isite,i)*orb(isite,j)*r_site(idim(1),isite)
                            dy=dy+orb(isite,i)*orb(isite,j)*r_site(idim(2),isite)
                        end do
                    
                        dmu_scf(ij,1)=dx
                        dmu_scf(ij,2)=dy
                        
			if(ndim==3)then
			  dz=0.d0
			  do isite=1,msite
		            dz=dz+orb(isite,i)*orb(isite,j)*r_site(idim(3),isite)
                          end do
			  dmu_scf(ij,3)=dz
                         end if
			
			if(ndim<=2)then 
                         if(dabs(dx) > xcut .OR. dabs(dy) > xcut) &
                         write(*,222)i,j,(eval(j)-eval(i))/scale,dx,dy
			else
			 if(dabs(dx) > xcut .OR. dabs(dy) > xcut.OR. dabs(dz) > xcut) &
                         write(*,223)i,j,(eval(j)-eval(i))/scale,dx,dy,dz
                        end if
                    
                    end do
                end do
            
            ! Absorption spectrum
            
                call spctrm_a(eval(1),eval(1+mocc_orb),mocc_orb,mvirt_orb &
                ,dmu_scf,omega1,omega2,domega,gamma, &
                ispecu,scale,ndim)
		
            elseif(nexcite(istate) == 1)then
                ij=0
                do  i=1,ntorb
                    do  j=1,i
                        ij=ij+1
                    
                        do  iaxis=1,ndim
                            dmu_scf(ij,iaxis)=0.d0
                        end do
                    
                        do  isite=1,msite
                            do  iaxis=1,ndim
                                dmu_scf(ij,iaxis)=dmu_scf(ij,iaxis)+ &
                                orb(isite,i)*orb(isite,j)*r_site(idim(iaxis),isite)
                            end do
                        end do

                    
                    
                    end do
                end do
                call spctrm_1ex(eval(1),nsite,ihole(istate), &
                ipart(istate),mocc_orb,dmu_scf,ndim,omega1, &
                omega2,domega,gamma,ispecu,scale)
            else
                print*,'Wrong value of nexcite=',nexcite(istate)

            end if
        end do
    end if

    if(isci == 0)return

    ij=0
    do i=1,ntorb
        do j=1,i
            ij=ij+1
        
            do iaxis=1,ndim
                dmu_scf(ij,iaxis)=0.d0
            end do
        
            do isite=1,msite
                do iaxis=1,ndim
                    dmu_scf(ij,iaxis)=dmu_scf(ij,iaxis)+orb(isite,i)* &
                    orb(isite,j)*r_site(idim(iaxis),isite)
                end do
            end do
        
        
        end do
    end do


    222 format(' i = ',i4.1,2x,'j = ',i4.1,2x,'energy =',1x,f8.4,2x, &
    'dx = ',f8.4,2x,'dy = ',f8.4)

    223 format(' i = ',i4.1,2x,'j = ',i4.1,2x,'energy =',1x,f8.4,2x, &
    'dx = ',f8.4,2x,'dy = ',f8.4,2x,'dz = ',f8.4)
    
! If SCI calculation was performed, compute absorption for the many
! -electron excited states from the states chosen

    write(*,*)' '
    write(*,43)xcut
    43 format(1x,'SCI dipole matrix elements with cutoff: ',1x,f5.2)
    write(*,*)' '
    call dipmat(dmu_csf,dmu_scf,mocc_orb,mvirt_orb,morb,mdimh,ndim)


! Loop over the states from which the absorption is to be computed

    do jabsorb=1,nabsorb
    
        istate=iabsorb(jabsorb)
    
        write(*,*)' '
        write(*,102)istate-1,e_ci(istate)
        write(*,*)' '
    
    ! Perform half transformation over the CSFs
    
        do iaxis=1,ndim
        
            do icsf=1,ndimh
            
                buff(icsf,iaxis)=0.d0
            
                do jcsf=1,icsf
                    ijcsf=(icsf*(icsf-1))/2+jcsf
                    buff(icsf,iaxis)=buff(icsf,iaxis)+ &
                    dmu_csf(ijcsf,iaxis)*ci_vec(jcsf,istate)
                end do
            
                do jcsf=icsf+1,ndimh
                
                    ijcsf=(jcsf*(jcsf-1))/2+icsf
                    buff(icsf,iaxis)=buff(icsf,iaxis)+ &
                    dmu_csf(ijcsf,iaxis)*ci_vec(jcsf,istate)
                
                end do
            
            end do
        
        end do
    
        write(*,*)'# '
        write(*,25)istate
        25 format(1x,'#SCI optical absorption spectrum from state: ',i4.1)
        write(*,*)'# '
        write(*,*)'# Omega   Probability'
        write(*,*)'# '
    
    ! Now use the array dmu_scf to store the dipole matrix elements
    ! w.r.t to the SCI states
    
        icount=0
        ncount=ndimh-istate
        dmu_scf=0.d0
        do jstate=istate+1,ndimh
        
            icount=icount+1
	    dxij=0.d0
            dyij=0.d0
	    dzij=0.d0
!                  
            do icsf=1,ndimh
             do iaxis=1,ndim
	       dmu_scf(icount,iaxis)=dmu_scf(icount,iaxis)+buff(icsf,iaxis)*ci_vec(icsf,jstate)
             end do
	      dxij=dxij+buff(icsf,1)*ci_vec(icsf,jstate)
              dyij=dyij+buff(icsf,2)*ci_vec(icsf,jstate)
              if(ndim==3)dzij=dzij+buff(icsf,3)*ci_vec(icsf,jstate)
	
            end do
        

           if(ndim<=2)then
            if(dabs(dxij) > xcut .OR. dabs(dyij) > xcut) &
            write(*,101)jstate-1,e_ci(jstate),dxij,dyij
	   else
	    if(dabs(dxij) > xcut .OR. dabs(dyij) > xcut.OR. dabs(dzij) > xcut) &
            write(*,103)jstate-1,e_ci(jstate),dxij,dyij,dzij
	   end if 
	    
        
        end do
    
    ! SCI absorption spectrum
    
      call spctrm_a(e_ci(istate),e_ci(istate+1),1,ncount, &
        dmu_scf,omega1,omega2,domega,gamma, &
        ispecu,scale,ndim)
    
    end do

    101 format('Target State:',1x,i4.1,2x,'energy = ',1x,f8.4,3x, &
    'x-dipole=',1x,f8.4,3x, 'y-dipole = ',1x,f8.4)

    102 format('Optical Absorption Data from State #',2x,i5.1,2x, &
    'with energy (eV):',1x,f10.4)
    
       103 format('Target State:',1x,i4.1,2x,'energy = ',1x,f8.4,3x, &
    'x-dipole=',1x,f8.4,3x, 'y-dipole = ',1x,f8.4,3x, 'z-dipole = ',1x,f8.4)


    return
    end subroutine optics

! cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

! ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

! This is the master routine for the general purpose Pariser-Parr-Pople model 
! Hamiltonian based electronic structure program written in Fortran 90 language.
! A Fortran 77 version of this program was developed by Alok Shukla during his 
! postdocs stints in the group of S. Mazumdar, at Physics Department, University 
! of Arizona, during 1998-99. That program lack the capability of doing unrestricted 
! Hartree-Fock (UHF) calculations, and thus, was limited only to closed-shell systems.
! Moreover, Fortran 77 language did not allow dynamic memory allocation, resulting in
! a dependence on compiled time dimensional parameters. 
!
! In order to go beyond these shortcomings, we, Priya Sony and Alok Shukla, undertook 
! a rewriting of the old code in Fortran 90 language, along with the implementation of 
! the dynamic memory allocation. Additionally, UHF module was also added. This task was
! carried out during May--August 2009, at Department of Physics, Indian Institute of 
! Technology Bombay, Mumbai,India, where Priya Sony was a postdoc in the group of Alok Shukla.

! copyright: P. Sony and A. Shukla (2009)

! contact details:
! email: psony11@gmail.com, shukla@phy.iitb.ac.in
! Postal Address:
! Alok Shukla
! Department of Physics
! Indian Institute of Technology Bombay
! Powai, Mumbai. 400076
! INDIA

! cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
program ppp_main
!
IMPLICIT NONE
!
integer::iham=0,ihampar=0,irhf=0,iuhf=0,isci=0,iciprp=0,&
        iorbdn=0,ioptics=0,nalphabeta=0,&
        nalpha=0,nbeta=0,inlo=0,inldel=0,ibands=0,&
	iprorb=0,icorfnc=0,ichngcl=0,iefield=0,&
	norbdo,norbdv,nread=0,icnto=0,j,j1,j2,icntv=0,&
	nabsorb=0,i,ndip,&
	ndelo1=0,ndelo2=0,ndelv1=0,ndelv2=0,idenfl=0,nosite,icharge,natom,nunit,&
	ibipart,maxiter,idamp=0,ndelatm=0,nelec,nsite,norb,igreen,&
	norbci=0,itran_ci=0,norbfr,itran_fr=0,idipint=0,ndir=0,ndirx,&
        natmi=0,natmj=0,ifcor=0,imin,imax,nni,ncount,ii,jmin,jmax,nnj,&
        natmic=0,natmjc=0,nij,itij,iihij,ivij,ieps,norb_o=0,norb_v=0,&
        ifkr,ifki,izkr,izki,ie,ifv1,ifv2,ifm1,nkval,n1n1,ntorb,ndimh,ndimd,ok,&
	maxitr,ifcorc=0,nsitex,nocc,nrem,nbas,nvrt,ifile,nn12,n1234,&
	ieocc,ievrt,idipmt,ir1int,ir2int,ibuf,nbuf,ndim=0,ialpbet,nfreez=0,ifrez,&
	ndelete=0,idel,idelmin,idelmax,msite,morb,morb_a,morb_b,icount,iorbu,len


 real(kind=8)::u=0.d0,r0=0.d0,diel=0.d0,gamma=0.d0,omega1=0.d0,&
               omega2=0.d0,domega=0.d0,scale=0.d0,charge=0.d0,conv=0.d0,xdamp=0.d0,&
             fkmin=0.d0,fkmax=0.d0,fkstp=0.d0,xcor=0.d0,unew=0.d0,r0new=0.d0,dielnew=0.d0,&
	     vv=0.d0,enuc=0.d0,escf=0.d0

 character(80)::hamltn,method,card,line,filenlo,file_orb1,file_orb2,fileopt
 character(32)::input_file

 integer,allocatable::iorbdo(:),iabsorb(:),nexcite(:),ihole(:),ipart(:),&
		 idip(:),iaxis(:),idelatm(:),idir(:),itype(:),jtype(:),itypec(:),&
		 jtypec(:),ihij(:),iosite(:),idim(:),ifreez(:),idelete(:)

 real(kind=8),allocatable::rr_atom(:,:),rorg(:),angle(:),tvec(:),rshft(:),&
                      efield(:),tij(:),e(:,:),hij(:),swape(:),swap(:,:),&
		      vij(:),eps(:),rr_site(:,:),orb(:,:),orb_alpha(:,:),orb_beta(:,:),&
		      eval(:),eval_alpha(:),eval_beta(:),e_ci(:),&
		      ci_vec(:,:),ep(:)

 complex*16,allocatable::fk(:)

!***********************************************************************************
!
! Read the title
!
!  call getarg(15,input_file)
!  print*,'input file name=',input_file
  
  read(*,*)line
  line=adjustl(line)
  print*,line
!
! Which Hamiltonian to use: Hubbard, Ext Hubbard, PPP?
!
  read(*,*)line
  read(*,*)hamltn
  hamltn=adjustl(hamltn)
!
  iham=0
  ihampar=0
  if(hamltn == 'PPP')then
        iham=1
        print*,'PPP Hamiltonian to be used'
  elseif(hamltn == 'EXTHUB')then
        iham=2
        print*,'Extended Hubbard Hamiltonian to be used'
  elseif(hamltn == 'HUBBARD')then
        iham=3
        print*,'Hubbard Hamiltonian to be used'
  elseif(hamltn == 'HUECKEL')then
        iham=4
        print*,'Hueckel Hamiltonian to be used'
  else
        print*,'Hamiltonian given by you: ',hamltn
        print*,'Valid options are: PPP, EXTHUB, OR HUBBARD'
        stop
  end if
! If PPP, is it Ohno, Mat(aga)-Nis(himoto), or exponential
! parameterization?

    if(iham == 1)then
       read(*,*)line
       read*,hamltn
       hamltn=adjustl(hamltn)
       if(hamltn(1:4) == 'OHNO')then
            ihampar=1
            print*,'Along with Ohno parameterization'
        elseif(hamltn(1:6) == 'MATNIS')then
            ihampar=2
            print*,'Along with Mataga-Nishimoto &
            parameterization'
        elseif(hamltn(1:3) == 'EXP')then
            ihampar=3
            print*,'Along with exponential parameterization'
        else
            print*,'PPP parameterization given by you: ',hamltn
            print*,'Allowed options: OHNO, MATNIS, OR EXP'
            stop
        end if
    end if

!
!parameters of the Hamiltonian
!
! Read the values of U, and, if needed, nrange (# of neighbors included
! in the extended Hubbard model) or r_0.
!
    if(iham == 1)then
        read(*,*)line
	if(ihampar == 1)then
	   read*,card
	   card=adjustl(card)
	   if(card(1:5)=='STAND')then
	      u=11.13d0
	      r0=1.2785884d0
	      diel=1.d0
              print*,'STANDARD OHNO PARAMETERS ARE BEING USED'
	   elseif(card(1:3)=='SCR')then
	      u=8.d0
	      r0=1.2785884d0
	      diel=2.d0   
	      print*,'SCREENED PARAMETERS ARE BEING USED'
	   elseif(card(1:4)=='PARA')then
	     read*,u,r0,diel
	   else
	     print*,'THE ALLOWED INPUT IS STAND,SCR, OR PARA'
	     print*,'BUT YOU GAVE',card
	     stop
	  end if     
	else
	   read*,u,r0
	end if
    elseif(iham == 2)then
         read*,u,vv
    elseif(iham==3)then
         read(*,*)line
	 read*,u
    end if	  
    
    write(*,101)u
    101 format(1x,'U=',1x,f8.4)       

    if(iham == 1)then
        write(*,102)r0
        write(*,202)diel
        102 format(1x,'r0=',1x,f8.4)
        202 format(1x,'Dielectric Constant =',1x,f8.4)
    elseif(iham == 2)then
        write(*,103)vv
        103 format(1x,'V=',1x,f8.4)
    end if
!    
! Read the ionicity (electric charge) of the system. For the neutral
! system the ionicity will be zero, and so on.
!
    icharge=0
    read(*,*)line
    read*,icharge
    charge=dble(icharge)

    If(icharge == 0)then
        print*,'The system is uncharged'
    else
        write(*,4)icharge
        4 format(1x,'Charge on the System:',1x,i3.1)
    end if

! Read the total number of atoms in the unit cell

    natom=0
    read(*,*)line
    read*,natom
    write(*,5)natom
  5 format(1x,'Total no. of atoms/cell:',1x,i3.1)
!  
! Read the coordinates of the atoms within the origin
! unit cell.
!
    allocate (rr_atom(3,natom),stat=ok)
      if(ok/=0)then
        print*,'NOT ENOUGH MEMORY IS AVAILABLE FOR ALLOCATING rr_atom'
	stop
      end if	
   allocate (rorg(3),stat=ok)
      if(ok/=0)then
        print*,'NOT ENOUGH MEMORY IS AVAILABLE FOR ALLOCATING rorg'
	stop
      end if      
   allocate (angle(3),stat=ok)
      if(ok/=0)then
        print*,'NOT ENOUGH MEMORY IS AVAILABLE FOR ALLOCATING angle'
	stop
      end if 
   allocate (iaxis(3),stat=ok)
      if(ok/=0)then
        print*,'NOT ENOUGH MEMORY IS AVAILABLE FOR ALLOCATING iaxis'
	stop
      end if 
!      
        call r_atom(rr_atom,angle,rorg,iaxis,natom)
    
!    
! What type of calculation? Possibilities should include: RHF, UHF, SCI,
! 

    irhf=0
    iuhf=0
    isci=0
    iciprp=0
    iorbdn=0
    ioptics=0
    inlo=0
    inldel=0
    ibands=0
    iprorb=0
    icorfnc=0
    ichngcl=0
    iefield=0
    
    if(iham == 4)irhf=1

    2 continue
    
    read(*,*)line
    read*,method(1:6)
    method=adjustl(method)

    if(method(1:3) == 'RHF')then
        irhf=1
        print*,'RHF Calculations will be performed'
        goto 2
    elseif(method(1:3) =='UHF')then
        iuhf=1
	print*,'UHF Calculations will be performed'
!       Give Number of electrons with up spin (nalpha) and down spins (nbeta)
        read(*,*)line
        read(*,*)nalpha,nbeta
        print*,'nalpha,nbeta',nalpha,nbeta
        goto 2	
    elseif(method(1:3) == 'SCI')then
        isci=1
        print*,'Singles CI Calculations will be performed'
        goto 2    
    elseif(method(1:6) == 'OPTICS')then
        ioptics=1
        print*,'Linear optical Properties will be computed'
!    read the components of dipole moments
        read(*,*)line
        read(*,*),ndim
 	allocate(idim(ndim),stat=ok)
	    if(ok/=0)then
              print*,'allocation failed for idim'
              stop
            end if

           read(*,*)(idim(i),i=1,ndim)
	
        fileopt='spec001.dat'

        print*,'Data for single-particle optics calculations will'
        print*,'put in the file',fileopt
 
    ! Read the line-width parameter
    
        read*,gamma,omega1,omega2,domega,scale
        if(scale <= 1.d-6)scale=1.d0
        write(*,171)gamma
        171 format(1x,'Wigner Line Width (eV): ',2x,f8.4)
    
    ! Read the state from which the absorption is to be computed
    
        if(isci > 0)then
	    read*,line
            read*,nabsorb
 	    if(nabsorb>0)then
	     allocate(iabsorb(nabsorb),stat=ok)
	       if(ok/=0)then
                 print*,'allocation failed for iabsorb'
                 stop
               end if

              read(*,*)(iabsorb(i),i=1,nabsorb)
        
              do i=1,nabsorb
                write(*,104)iabsorb(i)
              end do
            end if
        else
            read*,line
	    read*,nabsorb

	    if(nabsorb==0)then
	      print*,'nabsorb=0 is not allowed.' 
	      print*,'The minimum allowed value is 1.'
	      print*,'To compute ground state absorption set &
	            nabsorb=1, followed by 0 on the next line.'
	      stop
	    end if
	   
 	    if(nabsorb>0)then
              allocate(nexcite(nabsorb),ihole(nabsorb),ipart(nabsorb),stat=ok)
	      if(ok/=0)then
                print*,'allocation failed for nexcite or ihole or ipart'
                stop
              end if
	    	    
            do i=1,nabsorb
                read*,nexcite(i)
                if(nexcite(i) > 0)then
                    if(nexcite(i) > 1)then
                        print*,'sr input: nexcite > 1 not allowed'
                        stop
                    end if
                    read*,ihole(i),ipart(i)
                end if
            end do
            end if
        end if
    
        104 format(1x,'From the state #:',1x,i4.1)
        goto 2
    
    elseif(method(1:3) == 'NLO')then
    
        inlo=1
    print*,'Output data will be written for NLO calculations in file NLO001.DAT'

!read the components of dipole moments
        read(*,*)line
        read(*,*),ndip
	
	allocate(idip(ndip),stat=ok)
	    if(ok/=0)then
              print*,'allocation failed for idip'
              stop
            end if

           read(*,*)(idip(i),i=1,ndip)
	
        filenlo='NLO001.DAT'

        print*,'Data for single-particle NLO calculations will'
        print*,'put in the file',filenlo
!
!  Set operations like orbital deletion, orbital density analysis, etc.
!
        read(*,*)line
	read(*,*)card
	card=adjustl(card)
        if(card(1:6) == 'ORBDEL')then
            inldel=1
	    print*,'Orbitals will be deleted for NLO calculations'
            read*,ndelo1,ndelo2
            read*,ndelv1,ndelv2
            print*,'ndelo1,ndelo2',ndelo1,ndelo2
	    print*,'ndelv1,ndelv2',ndelv1,ndelv2
        end if
    
        goto 2
    
    elseif(method(1:6) == 'ORBDEN')then
    
    ! sitewise orbital density analysis
    
        iorbdn=1
	print*,'Orbital density analysis will be performed'
    
    ! idenfl=0 only the specified sites are to be analyzed
    ! idenfl=1 the sites given along with their periodic copies will be
    ! analyzed.

        read*,idenfl,nosite
    
        if(nosite > nsite)print*,'Error in Input', &
               'nosite > nsite'
        allocate(iosite(nosite),stat=ok)
	    if(ok/=0)then
              print*,'allocation failed for iosite'
              stop
            end if

        read(*,*)(iosite(i),i=1,nosite)
        print*,'orbital density analysis will be done'
        goto 2
    
    elseif(method(1:4) == 'BAND')then
        if(iham /= 4)print*,'Error in Input', &
        'Band structure calculation possible only with Hueckel model'
        print*,'Band structure to be computed with periodic boundary conditions'
    
        ibands=1
	print*,'Band structure will be calculated for the Hueckel model'
        read(*,*)fkmin,fkmax,fkstp
        goto 2
    
    elseif(method(1:4) == 'CIPR')then
    
      iciprp=1
      print*,'Files will be written for subsequent CI calculations'

      if(iuhf==1)then      ! If UHF method was used, find out whether alpha
       read*,ialpbet         ! or beta type orbitals will be used for 
       if(ialpbet==1)then    ! correlation treatment
         print*,'For correlation treatment Alpha-type orbitals will be used'
       elseif(ialpbet==2)then
         print*,'For correlation treatment Beta-type orbitals will be used'
       else
         print*,'For CIPRP with UHF method you gave ialpbet=',ialpbet
         print*,'Valid options are 1 and 2'
         print*,'Stopping the job'
         stop
        end if
       end if 
!
      read*,norbfr,ifrez           !if in CI calculation, some orbitals to be.
      !frozen. In that case total no. of orbitals to be
      !frozen will be read as input
      !
      !Allocating freezing oritals
      !
      nfreez=norbfr
      
      if(nfreez>0)then
        allocate(ifreez(nfreez),stat=ok)
        if(ok /= 0)then
           WRITE(*,*)'semi_emp:Allocation fails for ifreez'
           stop
        end if
        !
        if(ifrez==0)then
           do i=1,nfreez
              ifreez(i)=i     
           end do

        else                        
           read(*,*)(ifreez(i),i=1,nfreez)
           
        end if
      end if
      read*,ndelete,idel           
      !if in CI calculation, some orbitals need to be deleted, read the
      ! input here. 
      !
      !Take care of orbitals to be deleted
      !
      if(ndelete>0)then
        allocate(idelete(ndelete),stat=ok)
        if(ok /= 0)then
           WRITE(*,*)'semi_emp:Allocation fails for idelete'
           stop
        end if
        !
        if(idel==0)then
           read*,idelmin,idelmax
           if((idelmax-idelmin+1)/=ndelete)then
              print*,'semi_emp: idelmin, idelmax, not consistent with'
              print*,'the value of ndelete'
              print*,'idelmin,idelmax,ndelete',idelmin,idelmax,ndelete
              stop
           end if
           icount=0
           do i=idelmin,idelmax
              icount=icount+1
              idelete(icount)=i     
           end do
        else                        
           read*,(idelete(i),i=1,ndelete)
           !
        end if
      end if
    
    ! See if the DIPOLE matrix elements are needed
    
        read(*,*)line
        read(*,*)card
	card=adjustl(card)  
        if(card(1:6) == 'DIPINT')then
        
            idipint=1
	    print*,'dipole matrix elements file will be created'
   !read total number of directions     
            read(*,*)line
	    read(*,*)ndir
           allocate(idir(ndir),stat=ok)
	    if(ok/=0)then
              print*,'allocation failed for idip'
              stop
            end if	    
   
            read(*,*)(idir(i),i=1,ndir)

        
            goto 2
        
        else
        
            idipint=0

            goto 2
        
        end if
    
    elseif(method(1:5) == 'PRORB')then
        iprorb=1
	print*,'Orbitals will be printed in the output file'
        goto 2
    
    
    ! Whether to perform calculations in the presence of finite electric field
    
    elseif(method(1:6) == 'EFIELD')then
        iefield=1
	print*,'Scf calculations will be performed in the presence of electric field'
	allocate(efield(3),stat=ok)
	    if(ok/=0)then
              print*,'allocation failed for efield'
              stop
            end if
        read(*,*)(efield(i),i=1,3)
        print*,' '
        print*,'Calculations in the Presence of Finite Electric Field'
        write(*,191)efield(1),efield(2),efield(3)
        191 format(1x,'Ex, Ey, Ez',3x,f7.4,1x,f7.4,1x,f7.4)
        print*,' '
        goto 2
	
	
    elseif(method(1:3) == 'END' .OR. method(1:4) == 'ENDM' .OR. &
        method(1:5) == 'ENDME' .OR. method(1:6) == 'ENDMET' .OR. &
        method(1:6) == 'ENDMTH' .OR. method(1:6) == 'ENDMTD')then

        goto 3
    else
        print*,'The method of your choice:',method
        print*,'is not available presently'
        stop
    end if
    
    
   3 continue
 


!Skip all this if band structure calculations for Hueckel model is needed

  nunit=1
        
      allocate (tvec(3),stat=ok)
      if(ok/=0)then
        print*,'NOT ENOUGH MEMORY IS AVAILABLE FOR ALLOCATING tvec'
	stop
      end if 
    
      tvec=0.d0
      
  if (ibands > 0) go to 420
  
! Read the total number of unit cells
    
    nunit=0
    if(ibands==0)then 
      read(*,*)line
      read*,nunit
      write(*,6)nunit
      6 format(1x,'Total no. of unit cells:',1x,i3.1)
    end if
    
    
! Read the translation vector if there are more than one units
    ! For now only one dimensional lattices
     
       if(nunit > 1)then
        read(*,*)line
        read(*,*,err=222)(tvec(i),i=1,3)
	  goto 223
          222 print*,'Tvec','Error reading translational vectors'
	  stop
         223 continue	

        end if
	
    420 continue

! Check if any atomic operations (such as atom deletion, origin shift) need to
! be performed

    read(*,*)line
    read*,card
    card=adjustl(card)
    ndelatm=0
    
      if(card(1:7) == 'DELATOM')then
 
        read(*,*)ndelatm
	print*,'ndelatm=',ndelatm
	
	allocate (idelatm(ndelatm),stat=ok)
         if(ok/=0)then
           print*,'NOT ENOUGH MEMORY IS AVAILABLE FOR ALLOCATING idelat'
	   stop
         end if 
	read(*,*)(idelatm(i),i=1,ndelatm)
	
	
        
      end if

! Compute the total # of sites in the system 
    
        nsite=natom*nunit-ndelatm
    
        write(*,7)nsite
        7 format(1x,'Total no. of sites:',1x,i3.1)


! Compute the total # of electrons in the system

    nelec=nsite-icharge
      write(*,8)nelec
        8 format(1x,'Total no. of electrons:',1x,i3.1)

  !
  ! Compute the number of occupied orbitals
  ! norb: is the number of occupied orbitals for the RHF method
  ! nalpha/nbeta are the number of up/down spin orbitals for the UHF method
  !
 if(irhf==1)then
      if(mod(nelec,2)>0)then
        print*,'The system has odd number of electrons=',nelec
        print*,'So your choice of RHF method is not possible'
        stop
      end if
        norb=nelec/2
  elseif(iuhf==1)then
 !
      norb=nalpha+nbeta
      if(norb/=nelec)then
        !print*,'nalpha,nbeta',nalpha,nbeta
        print*,'For the UHF method, you gave nalpha, nbeta=',nalpha,nbeta
        print*,'Error: nalpha,nbeta donot add up to nelectron=',nelec
        stop
      end if
   else
      print*,"Neither RHF not UHF is given as input"
  end if

 if (ibands > 0) go to 421
! Read the convergence threshold on the total energy

    conv=1000.d0
    read(*,*)line
    read*,conv
    write(*,9)conv
  9 format(1x,'Convergence Threshold:',1x,f15.10)

! Read the maximum number of iterations

    maxitr=101
    read(*,*)line
    read*,maxitr
    write(*,10)maxitr
 10 format(1x,'Maximum number of iterations:',1x,i6.1)

! See if the Fock matrix needs to be damped

    idamp=0
    xdamp=1.d0
    read(*,*)line

    read*,card
    card=adjustl(card)
    if(card(1:4) == 'DAMP')then
        read*,xdamp
        write(*,152)xdamp
        idamp=1
        152 format(1x,'Damping on: Coeff of New Fock matrix',1x,f8.4)
	

    end if

 421 continue

    goto 121
  
! Fatal errors

    2001 print*,'Input','Error Reading the Title'
    2002 print*,'Input','Error Reading the File Extension'
    2003 print*,'Input','Error Reading the Method'
    2004 print*,'Input','Error Reading the Charge'
    2005 print*,'Input','Error Reading the no. of atoms in a unit cell'
    2006 print*,'Input','Error Reading the no. of unit cells'
    2007 print*,'Input','Error Reading the convergence threshold'
    2008 print*,'Input','Error Reading the maximum no. of iterations'
    2009 print*,'Input','Error Reading the Hamiltonian'
    2010 print*,'Input','Error Reading the PPP parameterization' 
    2011 print*,'Input','Error Reading Hamiltonian parameters' 
    2012 print*,'Input','Error Reading Bipartite lattice condition' 
    2013 print*,'Input','Error reading atomic operations card' 
    2014 print*,'Input','Error reading SCI related ORBDEL card'

    121 continue

! ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    print*,'done with input'
!
  nsitex=natom*nunit
!
   allocate( rr_site(3,nsitex),stat=ok)
      if(ok/=0)then
        print*,'NOT ENOUGH MEMORY IS AVAILABLE FOR ALLOCATING rr_site'
	stop
      end if   
!
    call r_site(rr_site,rr_atom,tvec,nsitex,natom,nunit)
!  
  allocate( swap(3,nsitex-ndelatm),stat=ok)
      if(ok/=0)then
        print*,'NOT ENOUGH MEMORY IS AVAILABLE FOR ALLOCATING swap'
	stop
      end if  
!
    if(ndelatm > 0)call delatm(rr_site,idelatm,swap, &
	    ndelatm,nsitex)
!
! Symmetrize the coordiates to place the origin at the center of
! mass

    call symsite(rr_site,nsite)
 
! Print out the coordinates of all the sites

    call printr(rr_site,nsite)
!
! itij--> Address for hopping elements
! ihij--> Address for complete one-electron matrix (hopping+on-site)
! ivij--> Address for the electron repulsion matrix elements
!
   allocate (tij((nsite*(nsite+1))/2),stat=ok)
      if(ok/=0)then
        print*,'NOT ENOUGH MEMORY IS AVAILABLE FOR ALLOCATING tij'
	stop
      end if
  allocate (hij((nsite*(nsite+1))/2),stat=ok)
      if(ok/=0)then
        print*,'NOT ENOUGH MEMORY IS AVAILABLE FOR ALLOCATING hij'
	stop
      end if 
  allocate (ihij((nsite*(nsite+1))/2),stat=ok)
      if(ok/=0)then
        print*,'NOT ENOUGH MEMORY IS AVAILABLE FOR ALLOCATING ihij'
	stop
      end if 
  allocate (vij((nsite*(nsite+1))/2),stat=ok)
      if(ok/=0)then
        print*,'NOT ENOUGH MEMORY IS AVAILABLE FOR ALLOCATING vij'
	stop
      end if 
  allocate (eps(nsite),stat=ok)
      if(ok/=0)then
        print*,'NOT ENOUGH MEMORY IS AVAILABLE FOR ALLOCATING eps'
	stop
      end if 
   allocate (fk((natom*(natom+1))/2),stat=ok)
      if(ok/=0)then
        print*,'NOT ENOUGH MEMORY IS AVAILABLE FOR ALLOCATING fk'
	stop
      end if               
    nij=(nsite*(nsite+1))/2

! Compute the one- and two-electron matrix elements in the AO representation.
    
    call matel(iham,ibands,enuc,tij,hij,ihij,vij,&
    rr_site,nunit,natom,nsite,icorfnc,xcor,natmi,natmj,natmic,natmjc, &
    itype,jtype,iefield,efield,ihampar,ichngcl,itypec,&
             jtypec,r0new,unew,dielnew,r0,u,vv,diel,ifcorc)    
!
! If the band structure within the Hueckel model is required, do it here
! 
    if(ibands > 0 .AND. iham == 4)then
   
        call bands(tij,hij,ihij,natom,fkmax,fkmin,fkstp)
    
        stop
    
    end if

    if(irhf == 0.AND.iuhf == 0)goto 1
 
    allocate (orb(nsite,nsite),stat=ok)
      if(ok/=0)then
        print*,'NOT ENOUGH MEMORY IS AVAILABLE FOR ALLOCATING orb'
        stop
      end if
      
    allocate (orb_alpha(nsite,nsite),orb_beta(nsite,nsite),stat=ok)
      if(ok/=0)then
        print*,'NOT ENOUGH MEMORY IS AVAILABLE FOR ALLOCATING orb_alpha or orb_beta'
        stop
      end if
      
    allocate (eval(nsite),stat=ok)
      if(ok/=0)then
        print*,'NOT ENOUGH MEMORY IS AVAILABLE FOR ALLOCATING eval'
        stop
      end if      
 
    allocate (eval_alpha(nsite),eval_beta(nsite),stat=ok)
      if(ok/=0)then
        print*,'NOT ENOUGH MEMORY IS AVAILABLE FOR ALLOCATING eval_alpha,eval_beta'
        stop
      end if
      
    write(*,*)'Entering routine SCF'
!
    if(irhf==1)then
      call scf_rhf(tij,hij,vij,nsite,idamp,xdamp,maxitr,norb,orb,iprorb,&
      conv,eval,enuc,escf,iham)
    elseif(iuhf==1)then
      call scf_uhf(tij,hij,vij,nsite,idamp,xdamp,maxitr,iprorb,conv,nalpha,&
                   nbeta,orb_alpha,orb_beta,eval_alpha,eval_beta,enuc)
    end if
! 
! norb_o: # of occupied orbitals
! norb_v: # of virtual orbitals
!
    write(*,*)'Done with SCF calculations'

    1 continue

    if(iorbdn > 0)then
!
      if(irhf==1)then
        call orbden(orb,iosite,nosite,nsite,nsite, &
        natom,nunit,idenfl)
      elseif(iuhf==1)then
       print*,'orbital density analysis for alpha spin'
       print*,' '
       call orbden(orb_alpha,iosite,nosite,nsite,nsite, &
        natom,nunit,idenfl)
       print*,'orbital density analysis for beta spin'
       print*,' '	
       call orbden(orb_beta,iosite,nosite,nsite,nsite, &
        natom,nunit,idenfl)	    
      end if
   end if

    
    if(irhf==1.AND.isci > 0)then
    
        write(*,*)' '
        write(*,*)'Entering routine SCI'
        write(*,*)' '
	
      norb_o=norb
      norb_v=nsite-norb_o
      ntorb=norb_o+norb_v
      ndimh=norb_o*norb_v+1
      
      allocate(e_ci(ndimh),ci_vec(ndimh,ndimh),stat=ok)
       if(ok/=0)then
        print*,'NOT ENOUGH MEMORY IS AVAILABLE FOR ALLOCATING e_ci or ci_vec'
 	stop
      end if 

      call sci(vij,hij,orb,eval,e_ci,ci_vec, &
             nsite,ntorb,norb_o,norb_v,ndimh,escf)

        write(*,*)'Done with SCI calculations'
    
    end if

! Optical absorption, if desired

    if(irhf==1.AND.ioptics > 0)then
    
      norb_o=norb
      norb_v=nsite-norb_o
      ntorb=norb_o+norb_v
    
        if(isci > 0)then
            ndimh=norb_v*norb_o+1
            ndimd=(ntorb*(ntorb+1))/2
        else
            ndimh=0
            if(nabsorb == 1.AND.nexcite(1) == 0)then
                ndimd=norb_v*norb_o
            else
                ndimd=(ntorb*(ntorb+1))/2
            end if
        end if
    
        write(*,*)' '
        write(*,*)'Entering routine OPTICS'
        write(*,*)' '
 	
        call optics(rr_site,orb,e_ci,eval,ci_vec, &
        nsite,norb_o,norb_v,ntorb,ndimh,ndimd,&
	isci,nabsorb,nexcite,iabsorb,ndim,idim,&
	omega1,omega2,domega,gamma,scale,ihole,ipart)    
    end if

! If desired, generate the output file for the data relevant for
! single-particle NLO calculations

    if(irhf==1.AND.inlo > 0)then
    
        nocc=nsite/2
        nrem=nsite-2*nocc
        nocc=nocc+nrem
        morb=nsite
        nbas=nsite
	
        if(inldel > 0)call delo_nlo(orb,eval,nbas, &
        morb,nocc,ndelo1,ndelo2,ndelv1,ndelv2)
      
        nvrt=nsite-nocc
!
        ifile=10
	
        call nlo(rr_site,orb,eval(1),eval(nocc+1), &
        nbas,nocc,nvrt,morb,ndip,idip,filenlo,ifile)    
    end if

! Output file generation for CI calculations, if desired.

    if(iciprp > 0)then        
!
        iorbu=11
      open(unit=iorbu,file='eval_orb.dat',form='unformatted',status='old')
	rewind(iorbu)

        if(irhf==1)then
	  read(iorbu)msite,morb
	  if(msite==nsite.AND.morb==nsite)then
           read(iorbu)(eval(i),i=1,msite)
           read(iorbu)((orb(j,i),j=1,msite),i=1,msite)
	  else
	   print*,'Error in ppp: msite/=nsite, stoping before ci_drv',msite,nsite
	   stop
	  end if
!	  
	elseif(iuhf==1)then
          read(iorbu)msite,morb_a,morb_b
	  if(msite==nsite.AND.morb_a==nsite.AND.morb_b==nsite)then
            read(iorbu)(eval(i),i=1,msite)
            read(iorbu)((orb(j,i),j=1,msite),i=1,msite)	   	  
	    if(ialpbet==2)then
	      read(iorbu)(eval(i),i=1,msite)
              read(iorbu)((orb(j,i),j=1,msite),i=1,msite)
	    end if
	  else
	     print*,'Error in ppp: msite/=nsite, stoping before ci_drv',msite,nsite
	     stop
	  end if
	 end if
	 close(unit=iorbu,status='keep')
!	
    call ci_drv(eval,orb,rr_site,hij,vij,nsite,nsite,&
       norbfr,nunit,nfreez,ifreez,idipint,ndelete,&
       idelete,enuc,idir,ndir,ifrez)
     
    end if

    end program

! ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

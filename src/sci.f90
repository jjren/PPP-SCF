! ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! This is the master routine meant for carrying out the single-CI
! using the HF occupied and virtual orbitals.

    subroutine sci(v,h,orb,eval,e_ci,ci_vec,&
    msite,morb,mocc_orb,mvirt_orb,ndimh,escf)

! cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

IMPLICIT NONE


    integer::mm,i,j,nmaxp,ok,info,norb_o
    integer,intent(in)::morb,msite,mocc_orb,mvirt_orb,ndimh
    real(kind=8)::x
    real(kind=8),intent(in)::v((msite*(msite+1))/2),h((msite*(msite+1))/2),&
                          eval(*),orb(msite,*),escf
    real(kind=8),intent(out)::e_ci(ndimh),ci_vec(ndimh,ndimh)
    
    real(kind=8),allocatable::ham(:),work(:)
     
    allocate(ham((ndimh*(ndimh+1))/2),stat=ok)
       if(ok/=0)then
        print*,'NOT ENOUGH MEMORY IS AVAILABLE FOR ALLOCATING ham'
 	stop
      end if 
          
   
      allocate(work(8*ndimh),stat=ok)
       if(ok/=0)then
        print*,"Allocation for work fails"
        stop
       end if
       
       work=0.d0  
       e_ci=0.d0
       ci_vec=0.d0
         
    write(*,31)morb
    31 format(1x,'SR sci: total no. of orbitals:',1x,i3.1)


! Now construct the single-CI Hamiltonian matrix

    write(*,*)' '
    write(*,*)'Starting SCI matrix construction'
    call ham_sci(orb,v,ham,eval,morb,mocc_orb,mvirt_orb,msite,ndimh)
    write(*,*)' '

! Diagonalize the single-CI matrix

    write(*,*)'Starting the diagonalization of the SCI matrix'
    write(*,*)' '
    
   call DSPEV('V','U',ndimh,ham,e_ci,ci_vec,ndimh,work,info)	     

!
       if(info/=0)then
        print*,"SCI: Hamiltonian diagonalization failed by DSPEV"
        print*,"info=",info
        stop
       end if
     !
! Labeling of CSFs
    norb_o=mocc_orb
    mm=0
    do i = 1,norb_o
        do j = norb_o+1,morb
            mm = mm+1
            write(*,1004)i,j,mm
        end do
    end do

    1004	format(' i = ',i3,1x,'j= ',i3,'mm = ',i4)


    write(*,51)escf
    51 format(1x,'E_SCF = ',f16.8)

! Print out the CI eigenvalues

    write(*,*)'CI Eigen values (relative to SCF energy) (eV)'
    write(*,*)' '

    nmaxp=min(ndimh,600)
    do i=2,nmaxp
        write(*,*)' '
        write(*,71)i-1,e_ci(i)
        write(*,*)' '
        do j=1,ndimh
            x=ci_vec(j,i)
            if(dabs(x) > 0.1d0)write(*,72)j-1,x
        end do
    end do

    71 format('Ex. State #: ',1x,i3.1,3x,'Energy ',1x,f12.6)
    72 format(' ',i7.1,2x,f8.4)



    deallocate(ham,work)
    return
    end subroutine sci

! ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

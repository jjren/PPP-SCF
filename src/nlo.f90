! ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    subroutine nlo(r_site,orb,eocc,evrt,nbas,nocc, &
    nvrt,norb,ndip,idip,filen,ifile)

! This routine generates the input file to be used for the single
! particle NLO calculations done by the program chi3.f
! cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


IMPLICIT NONE
  
    integer,intent(in)::norb,ndip,ifile,nocc,nvrt,nbas,idip(ndip)
    real(kind=8),intent(in)::r_site(3,nbas),orb(nbas,norb),&
                         eocc(nocc),evrt(nvrt)
    integer::ok,nxx,iflag,i,j,k,ij
    real(kind=8)::tol,xx
    real(kind=8),allocatable::dipmat(:,:)
    character*(*)filen

   allocate (dipmat((norb*(norb+1))/2,ndip),stat=ok)
      if(ok/=0)then
        print*,'NOT ENOUGH MEMORY IS AVAILABLE FOR ALLOCATING dipmat'
        stop
      end if 
      
    if(norb /= (nocc+nvrt))then
        write(*,*)'sr nlo: norb.ne.(nocc+nvrt)'
        write(*,*)'norb,nocc,nvrt',norb,nocc,nvrt
        stop
    end if

    open(file=filen,unit=ifile,status='unknown',access='sequential', &
    err=1)
    goto 2

    1 write(*,*)'sr nlo: error opening output file'
    write(*,*)'unit,file_name',ifile,filen
    stop

    2 continue

    write(ifile,*)'Chi3 calculation'
    write(ifile,*)'Specify the type of calculation (THG,TPA,...)'
    write(ifile,*)'No. of components of chi3 to be computed'
    write(ifile,*)'Write the total # of components here'
    write(ifile,*)'Components of Chi3 tensor to be computed'
    write(ifile,*)'Specify these here (1 1 1 1 etc.)'
    write(ifile,*)'It is an independent electron calculation'
    write(ifile,*)'ONE'
    write(ifile,*)'Total # of orbitals, # of occupied orbitals'
    write(ifile,*)norb,nocc

    write(ifile,*)'line width, lower frequency limit, upper frequency &
    limit, freq. step'
    write(ifile,*)'insert the above mentioned data here'

    write(ifile,*)'Energies of the occupied states'
    write(ifile,*)(eocc(i),i=1,nocc)

    write(ifile,*)'Energies of the virtual states'
    write(ifile,*)(evrt(i),i=1,nvrt)

! Calculate the dipole moments in the format required by the chi3
! program

!    call prblk(orb,nbas,nbas,norb,0,0,'bas','orb',6)



    call dipcal(r_site,orb,dipmat,idip,nbas,norb,ndip)

! Write out the dipole moments in the proper format


    write(ifile,*)'Here are the dipole moments'

    nxx=0

    tol=1.d-8
    ij=0

    do i=1,norb
        do j=1,i
            ij=ij+1
            iflag=0

            do k=1,ndip
                xx=abs(dipmat(ij,k))
                if(xx > tol)iflag=iflag+1
            end do
            if(iflag > 0)nxx=nxx+1
        end do
    end do

    write(ifile,*)nxx

    ij=0
    do i=1,norb
        do j=1,i
            ij=ij+1
            iflag=0
            do k=1,ndip
                xx=abs(dipmat(ij,k))
                if(xx > tol)iflag=iflag+1
            end do
            if(iflag > 0)write(ifile,*)i,j,(dipmat(ij,k),k=1,ndip)
        end do
    end do

    close(unit=ifile,status='keep')

    deallocate(dipmat)
    return
    end subroutine nlo

! ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

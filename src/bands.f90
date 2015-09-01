! cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    subroutine bands(tij,tij_rl,irij,matom,fkmax,fkmin,fkstp)

! This is the master routine that performs the band structure
! calculations by setting up the hopping matrix in the k-space.

! cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
IMPLICIT NONE

    integer,intent(in)::matom,irij((matom*(matom+1))/2) 
    integer::nkval,k,nwork,nrwork,iorbu,ibandu,info,i,j,ok,ij
    real(kind=8),intent(in)::tij((matom*(matom+1))/2),tij_rl((matom*(matom+1))/2),&
                             fkmax,fkmin,fkstp
    real(kind=8)::fkval
    complex*16,allocatable::fk(:),zk(:,:),work(:)
    character(72)::bndfil,orbfil,V,U   
    real(kind=8),allocatable::e(:,:),rwork(:)
!
  nwork=max(1,2*matom-1)
  nrwork=max(1, 3*matom-2)
  if(fkstp > 0.d0)nkval=nint((fkmax-fkmin)/fkstp)+1
!  print*,'nkval=',nkval
!    
    allocate (fk((matom*(matom+1))/2),zk(matom,matom),work(nwork),stat=ok)
      if(ok/=0)then
        print*,'NOT ENOUGH MEMORY IS AVAILABLE FOR ALLOCATING fk, zk, or work'
	stop
      end if 
    allocate (e(matom,nkval),rwork(nwork),stat=ok)
      if(ok/=0)then
        print*,'NOT ENOUGH MEMORY IS AVAILABLE FOR ALLOCATING e or rwork'
	stop
      end if

!
!     e=0.d0
! Open the file for the band structure data

    bndfil='bands.dat'
    ibandu=18
    open(unit=ibandu,file=bndfil,form='formatted',status='unknown',err=21)
    goto 22

    21 write(*,*)'sr band: error opening the band data file'
    stop
    
22 continue

    orbfil='bloch_orbitals.dat'
    iorbu=19
    open(unit=iorbu,file=orbfil,form='unformatted',status='unknown',err=23)
    goto 24

    23 write(*,*)'sr band: error opening the Bloech orbital output data file'
    stop
    
24 continue    

! Loop over the K vectors and get the band structure

    fkval=fkmin-fkstp
    do  k=1,nkval
    
        fkval=fkval+fkstp
    
    ! Fourier transform of the Hopping matrix

        call foutra_tb(fkval,tij,tij_rl,irij,fk,matom)
    
    
    ! Diagonalize the hopping

         call zhpev('V','U',matom,fk,e(1,k),zk,matom,work,rwork,info) 
	  
        if(info/= 0)then
        
            write(*,*)'k,info=',k,info
           print*,'Bands','Error diagonalizing hopping'
        
        end if
 
!    call prblk(e,matom,matom,matom,0,0,'site','site')
 
      write(iorbu)((zk(i,j),i=1,matom),j=1,matom)
 
    
    end do

! Print out the bands

    do  i=1,matom
        write(ibandu,*)'          '
        write(ibandu,*)'# BAND N0.',i
        fkval=fkmin-fkstp
        do  k=1,nkval
            fkval=fkval+fkstp
            write(ibandu,*)fkval,e(i,k)
         end do
    end do
 
   close(unit=iorbu,status='KEEP')   
   close(unit=ibandu,status='KEEP')
    deallocate(fk,zk,e,work,rwork)
    return
    end subroutine bands

! ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc




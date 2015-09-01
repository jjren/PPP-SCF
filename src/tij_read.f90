! cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! This subroutine reads the unique hopping matrix elements and their
! connectivities from the input file. All the reads are in free format.

    subroutine tij_read(tij,nsite,nunit,natom,ihopgen)

! cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
IMPLICIT NONE

!    implicit real*8(a-h,o-z)
!    implicit integer*4(i-n)

    integer::nij,ntij,nxij,itij,ixij,i,j,iunit,ii,jj,kl,ij,ijpk
    integer,intent(in)::nunit,natom,nsite
    integer,intent(out)::ihopgen
    real(kind=8),intent(out)::tij((nsite*(nsite+1))/2)
    real(kind=8)::x
    character(len=80)::card,line

!c

    ihopgen=0

    nij=(nsite*(nsite+1))/2

    read(*,*,err=101)line
    read(*,*)card
    card=adjustl(card)
    if(card(1:6) == 'HOPGEN')then
        ihopgen=1
        return
    end if

    read(*,*,err=101)ntij

    write(*,*)' '
    write(*,1)ntij
    1 format(1x,'Total no. of Unique hopping elements read=',1x,i6.1)
    write(*,*)' '

    do itij=1,ntij
          ! jjren read bondlink hopping integral from input file like DMRG-X
          read(*,*,err=103) i,j,x
          if(i <= nsite .AND. j <= nsite)then
   
              if(i >= j)then
                  ij=ijpk(i,j)
                  tij(ij)=x
  
              ! If more than one unit cell is there, generate rest of the hoppings
              ! here
              
                  do iunit=2,nunit
                      ii=i+(iunit-1)*natom
                      jj=j+(iunit-1)*natom
                      if(ii <= nsite .AND. jj <= nsite)then
                          kl=ijpk(ii,jj)
                          tij(kl)=x
                      end if

                  end do
                 
              else

                  write(*,*)'For itij,ixij=',itij,ixij
                  write(*,*)'i,j=',i,j
                  print*,'Error:Tij_read', &
                  'Hopping connectivity in wrong order'
              end if
          
          else
          
              write(*,*)'i,j',i,j
              print*,'Error:Tij_Read','i and/or j value(s) out of bounds'
          
          end if
    
    end do

!   call plblk(tij,nsite,0,'site',6)


! Successful, get out of here

    goto 105

    101 print*,'Tij_Read','Error reading ntij'

    102 write(*,*)'itij=',itij
    print*,'Error:Tij_Read','Error reading hopping element and nxij'
    103 write(*,*)'itij,ixij',itij,ixij
    print*,'Tij_Read','Error reading hopping connectivity'


    105 continue

 
! call PLBLK(tij,nsite,0,'fat',6)
    return
    end subroutine tij_read

! ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

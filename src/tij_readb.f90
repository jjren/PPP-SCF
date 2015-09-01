! cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! For the purpose of band structure calculationsThis subroutine reads the unique hopping matrix elements used in
! and their connectivities from the input file.
! All the reads are in free format.

    subroutine tij_readb(tij,tij_rl,irij,natom)

! tij:  intra unit cell hopping
! tij_rl: inter unit cell hopping
! irij: unit cells involved in the inter-cell hoppings.
! cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
IMPLICIT NONE

!    implicit real*8(a-h,o-z)
!    implicit integer*4(i-n)

    integer::nij,ntij1,ntij2,itij,ixij,nxij,i,j,k,ij,ijpk
    integer,intent(out)::irij((natom*(natom+1))/2)
    real(kind=8)::x
    real(kind=8),intent(out)::tij((natom*(natom+1))/2),tij_rl((natom*(natom+1))/2)
    integer,intent(in)::natom
    character(80)::line

    nij=(natom*(natom+1))/2

    do ij=1,nij
    
        tij(ij)=0.d0
        tij_rl(ij)=0.d0
    
    end do
    read(*,*)line
    read(*,*,err=101)ntij1,ntij2

    write(*,*)' '
    write(*,1)ntij1
    1 format(1x, &
    'No. of Unique intracell hopping elements read=',1x,i6.1)
    write(*,2)ntij2
    2 format(1x, &
    'No. of Unique intercell hopping elements read=',1x,i6.1)
    write(*,*)' '

! First read the intracell hoppings

    do itij=1,ntij1
        read(*,*,err=102)x,nxij
        do ixij=1,nxij
            read(*,*,err=103)i,j
            if(i <= natom .AND. j <= natom)then
            
                if(i >= j)then
                    ij=ijpk(i,j)
                    tij(ij)=x
                
                else
                    print*,'For itij,ixij=',itij,ixij
                    print*,'i,j=',i,j
                    print*,'Error:Tij_readb', &
                    'Hopping connectivity in wrong order for tij'
                end if
            
            else
            
                print*,'i,j',i,j
                print*,'Error:Tij_Readb','i and/or j value(s) out of bounds'
            
            end if
        
        end do
    
    end do

! Now read the intercell hoppings. Note that a linear tight-binding
! system with only nearest-neighbor coupling is assumed. Thus a
! given cell can only connect to at the most to two neighbors
! (one to the left and the other to the right). For anything else
! the code will have to be
! modified. Also one has to specify the unit cell to which it is
! connected. It is assumed to be in the format t_i(R_i)j(0)

    do itij=1,ntij2
        read(*,*,err=102)x,nxij
        do ixij=1,nxij
            read(*,*,err=103)i,j,k
            if(i <= natom .AND. j <= natom)then
            
                if(i >= j)then
                    ij=ijpk(i,j)
                    tij_rl(ij)=x
                    irij(ij)=k
                else
                
                    print*,'For itij,ixij=',itij,ixij
                    print*,'i,j=',i,j
                    print*,'Error:Tij_readb', &
                    'Hopping connectivity in wrong order for tij_rl'
                
                end if
            
            else
            
                print*,'i,j',i,j
                print*,'Error:Tij_Readb','i and/or j value(s) out of bounds'
            
            end if
        
        end do
    
    end do

! Successful, get out of here

    goto 104

    101 print*,'Tij_Read','Error reading ntij'

    102 write(*,*)'itij=',itij
    print*,'Error:Tij_Read','Error reading hopping element and nxij'
    103 write(*,*)'itij,ixij',itij,ixij
    print*,'Error:Tij_Read','Error reading hopping connectivity'
!c
    104 continue

    return
    end subroutine tij_readb

! ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


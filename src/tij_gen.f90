      subroutine tij_gen(tij,nsite,rr_site)

 implicit none
!
 integer,intent(in)::nsite
 integer::nbond,ok,i,j,ij,k,ii,jj
 real(kind=8),intent(out)::tij((nsite*(nsite+1))/2)
 real(kind=8),intent(in)::rr_site(3,nsite)
 real(kind=8),allocatable::bond_len(:),t_val(:)
 real(kind=8)::x1,y1,z1,dist,tol
!
 tol=0.01d0
 read(*,*)nbond
 print*,nbond
 allocate(bond_len(nbond),t_val(nbond),stat=ok)
 if(ok/=0)then
   print*,'sr tij_gen: error allocating bond_len, t_val'
   print*,'stopping the program'
   stop
 end if
 read(*,*)(bond_len(i),i=1,nbond)
 read(*,*)(t_val(i),i=1,nbond)
!
 tij=0.d0
 do i=1,nsite
   x1=rr_site(1,i)
   y1=rr_site(2,i)
   z1=rr_site(3,i)
   do j=1,i-1
      ij=(i*(i-1))/2+j
      dist=(x1-rr_site(1,j))**2+(y1-rr_site(2,j))**2+(z1-rr_site(3,j))**2
      dist=dsqrt(dist)
      do k=1,nbond
        if(dabs(dist-bond_len(k))<tol)tij(ij)=t_val(k)
      end do
   end do
 end do
!
 deallocate(bond_len,t_val)
 end  
     

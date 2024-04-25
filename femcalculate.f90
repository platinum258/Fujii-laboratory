program calculate
implicit none
! integer e1,e2!,nnode,node,mxb,mxne,ne,ndp,nepb,iof,iel
! real(8) mnx,mxx,mny,mxy
integer,parameter :: nd=3, mxe=4800, mxn=2501
! integer :: en0(mxe),en1(mxe),en2(mxe),en3(mxe),nodex(mxe,nd),eid(mxe)
real(8) :: xcoord(mxn),ycoord(mxn),rhs(mxe) 
real(8) :: nec(mxe,2)
integer i,j,iw,e1,e2
character*14 outfile

! write(*,*)'femout.txt solver'

!   IR=1
!   open(IR,file='MData.txt')
!   read(IR,*)NNODE
!   if(NNODE .gt. MXN) stop 'NNODE>MXN'
!   read(IR,*)NE
!   if(NE .gt. MXE) stop 'NE>MXE'
!   allocate(nec(ne,2))
!   read(IR,*)MXNE
!   read(IR,*)MXB
!   read(IR,*)NDP
!   read(IR,*)NEPB
!   read(IR,*)IOF
!   read(IR,*)MNX
!   read(IR,*)MXX
!   read(IR,*)MNY
!   read(IR,*)MXY
!   read(IR,*) (node,XCOORD(node),YCOORD(node),j=1,nnode)
!   read(ir,*) (iel,en0(iel),en1(iel),en2(iel),en3(iel),i=1,ne)
!   close(ir)
! write(*,*)'input over'  



   e1=0
   e2=0
   do i=0,60
     do j=1,40
       e2=i*80+j*2
       e1=e2-1
       nec(e2,1)=dble(i)-29.25d0
       nec(e2,2)=dble(j)-20.25d0
       nec(e1,1)=dble(i)-29.75d0
       nec(e1,2)=dble(j)-20.75d0
     end do
   end do
  
write(*,*) 'step1 complete'

do i=1,mxe
rhs(i)=0.5d0+nec(i,1)/60.0d0
end do

write(*,*)'calculated'
outfile='Treferance.txt'
iw=1
open(iw,file=outfile,status='new')
do j=1,mxe
write(iw,*) j,rhs(j)
end do
close(iw)
write(*,*)'output over'
end



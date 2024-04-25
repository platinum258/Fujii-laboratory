program fem4_2
implicit none
integer Number_Node,iel,iof,e1,e2,n1,n2,n3,num_threads
integer Number_EDesignArea, Number_SEDesignArea, Number_Element,Number_SElement,Number_Eins, Number_SEins
integer INFO
integer ID_Element_Objective_Function_1
integer,parameter :: nd=3,exx=1, eyy=1, exy=0, SE1=37616
integer :: rkr_id(SE1)
real(8) mnx,mxx,mny,mxy
real(8) :: x(nd),y(nd),a(nd),b(nd),sk(nd,nd)
integer,allocatable :: eid(:),nodex(:,:),nodex_s(:,:)
real(8),allocatable :: xcoord(:),ycoord(:),rhs(:),rk(:,:)
real(8),allocatable :: shape_function(:,:),nec(:,:),T_e(:)
real(8) :: rkr(SE1,SE1),rhsT(SE1),ipmkl(SE1)
real(8) :: rko(SE1,SE1),rko_id(SE1),rhso(SE1)
character*12 inpfile
character*11 outfile
integer,parameter :: ne=76800

integer ir,iw,i,j,j1,j2,k,l,l1,l2,n,node
real(8)m,p,p1,p2,e,q
      ir=1
 open(IR,file='MdataC_2.txt')
 read(IR,*)Number_Node
 allocate(xcoord(Number_Node),ycoord(Number_Node),rhs(Number_Node),rk(Number_Node,Number_Node))
 !read(IR,*)mepn
!  read(IR,*)npb
!  read(IR,*)npec
!  read(IR,*)epb
 read(IR,*)Number_Element
 allocate(nodex(Number_Element,3),eid(Number_Element),shape_function(Number_Element,3))
 allocate(nec(Number_Element,2),T_e(Number_Element))
 read(IR,*)Number_SElement
 allocate(nodex_s(Number_SElement,4))
 read(IR,*)Number_EDesignArea
 read(IR,*)Number_SEDesignArea
 read(IR,*)Number_Eins
 read(IR,*)Number_SEins
 read(IR,*)ID_Element_Objective_Function_1
 read(IR,*)MNX
 read(IR,*)MXX
 read(IR,*)MNY
 read(IR,*)MXY
 read(IR,*) (node,XCOORD(node),YCOORD(node),j=1,Number_Node)
 read(ir,*) (iel,eid(iel),(nodex(iel,j),j=1,nd),i=1,Number_Element)
 read(ir,*) (iel,eid(iel*2),(nodex_s(iel,k),k=1,4),i=1,Number_SElement)
 close(ir)
 write(*,*)'input over'

   e1=0
   e2=0
   do i=0,239
     do j=1,160
       e2=i*320+j*2
       e1=e2-1
       nec(e2,1)=dble(i)-119.25d0
       nec(e2,2)=dble(j)-80.25d0
       nec(e1,1)=dble(i)-119.75d0
       nec(e1,2)=dble(j)-80.75d0
     end do
   end do

call omp_set_num_threads(1)
call mkl_set_num_threads(1)

call gsm(nd,ne,Number_Node,eid,X,Y,A,B,sk,nodex,xcoord,ycoord,rk)
write(*,*)'local matrix calculated'

do i=1,Number_Node
rhs(i)=0.0d0
end do

do k=1,Number_Node
if(ycoord(k)==mny .or. ycoord(k)==mxy) then
rhs(k)=0.0d0
end if
end do

do k=1,Number_Node
if(xcoord(k) .eq. mnx) then
do j=1,Number_Node
rk(k,j)=0.0d0
end do
rk(k,k)=1.0d0
rhs(k)=0.0d0 
else if(xcoord(k) .eq. mxx) then
do l=1,Number_Node
rk(k,l)=0.0d0
end do
rk(k,k)=1.0d0
rhs(k)=1.0d0
end if
end do

l1=1
do i=1,Number_Node
  if(rk(i,i)/=0.0d0) then
  rko_id(l1)=i
  l1=l1+1
  end if
end do


do i=1,SE1
   do j=1,SE1
     rko(i,j)=rk(rko_id(i),rko_id(j))
   end do
  rhso(i)=rhs(rko_id(i))
end do   

    ! rkr(1:9720,1:9720)=rk(82:9801,82:9801)
    ! rhsT(1:9720)=rhs(82:9801)


call mkl_set_num_threads(num_threads)
call dgesv(SE1,1,rko,SE1,ipmkl,rhso,SE1,info)

! return rkr to rk
  !  rhs(82:9801)=rhsT(1:9720)
do i=1,SE1
rhs(rko_id(i))=rhso(i)
end do

   iw=203
open(iw,file='nodeoutput_2.txt',status='new')
do i=1,ne
if(eid(i)/=2) then
do j=1,nd
write(iw,*) xcoord(nodex(i,j)),ycoord(nodex(i,j)),rhs(nodex(i,j))
end do
end if
end do
close(iw)

call etemp(Number_Element,Number_Node,nd,eid,xcoord,ycoord,nec,x,y,a,b,n1,n2,n3,nodex,shape_function,rhs,T_e)
write(*,*) 'interpolated'

outfile='Tbare_2.txt'
iw=201
open(iw,file=outfile,status='new')
do i=1,ne
if(eid(i)/=2) then
write(iw,*) T_e(i)
end if
end do
close(iw)
write(*,*)'output over'


iw=202
open(iw,file='bareoutput_2.txt',status='new')
do i=1,ne
if(eid(i)/=2) then
write(iw,*) nec(i,1),nec(i,2),T_e(i)
end if
end do
close(iw)

deallocate(xcoord,ycoord,rhs,rk)
deallocate(nodex,eid,shape_function)
deallocate(nec,T_e)
deallocate(nodex_s)

end

SUBROUTINE GSM (ND,NE,Number_Node,eid,x,y,A,B,SK,NODEX,XCOORD,YCOORD,RK)
  IMPLICIT none
  integer nd,ne,Number_Node
  integer :: nodex(ne,nd),eid(ne)
  real(8) det
  real(8) :: x(nd),y(nd),a(nd),b(nd), sk(nd,nd), xcoord(Number_Node), ycoord(Number_Node), rk(Number_Node,Number_Node)
  integer i,j,k,l,ie,n

  DO I = 1 , Number_Node
  DO J = 1 , Number_Node
  RK ( I , J ) = 0.D0
  END DO
  END DO
  do ie=1,ne
  n=eid(ie)
  if(n/=2) then
  do i=1,nd
  x(i)=xcoord(nodex(ie,i))
  y(i)=ycoord(nodex(ie,i))
  end do
  B(1)=y(2)-y(3)
  B(2)=y(3)-y(1)
  B(3)=y(1)-y(2)
  A(1)=x(3)-x(2)
  A(2)=x(1)-x(3)
  A(3)=x(2)-x(1)
  det=A(3)*B(2)-B(3)*A(2)
  DO I = 1 , ND
  B(I)=B(I)/det
  A(I)=A(I)/det
  END DO
  DO I=1,ND
  DO J=1,ND
  SK(I,J)=(det/2.0d0)*(B(i)*B(j)+A(i)*A(j))
  END DO
  END DO
  DO K=1,ND
  I=NODEX(IE,K)
  DO L=1,ND
  J=NODEX(IE,L)
  RK(I,J)=RK(I,J)+SK(K,L)
  END DO
  END DO
  end if
  end do
  return
  end

  subroutine etemp(ne,Number_Node,nd,eid,xcoord,ycoord,nec,x,y,a,b,n1,n2,n3,nodex,shape_function,rhs,T_e)  
  implicit none
  integer i,j,n,ne,Number_Node,ND,n1,n2,n3
  integer :: nodex(ne,nd),eid(ne)
  real(8) det
  real(8) :: x(nd),y(nd),a(nd),b(nd), sk(nd,nd), xcoord(Number_Node), ycoord(Number_Node),nec(ne,2),rhs(Number_Node)
  real(8) :: shape_function(ne,3),T_e(ne)
  do n=1,ne
     if(eid(n)/=2) then
     do i=1,nd
       x(i)=xcoord(nodex(n,i))
       y(i)=ycoord(nodex(n,i))
     end do
      !  B(1)=y(2)-y(3)
      !  B(2)=y(3)-y(1)
      !  B(3)=y(1)-y(2)
      !  A(1)=x(3)-x(2)
      !  A(2)=x(1)-x(3)
      !  A(3)=x(2)-x(1)
       det=(x(2)-x(1))*(y(3)-y(1))-(y(1)-y(2))*(x(1)-x(3))
       shape_function(n,1)=(x(2)*y(3)-x(3)*y(2))/det+nec(n,1)*(y(2)-y(3))/det+nec(n,2)*(x(3)-x(2))/det 
       shape_function(n,2)=(x(3)*y(1)-x(1)*y(3))/det+nec(n,1)*(y(3)-y(1))/det+nec(n,2)*(x(1)-x(3))/det 
       shape_function(n,3)=(x(1)*y(2)-x(2)*y(1))/det+nec(n,1)*(y(1)-y(2))/det+nec(n,2)*(x(2)-x(1))/det 
       n1=nodex(n,1)
       n2=nodex(n,2)
       n3=nodex(n,3)
       T_e(n)=shape_function(n,1)*rhs(n1)+shape_function(n,2)*rhs(n2)+shape_function(n,3)*rhs(n3)
     end if
  end do
  end  
 

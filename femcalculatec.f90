implicit none

integer ir,Number_Element,Number_Node,Number_EDesignArea,Number_SEDesignArea,Number_SEins,Number_Eins
integer Number_SElement
integer,parameter :: nd=3
integer ID_Element_Objective_Function
real(8) MNX,MXX,MNY,MXY
integer i,j,k,l,m,n,iel,node
integer,allocatable :: eid(:),nodex(:,:),nodex_s(:,:)
real(8),allocatable :: XCOORD(:),YCOORD(:),T(:),T_e(:)
character*12 inpfile
character*16 outfile


   ir=1
   open(IR,file='MdataC_2.txt')
   read(IR,*)Number_Node
   allocate(T(Number_Node),XCOORD(Number_Node),YCOORD(Number_Node))
   !read(IR,*)mepn
  !  read(IR,*)npb
  !  read(IR,*)npec
  !  read(IR,*)epb
   read(IR,*)Number_Element
   allocate(nodex(Number_Element,3),eid(Number_Element),T_e(Number_Element))
   read(IR,*)Number_SElement
   allocate(nodex_s(Number_SElement,4))
   read(IR,*)Number_EDesignArea
   read(IR,*)Number_SEDesignArea
   read(IR,*)Number_Eins
   read(IR,*)Number_SEins
   read(IR,*)ID_Element_Objective_Function
   read(IR,*)MNX
   read(IR,*)MXX
   read(IR,*)MNY
   read(IR,*)MXY
   read(IR,*) (node,XCOORD(node),YCOORD(node),j=1,Number_Node)
   read(ir,*) (iel,eid(iel),(nodex(iel,j),j=1,nd),i=1,Number_Element)
   read(ir,*) (iel,eid(iel*2),(nodex_s(iel,k),k=1,4),m=1,Number_SElement)
   close(ir)
   write(*,*)'input over'

   do i=1,Number_Node
     T(i)=0.0d0
   end do  

   do i=0,239
     do j=1,160
       k=(i*160+J)*2
       l=k-1
       T_e(k)=0.0d0+(dble(i)+0.75d0)/240.0d0
       T_e(l)=0.0d0+(dble(i)+0.25d0)/240.0d0
     end do
   end do    

write(*,*)'calculated'
outfile='Treferance_2.txt'
ir=2
open(IR,file=outfile,status='new')
do j=1,Number_Element
write(IR,*) j,T_e(j)
end do
close(IR)
write(*,*)'output over'

deallocate(T)
deallocate(XCOORD)
deallocate(YCOORD)
deallocate(eid)
deallocate(nodex)
deallocate(nodex_s)

end program
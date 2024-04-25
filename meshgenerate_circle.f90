program generate
implicit none
integer,parameter :: nes=30000,ne=60000,nnode=30351,nd=4
real(8),parameter :: mnx=-100.0d0,mxx=100.0d0,mny=-75.0d0,mxy=75.0d0
real(8) :: xcoord(nnode),ycoord(nnode),nesc(nes,2)
integer i,iw,j,k,l,m,n,pe,x,y
integer Number_EDesignArea,Number_Eins,Number_SEDesignArea,Number_SEins
real(8) rd,rins,len,p,q
integer :: nodex_s(nes,nd),nodex(ne,3),eid(ne)
character*12 outfile

do i=0,200
do j=1,151
k=j+i*151
xcoord(k)=-100.0d0+dble(i)*1.0d0
end do
end do

do i=0,200
do j=1,151
l=j+i*151
ycoord(l)=-76.0d0+dble(j)*1.0d0
end do
end do

do i=0,199
do j=1,150
m=j+i*150

nodex_s(m,1)=j+i*151
nodex_s(m,2)=j+I*151+151
nodex_s(m,3)=j+i*151+152
nodex_s(m,4)=j+i*151+1

end do
write(*,*) "loop" ,i
end do

do i=1,nes
nodex(i*2,1)=nodex_s(i,4)
nodex(i*2,2)=nodex_s(i,2)
nodex(i*2,3)=nodex_s(i,3)
nodex(i*2-1,1)=nodex_s(i,4)
nodex(i*2-1,2)=nodex_s(i,1)
nodex(i*2-1,3)=nodex_s(i,2)
write(*,*) "loop2" ,i
end do

do i=0,199
do j=1,150
n=j+i*150
nesc(n,1)=i-99.5d0
nesc(n,2)=j-75.5d0
end do
end do

do i=1,ne
eid(i)=5
end do

do i=1,nes
p=nesc(i,1)*nesc(i,1)
q=nesc(i,2)*nesc(i,2)
len=p+q
if(len<=2500.0d0)then
eid(i*2)=4
eid(i*2-1)=4
end if
end do

p=0.0d0
q=0.0d0
len=0.0d0
do i=1,nes
p=nesc(i,1)*nesc(i,1)
q=nesc(i,2)*nesc(i,2)
len=p+q
if(len<=400.0d0)then
eid(i*2)=2
eid(i*2-1)=2
end if
end do

Number_EDesignArea=0
do i=1,ne
if(eid(i)==4)then
Number_EDesignArea=Number_EDesignArea+1
end if 
end do
Number_SEDesignArea=Number_EDesignArea/2

Number_Eins=0
do i=1,ne
if(eid(i)==2)then
Number_Eins=Number_Eins+1
end if
end do
Number_SEins=Number_Eins/2

outfile='MdataC_2.txt'
iw=1
open(iw,file=outfile,status='new')
write(iw,*) 38801
write(iw,*) 76800
write(iw,*) 38400
write(iw,*) Number_EDesignArea
write(iw,*) Number_SEDesignArea
write(iw,*) Number_Eins
write(iw,*) Number_SEins
write(iw,*) -120.0d0
write(iw,*) 120.0d0
write(iw,*) -80.0d0
write(iw,*) 80.0d0
do j=1,nnode
write(iw,*) j,xcoord(j),ycoord(j)
end do
do k=1,ne
write(iw,*) k,eid(k),nodex(k,1),nodex(k,2),nodex(k,3)
end do
do l=1,nes
write(iw,*) l,eid(l*2),nodex_s(l,1),nodex_s(l,2),nodex_s(l,3),nodex_s(l,4)
end do
close(iw)
write(*,*)'output over'
end


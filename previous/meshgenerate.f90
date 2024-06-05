program generate
implicit none
integer,parameter :: ne=19200,nnode=9801,nd=3
real(8),parameter :: mnx=-60.0d0,mxx=60.0d0,mny=-40.0d0,mxy=40.0d0,mnxd=-25.0d0,mxxd=25.0d0
real(8),parameter :: mnyd=-25.0d0,mxyd=25.0d0,mnxp=-10.0d0,mxxp=10.0d0,mnyp=-10.0d0,mxyp=10.0d0
real(8) :: xcoord(nnode),ycoord(nnode)
integer i,iw,j,k,l,m,n,p,q,pe,qe
integer :: nodex(ne,nd),eid(ne)
character*9 outfile

do i=0,120
do j=1,81
k=j+i*81
xcoord(k)=-60.0d0+dble(i)*1.0d0
end do
end do

do i=0,120
do j=1,81
l=j+i*81
ycoord(l)=-41.0d0+dble(j)*1.0d0
end do
end do

do i=0,119
do j=1,80
m=-1+j*2+i*160
n=m+1
nodex(m,1)=j+i*81+1
nodex(m,2)=j+I*81
nodex(m,3)=j+i*81+81

nodex(n,1)=nodex(m,3)
nodex(n,2)=nodex(m,3)+1
nodex(n,3)=nodex(m,1)

end do
write(*,*) "loop" ,i
end do

do i=1,ne
eid(i)=5
end do

do i=35,84
do j=16,65
pe=j*2+i*160
qe=pe-1
eid(pe)=4
eid(qe)=4
end do
end do

do i=50,69
do j=31,50
p=j*2+i*160
q=p-1
eid(p)=2
eid(q)=2
end do
end do

outfile='Mdata.txt'
iw=1
open(iw,file=outfile,status='new')
write(iw,*) -60.0d0
write(iw,*) 60.0d0
write(iw,*) -40.0d0
write(iw,*) 40.0d0
do j=1,nnode
write(iw,*) j,xcoord(j),ycoord(j)
end do
do k=1,ne
write(iw,*) k,eid(k),nodex(k,1),nodex(k,2),nodex(k,3)
end do
close(iw)
write(*,*)'output over'
end


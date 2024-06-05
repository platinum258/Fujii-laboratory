program generate
implicit none
integer,parameter :: ne=4800,nnode=2501,nd=3
real(8),parameter :: mnx=-30.0d0,mxx=30.0d0,mny=-20.0d0,mxy=20.0d0,mnxd=-12.0d0,mxxd=12.0d0
real(8),parameter :: mnyd=-12.0d0,mxyd=12.0d0,mnxp=-5.0d0,mxxp=5.0d0,mnyp=-5.0d0,mxyp=5.0d0
real(8) :: xcoord(nnode),ycoord(nnode)
integer i,iw,j,k,l,m,n,p,q,pe,qe
integer :: nodex(ne,nd),eid(ne)
character*10 outfile

do i=0,60
do j=1,41
k=j+i*41
xcoord(k)=-30.0d0+dble(i)*1.0d0
end do
end do

do i=0,60
do j=1,41
l=j+i*41
ycoord(l)=-21.0d0+dble(j)*1.0d0
end do
end do

do i=0,59
do j=1,40
m=-1+j*2+i*80
n=m+1
nodex(m,1)=j+i*41+1
nodex(m,2)=j+I*41
nodex(m,3)=j+i*41+41

nodex(n,1)=nodex(m,3)
nodex(n,2)=nodex(m,3)+1
nodex(n,3)=nodex(m,1)

end do
write(*,*) "loop" ,i
end do

do i=1,ne
eid(i)=5
end do

do i=18,41
do j=9,32
pe=j*2+i*80
qe=pe-1
eid(pe)=4
eid(qe)=4
end do
end do

do i=25,34
do j=16,25
p=j*2+i*80
q=p-1
eid(p)=2
eid(q)=2
end do
end do

outfile='Mdata1.txt'
iw=1
open(iw,file=outfile,status='new')
write(iw,*) -30.0d0
write(iw,*) 30.0d0
write(iw,*) -20.0d0
write(iw,*) 20.0d0
do j=1,nnode
write(iw,*) j,xcoord(j),ycoord(j)
end do
do k=1,ne
write(iw,*) k,eid(k),nodex(k,1),nodex(k,2),nodex(k,3)
end do
close(iw)
write(*,*)'output over'
end


implicit none
integer Number_Step
integer,parameter :: Number_Objective_Function=29200,Number_average=292
real(8),allocatable :: Objective_Function(:), Objective_Function_average(:)
integer i,j,k,m,n 
real(8) l
character filename*128

allocate(Objective_Function(Number_Objective_Function),Objective_Function_average(Number_average))

n=1
open(n,file='Objective_Function.dat')
do i=1,Number_Objective_Function
read(n,*) Number_Step,Objective_Function(i)
end do
close(n)
write(*,*)'input over'

k=0
do i=0,Number_average-1
do j=1,100
k=i*100+j
l=l+Objective_Function(k)
end do
Objective_Function_average(i+1)=l/100.0d0
l=0.0d0
end do
write(*,*) 'calculated'

n=2
open(n,file='Objective_Function_average.dat',position='append')
do i=1,Number_average
write(n,*) i,Objective_Function_average(i)
end do
close(n)

deallocate(Objective_Function)
deallocate(Objective_Function_average)

end program

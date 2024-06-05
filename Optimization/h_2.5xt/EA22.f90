
! CODE TO ILLUSTRATE THE USE OF EA22 integer NX,NMAXE,LWMAX
subroutine EA22( N, EPS, NZ, IRN, ICN, A, D, X )

integer, intent(in) :: N, NZ
integer, intent(in) :: IRN(NZ), ICN(NZ)
double precision, intent(in) :: EPS, A(NZ)

double precision, intent(out) :: D(N), X(N,N) 

integer i, j, k, l, cp

integer NX,NUMEIG,NUMCOL,LW,KMULT,IPOS 
integer ICNTL(10),KEEP(30) 

doUBLE PRECISION RKEEP(20)

integer, allocatable, dimension(:) :: IW
double precision, allocatable, dimension(:) :: U, W

!write(*,*)'==========================================='
write(*,*)'   call EA22'
!write(*,*)'==========================================='

! INITIALISE
call EA22ID(ICNTL,KEEP,RKEEP) 

NX = N ! NX >= N
NUMEIG = N-1 ! N
NUMCOL = N ! NUMEIG <= NUMCOL <= N 
LW = ( 2*N + 2*NUMCOL*NUMCOL )*5 
KMULT = ( N+1 )*( N+1 ) 
IPOS = 1

!write(*,*)'      N=', N
!write(*,*)'      EPS=', EPS
!write(*,*)'      NZ=', NZ
!write(*,*)'      NUMEIG=', NUMEIG
!write(*,*)'      NUMCOL=', NUMCOL
!write(*,*)'      Restriction :: NUMCOL ≤ N and NUMCOL ≥ NUMEIG'
!write(*,*)'      Restriction :: ', NUMCOL, ' ≤ ', N, ' and ', NUMCOL, ' ≥ ', NUMEIG

do i= 1, NZ
   if( IRN(i) <= 0 .or. ICN(i) <= 0 .or. N < IRN(i) .or. N < ICN(i) )then
      write(*,*)'IRN(i) <= 0 .or. ICN(i) <= 0 .or. N < IRN(i) .or. N < ICN(i)'
      write(*,*)'i=', i, 'NZ=', NZ
      write(*,*)'IRN(i)=', IRN(i)
      write(*,*)'ICN(i)=', ICN(i)
      write(*,*)'EA22.f90 57'
      stop
   end if
end do

allocate( IW( NUMCOL ) )
allocate( U( N ), W( LW ) )

if( NUMCOL > N .or. NUMCOL < NUMEIG )then
   write(*,*)'EA22.f90 26'
   stop
end if

do k=1,KMULT

   if( k== KMULT/100*cp )then
     !write(*,'(a,i3,a)', advance='no')'[', cp, '%]'
      write(*,'(a,i3,a,$)')'[', cp, '%]'
      cp= cp +1
   end if 

   call EA22AD(N,NUMEIG,NUMCOL,EPS,X,NX,D,U,W,LW,IW,KMULT,IPOS,ICNTL,KEEP,RKEEP) 

   if(IPOS.EQ.1) go to 40
   if(IPOS.LT.0) go to 100 
   do i=1,N
      U(i) = 0.0D0 
   end do
   do l=1,NZ
      i=IRN(l)
      j=ICN(l)
      U(i) = U(i) + A(l)*W(j)
   end do
end do

40 continue
100 continue

deallocate( IW )
deallocate( U, W )

if( IPOS/=1 )then
   write(*,*)'EA22.f90 61'
   stop
end if

!write(*,*)'==========================================='
write(*,*)'   end EA22'
!write(*,*)'==========================================='

return
end subroutine EA22 


subroutine Check_Matrix_EA22(  N, NZ, IRN, ICN, A )
   implicit none
   integer, intent(in) :: N, NZ, IRN( NZ ), ICN( NZ )
   double precision, intent(in) :: A( NZ )

   integer i, j, Flag( NZ )
   integer cu, cl, cd

   double precision, allocatable, dimension(:) :: upper, lower, diag 
   double precision, allocatable, dimension(:) :: upper_tmp, lower_tmp, diag_tmp
   integer, allocatable, dimension(:) :: iu_tmp, ju_tmp, il_tmp, jl_tmp
   integer, allocatable, dimension(:) :: iu, ju, il, jl

   write(*,*)'   call Check_Matrix_EA22'

   cu = 0
   cl = 0
   cd = 0

   allocate( upper_tmp( NZ ), lower_tmp( NZ ), diag_tmp( NZ ) )
   allocate( iu_tmp( NZ ), ju_tmp( NZ ), il_tmp( NZ ), jl_tmp( NZ ) )
   upper_tmp = 0.0d0
   lower_tmp = 0.0d0
   diag_tmp = 0.0d0

   do i = 1, NZ
      if( IRN( i ) < ICN( i ) )then
         cu = cu+1
         upper_tmp( cu ) = A( i )
         iu_tmp( cu ) = IRN( i )
         ju_tmp( cu ) = ICN( i )
      else if( IRN( i ) > ICN( i ) )then
         cl = cl+1
         lower_tmp( cl ) = A( i )
         il_tmp( cl ) = irn( i )
         jl_tmp( cl ) = icn( i )
      else
         cd = cd+1
         diag_tmp( cd ) = A( i )
      end if
   end do

   allocate( upper( cu ), lower( cl ), diag( cd ) )
   allocate( iu( cu ), ju( cu ) )
   allocate( il( cl ), jl( cl ) )

   upper( 1:cu ) = upper_tmp( 1:cu )
   lower( 1:cl ) = lower_tmp( 1:cl )
   diag( 1:cd ) = diag_tmp( 1:cd )

   
   iu( 1:cu ) = iu_tmp( 1:cu )
   ju( 1:cu ) = ju_tmp( 1:cu )
   il( 1:cl ) = il_tmp( 1:cl )
   jl( 1:cl ) = jl_tmp( 1:cl )

   deallocate( upper_tmp, lower_tmp, diag_tmp )
   deallocate( iu_tmp, ju_tmp )
   deallocate( il_tmp, jl_tmp )

   do i= 1, cu
      do j= 1, cl
         if( iu( i )==jl( j ) .and. ju( i )==il( j ) )then          
            if( abs( upper( i )-lower( j ) ) <=  1d-8 )then
               go to 67
            else       
               write(*,*)'=========================='
               write(*,*)'ERROR'
               write(*,*)'=========================='
               write(*,*)'iu(i)=', iu( i )
               write(*,*)'ju(i)=', ju( i )
               write(*,*)'upper(i)=', upper( i )
               write(*,*)'il(i)=', il( i )
               write(*,*)'jl(i)=', jl( i )
               write(*,*)'lower(i)=', lower( i )
               write(*,*)'Check_Matrix_EA22.f90  70'
               stop
            end if
         end if 
      end do
      write(*,*)'=========================='
      write(*,*)'ERROR'
      write(*,*)'=========================='
      write(*,*)'iu(i)=', iu( i )
      write(*,*)'ju(i)=', ju( i )
      write(*,*)'upper(i)=', upper( i )
      write(*,*)'Check_Matrix_EA22.f90  78'
      stop
     
      67 continue
   end do

   deallocate( upper, lower, diag )
   deallocate( iu, ju, il, jl )

   write(*,*)'   end Check_Matrix_EA22'
   return
end subroutine Check_Matrix_EA22




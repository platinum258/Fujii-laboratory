
subroutine ME57( N, NE, A, IRN, JCN, RHS, X )

   implicit none   

   integer, intent(in) :: N, NE
   complex( kind( 0d0 ) ), intent(in) :: A( NE )
   integer, intent(in) :: IRN( NE ), JCN( NE )
   complex( kind( 0d0 ) ), intent(in) :: RHS( N )
   complex( kind( 0d0 ) ), intent(out) :: X( N )

   integer :: LFACT, LIFACT, JOB, KIND_ME57, LKEEP, LWORK, NRHS, LRHS
   integer, allocatable, dimension(:) :: IFACT, IWORK, ICNTL, INFO, KEEP
   double precision, allocatable, dimension(:) :: CNTL, RINFO, RWORK 
   complex( kind( 0d0 ) ), allocatable, dimension(:) :: WORK, FACT 

   allocate( ICNTL( 20 ) )
   allocate( CNTL( 5 ) )
 
   CALL ME57ID(CNTL,ICNTL)
   !ICNTL(5) = 4
   !ICNTL(3) = 6

   !Matrix is complex symmetric (KIND_ME57 = 1) 
   !Matrix is Hermitian (KIND_ME57 = 2)
   KIND_ME57 = 1
         
   LKEEP= 5*N+ NE+ max( N,NE ) +42 !OK
   allocate( KEEP( LKEEP ) ) !OK

   allocate( IWORK( 5*N ) ) !OK
   allocate( INFO( 40 ) ) !OK
   allocate( RINFO( 20 ) ) !OK
        
   CALL ME57AD(KIND_ME57, N, NE, IRN, JCN, LKEEP,KEEP,IWORK,ICNTL,INFO,RINFO)

   if( INFO(1) < 0 )then
      write(*,*)'INFO(1)=', INFO(1)
      write(*,*)'ME57.f90 38'
      stop
   end if

   LIFACT=INFO(10)
   !LIFACT=100000 !INFO(10)
   allocate( IFACT( LIFACT ) )

   LFACT=INFO(9)
   !LFACT=100000 !INFO(9)
   allocate( FACT( LFACT ) )

   allocate( RWORK( 2*NE+3*N ) )
         
   CALL ME57BD(KIND_ME57, N, NE, A,FACT,LFACT,IFACT,LIFACT,LKEEP,KEEP,IWORK,RWORK,ICNTL,CNTL,INFO,RINFO)

   deallocate( RWORK )

   NRHS= 1 !OK
   LRHS= N !OK

   LWORK=N*NRHS !OK 
   allocate( WORK( LWORK ) )

   JOB = 1

   CALL ME57CD(KIND_ME57,JOB, N,FACT,LFACT,IFACT,LIFACT, NRHS, RHS, LRHS,WORK,LWORK,IWORK,ICNTL,INFO)

   X( : )=RHS(:)
         
   deallocate( KEEP )
   deallocate( ICNTL )
   deallocate( CNTL )
   deallocate( FACT )
   deallocate( IFACT )
   deallocate( INFO )
   deallocate( RINFO )
   deallocate( IWORK )
   deallocate( WORK )

   return
end subroutine ME57


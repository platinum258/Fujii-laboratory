

subroutine Judge_Cross_Product&
       ( Number_Vector, Vector_1, Vector_2, Flag_Error, Flag_Direction_Vector )

   !$ use omp_lib
   implicit none

   integer, intent(in) :: Number_Vector, Flag_Error 
   double precision, intent(in) :: Vector_1( 2, Number_Vector ), Vector_2( 2, Number_Vector ) 

   integer, intent(out) :: Flag_Direction_Vector( Number_Vector ) 
   integer :: i
   double precision :: error
   
   double precision, allocatable, dimension(:) :: Outer_Product 
   integer, allocatable, dimension(:) :: Flag_Direction_Vector_tmp 

   !================================================================================================================
   !write(*,*)'            call Judge_Cross_Product '
   !================================================================================================================

   if( Flag_Error==0 )then
      error= 0.0d0
   else if( Flag_Error==1 )then
      error= -1.0d-8 
   else

   end if

   allocate( Outer_Product( Number_Vector ) )
   allocate( Flag_Direction_Vector_tmp( Number_Vector ) )

   do i= 1, Number_Vector
      Flag_Direction_Vector_tmp( i )= 1

      Outer_Product( i ) &
      = Vector_1( 1, i ) *Vector_2( 2, i ) -Vector_1( 2, i ) *Vector_2( 1, i ) 

      if( Outer_Product( 1 )*Outer_Product( i ) < error )then
         Flag_Direction_Vector_tmp( i )= -1
      end if 
   end do

   deallocate( Outer_Product )

   do i= 1, Number_Vector
      Flag_Direction_Vector( i )= Flag_Direction_Vector_tmp( i )
   end do

   deallocate( Flag_Direction_Vector_tmp )

   !================================================================================================================
   !write(*,*)'            end Judge_Cross_Product '
   !================================================================================================================

   return
end subroutine Judge_Cross_Product


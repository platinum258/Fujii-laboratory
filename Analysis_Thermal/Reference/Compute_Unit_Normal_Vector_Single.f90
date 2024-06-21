
subroutine Compute_Unit_Normal_Vector_Single &
       ( Position_Node_on_Line, Unit_Normal_Vector )

   !use omp_lib
   implicit none   
 
   double precision, intent(in) :: Position_Node_on_Line( 2, 2 )

   double precision, intent(out) :: Unit_Normal_Vector( 2 )

   double precision, allocatable, dimension(:) :: Original_Vector
   double precision :: Length_Original_Vector, Length_Unit_Normal_Vector
   double precision, parameter :: Tolerance= 1d-8
   integer :: i

   !====================================================================================================
   !write(*,*)'       call Compute_Unit_Normal_Vector_Single'
   !====================================================================================================

   allocate( Original_Vector( 2 ) )

   do i= 1, 2
      Original_Vector( i )= Position_Node_on_Line( i, 2 ) -Position_Node_on_Line( i, 1 )
   end do

   Length_Original_Vector= dsqrt( Original_Vector( 1 )*Original_Vector( 1 ) +Original_Vector( 2 )*Original_Vector( 2 ) )

   Unit_Normal_Vector( 1 )= -Original_Vector( 2 )/Length_Original_Vector 
   Unit_Normal_Vector( 2 )=  Original_Vector( 1 )/Length_Original_Vector

   ! Compute Rotation of original and unit normal vectors
   if( Unit_Normal_Vector( 2 )*Original_Vector( 1 ) -Unit_Normal_Vector( 1 )*Original_Vector( 2 ) <= 0.0d0 )then
      Unit_Normal_Vector( 1 )=  Original_Vector( 2 )/Length_Original_Vector 
      Unit_Normal_Vector( 2 )= -Original_Vector( 1 )/Length_Original_Vector
   end if

   Length_Unit_Normal_Vector= dsqrt( Unit_Normal_Vector( 1 )*Unit_Normal_Vector( 1 ) +Unit_Normal_Vector( 2 )*Unit_Normal_Vector( 2 ) )

   if( Length_Unit_Normal_Vector < 1d0*( 1d0 -Tolerance ) .or. &
       Length_Unit_Normal_Vector > 1d0*( 1d0 +Tolerance ) )then
      write(*,*)'Length_Unit_Normal_Vector=', Length_Unit_Normal_Vector
      call Output_Error( 'Compute_Unit_Normal_Vector_Single', 42 )
   end if

   deallocate( Original_Vector )

   return
end subroutine Compute_Unit_Normal_Vector_Single




subroutine Compute_Unit_Normal_Vector_All &
       ( Number_Edge, Position_Node_on_Boundary, Unit_Normal_Vector )

   implicit none

   integer, intent(in) :: Number_Edge 
   double precision, intent(in) :: Position_Node_on_Boundary( 2, 2, Number_Edge )

   double precision, intent(out) :: Unit_Normal_Vector( 2, Number_Edge )

   integer :: i, j, k
   double precision, allocatable :: Position_Node_tmp(:,:), Unit_Normal_Vector_tmp(:) 

   !====================================================================================================
   write(*,*)'    call Compute_Unit_Normal_Vector_All'
   !====================================================================================================

   do i= 1, Number_Edge
      allocate( Position_Node_tmp( 2, 2 ) )

      do j= 1, 2
         do k= 1, 2
            Position_Node_tmp( k, j )= Position_Node_on_Boundary( k, j, i ) 
         end do
      end do

      allocate( Unit_Normal_Vector_tmp( 2 ) )

      call Compute_Unit_Normal_Vector_Single &
         ( Position_Node_tmp, Unit_Normal_Vector_tmp )

      do j= 1, 2
         Unit_Normal_Vector( j, i )= Unit_Normal_Vector_tmp( j )
      end do

      deallocate( Position_Node_tmp )
      deallocate( Unit_Normal_Vector_tmp )

   end do

   return
end subroutine Compute_Unit_Normal_Vector_All 



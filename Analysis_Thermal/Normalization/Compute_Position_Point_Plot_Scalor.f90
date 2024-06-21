
subroutine Compute_Position_Point_Plot_Scalor &
       ( Value_Boundary_Level_Line, &
         Position_Node_Line_1, Position_Node_Line_2, &
         Value_Node_Line_1, Value_Node_Line_2, &
         !======================================================================
         Position_Interpolated )

   use Parameters
   implicit none
  
   double precision, intent(in) :: Position_Node_Line_1( 2 ), Position_Node_Line_2( 2 )
   double precision, intent(in) :: Value_Node_Line_1, Value_Node_Line_2
   double precision, intent(in) :: Value_Boundary_Level_Line 

   double precision, intent(out) :: Position_Interpolated( 2 )
  
   integer :: i, j
  
   !==============================================================================
   !write(*,*)'call Compute_Position_Point_Plot_Scalor'
   !==============================================================================
  
   do i= 1, 2 !x, y
  
          Position_Interpolated( i ) &
         = Position_Node_Line_1( i ) &
         +( Position_Node_Line_2( i ) -Position_Node_Line_1( i ) ) &
         *( Value_Boundary_Level_Line -Value_Node_Line_1 ) &
         /( Value_Node_Line_2 -Value_Node_Line_1 ) 
  
   end do

   return
end subroutine Compute_Position_Point_Plot_Scalor



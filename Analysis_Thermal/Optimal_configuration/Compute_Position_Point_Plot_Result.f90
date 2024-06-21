
subroutine Compute_Position_Point_Plot_Result &
       ( Value_Plotted, Value_Boundary_Level_Line, &
         Local_Node_Number_Edge_Line, &
         Position_Node, Number_Node, Number_Node_Difference_Level, &
         Index_Element_2_Node, Number_Element, Element_Number, Max_Number_Level_Element, &
         !======================================================================
         Position_Node_Boundary_Level_Plot )

   use Parameters
   implicit none
  
   integer, intent(in) :: Number_Node, Number_Element 
   integer, intent(in) :: Element_Number, Number_Node_Difference_Level
   double precision, intent(in) :: Position_Node( 2, Number_Node )
   integer, intent(in) :: Index_Element_2_Node( 4, Number_Element )
   double precision, intent(in) :: Value_Plotted( Number_Node )
   integer, intent(in) :: Local_Node_Number_Edge_Line( 2, Number_Element )
   integer, intent(in) :: Max_Number_Level_Element
 
   double precision, intent(in) :: Value_Boundary_Level_Line( Max_Number_Level_Element, Number_Element ) 

   double precision, intent(out) :: Position_Node_Boundary_Level_Plot( 2, Max_Number_Level_Element, Number_Element )
  
   integer :: i, j
  
   !==============================================================================
   !write(*,*)'call Compute_Position_Point_Plot_Result'
   !==============================================================================
  
  
   do j= 1, Number_Node_Difference_Level
       do i= 1, 2 !x, y
  
          Position_Node_Boundary_Level_Plot( i, j, Element_Number ) &
         = Position_Node( i, Index_Element_2_Node( Local_Node_Number_Edge_Line( 1, Element_Number ), Element_Number ) ) &
         +( Position_Node( i, Index_Element_2_Node( Local_Node_Number_Edge_Line( 2, Element_Number ), Element_Number ) ) &
           -Position_Node( i, Index_Element_2_Node( Local_Node_Number_Edge_Line( 1, Element_Number ), Element_Number ) ) ) &
         *( Value_Boundary_Level_Line( j, Element_Number ) &
           -Value_Plotted( Index_Element_2_Node( Local_Node_Number_Edge_Line( 1, Element_Number ), Element_Number ) ) ) &
         /( Value_Plotted( Index_Element_2_Node( Local_Node_Number_Edge_Line( 2, Element_Number ), Element_Number ) ) & 
           -Value_Plotted( Index_Element_2_Node( Local_Node_Number_Edge_Line( 1, Element_Number ), Element_Number ) ) ) 
  
       end do
   end do

   return
end subroutine Compute_Position_Point_Plot_Result



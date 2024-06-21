
subroutine Create_PlotData_Element_Homogeneous &
           ( Position_Node_Plot_Mesh, Index_Element_2_Node, Level_Node, & 
             Element_Number, Number_Node, Number_Element, Number_Level, Max_Number_Level_Element, &
             !================================================================================================= 
             Number_Level_Element, Number_Node_Level_Element, Position_Node_Plot_Mesh_Element, Level_Element )

     !use omp_lib
     use Parameters
     implicit none

     integer, intent(in) :: Element_Number, Number_Node, Number_Element
     integer, intent(in) :: Index_Element_2_Node( 4, Number_Element )
     double precision, intent(in) :: Position_Node_Plot_Mesh( 2, Number_Node )
     integer, intent(in) :: Level_Node( Number_Node )
     integer, intent(in) :: Number_Level, Max_Number_Level_Element

     integer, intent(out) :: Number_Level_Element( Number_Element )
     integer, intent(out) :: Number_Node_Level_Element( Max_Number_Level_Element, Number_Element )
     double precision, intent(out) :: Position_Node_Plot_Mesh_Element( 2, 4, Max_Number_Level_Element, Number_Element )
     integer, intent(out) :: Level_Element( Max_Number_Level_Element, Number_Element )

     integer :: i, j 

     Number_Level_Element( Element_Number )= 1

     if( Index_Element_2_Node( 3, Element_Number )==Index_Element_2_Node( 4, Element_Number ) )then

          Number_Node_Level_Element( 1, Element_Number )= 3 
     
          do i= 1, Number_Node_Level_Element( 1, Element_Number )
               do j= 1, 2 
                    Position_Node_Plot_Mesh_Element( j, i, 1, Element_Number )= Position_Node_Plot_Mesh( j, Index_Element_2_Node( i, Element_Number ) )
               end do
          end do

     else
  
          Number_Node_Level_Element( 1, Element_Number )= 4 
     
          do i= 1, Number_Node_Level_Element( 1, Element_Number )
               do j= 1, 2 
                    Position_Node_Plot_Mesh_Element( j, i, 1, Element_Number )= Position_Node_Plot_Mesh( j, Index_Element_2_Node( i, Element_Number ) )
               end do
          end do

     end if


     if( ( Level_Node( Index_Element_2_Node( 1, j ) ) <= 0 .and. Level_Node( Index_Element_2_Node( 2, j ) ) <= 0 ) .or. & 
         ( Level_Node( Index_Element_2_Node( 2, j ) ) <= 0 .and. Level_Node( Index_Element_2_Node( 3, j ) ) <= 0 ) .or. & 
         ( Level_Node( Index_Element_2_Node( 3, j ) ) <= 0 .and. Level_Node( Index_Element_2_Node( 1, j ) ) <= 0 ) )then 

          Level_Element( 1, Element_Number )= 0 

     else if( ( Level_Node( Index_Element_2_Node( 1, j ) ) > Number_Level .and. &
                Level_Node( Index_Element_2_Node( 2, j ) ) > Number_Level ) .or. & 
              ( Level_Node( Index_Element_2_Node( 2, j ) ) > Number_Level .and. &
                Level_Node( Index_Element_2_Node( 3, j ) ) > Number_Level ) .or. & 
              ( Level_Node( Index_Element_2_Node( 3, j ) ) > Number_Level .and. &
                Level_Node( Index_Element_2_Node( 1, j ) ) > Number_Level ) )then

          Level_Element( 1, Element_Number )= Number_Level +1
   
     else if( Level_Node( Index_Element_2_Node( 1, j ) ) <= 0 .or. & 
              Level_Node( Index_Element_2_Node( 2, j ) ) <= 0 .or. & 
              Level_Node( Index_Element_2_Node( 3, j ) ) <= 0 )then

          Level_Element( 1, Element_Number )= 0 

     else if( Level_Node( Index_Element_2_Node( 1, j ) ) > Number_Level .or. & 
              Level_Node( Index_Element_2_Node( 2, j ) ) > Number_Level .or. & 
              Level_Node( Index_Element_2_Node( 3, j ) ) > Number_Level )then

          Level_Element( 1, Element_Number )= Number_Level +1
     else  
          Level_Element( 1, Element_Number )= Level_Node( Index_Element_2_Node( 1, Element_Number ) )
     end if

     return
end subroutine Create_PlotData_Element_Homogeneous


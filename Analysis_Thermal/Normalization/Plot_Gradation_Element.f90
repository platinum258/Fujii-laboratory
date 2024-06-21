

subroutine Plot_Gradation_Element &
       ( File_Number, Element_Number, Number_Element, Number_Level, Max_Number_Level_Element, & 
         Number_Level_Element, Number_Node_Level_Element, Position_Node_Plot_Mesh_Element, RGB, Level_Element )

   use Parameters
   implicit none

   integer, intent(in) :: File_Number, Element_Number, Number_Element, Number_Level, Max_Number_Level_Element

   integer, intent(in) :: Number_Level_Element( Number_Element )
   integer, intent(in) :: Number_Node_Level_Element( Max_Number_Level_Element, Number_Element )
   double precision, intent(in) :: Position_Node_Plot_Mesh_Element( 2, 4, Max_Number_Level_Element, Number_Element )
   double precision, intent(in) :: RGB( 3, Number_Level )
   integer, intent(in) :: Level_Element( Max_Number_Level_Element, Number_Element )

   integer :: i, j

   !===============================================================================================================
   !write(*,*)'      Plot_Gradation_Element'
   !write(*,*)'      Number_Level_Element( Element_Number )=', Number_Level_Element( Element_Number )
   !write(*,*)'      Number_Node_Level_Element( 1, Element_Number )=', Number_Node_Level_Element( 1, Element_Number )
   !===============================================================================================================
   
   do i= 1, Number_Level_Element( Element_Number ) 

      write(File_Number,*)'gsave'
      write(File_Number,*) Position_Node_Plot_Mesh_Element( 1, 1, i, Element_Number ), &
                   Position_Node_Plot_Mesh_Element( 2, 1, i, Element_Number ), &
                  'moveto'

      do j= 2, Number_Node_Level_Element( i, Element_Number )
         write(File_Number,*) Position_Node_Plot_Mesh_Element( 1, j, i, Element_Number ), &
                      Position_Node_Plot_Mesh_Element( 2, j, i, Element_Number ), &
                     'lineto'
      end do

      write(File_Number,*) Position_Node_Plot_Mesh_Element( 1, 1, i, Element_Number ), &
                   Position_Node_Plot_Mesh_Element( 2, 1, i, Element_Number ), &
                  'lineto'
 
      write(File_Number,*)'closepath'
      write(File_Number,*) RGB( 1, Level_Element( i, Element_Number ) ), & 
                   RGB( 2, Level_Element( i, Element_Number ) ), & 
                   RGB( 3, Level_Element( i, Element_Number ) ), & 
                   'setrgbcolor' 
      write(File_Number,*)'fill'
      write(File_Number,*)'grestore'
      write(File_Number,*)'newpath'

   end do

   return 
end subroutine Plot_Gradation_Element



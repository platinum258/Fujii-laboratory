
subroutine Plot_Element_InHomogeneous &
       ( Position_Node_Plot_Mesh, Index_Element_2_Node, Level_Node, RGB, & 
         File_Number, Element_Number, Number_Node, Number_Element, Number_Level, Max_Number_Level_Element, & 
         Number_Node_Difference_Level, Local_Node_Number_Different_Level, Local_Node_Nubmer_Edge_Line, &
         Position_Node_Boundary_Level_Plot_1, Position_Node_Boundary_Level_Plot_2 )

   !use omp_lib
   use Parameters
   implicit none

   integer, intent(in) :: File_Number, Element_Number, Number_Node, Number_Element, Number_Level, Max_Number_Level_Element
   integer, intent(in) :: Number_Node_Difference_Level, Local_Node_Number_Different_Level
   integer, intent(in) :: Index_Element_2_Node( 4, Number_Element )
   double precision, intent(in) :: Position_Node_Plot_Mesh( 2, Number_Node )
   double precision, intent(in) :: RGB( 3, Number_Level )
   integer, intent(in) :: Level_Node( Number_Node )
   integer, intent(in) :: Local_Node_Nubmer_Edge_Line( 2, 2, Number_Element )

   double precision, intent(in) :: Position_Node_Boundary_Level_Plot_1( 2, Max_Number_Level_Element, Number_Element ) 
   double precision, intent(in) :: Position_Node_Boundary_Level_Plot_2( 2, Max_Number_Level_Element, Number_Element ) 

   integer, allocatable, dimension(:) :: Level_Element

   integer :: i 
   integer :: Coefficient_Plot 

   if( Level_Node( Index_Element_2_Node( Local_Node_Nubmer_Edge_Line( 1, 1, Element_Number ), Element_Number ) )  & 
     > Level_Node( Index_Element_2_Node( Local_Node_Nubmer_Edge_Line( 2, 1, Element_Number ), Element_Number ) ) )then

      Coefficient_Plot= -1

   else if( Level_Node( Index_Element_2_Node( Local_Node_Nubmer_Edge_Line( 1, 1, Element_Number ), Element_Number ) ) & 
        < Level_Node( Index_Element_2_Node( Local_Node_Nubmer_Edge_Line( 2, 1, Element_Number ), Element_Number ) ) )then

      Coefficient_Plot= 1

   else
      write(*,*)'Plot_Element_InHomogeneous.f90 39'
      stop
   end if



   if( Index_Element_2_Node( 3, Element_Number )==Index_Element_2_Node( 4, Element_Number ) )then

      !==============================================================================================================================
      write(File_Number,*)'gsave'
      write(File_Number,*) Position_Node_Plot_Mesh( 1, Index_Element_2_Node( Local_Node_Number_Different_Level, Element_Number ) ), &
                   Position_Node_Plot_Mesh( 2, Index_Element_2_Node( Local_Node_Number_Different_Level, Element_Number ) ), &
                   'moveto'
       
      write(File_Number,*) Position_Node_Boundary_Level_Plot_1( 1, 1, Element_Number ), &
                   Position_Node_Boundary_Level_Plot_1( 2, 1, Element_Number ), &
                   'lineto'

      write(File_Number,*) Position_Node_Boundary_Level_Plot_2( 1, 1, Element_Number ), &
                   Position_Node_Boundary_Level_Plot_2( 2, 1, Element_Number ), &
                   'lineto'
      
      allocate( Level_Element( Number_Element ) )
   
      Level_Element( Element_Number )= Level_Node( Index_Element_2_Node( Local_Node_Number_Different_Level, Element_Number ) )
    
      write(File_Number,*)'closepath'
      write(File_Number,*) RGB( 1, Level_Element( Element_Number ) ), & 
                   RGB( 2, Level_Element( Element_Number ) ), & 
                   RGB( 3, Level_Element( Element_Number ) ), & 
                   'setrgbcolor' 
      !write(File_Number,*) 0.0d0, 0.0d0, 0.0d0, 'setrgbcolor' 
      write(File_Number,*)'fill'
      write(File_Number,*)'grestore'
      write(File_Number,*)'newpath'
 
      deallocate( Level_Element )
      !==============================================================================================================================

      do i= 1, Number_Node_Difference_Level -1 
         !==============================================================================================================================
         write(File_Number,*)'gsave'
         write(File_Number,*) Position_Node_Boundary_Level_Plot_1( 1, i, Element_Number ), &
                      Position_Node_Boundary_Level_Plot_1( 2, i, Element_Number ), &
                      'moveto' 
   
         write(File_Number,*) Position_Node_Boundary_Level_Plot_1( 1, i+1, Element_Number ), &
                      Position_Node_Boundary_Level_Plot_1( 2, i+1, Element_Number ), &
                      'lineto'
         
         write(File_Number,*) Position_Node_Boundary_Level_Plot_2( 1, i+1, Element_Number ), &
                      Position_Node_Boundary_Level_Plot_2( 2, i+1, Element_Number ), &
                      'lineto'
         
         write(File_Number,*) Position_Node_Boundary_Level_Plot_2( 1, i, Element_Number ), &
                      Position_Node_Boundary_Level_Plot_2( 2, i, Element_Number ), &
                      'lineto'
         
         allocate( Level_Element( Number_Element ) )
      
         Level_Element( Element_Number ) &
         = Level_Node( Index_Element_2_Node( Local_Node_Number_Different_Level, Element_Number ) ) +Coefficient_Plot*i 
       
         write(File_Number,*)'closepath'
         write(File_Number,*) RGB( 1, Level_Element( Element_Number ) ), & 
                      RGB( 2, Level_Element( Element_Number ) ), & 
                      RGB( 3, Level_Element( Element_Number ) ), & 
                      'setrgbcolor' 
         !write(File_Number,*) 0.0d0, 0.0d0, 0.0d0, 'setrgbcolor' 
         write(File_Number,*)'fill'
         write(File_Number,*)'grestore'
         write(File_Number,*)'newpath'
      
         deallocate( Level_Element )
         !==============================================================================================================================
      end do

      !==============================================================================================================================
      write(File_Number,*)'gsave'
      write(File_Number,*) Position_Node_Boundary_Level_Plot_1( 1, Number_Node_Difference_Level, Element_Number ), &
                   Position_Node_Boundary_Level_Plot_1( 2, Number_Node_Difference_Level, Element_Number ), &
                   'moveto' 
   
      write(File_Number,*) Position_Node_Plot_Mesh( 1, Index_Element_2_Node( Local_Node_Nubmer_Edge_Line( 2, 1, Element_Number ), Element_Number ) ), &
                   Position_Node_Plot_Mesh( 2, Index_Element_2_Node( Local_Node_Nubmer_Edge_Line( 2, 1, Element_Number ), Element_Number ) ), &
                   'lineto'
         
      write(File_Number,*) Position_Node_Plot_Mesh( 1, Index_Element_2_Node( Local_Node_Nubmer_Edge_Line( 2, 2, Element_Number ), Element_Number ) ), &
                   Position_Node_Plot_Mesh( 2, Index_Element_2_Node( Local_Node_Nubmer_Edge_Line( 2, 2, Element_Number ), Element_Number ) ), &
                   'lineto'
         
      write(File_Number,*) Position_Node_Boundary_Level_Plot_2( 1, Number_Node_Difference_Level, Element_Number ), &
                   Position_Node_Boundary_Level_Plot_2( 2, Number_Node_Difference_Level, Element_Number ), &
                   'lineto'
         
      allocate( Level_Element( Number_Element ) )
      
      Level_Element( Element_Number ) &
      =  Level_Node( Index_Element_2_Node( Local_Node_Number_Different_Level, Element_Number ) ) & 
        +Coefficient_Plot*Number_Node_Difference_Level 
       
      write(File_Number,*)'closepath'
      write(File_Number,*) RGB( 1, Level_Element( Element_Number ) ), & 
                   RGB( 2, Level_Element( Element_Number ) ), & 
                   RGB( 3, Level_Element( Element_Number ) ), & 
                   'setrgbcolor' 
      !write(File_Number,*) 0.0d0, 0.0d0, 0.0d0, 'setrgbcolor' 
      write(File_Number,*)'fill'
      write(File_Number,*)'grestore'
      write(File_Number,*)'newpath'
     
      deallocate( Level_Element )

    
   else
      write(*,*)'Plot_Element_InHomogeneous.f90 64 '
      stop

   end if
 
   return
end subroutine Plot_Element_InHomogeneous


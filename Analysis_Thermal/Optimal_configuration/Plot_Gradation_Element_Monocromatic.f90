

subroutine Plot_Gradation_Element_Monocromatic &
       ( File_Number, Element_Number, Number_Element, Max_Number_Level_Element, & 
         Number_Level_Element, Number_Node_Level_Element, Position_Node_Plot_Mesh_Element, Flag_RGB_Monocromatic, &
         RGB_Lowest, RGB_Highest )

   use Parameters
   implicit none

   integer, intent(in) :: File_Number, Element_Number, Number_Element, Max_Number_Level_Element 

   integer, intent(in) :: Number_Level_Element( Number_Element )
   integer, intent(in) :: Number_Node_Level_Element( Max_Number_Level_Element, Number_Element )
   double precision, intent(in) :: Position_Node_Plot_Mesh_Element( 2, 4, Max_Number_Level_Element, Number_Element )
   integer, intent(in) :: Flag_RGB_Monocromatic
   double precision, intent(in) :: RGB_Lowest( 3 ), RGB_Highest( 3 ) 

   integer :: i, j
   double precision, allocatable, dimension(:) :: RGB_Monocromatic

   !===============================================================================================================
   !write(*,*)'      Plot_Gradation_Element_Monocromatic'
   !write(*,*)'      Number_Level_Element( Element_Number )=', Number_Level_Element( Element_Number )
   !write(*,*)'      Number_Node_Level_Element( 1, Element_Number )=', Number_Node_Level_Element( 1, Element_Number )
   !===============================================================================================================
    
   allocate( RGB_Monocromatic( 3 ) )

   if( Flag_RGB_Monocromatic==0 )then
      RGB_Monocromatic( 1 )= RGB_Lowest( 1 ) 
      RGB_Monocromatic( 2 )= RGB_Lowest( 2 )
      RGB_Monocromatic( 3 )= RGB_Lowest( 3 )
   else if( Flag_RGB_Monocromatic==1 )then
      RGB_Monocromatic( 1 )= RGB_Highest( 1 )
      RGB_Monocromatic( 2 )= RGB_Highest( 2 )
      RGB_Monocromatic( 3 )= RGB_Highest( 3 )
   end if
 
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
      write(File_Number,*) RGB_Monocromatic( 1 ), RGB_Monocromatic( 2 ), RGB_Monocromatic( 3 ), 'setrgbcolor' 
      write(File_Number,*)'fill'
      write(File_Number,*)'grestore'
      write(File_Number,*)'newpath'

   end do

   deallocate( RGB_Monocromatic )

   return 
end subroutine Plot_Gradation_Element_Monocromatic



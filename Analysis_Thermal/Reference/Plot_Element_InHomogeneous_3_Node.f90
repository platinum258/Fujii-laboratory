
subroutine Plot_Element_InHomogeneous_3_Node &
       ( Position_Node_Plot_Mesh, Index_Element_2_Node, Level_Node, RGB, & 
         File_Number, Element_Number, Number_Node, Number_Element, Number_Level, Max_Number_Level_Element, & 
         Number_Node_Difference_Level_3_Node, Local_Node_Number_Level_Max_2_Min, &
         Position_Node_Boundary_Level_Plot )

   !use omp_lib
   use Parameters
   implicit none

   integer, intent(in) :: File_Number, Element_Number, Number_Node, Number_Element, Number_Level, Max_Number_Level_Element
   integer, intent(in) :: Number_Node_Difference_Level_3_Node( 3, Number_Element )
   integer, intent(in) :: Local_Node_Number_Level_Max_2_Min( 3, Number_Element )
   integer, intent(in) :: Index_Element_2_Node( 4, Number_Element )
   double precision, intent(in) :: Position_Node_Plot_Mesh( 2, Number_Node )
   double precision, intent(in) :: RGB( 3, Number_Level )
   integer, intent(in) :: Level_Node( Number_Node )

   double precision, intent(in) :: Position_Node_Boundary_Level_Plot( 2, Max_Number_Level_Element, 3, Number_Element ) 

   integer, allocatable, dimension(:) :: Level_Element

   integer :: i 

      !==============================================================================================================================
      write(File_Number,*)'gsave'
      write(File_Number,*) Position_Node_Plot_Mesh( 1, Index_Element_2_Node( Local_Node_Number_Level_Max_2_Min( 1, Element_Number ), Element_Number ) ), &
                   Position_Node_Plot_Mesh( 2, Index_Element_2_Node( Local_Node_Number_Level_Max_2_Min( 1, Element_Number ), Element_Number ) ), &
                   'moveto'
       
      write(File_Number,*) Position_Node_Boundary_Level_Plot( 1, 1, 1, Element_Number ), &
                   Position_Node_Boundary_Level_Plot( 2, 1, 1, Element_Number ), &
                   'lineto'

      write(File_Number,*) Position_Node_Boundary_Level_Plot( 1, 1, 2, Element_Number ), &
                   Position_Node_Boundary_Level_Plot( 2, 1, 2, Element_Number ), &
                   'lineto'
      
      allocate( Level_Element( Number_Element ) )
   
      Level_Element( Element_Number )= Level_Node( Index_Element_2_Node( Local_Node_Number_Level_Max_2_Min( 1, Element_Number ), Element_Number ) )
    
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

      do i= 1, Number_Node_Difference_Level_3_Node( 1, Element_Number ) -1 
         !==============================================================================================================================
         write(File_Number,*)'gsave'
         write(File_Number,*) Position_Node_Boundary_Level_Plot( 1, i, 1, Element_Number ), &
                      Position_Node_Boundary_Level_Plot( 2, i, 1, Element_Number ), &
                      'moveto' 
   
         write(File_Number,*) Position_Node_Boundary_Level_Plot( 1, i+1, 1, Element_Number ), &
                      Position_Node_Boundary_Level_Plot( 2, i+1, 1, Element_Number ), &
                      'lineto'
         
         write(File_Number,*) Position_Node_Boundary_Level_Plot( 1, i+1, 2, Element_Number ), &
                      Position_Node_Boundary_Level_Plot( 2, i+1, 2, Element_Number ), &
                      'lineto'
         
         write(File_Number,*) Position_Node_Boundary_Level_Plot( 1, i, 2, Element_Number ), &
                      Position_Node_Boundary_Level_Plot( 2, i, 2, Element_Number ), &
                      'lineto'
         
         allocate( Level_Element( Number_Element ) )
      
         Level_Element( Element_Number ) &
         = Level_Node( Index_Element_2_Node( Local_Node_Number_Level_Max_2_Min( 1, Element_Number ), Element_Number ) ) -i 
       
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
      write(File_Number,*) Position_Node_Plot_Mesh( 1, Index_Element_2_Node( Local_Node_Number_Level_Max_2_Min( 2, Element_Number ), Element_Number ) ), &
                   Position_Node_Plot_Mesh( 2, Index_Element_2_Node( Local_Node_Number_Level_Max_2_Min( 2, Element_Number ), Element_Number ) ), &
                   'moveto' 
         
      write(File_Number,*) Position_Node_Boundary_Level_Plot( 1, Number_Node_Difference_Level_3_Node( 1, Element_Number ), 1, Element_Number ), &
                   Position_Node_Boundary_Level_Plot( 2, Number_Node_Difference_Level_3_Node( 1, Element_Number ), 1, Element_Number ), &
                   'lineto'
   
      write(File_Number,*) Position_Node_Boundary_Level_Plot( 1, Number_Node_Difference_Level_3_Node( 1, Element_Number ), 2, Element_Number ), &
                   Position_Node_Boundary_Level_Plot( 2, Number_Node_Difference_Level_3_Node( 1, Element_Number ), 2, Element_Number ), &
                   'lineto'
   
      write(File_Number,*) Position_Node_Boundary_Level_Plot( 1, Number_Node_Difference_Level_3_Node( 1, Element_Number )+1, 2, Element_Number ), &
                   Position_Node_Boundary_Level_Plot( 2, Number_Node_Difference_Level_3_Node( 1, Element_Number )+1, 2, Element_Number ), &
                   'lineto'
   
      write(File_Number,*) Position_Node_Boundary_Level_Plot( 1, 1, 3, Element_Number ), &
                   Position_Node_Boundary_Level_Plot( 2, 1, 3, Element_Number ), &
                   'lineto'
      
      allocate( Level_Element( Number_Element ) )
      
      Level_Element( Element_Number ) &
      =  Level_Node( Index_Element_2_Node( Local_Node_Number_Level_Max_2_Min( 1, Element_Number ), Element_Number ) ) & 
        -Number_Node_Difference_Level_3_Node( 1, Element_Number ) 
       
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

      !===========================================================================================================

      do i= 1, Number_Node_Difference_Level_3_Node( 3, Element_Number ) -1 

         write(File_Number,*)'gsave'
         write(File_Number,*) Position_Node_Boundary_Level_Plot( 1, i, 3, Element_Number ), &
                      Position_Node_Boundary_Level_Plot( 2, i, 3, Element_Number ), &
                      'moveto' 
   
         write(File_Number,*) Position_Node_Boundary_Level_Plot( 1, i+1, 3, Element_Number ), &
                      Position_Node_Boundary_Level_Plot( 2, i+1, 3, Element_Number ), &
                      'lineto'
         
         write(File_Number,*) Position_Node_Boundary_Level_Plot( 1, i +1 +Number_Node_Difference_Level_3_Node( 1, Element_Number ), 2, Element_Number ), &
                      Position_Node_Boundary_Level_Plot( 2, i +1 +Number_Node_Difference_Level_3_Node( 1, Element_Number ), 2, Element_Number ), &
                      'lineto'
         
         write(File_Number,*) Position_Node_Boundary_Level_Plot( 1, i +Number_Node_Difference_Level_3_Node( 1, Element_Number ), 2, Element_Number ), &
                      Position_Node_Boundary_Level_Plot( 2, i +Number_Node_Difference_Level_3_Node( 1, Element_Number ), 2, Element_Number ), &
                      'lineto'
 
         allocate( Level_Element( Number_Element ) )
         
         Level_Element( Element_Number ) &
         =  Level_Node( Index_Element_2_Node( Local_Node_Number_Level_Max_2_Min( 1, Element_Number ), Element_Number ) ) & 
            -Number_Node_Difference_Level_3_Node( 1, Element_Number ) -i  
          
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

      end do

      !=================================================================================================================

      write(File_Number,*)'gsave'
      write(File_Number,*) Position_Node_Boundary_Level_Plot( 1, Number_Node_Difference_Level_3_Node( 3, Element_Number ), 3, Element_Number ), &
                   Position_Node_Boundary_Level_Plot( 2, Number_Node_Difference_Level_3_Node( 3, Element_Number ), 3, Element_Number ), &
                   'moveto' 
   
      write(File_Number,*) Position_Node_Boundary_Level_Plot( 1, Number_Node_Difference_Level_3_Node( 2, Element_Number ), 2, Element_Number ), &
                   Position_Node_Boundary_Level_Plot( 2, Number_Node_Difference_Level_3_Node( 2, Element_Number ), 2, Element_Number ), &
                   'lineto'
         
      write(File_Number,*) Position_Node_Plot_Mesh( 1, Index_Element_2_Node( Local_Node_Number_Level_Max_2_Min( 3, Element_Number ), Element_Number ) ), &
                   Position_Node_Plot_Mesh( 2, Index_Element_2_Node( Local_Node_Number_Level_Max_2_Min( 3, Element_Number ), Element_Number ) ), &
                   'lineto'
 
      allocate( Level_Element( Number_Element ) )
         
      Level_Element( Element_Number ) &
      =  Level_Node( Index_Element_2_Node( Local_Node_Number_Level_Max_2_Min( 3, Element_Number ), Element_Number ) ) 
          
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



      return
end subroutine Plot_Element_InHomogeneous_3_Node


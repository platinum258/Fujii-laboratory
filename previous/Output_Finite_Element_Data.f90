

subroutine Output_Finite_Element_Data &
       ( Number_Node, Position_Node, &
         Max_Number_Element_Share_OneNode, Max_Position_Node, &
         Number_Element_Triangle, Index_Element_2_Node_Triangle, Class_Element_Triangle, &
         Element_and_LocalNode_Number_on_PEC_BC, Number_Node_on_PEC_Boundary, Number_Node_in_PEC, Number_Element_PEC_Boundary, &
         Flag_Output, Optimization_Step )!, Value_ObjectiveFunction_No_Cloak_TM, Value_ObjectiveFunction_No_Cloak_TE )

   !use omp_lib
   use Parameters
   implicit none

   integer, intent(in) :: Number_Node, Number_Element_Triangle 
   double precision, intent(in) :: Position_Node( 2, Number_Node ) 

   integer, intent(in) :: Index_Element_2_Node_Triangle( 3, Number_Element_Triangle ) 
   integer, intent(in) :: Class_Element_Triangle( Number_Element_Triangle ) 

   integer, intent(in) :: Number_Node_on_PEC_Boundary, Number_Node_in_PEC, Number_Element_PEC_Boundary
   integer, intent(in) :: Element_and_LocalNode_Number_on_PEC_BC( 3, Number_Element_PEC_Boundary ) 

   integer, intent(in) :: Max_Number_Element_Share_OneNode
   double precision, intent(in) :: Max_Position_Node( 2, 2 )

   integer, intent(in) :: Optimization_Step
   integer, intent(in) :: Flag_Output
 
   !double precision, intent(in) :: Value_ObjectiveFunction_No_Cloak_TM, Value_ObjectiveFunction_No_Cloak_TE
 
   character(len=8) :: Format_Filenumber
   character(len=Length_Character_Optimization_Step) :: Optimization_Step_Character
   character(len=256) :: Filename_FEM_Data
   integer :: i

   !================================================
   write(*,*)'   call Output_Finite_Element_Data'
   !================================================

   write( Format_Filenumber, '(A2,I1,A1,I1,A1)' ) & 
   '(i', Length_Character_Optimization_Step,'.',Length_Character_Optimization_Step,')'
  
   write( Optimization_Step_Character, Format_Filenumber ) Optimization_Step
 
   if( Flag_Output==0 )then 
      Filename_FEM_Data= trim( "FEM_Data_"//trim(Optimization_Step_Character) )
   else if( Flag_Output==1 )then 
      Filename_FEM_Data= trim( "FEM_Data_Convergence_"//trim(Optimization_Step_Character) )
   else if( Flag_Output==2 )then 
      Filename_FEM_Data= trim( "FEM_Data_Interval_"//trim(Optimization_Step_Character) )
   else if( Flag_Output==3 )then 
      Filename_FEM_Data= trim( "FEM_Data_Avr_Convergence_"//trim(Optimization_Step_Character) )
   else if( Flag_Output==4 )then 
      Filename_FEM_Data= trim( "FEM_Data_xmean_Convergence_"//trim(Optimization_Step_Character) )
   else
      write(*,*)'Flag_Output=', Flag_Output, '/= 0, 1, 2'
      write(*,*)'Error : Output_Finite_Element_Data 58'
      stop
   end if

   open( 10, file=Filename_FEM_Data )
      write(10,*) Flag_Electrical_Insulation_BC, '# Flag_Electrical_Insulation_BC '
      write(10,*) Number_Node, '# The Number of Nodes'
      write(10,*) Number_Element_Triangle, '# The Number of Triangle Elements'
      write(10,*) Max_Number_Element_Share_OneNode, '# Maximum Number of Elements sharing one Nodes '
      write(10,*) Number_Node_on_PEC_Boundary, '# The Number of Nodes on PEC Boundary'
      write(10,*) Number_Node_in_PEC, '# The Number of Nodes in PEC'
      write(10,*) Number_Element_PEC_Boundary, '# The Number of Element on PEC Boundary'

      write(10,*) ID_Element_Obj_Func_1, '# ID_Element_Obj_Func_1' 
      write(10,*) ID_Element_Material, ID_Element_OuterDomain, ID_Element_FixedDomain, ID_Element_Base_Material, '# ID_Element'

      write(10,*) ID_Element_Material_Exterior, ID_Element_Base_Material_Exterior, '# ID_Element_Exterior'

      write(10,*) Max_Position_Node( 1, 1 ), Max_Position_Node( 1, 2 ), Max_Position_Node( 2, 1 ), Max_Position_Node( 2, 2 )

      do i= 1, Number_Node
         write(10, 78 ) i, Position_Node( 1, i ), 'd0', Position_Node( 2, i ), 'd0'
      end do
      78 format( 1x, i10, 5x, f15.12, a2, 5x, f15.12, a2 )

      do i= 1, Number_Element_Triangle
         write(10,*) i, Class_Element_Triangle( i ), &
                     Index_Element_2_Node_Triangle( 1, i ), Index_Element_2_Node_Triangle( 2, i ), Index_Element_2_Node_Triangle( 3, i )
      end do

      !do i= 1, Number_Element_PEC_Boundary
      !   write(10,*) i, Element_and_LocalNode_Number_on_PEC_BC( 1, i ), &
      !           Element_and_LocalNode_Number_on_PEC_BC( 2, i ), &
      !           Element_and_LocalNode_Number_on_PEC_BC( 3, i )
      !end do

   close( 10 )

   !8743    FORMAT(I7, 2x, e15.8, 2x, e15.8, 2x, e15.8)
 
   return
end subroutine Output_Finite_Element_Data



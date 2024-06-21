
subroutine Counter_Number_Node_in_Element &
       ( Number_Node, Number_Element, Index_Element_2_Node, Class_Element, ID_Element, &
         Number_Node_in_Element  )

   !$use omp_lib
   use Parameters
   implicit none

   integer, intent(in) :: Number_Node, Number_Element, ID_Element
   integer, intent(in) :: Index_Element_2_Node( 3, Number_Element ), Class_Element( Number_Element ) 
   integer, intent(out) :: Number_Node_in_Element

   integer :: e, i 
   integer, allocatable, dimension(:) :: Flag_Node_Objective_Function
   integer :: Counter_Node

   !====================================================================================
   write(*,*)'    call Counter_Number_Node_in_Element'
   !====================================================================================

   allocate( Flag_Node_Objective_Function( Number_Node ) )

   do i= 1, Number_Node
      Flag_Node_Objective_Function( i )= 0 
   end do

   do e= 1, Number_Element
      if( Class_Element( e )==ID_Element )then
         do i= 1, 3
            Flag_Node_Objective_Function( Index_Element_2_Node( i, e ) )= 1 
         end do
      end if
   end do

   Counter_Node=0

    do i= 1, Number_Node
      if( Flag_Node_Objective_Function( i )==1 )then
         Counter_Node= Counter_Node +1
      end if
   end do

   Number_Node_in_Element= Counter_Node

   deallocate( Flag_Node_Objective_Function )

   return
end subroutine Counter_Number_Node_in_Element


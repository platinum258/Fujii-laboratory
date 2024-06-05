
subroutine Implement_Neumann_BC_complex_CSRF &
       ( Number_Element, Flag_Element_Neumann_BC, Index_Element_2_Node_Triangle, &
         Number_Node, Width_Matrix_LHS, &
         !================================================================
         GlobalMatrix, J_GlobalMatrix, Global_Vector_RHS )
    
   !$use omp_lib
   implicit none

   integer, intent(in) :: Number_Element, Number_Node, Width_Matrix_LHS

   complex( kind( 0d0 ) ), intent(inout) :: GlobalMatrix( Number_Node, Width_Matrix_LHS )
   integer, intent(inout) :: J_GlobalMatrix( Number_Node, Width_Matrix_LHS )
   complex( kind( 0d0 ) ), intent(inout) :: Global_Vector_RHS( Number_Node )
   integer, intent(in) :: Flag_Element_Neumann_BC( Number_Element )
   integer, intent(in) :: Index_Element_2_Node_Triangle( 3, Number_Element )

   complex( kind( 0d0 ) ), parameter :: Zero= dcmplx( 0.0d0, 0.0d0 )
   complex( kind( 0d0 ) ), parameter :: One= dcmplx( 1.0d0, 0.0d0 )

   integer :: i, j, k

   complex( kind( 0d0 ) ), allocatable, dimension(:) :: Global_Vector_RHS_PEC

   integer :: Number_Group_Node_Neumann_BC
   integer, allocatable, dimension(:,:) :: NodeNumber_Group, NodeNumber_Group_tmp

   integer :: Number_Element_Neumann_BC
   integer, allocatable, dimension(:,:) :: Index_Element_2_Node_Neumann_BC, Index_Element_2_Node_Neumann_BC_tmp

   integer, allocatable, dimension(:) :: Number_Node_Group_Neumann_BC 

   integer, parameter :: Flag_On= 1
   integer, parameter :: Flag_Off= 0

   !================================================================================================================
   write(*,*)'            call Implement_Neumann_BC_complex_CSRF'
   !================================================================================================================

      !===============================================
      write(*,*) '            Check Matrix'
      !===============================================

      !$omp parallel do default( none ) &
      !$omp private( i ) &
      !$omp shared( Number_Node, Flag_Element_Neumann_BC ) 
      do i= 1, Number_Node
         if( Flag_Element_Neumann_BC( i )/=Flag_On .and. &
             Flag_Element_Neumann_BC( i )/=Flag_Off )then
            call Output_Error( 'Implement_Neumann_BC_complex_CSRF' , 46 )
         end if
      end do

      !==========================================================
      write(*,*) '            Create Neumann B.C. Data'
      !==========================================================

      allocate( Index_Element_2_Node_Neumann_BC_tmp( 3, Number_Node ) )

      Number_Element_Neumann_BC= 0

      ! No Parallel
      do i= 1, Number_Element
         if( Flag_Element_Neumann_BC( i )==Flag_On )then
            Number_Element_Neumann_BC= Number_Element_Neumann_BC +1
            do j= 1, 3
               Index_Element_2_Node_Neumann_BC_tmp( j, Number_Element_Neumann_BC )= Index_Element_2_Node_Triangle( j, i )
            end do
         end if
      end do

      allocate( Index_Element_2_Node_Neumann_BC( 3, Number_Element_Neumann_BC ) )

      do i= 1, Number_Element_Neumann_BC
         do j= 1, 3
            Index_Element_2_Node_Neumann_BC( j, i )= Index_Element_2_Node_Neumann_BC_tmp( j, i )
         end do
      end do

      deallocate( Index_Element_2_Node_Neumann_BC_tmp )

      !==========================================================
      write(*,*) '            Implement Neumann B.C. to GlobalMatrix '
      !==========================================================















      deallocate( Index_Element_2_Node_Neumann_BC )

      !==========================================================================
      write(*,*)'            NodeNumber_Group_tmp --> NodeNumber_Group '
      !==========================================================================

   return
end subroutine Implement_Neumann_BC_complex_CSRF 



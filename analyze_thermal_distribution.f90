subroutine analyze_thermal_distribution(LocalMatrix_Triangle, Index_Element_2_Node_Triangle, Class_Element_Triangle, &
           Flag_Thermal_Insulation_BC, Number_Node, Number_Element_Triangle, Width_Matrix_LHS, &
           position_node_x, position_node_y, t_boundary_higher_side, t_boundary_lower_side, &
           max_node_position_x,max_node_position_y,min_node_position_x,min_node_position_y, &
           !=======================================================================================================
           Temperature_Solution)
   implicit none

   integer, intent(in) :: Number_Node
   integer, intent(in) :: Number_Element_Triangle 
   integer, intent(in) :: Width_Matrix_LHS
   integer, intent(in) :: Flag_Thermal_Insulation_BC
   real(8), intent(in) :: t_boundary_higher_side, t_boundary_lower_side
   real(8), intent(in) :: max_node_position_x,max_node_position_y,min_node_position_x,min_node_position_y   

   integer, intent(in) :: Index_Element_2_Node_Triangle( 3, Number_Element_Triangle ) 
   integer, intent(in) :: Class_Element_Triangle( Number_Element_Triangle )
   real(8), intent(in) :: LocalMatrix_Triangle(3,3,Number_Element_Triangle)
   real(8), intent(in) :: Position_Node_X(number_node), Position_Node_Y(number_node)

   real(8), intent(out) :: Temperature_Solution(Number_Node)

   integer, parameter :: Flag_On= 1
   integer, parameter :: Flag_Off= 0
   
   integer, parameter :: Flag_Symmetric= Flag_On
   integer, parameter :: Flag_Matrix_Upper= Flag_On
   integer, parameter :: ID_Element_Thermal_Insulator= 2

   integer :: Flag_Remove_Element
   integer :: Number_nonzero

   integer :: i, j, k, e
   integer :: counter
   character filename*128

   integer, allocatable :: Flag_Element_Triangle_Removed(:),Flag_Element_Neumann_BC(:),Flag_Node_Dirichlet_BC(:)

   real(8), allocatable :: GlobalMatrix(:,:), J_GlobalMatrix(:,:)
   real(8), allocatable :: Global_Vector_RHS(:),Nodal_Value_on_Dirichlet_BC(:),Nodal_Value_Fixed_in_PEC(:)
   real(8), allocatable :: aa_zmumps(:), ia_zmumps(:), ja_zmumps(:)

   Flag_Remove_Element= Flag_Thermal_Insulation_BC!中心断熱領域を解析に加わる場合（リバーサル、コンセントレーターなど）は0,他は1   

   allocate( Flag_Element_Triangle_Removed( Number_Element_Triangle ) ) 

   Flag_Element_Triangle_Removed= Flag_Off

   allocate( GlobalMatrix( Number_Node, Width_Matrix_LHS ) ) 
   allocate( J_GlobalMatrix( Number_Node, Width_Matrix_LHS ) )

   call Format_Local_Global_DMUMPS&
   ( LocalMatrix_Triangle, Index_Element_2_Node_Triangle, &
     3, Number_Node, Number_Element_Triangle, &
     1, Width_Matrix_LHS, Flag_Symmetric, &
     Class_Element_Triangle, Flag_Element_Triangle_Removed, Flag_Remove_Element, &
     !=====================================================
     GlobalMatrix, J_GlobalMatrix, Number_NonZero )

   deallocate( Flag_Element_Triangle_Removed ) 

   !check GlobalMatrix
   open(unit=10, file='G_Matrix_MA57.txt', status='replace')
   ! Write the array to the file
   do i = 1, Number_Node
      write(10, '(F12.5)', advance='no') GlobalMatrix(i, 1)
      do j = 2, Width_Matrix_LHS
         write(10, '(F12.5)', advance='no') GlobalMatrix(i, j)
      end do
      write(10, *)
   end do 
   ! Close the file
   close(10)

   !check J_GlobalMatrix
   open(unit=11, file='J_Matrix.txt', status='replace')
   ! Write the array to the file
   do i = 1, Number_Node
      write(11, '(I10)', advance='no') J_GlobalMatrix(i, 1)
      do j = 2, Width_Matrix_LHS
       write(11, '(I10)', advance='no') J_GlobalMatrix(i, j)
      end do
      write(11, *)
   end do 
   ! Close the file
   close(11)

   allocate( Global_Vector_RHS( Number_Node ) )
   Global_Vector_RHS= 0.0d0

   !set Neumann BC
   allocate(Flag_Element_Neumann_BC( number_node ))
   Flag_Element_Neumann_BC=Flag_Off

   do i= 1, number_node
      if( Position_Node_Y( i ) == min_node_position_y )then
         Flag_Element_Neumann_BC(i)= Flag_On
      else if( Position_Node_Y( i ) == max_node_position_y )then 
         Flag_Element_Neumann_BC(i)= Flag_On 
      end if     
   end do

   call Implement_Neumann_BC_complex_CSRF &
          ( Number_Element_Triangle, Flag_Element_Neumann_BC, Index_Element_2_Node_Triangle, &
            Number_Node, Width_Matrix_LHS, &
            !================================================================
            GlobalMatrix, J_GlobalMatrix, Global_Vector_RHS )

   deallocate(Flag_Element_Neumann_BC)
   
   !set Dirichlet BC
   allocate( Nodal_Value_on_Dirichlet_BC( Number_Node ) )
   allocate( Flag_Node_Dirichlet_BC( Number_Node ) )

   Nodal_Value_on_Dirichlet_BC= 0.0d0
   Flag_Node_Dirichlet_BC= Flag_Off
   
   do i= 1, number_node
      if( Position_Node_X(i) == min_node_position_x)then
         Flag_Node_Dirichlet_BC(i)= Flag_On
         Nodal_Value_on_Dirichlet_BC(i)= t_boundary_lower_side
      else if( Position_Node_X(i) == max_node_position_x )then 
         Flag_Node_Dirichlet_BC(i)= Flag_On
         Nodal_Value_on_Dirichlet_BC(i)= t_boundary_higher_side
      end if     
   end do

   call Implement_Dirichlet_BC_real_CSRF &
      ( Nodal_Value_on_Dirichlet_BC, Flag_Node_Dirichlet_BC, &
        Number_Node, Width_Matrix_LHS, &
        !======================================================
        GlobalMatrix, J_GlobalMatrix, Global_Vector_RHS )

   deallocate( Flag_Node_Dirichlet_BC )
   deallocate( Nodal_Value_on_Dirichlet_BC )

   !==============================================================================================================
   if( Flag_Thermal_Insulation_BC==Flag_On )then
   !==============================================================================================================
   
      allocate( Flag_Node_Dirichlet_BC( Number_Node ) )
   
      do i= 1, Number_Node
         Flag_Node_Dirichlet_BC( i )= Flag_Off
      end do
   
      do e= 1, Number_Element_Triangle
         if( Class_Element_Triangle( e ) == ID_Element_Thermal_Insulator )then
            do i= 1, 3
               Flag_Node_Dirichlet_BC( Index_Element_2_Node_Triangle( i,e ) )= Flag_On
            end do
         end if
      end do
 
      ! Flag_MultiPhysics_Device==20/30/40 : No insulatorの場合下記do loopは省略
         do e= 1, Number_Element_Triangle
            if( Class_Element_Triangle( e ) /= ID_Element_Thermal_Insulator )then
               do i= 1, 3
                  Flag_Node_Dirichlet_BC( Index_Element_2_Node_Triangle( i,e ) )= Flag_Off
               end do
            end if
         end do
   
      Counter= 0
      do i= 1, Number_Node
         if( Flag_Node_Dirichlet_BC( i )==Flag_On )then
            Counter= Counter +1
         end if
      end do
   
      allocate( Nodal_Value_Fixed_in_PEC( Number_Node ) )
    
      do i= 1, Number_Node 
         Nodal_Value_Fixed_in_PEC( i )= 0d0 
      end do
   
      if( Counter > 0 )then
         !call Implement_Dirichlet_BC_complex_CSRF &
         call Implement_Dirichlet_BC_real_CSRF &
            ( Nodal_Value_Fixed_in_PEC, Flag_Node_Dirichlet_BC, &
              Number_Node, Width_Matrix_LHS, &
             !========================================================
              GlobalMatrix, J_GlobalMatrix, Global_Vector_RHS )
      end if
   
      deallocate( Nodal_Value_Fixed_in_PEC )
      deallocate( Flag_Node_Dirichlet_BC )

   end if   

   !count Numbver_nonzero
   Number_NonZero= 0

   do j= 1, Width_Matrix_LHS
      do i= 1, Number_Node
         if( Flag_Symmetric==Flag_On )then
            if( i >= J_GlobalMatrix( i, j ) .and. J_GlobalMatrix( i, j ) > 0 )then
               Number_NonZero= Number_NonZero +1
            end if
         else if( Flag_Symmetric==Flag_Off )then
            if( 0 < J_GlobalMatrix( i, j ) .and. J_GlobalMatrix( i, j ) <= Number_Node  )then
               Number_NonZero= Number_NonZero +1
            end if
         else
            write(*,*)'Flag_Symmetric=', Flag_Symmetric
            call Output_Error( 'Analyze_Steady_State_Heat_Conduction', 1118 )
         end if
      end do
   end do

   write(*,*)'         Number_NonZero=', Number_NonZero
   write(*,*)'         Number_Node=',  Number_Node 
   write(*,*)'         Width_Matrix_LHS=', Width_Matrix_LHS

   !==================================================================================================
   write(*,*)'         Global Matrix --> MUMPS Format '
   !==================================================================================================


   allocate( aa_zmumps( Number_NonZero ) )
   allocate( ia_zmumps( Number_NonZero ) )
   allocate( ja_zmumps( Number_NonZero ) )

   !call Format_Global_aa_ZMUMPS&
   call Format_Global_aa_DMUMPS&
      ( GlobalMatrix, J_GlobalMatrix, & 
        Number_Node, Width_Matrix_LHS, Number_NonZero, Flag_Symmetric, Flag_Matrix_Upper, & 
        !==================================================== 
        aa_zmumps, ia_zmumps, ja_zmumps )    

   deallocate( GlobalMatrix )
   deallocate( J_GlobalMatrix )

   !==================================================================================================
   write(*,*)'         Solve Linear Equation'
   write(*,*)'            Number_Node=', Number_Node
   write(*,*)'            Number_NonZero=', Number_NonZero
   !==================================================================================================

   call MA57( Number_Node, Number_NonZero, aa_zmumps, ia_zmumps, ja_zmumps, Global_Vector_RHS, Temperature_Solution)   

   deallocate( aa_zmumps )
   deallocate( ia_zmumps )
   deallocate( ja_zmumps )

   !==================================================================================================
   write(*,*)'         end subroutine analyze_thermal_distribution'
   !==================================================================================================
   return
end subroutine analyze_thermal_distribution

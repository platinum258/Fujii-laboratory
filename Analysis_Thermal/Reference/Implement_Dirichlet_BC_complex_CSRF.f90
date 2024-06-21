
subroutine Implement_Dirichlet_BC_complex_CSRF &
       ( Nodal_Value_on_Dirichlet_BC_All, Flag_Node_Dirichlet_BC, &
         Number_Node, Width_Matrix_LHS, &
         !================================================================
         GlobalMatrix, J_GlobalMatrix, Global_Vector_RHS )
    
   !$use omp_lib
   implicit none

   integer, intent(in) :: Number_Node, Width_Matrix_LHS

   complex( kind( 0d0 ) ), intent(inout) :: GlobalMatrix( Number_Node, Width_Matrix_LHS )
   integer, intent(inout) :: J_GlobalMatrix( Number_Node, Width_Matrix_LHS )
   complex( kind( 0d0 ) ), intent(inout) :: Global_Vector_RHS( Number_Node )
   integer, intent(in) :: Flag_Node_Dirichlet_BC( Number_Node )
   complex( kind( 0d0 ) ), intent(in) :: Nodal_Value_on_Dirichlet_BC_All( Number_Node )

   complex( kind( 0d0 ) ), parameter :: Zero= dcmplx( 0.0d0, 0.0d0 )
   complex( kind( 0d0 ) ), parameter :: One= dcmplx( 1.0d0, 0.0d0 )

   integer :: i, j, k, Counter_Percent

   complex( kind( 0d0 ) ), allocatable, dimension(:) :: Global_Vector_RHS_PEC

   integer :: Number_Group_Node_Dirichlet_BC
   integer, allocatable, dimension(:,:) :: NodeNumber_Group, NodeNumber_Group_tmp

   integer :: Number_Node_Dirichlet_BC
   integer, allocatable, dimension(:) :: NodeNumber_Dirichlet_BC, NodeNumber_Dirichlet_BC_tmp
   complex( kind( 0d0 ) ), allocatable, dimension(:) :: Nodal_Value_on_Dirichlet_BC, Nodal_Value_on_Dirichlet_BC_tmp

   integer, allocatable, dimension(:) :: Number_Node_Group_Dirichlet_BC 

   integer, parameter :: Flag_On= 1
   integer, parameter :: Flag_Off= 0

   !================================================================================================================
   write(*,*)'            call Implement_Dirichlet_BC_complex_CSRF'
   !================================================================================================================

      !===============================================
      write(*,*) '            Check Matrix'
      !===============================================

      !$omp parallel do default( none ) &
      !$omp private( i ) &
      !$omp shared( Number_Node, Flag_Node_Dirichlet_BC ) 
      do i= 1, Number_Node
         if( Flag_Node_Dirichlet_BC( i )/=Flag_On .and. &
             Flag_Node_Dirichlet_BC( i )/=Flag_Off )then
            call Output_Error( 'Implement_Dirichlet_BC_complex_CSRF' , 46 )
         end if
      end do

      !==========================================================
      write(*,*) '            Create Dirichlet B.C. Data'
      !==========================================================


      allocate( Nodal_Value_on_Dirichlet_BC_tmp( Number_Node ) )
      allocate( NodeNumber_Dirichlet_BC_tmp( Number_Node ) )


      Number_Node_Dirichlet_BC= 0

      ! No Parallel
      do i= 1, Number_Node
         if( Flag_Node_Dirichlet_BC( i )==Flag_On )then
            Number_Node_Dirichlet_BC= Number_Node_Dirichlet_BC +1
            Nodal_Value_on_Dirichlet_BC_tmp( Number_Node_Dirichlet_BC )= Nodal_Value_on_Dirichlet_BC_All( i )
            NodeNumber_Dirichlet_BC_tmp( Number_Node_Dirichlet_BC )= i
         end if
      end do

      allocate( Nodal_Value_on_Dirichlet_BC( Number_Node_Dirichlet_BC ) )
      allocate( NodeNumber_Dirichlet_BC( Number_Node_Dirichlet_BC ) )

      do i= 1, Number_Node_Dirichlet_BC
         Nodal_Value_on_Dirichlet_BC( i )= Nodal_Value_on_Dirichlet_BC_tmp( i )
         NodeNumber_Dirichlet_BC( i )= NodeNumber_Dirichlet_BC_tmp( i )
      end do

      deallocate( Nodal_Value_on_Dirichlet_BC_tmp )
      deallocate( NodeNumber_Dirichlet_BC_tmp )


      !==========================================================
      write(*,*) '            Create Group of Node Data on Dirichlet B.C.'
      !==========================================================

      allocate( NodeNumber_Group_tmp( 2, Number_Node ) )

      Number_Group_Node_Dirichlet_BC= 1
      NodeNumber_Group_tmp( 1, Number_Group_Node_Dirichlet_BC )= NodeNumber_Dirichlet_BC( 1 )

      do i= 1, Number_Node_Dirichlet_BC-1
         if( NodeNumber_Dirichlet_BC( i+1 )-NodeNumber_Dirichlet_BC( i )/=1 )then

            NodeNumber_Group_tmp( 2, Number_Group_Node_Dirichlet_BC ) &
            = NodeNumber_Dirichlet_BC( i )

            Number_Group_Node_Dirichlet_BC= Number_Group_Node_Dirichlet_BC +1

            NodeNumber_Group_tmp( 1, Number_Group_Node_Dirichlet_BC ) &
            = NodeNumber_Dirichlet_BC( i+1 )
         end if
      end do

      NodeNumber_Group_tmp( 2, Number_Group_Node_Dirichlet_BC ) &
      = NodeNumber_Dirichlet_BC( Number_Node_Dirichlet_BC )

      deallocate( NodeNumber_Dirichlet_BC )

      !==========================================================================
      write(*,*)'            NodeNumber_Group_tmp --> NodeNumber_Group '
      !==========================================================================

      allocate( NodeNumber_Group( 2, Number_Group_Node_Dirichlet_BC ) )

      do i= 1, Number_Group_Node_Dirichlet_BC
         do j= 1, 2
            NodeNumber_Group( j, i )= NodeNumber_Group_tmp( j, i )
         end do
      end do

      deallocate( NodeNumber_Group_tmp )

      allocate( Number_Node_Group_Dirichlet_BC( Number_Group_Node_Dirichlet_BC ) )

      do i= 1, Number_Group_Node_Dirichlet_BC
         Number_Node_Group_Dirichlet_BC( i ) &
         = NodeNumber_Group( 2, i ) -NodeNumber_Group( 1, i ) +1
      end do

      deallocate( Nodal_Value_on_Dirichlet_BC )
      deallocate( Number_Node_Group_Dirichlet_BC )

      Counter_Percent= 1
      !===============================================
      write(*,*) '            ================================================='
      write(*,*) '            Number_Group_Node_Dirichlet_BC=', Number_Group_Node_Dirichlet_BC
      write(*,*) '            ================================================='
      !===============================================
      do k= 1, Number_Group_Node_Dirichlet_BC
    
         if( k== Int(Number_Group_Node_Dirichlet_BC*Counter_Percent/100d0) )then
            write(*,'(a,i3,a,$)')'[', Counter_Percent, '%]'
            !write(*,'(a,i3,a,$,i3)')'[', Counter_Percent, '%]', k 
            Counter_Percent= Counter_Percent +1
         end if
 
         !===============================================
         !write(*,*) '            Right Hand Side'
         !===============================================
         do i= NodeNumber_Group( 1, k ), NodeNumber_Group( 2, k ), 1
            Global_Vector_RHS( i )= Nodal_Value_on_Dirichlet_BC_All( i ) 
         end do
   
         allocate( Global_Vector_RHS_PEC( Number_Node ) )
   
         do i= 1, Number_Node, 1 
            Global_Vector_RHS_PEC( i )= Zero
         end do

         !$omp parallel do default( none ) &
         !$omp private( i, j ) &
         !$omp shared( k, NodeNumber_Group, Width_Matrix_LHS ) &
         !$omp shared( GlobalMatrix, J_GlobalMatrix, Global_Vector_RHS_PEC ) &
         !$omp shared( Nodal_Value_on_Dirichlet_BC_All )
         do i= 1, NodeNumber_Group( 1, k ) -1, 1
            do j= 1, Width_Matrix_LHS 
               if( J_GlobalMatrix( i, j ) >= NodeNumber_Group( 1, k ) .and. & 
                   J_GlobalMatrix( i, j ) <= NodeNumber_Group( 2, k ) )then
                  Global_Vector_RHS_PEC( i ) &
                  = Global_Vector_RHS_PEC( i ) &
                   +GlobalMatrix( i, j )*( Nodal_Value_on_Dirichlet_BC_All( J_GlobalMatrix( i, j ) ) )
               end if 
            end do
         end do
    
         !$omp parallel do default( none ) &
         !$omp private( i, j ) &
         !$omp shared( k, NodeNumber_Group, Number_Node, Width_Matrix_LHS ) &
         !$omp shared( GlobalMatrix, J_GlobalMatrix, Global_Vector_RHS_PEC ) &
         !$omp shared( Nodal_Value_on_Dirichlet_BC_All )
         do i= NodeNumber_Group( 2, k ) +1, Number_Node, 1 
            do j= 1, Width_Matrix_LHS 
               if( J_GlobalMatrix( i, j ) >= NodeNumber_Group( 1, k ) .and. & 
                   J_GlobalMatrix( i, j ) <= NodeNumber_Group( 2, k ) )then
                  Global_Vector_RHS_PEC( i ) &
                  = Global_Vector_RHS_PEC( i ) &
                   +GlobalMatrix( i, j )*( Nodal_Value_on_Dirichlet_BC_All( J_GlobalMatrix( i, j ) ) )
               end if 
            end do
         end do
   
         do i= 1, Number_Node, 1 
            Global_Vector_RHS( i )= Global_Vector_RHS( i ) -Global_Vector_RHS_PEC( i )
         end do
   
         deallocate( Global_Vector_RHS_PEC )
   
         !===============================================
         !write(*,*) '            Left Hand Side'
         !===============================================
   
         !$omp parallel do default( none ) &
         !$omp private( i, j ) &
         !$omp shared( NodeNumber_Group, Width_Matrix_LHS, k ) & 
         !$omp shared( GlobalMatrix, J_GlobalMatrix )
         do j= 1, Width_Matrix_LHS, 1
            do i= NodeNumber_Group( 1, k ), NodeNumber_Group( 2, k ) 
               GlobalMatrix( i, j )= Zero 
               J_GlobalMatrix( i, j )= 0 
            end do
         end do
   
         !$omp parallel do default( none ) &
         !$omp private( i, j ) &
         !$omp shared( NodeNumber_Group, Width_Matrix_LHS, k ) & 
         !$omp shared( GlobalMatrix, J_GlobalMatrix, Number_Node )
         do j= 1, Width_Matrix_LHS, 1
            do i= 1, NodeNumber_Group( 1, k ) -1, 1  
               if( J_GlobalMatrix( i, j ) >= NodeNumber_Group( 1, k ) .and. & 
                   J_GlobalMatrix( i, j ) <= NodeNumber_Group( 2, k ) )then
   
                  GlobalMatrix( i, j )= Zero 
                  J_GlobalMatrix( i, j )= 0 
               end if
            end do

            do i= NodeNumber_Group( 2, k ) +1, Number_Node 
               if( J_GlobalMatrix( i, j ) >= NodeNumber_Group( 1, k ) .and. & 
                   J_GlobalMatrix( i, j ) <= NodeNumber_Group( 2, k ) )then
   
                  GlobalMatrix( i, j )= Zero 
                  J_GlobalMatrix( i, j )= 0 
               end if
            end do
         end do
    
         do i= NodeNumber_Group( 1, k ), NodeNumber_Group( 2, k ) 
            GlobalMatrix( i, 1 )= One 
            J_GlobalMatrix( i, 1 )= i 
         end do
   
         !====================================================
         !write(*,*) '            Check Left Hand Side'
         !====================================================
   
         !$omp parallel do default( none ) &
         !$omp private( i, j ) &
         !$omp shared( NodeNumber_Group, Width_Matrix_LHS, k ) & 
         !$omp shared( GlobalMatrix, J_GlobalMatrix, Number_Node )
         do j= 1, Width_Matrix_LHS
            do i= 1, Number_Node
 
               if( ( J_GlobalMatrix( i, j ) >= NodeNumber_Group( 1, k ) .and. &
                   J_GlobalMatrix( i, j ) <= NodeNumber_Group( 2, k ) ) & 
                   .or. &
                   ( i >= NodeNumber_Group( 1, k ) .and. &
                   i <= NodeNumber_Group( 2, k ) ) )then

                  if( i==J_GlobalMatrix( i, j ) )then 
                     if( GlobalMatrix( i, j )/= One )then
                        write(*,*)'GlobalMatrix( i, j )=', GlobalMatrix( i, j )
                        write(*,*)'i=', i
                        write(*,*)'J_GlobalMatrix( i, j )=', J_GlobalMatrix( i, j )
                        write(*,*)'NodeNumber_Group( 1, k )=', NodeNumber_Group( 1, k )
                        write(*,*)'NodeNumber_Group( 2, k )=', NodeNumber_Group( 2, k )
                        call Output_Error( 'Implement_Dirichlet_BC_complex_CSRF' , 248 )
                     end if

                  else if(  GlobalMatrix( i, j )/= Zero )then
                     write(*,*)'GlobalMatrix( i, j )=', GlobalMatrix( i, j )
                     write(*,*)'i=', i
                     write(*,*)'J_GlobalMatrix( i, j )=', J_GlobalMatrix( i, j )
                     write(*,*)'NodeNumber_Group( 1, k )=', NodeNumber_Group( 1, k )
                     write(*,*)'NodeNumber_Group( 2, k )=', NodeNumber_Group( 2, k )
                     call Output_Error( 'Implement_Dirichlet_BC_complex_CSRF' , 254 )
                  end if
                  
               end if
            end do
         end do

      end do ! k

      deallocate( NodeNumber_Group )
   
   return
end subroutine Implement_Dirichlet_BC_complex_CSRF 



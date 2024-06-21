
subroutine Change_Direction_Vector_for_Unit_Normal_Vector &
       ( Number_Edge, Position_Node_Neumann_BC, Flag_Direction_Vector )

   !use omp_lib
   implicit none   

   integer, intent(in) :: Number_Edge, Flag_Direction_Vector 
   double precision, intent(inout) :: Position_Node_Neumann_BC( 2, 2, Number_Edge )

   integer :: i, j 
   double precision, allocatable, dimension(:,:) :: Center_Position_Edge_Neumann_BC
   double precision, allocatable, dimension(:,:) :: Double_Precision_Vector_tmp 
   double precision, allocatable, dimension(:) :: Outer_Product

   !================================================================================
   write(*,*)'         call Change_Direction_Vector_for_Unit_Normal_Vector'
   !================================================================================

   allocate( Center_Position_Edge_Neumann_BC( 2, Number_Edge ) )

   do i= 1, Number_Edge 
      do j= 1, 2
         Center_Position_Edge_Neumann_BC( j, i )= ( Position_Node_Neumann_BC( j, 1, i ) +Position_Node_Neumann_BC( j, 2, i ) )/2d0
      end do
   end do

   !================================================================================
   write(*,*)'            Check Original Vector'
   !================================================================================
   allocate( Double_Precision_Vector_tmp( 2, Number_Edge ) )
   allocate( Outer_Product( Number_Edge ) )

   do i= 1, Number_Edge 

      Outer_Product( i ) &
      = ( Position_Node_Neumann_BC( 1, 2, i ) -Position_Node_Neumann_BC( 1, 1, i ) ) &
       *Center_Position_Edge_Neumann_BC( 2, i ) & 
       -( Position_Node_Neumann_BC( 2, 2, i ) -Position_Node_Neumann_BC( 2, 1, i ) ) &
       *Center_Position_Edge_Neumann_BC( 1, i ) 

      if( Outer_Product( 1 )*Outer_Product( i ) < 0d0 )then
         do j= 1, 2
            Double_Precision_Vector_tmp( j, i ) &
            = Position_Node_Neumann_BC( j, 1, i ) 
         end do
         
         do j= 1, 2
            Position_Node_Neumann_BC( j, 1, i ) &
            = Position_Node_Neumann_BC( j, 2, i ) 
         end do
         
         do j= 1, 2
            Position_Node_Neumann_BC( j, 2, i ) &
            = Double_Precision_Vector_tmp( j, i ) 
         end do
      end if
   end do

   if( Flag_Direction_Vector==-1 )then
      do i= 1, Number_Edge 
         do j= 1, 2
            Double_Precision_Vector_tmp( j, i ) &
            = Position_Node_Neumann_BC( j, 1, i ) 
         end do
         
         do j= 1, 2
            Position_Node_Neumann_BC( j, 1, i ) &
            = Position_Node_Neumann_BC( j, 2, i ) 
         end do
         
         do j= 1, 2
            Position_Node_Neumann_BC( j, 2, i ) &
            = Double_Precision_Vector_tmp( j, i ) 
         end do
      end do
   end if

   deallocate( Double_Precision_Vector_tmp )
   deallocate( Outer_Product )
   deallocate( Center_Position_Edge_Neumann_BC )
   !================================================================================
   write(*,*)'            end Change_Direction_Vector_for_Unit_Normal_Vector'
   !================================================================================

   return
end subroutine Change_Direction_Vector_for_Unit_Normal_Vector


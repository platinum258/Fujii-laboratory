
subroutine Compute_Unit_Normal_Vector &
       ( Number_Edge, Position_Node_Neumann_BC, Unit_Normal_Vector, Flag_Output_Unit_Normal_Vector )

   !use omp_lib
   implicit none   

   integer, intent(in) :: Number_Edge, Flag_Output_Unit_Normal_Vector
   double precision, intent(in) :: Position_Node_Neumann_BC( 2, 2, Number_Edge )

   double precision, intent(out) :: Unit_Normal_Vector( 2, Number_Edge )

   integer :: i, j 
   double precision, allocatable, dimension(:,:) :: Original_Vector, Center_Position_Edge_Neumann_BC
   double precision, allocatable, dimension(:,:) :: Double_Precision_Vector_tmp 
   double precision, allocatable, dimension(:) :: Length, Length_Evaluate, Outer_Product

   integer, parameter :: Flag_On= 1
   integer, parameter :: Flag_Off= 0
   double precision, parameter :: Tolerance_Length= 1d-8

   !================================================================================
   write(*,*)'         call Compute_Unit_Normal_Vector'
   !================================================================================

   allocate( Center_Position_Edge_Neumann_BC( 2, Number_Edge ) )
   allocate( Original_Vector( 2, Number_Edge ) )

   do i= 1, Number_Edge 
      do j= 1, 2
         Center_Position_Edge_Neumann_BC( j, i )= ( Position_Node_Neumann_BC( j, 1, i ) +Position_Node_Neumann_BC( j, 2, i ) )/2d0
      end do
      do j= 1, 2
         Original_Vector( j, i )= Position_Node_Neumann_BC( j, 2, i ) -Position_Node_Neumann_BC( j, 1, i )
      end do
   end do

   !================================================================================
   write(*,*)'            Compute Unit Normal Vector'
   !================================================================================
   allocate( Length( Number_Edge ) )
   allocate( Length_Evaluate( Number_Edge ) )

   do i= 1, Number_Edge
   
      Length( i )= dsqrt( Original_Vector( 1, i )*Original_Vector( 1, i ) +Original_Vector( 2, i )*Original_Vector( 2, i ) )
   
      Unit_Normal_Vector( 1, i )=  Original_Vector( 2, i )/Length( i ) 
      Unit_Normal_Vector( 2, i )= -Original_Vector( 1, i )/Length( i )
   
      Length_Evaluate( i )= sqrt( Unit_Normal_Vector( 1, i )*Unit_Normal_Vector( 1, i ) +Unit_Normal_Vector( 2, i )*Unit_Normal_Vector( 2, i ) )
   
   end do

   if( Flag_Output_Unit_Normal_Vector==Flag_On )then
   !================================================================================
   write(*,*)'            Ouput Original Vector & Unit Normal Vector'
   !================================================================================
      open( 104, file='Unit_Normal_Vector.dat', status='replace' ) 
         do i= 1, Number_Edge
            write(104,112) & 
            Center_Position_Edge_Neumann_BC( 1, i ), Center_Position_Edge_Neumann_BC( 2, i ), & 
            Unit_Normal_Vector( 1, i ), Unit_Normal_Vector( 2, i )
         end do
      112 format(es15.8,1x,es15.8,1x,es15.8,1x,es15.8)
      close( 104 )

      open( 105, file='Original_Vector_Neumann_BC.dat', status='replace' ) 
         do i= 1, Number_Edge
            write(105,113) & 
            Position_Node_Neumann_BC( 1, 1, i ), Position_Node_Neumann_BC( 2, 1, i ), &
            Original_Vector( 1, i ), Original_Vector( 2, i )
         end do
      113 format(es15.8,1x,es15.8,1x,es15.8,1x,es15.8)
      close( 105 )
   end if

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
         !do j= 1, 2
         !   Double_Precision_Vector_tmp( j, i ) &
         !   = Position_Node_Neumann_BC( j, 1, i ) 
         !end do
         !
         !do j= 1, 2
         !   Position_Node_Neumann_BC( j, 1, i ) &
         !   = Position_Node_Neumann_BC( j, 2, i ) 
         !end do
         !
         !do j= 1, 2
         !   Position_Node_Neumann_BC( j, 2, i ) &
         !   = Double_Precision_Vector_tmp( j, i ) 
         !end do
         call Output_Error( 'Compute_Unit_Normal_Vector', 188 )
      end if
   end do

   deallocate( Double_Precision_Vector_tmp )
   deallocate( Outer_Product )
   deallocate( Center_Position_Edge_Neumann_BC )
   deallocate( Original_Vector )
   deallocate( Length )

   do i= 1, Number_Edge
      if( 1d0 -Tolerance_Length > Length_Evaluate( i ) .or. Length_Evaluate( i ) > 1d0 +Tolerance_Length )then
          write(*,*)'Error : Compute_Unit_Normal_Vector.f90 38'
          write(*,*)'Length/=1d0'
          write(*,*)'Length=', Length_Evaluate( i )
          write(*,*)'Range : ', 1d0 -Tolerance_Length, '--', 1d0 +Tolerance_Length
          stop
      end if
   end do

   deallocate( Length_Evaluate )

   !================================================================================
   write(*,*)'            end Compute_Unit_Normal_Vector'
   !================================================================================

   return
end subroutine Compute_Unit_Normal_Vector


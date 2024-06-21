
subroutine Compute_Line_Integral_Poynting_Vector &
       ( Number_Node, Position_Node, &
         Number_Element, Shape_Element, Index_Element_2_Node, &
         Number_Element_on_Boundary_1, Element_Number_Boundary_1, Local_Node_Number_1, &
         Number_Element_on_Boundary_2, Element_Number_Boundary_2, Local_Node_Number_2, &
         Value_Node_Integrated, &
       !====================================================================================================
         Value_Integrated_Poynting_Vector )

   implicit none

   integer, intent(in) :: Number_Node, Number_Element, Shape_Element 
   integer, intent(in) :: Number_Element_on_Boundary_1, Number_Element_on_Boundary_2

   double precision, intent(in) :: Position_Node( 2, Number_Node ) 
   integer, intent(in) :: Index_Element_2_Node( Shape_Element, Number_Element )

   integer, intent(in) :: Element_Number_Boundary_1( Number_Element_on_Boundary_1 ), &
                  Element_Number_Boundary_2( Number_Element_on_Boundary_2 )
   integer, intent(in) :: Local_Node_Number_1( 2, Number_Element_on_Boundary_1 ), &
                  Local_Node_Number_2( 2, Number_Element_on_Boundary_1 )

   complex( kind(0d0) ), intent(in) :: Value_Node_Integrated( Number_Node )

   double precision, intent(out) :: Value_Integrated_Poynting_Vector 

   integer :: e, i, j 
   double precision, allocatable, dimension(:,:,:) :: Position_Node_on_Boundary_1, Position_Node_on_Boundary_2 
   double precision, allocatable, dimension(:) :: Edge_Length_1, Edge_Length_2
   complex( kind(0d0) ), allocatable, dimension(:,:) :: Value_Node_in_Element_on_Boundary_1, Value_Node_in_Element_on_Boundary_2 

   double precision, allocatable, dimension(:,:) :: Unit_Normal_Vector_1, Unit_Normal_Vector_2

   double precision, allocatable, dimension(:,:) :: BasisFunction_A, BasisFunction_B, BasisFunction_C
   double precision, allocatable, dimension(:,:) :: BasisFunction_A_Boundary_1, BasisFunction_B_Boundary_1, BasisFunction_C_Boundary_1 
   double precision, allocatable, dimension(:,:) :: BasisFunction_A_Boundary_2, BasisFunction_B_Boundary_2, BasisFunction_C_Boundary_2 

   complex( kind(0d0) ), allocatable, dimension(:,:) :: Term_Unit_Normal_Vector_1, Term_Unit_Normal_Vector_2 
   complex( kind(0d0) ), allocatable, dimension(:,:) :: Term_Line_Integral_1, Term_Line_Integral_2 
   complex( kind(0d0) ), allocatable, dimension(:) :: Term_Unit_Normal_Vector_1_Summed, Term_Unit_Normal_Vector_2_Summed
   complex( kind(0d0) ), allocatable, dimension(:) :: Term_Line_Integral_1_Summed, Term_Line_Integral_2_Summed

   double precision :: Value_Integrated_Poynting_Vector_1, Value_Integrated_Poynting_Vector_2 

   complex( kind(0d0) ), allocatable, dimension(:,:) :: Poynting_Vector_Nabla_Factor_1, Poynting_Vector_Nabla_Factor_2
   complex( kind(0d0) ), allocatable, dimension(:) :: Poynting_Vector_Factor_1, Poynting_Vector_Factor_2
   double precision, allocatable, dimension(:,:) :: Poynting_Vector, Poynting_Vector_1, Poynting_Vector_2

   integer, parameter :: Flag_On= 1
   integer, parameter :: Flag_Off= 0
   complex( kind(0d0) ), parameter :: Zero= dcmplx( 0d0, 0d0 )
   complex( kind(0d0) ), parameter :: Complex_Unit= dcmplx( 0d0, 1d0 )

   !====================================================================================================
   write(*,*)' call Compute_Line_Integral_Poynting_Vector '
   !====================================================================================================
 
   !====================================================================================================
   write(*,*)'    Imput_Data --> Computing_Data '
   !====================================================================================================

   allocate( Position_Node_on_Boundary_1( 2, 2, Number_Element_on_Boundary_1 ) )
   allocate( Position_Node_on_Boundary_2( 2, 2, Number_Element_on_Boundary_2 ) )

   do e= 1, Number_Element_on_Boundary_1
      do i= 1, 2
         do j= 1, 2
            Position_Node_on_Boundary_1( j, i, e )= Position_Node( j, Index_Element_2_Node( Local_Node_Number_1( i, e ), Element_Number_Boundary_1( e ) ) )
         end do
      end do
   end do

   do e= 1, Number_Element_on_Boundary_2
      do i= 1, 2
         do j= 1, 2
            Position_Node_on_Boundary_2( j, i, e )= Position_Node( j, Index_Element_2_Node( Local_Node_Number_2( i, e ), Element_Number_Boundary_2( e ) ) )
         end do
      end do
   end do
 
   allocate( Edge_Length_1( Number_Element_on_Boundary_1 ) )
   allocate( Edge_Length_2( Number_Element_on_Boundary_2 ) )

   do e= 1, Number_Element_on_Boundary_1
      Edge_Length_1( e )= sqrt( ( Position_Node_on_Boundary_1( 1, 1, e ) -Position_Node_on_Boundary_1( 1, 2, e ) )**2 & 
                     +( Position_Node_on_Boundary_1( 2, 1, e ) -Position_Node_on_Boundary_1( 2, 2, e ) )**2 ) 
   end do

   do e= 1, Number_Element_on_Boundary_2
      Edge_Length_2( e )= sqrt( ( Position_Node_on_Boundary_2( 1, 1, e ) -Position_Node_on_Boundary_2( 1, 2, e ) )**2 & 
                     +( Position_Node_on_Boundary_2( 2, 1, e ) -Position_Node_on_Boundary_2( 2, 2, e ) )**2 ) 
   end do

   allocate( Value_Node_in_Element_on_Boundary_1( 3, Number_Element_on_Boundary_1 ) )
   allocate( Value_Node_in_Element_on_Boundary_2( 3, Number_Element_on_Boundary_2 ) )

   do e= 1, Number_Element_on_Boundary_1
      do i= 1, 3
         Value_Node_in_Element_on_Boundary_1( i, e )= Value_Node_Integrated( Index_Element_2_Node( i, Element_Number_Boundary_1( e ) ) )
      end do
   end do

   do e= 1, Number_Element_on_Boundary_2
      do i= 1, 3
         Value_Node_in_Element_on_Boundary_2( i, e )= Value_Node_Integrated( Index_Element_2_Node( i, Element_Number_Boundary_2( e ) ) )
      end do
   end do

   !====================================================================================================
   write(*,*)'    Compute Unit Normal Vector '
   !====================================================================================================

   allocate( Unit_Normal_Vector_1( 2, Number_Element_on_Boundary_1 ) )
   allocate( Unit_Normal_Vector_2( 2, Number_Element_on_Boundary_2 ) )

   call Compute_Unit_Normal_Vector_All &
      ( Number_Element_on_Boundary_1, Position_Node_on_Boundary_1, Unit_Normal_Vector_1 )

   call Compute_Unit_Normal_Vector_All &
      ( Number_Element_on_Boundary_2, Position_Node_on_Boundary_2, Unit_Normal_Vector_2 )

   !====================================================================================================
   write(*,*)'    Compute Basis Finctions'
   !====================================================================================================

   if( Shape_Element==3 )then

      allocate( BasisFunction_A( Shape_Element, Number_Element ) )
      allocate( BasisFunction_B( Shape_Element, Number_Element ) )
      allocate( BasisFunction_C( Shape_Element, Number_Element ) )

      call Compute_Basis_Function_Triangle_Element &
         ( Number_Node, Position_Node, &
           Number_Element, Index_Element_2_Node, &
         !====================================================================================================
           BasisFunction_A, BasisFunction_B, BasisFunction_C )

   else
      call Output_Error( 'Compute_Line_Integral_Poynting_Vector', 127 )
   end if

   allocate( BasisFunction_A_Boundary_1( Shape_Element, Number_Element_on_Boundary_1 ) )
   allocate( BasisFunction_B_Boundary_1( Shape_Element, Number_Element_on_Boundary_1 ) )
   allocate( BasisFunction_C_Boundary_1( Shape_Element, Number_Element_on_Boundary_1 ) )

   allocate( BasisFunction_A_Boundary_2( Shape_Element, Number_Element_on_Boundary_2 ) )
   allocate( BasisFunction_B_Boundary_2( Shape_Element, Number_Element_on_Boundary_2 ) )
   allocate( BasisFunction_C_Boundary_2( Shape_Element, Number_Element_on_Boundary_2 ) )

   do e= 1, Number_Element_on_Boundary_1
      do i= 1, Shape_Element
         BasisFunction_A_Boundary_1( i, e )= BasisFunction_A( i, Element_Number_Boundary_1( e ) )
         BasisFunction_B_Boundary_1( i, e )= BasisFunction_B( i, Element_Number_Boundary_1( e ) )
         BasisFunction_C_Boundary_1( i, e )= BasisFunction_C( i, Element_Number_Boundary_1( e ) )
      end do
   end do

   do e= 1, Number_Element_on_Boundary_2
      do i= 1, Shape_Element
         BasisFunction_A_Boundary_2( i, e )= BasisFunction_A( i, Element_Number_Boundary_2( e ) )
         BasisFunction_B_Boundary_2( i, e )= BasisFunction_B( i, Element_Number_Boundary_2( e ) )
         BasisFunction_C_Boundary_2( i, e )= BasisFunction_C( i, Element_Number_Boundary_2( e ) )
      end do
   end do

   deallocate( BasisFunction_A )
   deallocate( BasisFunction_B )
   deallocate( BasisFunction_C )

   !====================================================================================================
   write(*,*)'    Compute Magnetic Field from Obtained Electric Field '
   !====================================================================================================

   allocate( Term_Unit_Normal_Vector_1( Shape_Element, Number_Element_on_Boundary_1 ) )
   allocate( Term_Unit_Normal_Vector_2( Shape_Element, Number_Element_on_Boundary_2 ) )

   do e= 1, Number_Element_on_Boundary_1
      do i= 1, Shape_Element
         Term_Unit_Normal_Vector_1( i, e ) &
         = ( Unit_Normal_Vector_1( 1, e )*BasisFunction_B_Boundary_1( i, e ) +Unit_Normal_Vector_1( 2, e )*BasisFunction_C_Boundary_1( i, e ) ) &
          *conjg( Value_Node_in_Element_on_Boundary_1( i, e ) ) 
      end do
   end do
 
   do e= 1, Number_Element_on_Boundary_2
      do i= 1, Shape_Element
         Term_Unit_Normal_Vector_2( i, e ) &
         = ( Unit_Normal_Vector_2( 1, e )*BasisFunction_B_Boundary_2( i, e ) +Unit_Normal_Vector_2( 2, e )*BasisFunction_C_Boundary_2( i, e ) ) &
          *conjg( Value_Node_in_Element_on_Boundary_2( i, e ) ) 
      end do
   end do
 
   deallocate( Unit_Normal_Vector_1 )
   deallocate( Unit_Normal_Vector_2 )

   allocate( Term_Line_Integral_1( Shape_Element, Number_Element_on_Boundary_1 ) )
   allocate( Term_Line_Integral_2( Shape_Element, Number_Element_on_Boundary_2 ) )

   do e= 1, Number_Element_on_Boundary_1
      do i= 1, Shape_Element
         Term_Line_Integral_1( i, e ) &
         = Value_Node_in_Element_on_Boundary_1( i, e ) *Edge_Length_1( e ) &
          *( BasisFunction_A_Boundary_1( i, e ) &
            +( BasisFunction_B_Boundary_1( i, e )*( Position_Node_on_Boundary_1( 1, 1, e ) +Position_Node_on_Boundary_1( 1, 2, e ) ) )/2d0 & 
            +( BasisFunction_C_Boundary_1( i, e )*( Position_Node_on_Boundary_1( 2, 1, e ) +Position_Node_on_Boundary_1( 2, 2, e ) ) )/2d0 ) 
      end do
   end do

   do e= 1, Number_Element_on_Boundary_2
      do i= 1, Shape_Element
         Term_Line_Integral_2( i, e ) &
         = Value_Node_in_Element_on_Boundary_2( i, e ) *Edge_Length_2( e ) &
          *( BasisFunction_A_Boundary_2( i, e ) &
            +( BasisFunction_B_Boundary_2( i, e )*( Position_Node_on_Boundary_2( 1, 1, e ) +Position_Node_on_Boundary_2( 1, 2, e ) ) )/2d0 & 
            +( BasisFunction_C_Boundary_2( i, e )*( Position_Node_on_Boundary_2( 2, 1, e ) +Position_Node_on_Boundary_2( 2, 2, e ) ) )/2d0 )  
      end do
   end do

   
   !====================================================================================================
   write(*,*)'    Investigate Poynting Vector '
   !====================================================================================================

   allocate( Poynting_Vector_Nabla_Factor_1( 2, Number_Element_on_Boundary_1 ) )
   allocate( Poynting_Vector_Nabla_Factor_2( 2, Number_Element_on_Boundary_2 ) )

   do e= 1, Number_Element_on_Boundary_1
      do i= 1, 2 ! x, y
         Poynting_Vector_Nabla_Factor_1( i, e )= Zero
      end do
   end do

   do e= 1, Number_Element_on_Boundary_2
      do i= 1, 2 ! x, y
         Poynting_Vector_Nabla_Factor_2( i, e )= Zero
      end do
   end do

   do e= 1, Number_Element_on_Boundary_1
      do i= 1, Shape_Element ! Local Node
         Poynting_Vector_Nabla_Factor_1( 1, e )= Poynting_Vector_Nabla_Factor_1( 1, e ) +BasisFunction_B_Boundary_1( i, e ) *conjg( Value_Node_in_Element_on_Boundary_1( i, e ) )
         Poynting_Vector_Nabla_Factor_1( 2, e )= Poynting_Vector_Nabla_Factor_1( 2, e ) +BasisFunction_C_Boundary_1( i, e ) *conjg( Value_Node_in_Element_on_Boundary_1( i, e ) )
      end do
   end do

   do e= 1, Number_Element_on_Boundary_2
      do i= 1, Shape_Element ! Local Node
         Poynting_Vector_Nabla_Factor_2( 1, e )= Poynting_Vector_Nabla_Factor_2( 1, e ) +BasisFunction_B_Boundary_2( i, e ) *conjg( Value_Node_in_Element_on_Boundary_2( i, e ) )
         Poynting_Vector_Nabla_Factor_2( 2, e )= Poynting_Vector_Nabla_Factor_2( 2, e ) +BasisFunction_C_Boundary_2( i, e ) *conjg( Value_Node_in_Element_on_Boundary_2( i, e ) )
      end do
   end do

   allocate( Poynting_Vector_Factor_1( Number_Element_on_Boundary_1 ) )
   allocate( Poynting_Vector_Factor_2( Number_Element_on_Boundary_2 ) )

   do e= 1, Number_Element_on_Boundary_1
      Poynting_Vector_Factor_1( e )= Zero
   end do

   do e= 1, Number_Element_on_Boundary_2
      Poynting_Vector_Factor_2( e )= Zero
   end do

   do e= 1, Number_Element_on_Boundary_1
      do i= 1, Shape_Element ! Local Node
         Poynting_Vector_Factor_1( e )= Poynting_Vector_Factor_1( e ) +Term_Line_Integral_1( i, e )
      end do
   end do

   do e= 1, Number_Element_on_Boundary_2
      do i= 1, Shape_Element ! Local Node
         Poynting_Vector_Factor_2( e )= Poynting_Vector_Factor_2( e ) +Term_Line_Integral_2( i, e )
      end do
   end do

   allocate( Poynting_Vector_1( 2, Number_Element_on_Boundary_1 ) )
   allocate( Poynting_Vector_2( 2, Number_Element_on_Boundary_2 ) )

   do e= 1, Number_Element_on_Boundary_1
      do i= 1, 2 ! x, y
         Poynting_Vector_1( i, e )= real( Complex_Unit *Poynting_Vector_Factor_1( e ) *Poynting_Vector_Nabla_Factor_1( i, e ) )
      end do
   end do

   do e= 1, Number_Element_on_Boundary_2
      do i= 1, 2 ! x, y
         Poynting_Vector_2( i, e )= real( Complex_Unit *Poynting_Vector_Factor_2( e ) *Poynting_Vector_Nabla_Factor_2( i, e ) )
      end do
   end do

   allocate( Poynting_Vector( 2, Number_Element_on_Boundary_1 ) )

   do e= 1, Number_Element_on_Boundary_1
      do i= 1, 2 ! x, y
         Poynting_Vector( i, e )= Poynting_Vector_1( i, e ) !+Poynting_Vector_2( i, e )
         ! e_1 /= e_2
         !Poynting_Vector( i, e )= ( Poynting_Vector_1( i, e ) +Poynting_Vector_2( i, e ) )/2d0
      end do
   end do

!   do e= 1, Number_Element_on_Boundary_1
!      write(100000,'(2x,es15.8,2x,es15.8,2x,es15.8,2x,es15.8)') &
!      ( Position_Node_on_Boundary_1( 1, 1, e ) +Position_Node_on_Boundary_1( 1, 2, e ) )/2d0, & 
!      ( Position_Node_on_Boundary_1( 2, 1, e ) +Position_Node_on_Boundary_1( 2, 2, e ) )/2d0, & 
!      Poynting_Vector( 1, e ), Poynting_Vector( 2, e )
!   end do
!   write(100000,*)' '
!
!   do e= 1, Number_Element_on_Boundary_2
!      write(100001,'(2x,es15.8,2x,es15.8,2x,es15.8,2x,es15.8)') &
!      ( Position_Node_on_Boundary_2( 1, 1, e ) +Position_Node_on_Boundary_2( 1, 2, e ) )/2d0, & 
!      ( Position_Node_on_Boundary_2( 2, 1, e ) +Position_Node_on_Boundary_2( 2, 2, e ) )/2d0, & 
!      Poynting_Vector_2( 1, e ), Poynting_Vector_2( 2, e )
!   end do
!   write(100001,*)' '


   deallocate( Poynting_Vector )

   deallocate( Poynting_Vector_1 )
   deallocate( Poynting_Vector_2 )

   deallocate( Poynting_Vector_Factor_1 )
   deallocate( Poynting_Vector_Factor_2 )

   deallocate( Poynting_Vector_Nabla_Factor_1 )
   deallocate( Poynting_Vector_Nabla_Factor_2 )

   !====================================================================================================

   deallocate( Value_Node_in_Element_on_Boundary_1 )
   deallocate( Value_Node_in_Element_on_Boundary_2 )

   deallocate( Edge_Length_1 )
   deallocate( Edge_Length_2 )

   deallocate( Position_Node_on_Boundary_1 )
   deallocate( Position_Node_on_Boundary_2 )

   deallocate( BasisFunction_A_Boundary_1 )
   deallocate( BasisFunction_B_Boundary_1 )
   deallocate( BasisFunction_C_Boundary_1 )

   deallocate( BasisFunction_A_Boundary_2 )
   deallocate( BasisFunction_B_Boundary_2 )
   deallocate( BasisFunction_C_Boundary_2 )

   allocate( Term_Unit_Normal_Vector_1_Summed( Number_Element_on_Boundary_1 ) )
   allocate( Term_Unit_Normal_Vector_2_Summed( Number_Element_on_Boundary_2 ) )

   allocate( Term_Line_Integral_1_Summed( Number_Element_on_Boundary_1 ) )
   allocate( Term_Line_Integral_2_Summed( Number_Element_on_Boundary_2 ) )

   do e= 1, Number_Element_on_Boundary_1
      Term_Unit_Normal_Vector_1_Summed( e )= Zero
      Term_Line_Integral_1_Summed( e )= Zero
   end do

   do e= 1, Number_Element_on_Boundary_2
      Term_Unit_Normal_Vector_2_Summed( e )= Zero
      Term_Line_Integral_2_Summed( e )= Zero
   end do

   do e= 1, Number_Element_on_Boundary_1
      do i= 1, 3
         Term_Unit_Normal_Vector_1_Summed( e )= Term_Unit_Normal_Vector_1_Summed( e ) +Term_Unit_Normal_Vector_1( i, e ) 
      end do
      do i= 1, 3
         Term_Line_Integral_1_Summed( e )= Term_Line_Integral_1_Summed( e ) +Term_Line_Integral_1( i, e ) 
      end do
   end do

   do e= 1, Number_Element_on_Boundary_2
      do i= 1, 3
         Term_Unit_Normal_Vector_2_Summed( e )= Term_Unit_Normal_Vector_2_Summed( e ) +Term_Unit_Normal_Vector_2( i, e ) 
      end do
      do i= 1, 3
         Term_Line_Integral_2_Summed( e )= Term_Line_Integral_2_Summed( e ) +Term_Line_Integral_2( i, e ) 
      end do
   end do

   deallocate( Term_Unit_Normal_Vector_1 )
   deallocate( Term_Unit_Normal_Vector_2 )
   deallocate( Term_Line_Integral_1 )
   deallocate( Term_Line_Integral_2 )

   !if( Number_Element_on_Boundary_1/=Number_Element_on_Boundary_2 )then
   !   call Output_Error( 'Compute_Line_Integral_Poynting_Vector', 218 )
   !end if

   Value_Integrated_Poynting_Vector_1= 0d0
   Value_Integrated_Poynting_Vector_2= 0d0

   do e= 1, Number_Element_on_Boundary_1
      Value_Integrated_Poynting_Vector_1 &
      = Value_Integrated_Poynting_Vector_1 &
       +real( Complex_Unit *Term_Unit_Normal_Vector_1_Summed( e )* Term_Line_Integral_1_Summed( e ) )
   end do

   do e= 1, Number_Element_on_Boundary_2
      Value_Integrated_Poynting_Vector_2 &
      = Value_Integrated_Poynting_Vector_2 &
       +real( Complex_Unit *Term_Unit_Normal_Vector_2_Summed( e )* Term_Line_Integral_2_Summed( e ) )
   end do

   Value_Integrated_Poynting_Vector= ( Value_Integrated_Poynting_Vector_1 +Value_Integrated_Poynting_Vector_2 )/2d0

   !do e= 1, Number_Element_on_Boundary_1
   !   write(99999,*) 'Boundary_1', e, real( Complex_Unit *Term_Unit_Normal_Vector_1_Summed( e ) *Term_Line_Integral_1_Summed( e ) ) 
   !end do

   !do e= 1, Number_Element_on_Boundary_2
   !   write(99999,*) 'Boundary_2', e, real( Complex_Unit *Term_Unit_Normal_Vector_2_Summed( e ) *Term_Line_Integral_2_Summed( e ) ) 
   !end do

   deallocate( Term_Line_Integral_1_Summed )
   deallocate( Term_Line_Integral_2_Summed )

   deallocate( Term_Unit_Normal_Vector_1_Summed )
   deallocate( Term_Unit_Normal_Vector_2_Summed )

   return
end subroutine Compute_Line_Integral_Poynting_Vector 



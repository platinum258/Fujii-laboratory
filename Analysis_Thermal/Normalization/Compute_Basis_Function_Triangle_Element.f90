
subroutine Compute_Basis_Function_Triangle_Element &
       ( Number_Node, Position_Node, &
         Number_Element, Index_Element_2_Node, &
       !====================================================================================================
         BasisFunction_A, BasisFunction_B, BasisFunction_C )

   implicit none
   integer, intent(in) :: Number_Node, Number_Element

   double precision, intent(in) :: Position_Node( 2, Number_Node ) 
   integer, intent(in) :: Index_Element_2_Node( 3, Number_Element )

   double precision, intent(out) :: BasisFunction_A( 3, Number_Element ), BasisFunction_B( 3, Number_Element ), BasisFunction_C( 3, Number_Element )

   integer :: e, i, j, k

   double precision, allocatable, dimension(:,:,:,:) :: Difference_Position_Node_Element 
   double precision, allocatable, dimension(:,:,:) :: Position_Node_Element 
   double precision, allocatable, dimension(:) :: Area_Element_Triangle

   integer, parameter :: Shape_Element= 3

   !====================================================================================================
   write(*,*)'    call Compute_Basis_Function_Triangle_Element'
   !====================================================================================================

   !====================================================================================================
   write(*,*)'    Check Input Data'
   !====================================================================================================

   do e= 1, Number_Element 
      do i= 1, Shape_Element 
         if( Index_Element_2_Node( i, e )==0 )then
            write(*,*) e, i, Index_Element_2_Node( i, e )
            call Output_Error( 'Compute_Basis_Function_Triangle_Element', 36 )
         end if
      end do
   end do

   allocate( Position_Node_Element( 2, Shape_Element, Number_Element ) )
 
   !$omp parallel do default( none ) &
   !$omp private( e, i, j ) & 
   !$omp shared( Number_Element, Position_Node, Position_Node_Element, Index_Element_2_Node ) 
   do e= 1, Number_Element
      do i= 1, Shape_Element
         do j= 1, 2
            Position_Node_Element( j, i, e )= Position_Node( j, Index_Element_2_Node( i, e ) )
         end do
      end do
   end do

   allocate( Difference_Position_Node_Element( 2, Shape_Element, Shape_Element, Number_Element ) )

   !$omp parallel do default( none ) &
   !$omp private( e, i, j ) & 
   !$omp shared( Number_Element, Difference_Position_Node_Element, Position_Node_Element ) 
   do e= 1, Number_Element
      do i= 1, Shape_Element
         do j= 1, Shape_Element
            do k= 1, 2
               Difference_Position_Node_Element( k, j, i, e ) &
               =Position_Node_Element( k, j, e ) -Position_Node_Element( k, i, e )
            end do
         end do
      end do
   end do

   allocate( Area_Element_Triangle( Number_Element ) )

   !$omp parallel do default( none ) &
   !$omp private( e ) & 
   !$omp shared( Number_Element, Area_Element_Triangle, Difference_Position_Node_Element ) 
   do e= 1, Number_Element
      Area_Element_Triangle( e ) &
      =abs( Difference_Position_Node_Element( 1, 1, 3, e ) & 
         *Difference_Position_Node_Element( 2, 2, 3, e ) &
         -Difference_Position_Node_Element( 1, 2, 3, e ) & 
         *Difference_Position_Node_Element( 2, 1, 3, e ) )/2.0D0
   end do 

   !$omp parallel do default( none ) &
   !$omp private( e ) & 
   !$omp shared( Number_Element, Area_Element_Triangle, Position_Node_Element ) 
   do e= 1, Number_Element
      if( Area_Element_Triangle( e )==0d0 )then
         write(*,*)'Area_Element_Triangle( e )==0.0d0'
         write(*,*)'Element Number=', e
         call Output_Error( 'Compute_Basis_Function_Triangle_Element', 166 )
      end if
   end do 

   !$omp parallel do default( none )   &
   !$omp private( e ) &
   !$omp shared( Number_Element, Position_Node_Element ) &
   !$omp shared( Area_Element_Triangle, BasisFunction_A )
   do e= 1, Number_Element
  
      BasisFunction_A( 1, e ) &
      = ( Position_Node_Element( 1, 2, e ) *Position_Node_Element( 2, 3, e ) &
         -Position_Node_Element( 1, 3, e ) *Position_Node_Element( 2, 2, e ) ) &
       /( 2.0D0 *Area_Element_Triangle( e ) )

      BasisFunction_A( 2, e ) &
      = ( Position_Node_Element( 1, 3, e ) *Position_Node_Element( 2, 1, e ) & 
         -Position_Node_Element( 1, 1, e ) *Position_Node_Element( 2, 3, e ) ) &
       /( 2.0D0 *Area_Element_Triangle( e ) )

      BasisFunction_A( 3, e ) & 
      = ( Position_Node_Element( 1, 1, e ) *Position_Node_Element( 2, 2, e ) & 
         -Position_Node_Element( 1, 2, e ) *Position_Node_Element( 2, 1, e ) ) & 
       /( 2.0D0 *Area_Element_Triangle( e ) )
      
   end do
 
   !$omp parallel do default( none )   &
   !$omp private( e ) &
   !$omp shared( Number_Element, Difference_Position_Node_Element ) &
   !$omp shared( Area_Element_Triangle, BasisFunction_B )
   do e= 1, Number_Element

      BasisFunction_B( 1, e ) &
      = Difference_Position_Node_Element( 2, 2, 3, e )/( 2.0D0 *Area_Element_Triangle( e ) )

      BasisFunction_B( 2, e ) & 
      = Difference_Position_Node_Element( 2, 3, 1, e )/( 2.0D0 *Area_Element_Triangle( e ) )

      BasisFunction_B( 3, e ) & 
      = Difference_Position_Node_Element( 2, 1, 2, e )/( 2.0D0 *Area_Element_Triangle( e ) )
      
   end do


   !$omp parallel do default( none )   &
   !$omp private( e ) &
   !$omp shared( Number_Element, Difference_Position_Node_Element ) &
   !$omp shared( Area_Element_Triangle, BasisFunction_C )
   do e= 1, Number_Element

      BasisFunction_C( 1, e ) & 
      = Difference_Position_Node_Element( 1, 3, 2, e )/( 2.0D0 *Area_Element_Triangle( e ) )

      BasisFunction_C( 2, e ) & 
      = Difference_Position_Node_Element( 1, 1, 3, e )/( 2.0D0 *Area_Element_Triangle( e ) )

      BasisFunction_C( 3, e ) & 
      = Difference_Position_Node_Element( 1, 2, 1, e )/( 2.0D0 *Area_Element_Triangle( e ) )
      
   end do
 
   deallocate( Position_Node_Element )
   deallocate( Difference_Position_Node_Element )
   deallocate( Area_Element_Triangle )


   return
end subroutine Compute_Basis_Function_Triangle_Element 


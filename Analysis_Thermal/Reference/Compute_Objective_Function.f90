
subroutine Compute_Objective_Function &
       ( Flag_No_Structure, ID_Element, Core_Number, &
         Value_Objective_Function, &
         Position_Node, Number_Node, &
         Index_Element_2_Node, Class_Element, Number_Element, &
         Value_ObjectiveFunction_Normalize,  &
         !=================================================================
         Objective_Function, Objective_Function_Value )

   !$use onp_lib
   use Parameters
   implicit none

   integer, intent(in) :: Flag_No_Structure
   integer, intent(in) :: ID_Element 
   integer, intent(in) :: Core_Number 
   integer, intent(in) :: Number_Node, Number_Element

   double precision, intent(in) :: Value_Objective_Function( Core_Number, Number_Node )

   double precision, intent(in) :: Position_Node( 2, Number_Node )
   integer, intent(in) :: Index_Element_2_Node( 3, Number_Element ) 
   integer, intent(in) :: Class_Element( Number_Element ) 

   double precision, intent(in) :: Value_ObjectiveFunction_Normalize

   double precision, intent(out) :: Objective_Function
   double precision, intent(out) :: Objective_Function_Value

   integer :: e, i, j, k
   double precision, allocatable, dimension(:) :: Objective_Function_Element
   double precision, allocatable, dimension(:) :: Objective_Function_Element_tmp
   double precision, allocatable, dimension(:,:,:) :: Position_Node_Element_Triangle
   double precision, allocatable, dimension(:,:,:,:) :: Difference_Position_Node_Element_Triangle
   double precision, allocatable, dimension(:) :: Area_Element  
   double precision, allocatable, dimension(:,:) :: BasisFunction_A, BasisFunction_B, BasisFunction_C
   double precision, allocatable, dimension(:,:,:) :: BasisFunctionBC_Times_Difference_Position_Node_TE
   double precision, allocatable, dimension(:,:) :: Basis_Function_Times_Position_Node_TE
   !complex( kind( 0d0 ) ), allocatable, dimension(:,:,:) :: Coefficient_C, Coefficient_xi, Coefficient_eta 
   !complex( kind( 0d0 ) ), allocatable, dimension(:,:,:) :: Coefficient_xi_eta, Coefficient_xi_xi, Coefficient_eta_eta
   !complex( kind( 0d0 ) ), allocatable, dimension(:,:,:) :: Integral_Element_Triangle
   double precision, allocatable, dimension(:,:,:) :: Coefficient_C, Coefficient_xi, Coefficient_eta 
   double precision, allocatable, dimension(:,:,:) :: Coefficient_xi_eta, Coefficient_xi_xi, Coefficient_eta_eta
   double precision, allocatable, dimension(:,:,:) :: Integral_Element_Triangle

   integer :: Counter_Element_OnjectiveFunction, Number_Element_ObjectiveFunction
   integer, allocatable, dimension(:,:) :: Index_Element_2_Node_ObjectiveFunction, Index_Element_2_Node_ObjectiveFunction_tmp

   double precision, allocatable, dimension(:,:,:) :: Intensity_Value_Objective_Function

   !=======================================================================================================
   write(*,*)'   call Compute_Objective_Function ' 
   !=======================================================================================================

   allocate( Index_Element_2_Node_ObjectiveFunction_tmp( 3, Number_Element ) ) 

   Counter_Element_OnjectiveFunction= 0

   ! No Parallel, Vectorize
   do e= 1, Number_Element
      if( Class_Element( e )==ID_Element )then

         Counter_Element_OnjectiveFunction= Counter_Element_OnjectiveFunction +1

         do i= 1, 3
            Index_Element_2_Node_ObjectiveFunction_tmp( i, Counter_Element_OnjectiveFunction ) & 
            = Index_Element_2_Node( i, e ) 
         end do
      end if
   end do

   Number_Element_ObjectiveFunction= Counter_Element_OnjectiveFunction

   allocate( Index_Element_2_Node_ObjectiveFunction( 3, Number_Element_ObjectiveFunction ) ) 
   Index_Element_2_Node_ObjectiveFunction( 1:3, 1:Number_Element_ObjectiveFunction )&
   = Index_Element_2_Node_ObjectiveFunction_tmp( 1:3, 1:Number_Element_ObjectiveFunction )

   deallocate( Index_Element_2_Node_ObjectiveFunction_tmp ) 

   !===================================================================================
   write(*,*)'      Position_Node --> Position_Node_Element_Triangle '
   !===================================================================================
   allocate( Position_Node_Element_Triangle( 2, 3, Number_Element_ObjectiveFunction ) )
   
   !$omp parallel do default( none )   &
   !$omp private( e, i, j ) &
   !$omp shared( Number_Element_ObjectiveFunction, Position_Node_Element_Triangle, Position_Node, Index_Element_2_Node_ObjectiveFunction ) 
   do e= 1, Number_Element_ObjectiveFunction
      do i= 1, 3
         do j= 1, 2
            Position_Node_Element_Triangle( j, i, e )= Position_Node( j, Index_Element_2_Node_ObjectiveFunction( i, e ) )
         end do
      end do
   end do

   allocate( Difference_Position_Node_Element_Triangle( 2, 3, 3, Number_Element_ObjectiveFunction ) )

   do e= 1, Number_Element_ObjectiveFunction
      do i= 1, 3
         do j= 1, 3
            do k= 1, 2
               Difference_Position_Node_Element_Triangle( k, j, i, e ) &
               =Position_Node_Element_Triangle( k, j, e ) -Position_Node_Element_Triangle( k, i, e )
            end do
         end do
      end do
   end do

   allocate( Area_Element( Number_Element_ObjectiveFunction ) )
   Area_Element( : ) &
   =abs( Difference_Position_Node_Element_Triangle( 1, 1, 3, : ) & 
      *Difference_Position_Node_Element_Triangle( 2, 2, 3, : ) &
      -Difference_Position_Node_Element_Triangle( 1, 2, 3, : ) & 
      *Difference_Position_Node_Element_Triangle( 2, 1, 3, : ) )/2.0D0

   !$omp parallel do default( none ) &
   !$omp private( e ) & 
   !$omp shared( Number_Element_ObjectiveFunction, Area_Element, Position_Node_Element_Triangle ) 
   do e= 1, Number_Element_ObjectiveFunction
      if( Area_Element( e )==0.0d0 )then
         write(*,*)'Area_Element( e )==0.0d0'
         write(*,*)'Element Number=', e
         write(*,*)'Number_Element_ObjectiveFunction=', Number_Element_ObjectiveFunction
         call Output_Error( 'Compute_Objective_Function', 166 )
      end if
   end do 
    
   allocate( BasisFunction_A( 3, Number_Element_ObjectiveFunction ) )
   allocate( BasisFunction_B( 3, Number_Element_ObjectiveFunction ) )
   allocate( BasisFunction_C( 3, Number_Element_ObjectiveFunction ) )

   BasisFunction_A( 1, : ) &
   = ( Position_Node_Element_Triangle( 1, 2, : ) *Position_Node_Element_Triangle( 2, 3, : ) &
      -Position_Node_Element_Triangle( 1, 3, : ) *Position_Node_Element_Triangle( 2, 2, : ) ) &
    /( 2.0D0 *Area_Element( : ) )

   BasisFunction_A( 2, : ) &
   = ( Position_Node_Element_Triangle( 1, 3, : ) *Position_Node_Element_Triangle( 2, 1, : ) & 
      -Position_Node_Element_Triangle( 1, 1, : ) *Position_Node_Element_Triangle( 2, 3, : ) ) &
    /( 2.0D0 *Area_Element( : ) )

   BasisFunction_A( 3, : ) & 
   = ( Position_Node_Element_Triangle( 1, 1, : ) *Position_Node_Element_Triangle( 2, 2, : ) & 
      -Position_Node_Element_Triangle( 1, 2, : ) *Position_Node_Element_Triangle( 2, 1, : ) ) & 
    /( 2.0D0 *Area_Element( : ) )

   BasisFunction_B( 1, : ) &
   = Difference_Position_Node_Element_Triangle( 2, 2, 3, : )/( 2.0D0 *Area_Element( : ) )

   BasisFunction_B( 2, : ) & 
   = Difference_Position_Node_Element_Triangle( 2, 3, 1, : )/( 2.0D0 *Area_Element( : ) )

   BasisFunction_B( 3, : ) & 
   = Difference_Position_Node_Element_Triangle( 2, 1, 2, : )/( 2.0D0 *Area_Element( : ) )
      
   BasisFunction_C( 1, : ) & 
   = Difference_Position_Node_Element_Triangle( 1, 3, 2, : )/( 2.0D0 *Area_Element( : ) )

   BasisFunction_C( 2, : ) & 
   = Difference_Position_Node_Element_Triangle( 1, 1, 3, : )/( 2.0D0 *Area_Element( : ) )

   BasisFunction_C( 3, : ) & 
   = Difference_Position_Node_Element_Triangle( 1, 2, 1, : )/( 2.0D0 *Area_Element( : ) )

   allocate( Basis_Function_Times_Position_Node_TE( 3, Number_Element_ObjectiveFunction ) ) 

   !$omp parallel do default( none )   &
   !$omp private( e, i ) &
   !$omp shared( Number_Element_ObjectiveFunction, Basis_Function_Times_Position_Node_TE ) &
   !$omp shared( BasisFunction_A, BasisFunction_B, BasisFunction_C, Position_Node_Element_Triangle )
   do e= 1, Number_Element_ObjectiveFunction
      do i= 1, 3
         Basis_Function_Times_Position_Node_TE( i, e ) & 
         = BasisFunction_A( i, e ) & 
          +BasisFunction_B( i, e ) *Position_Node_Element_Triangle( 1, 3, e ) & 
          +BasisFunction_C( i, e ) *Position_Node_Element_Triangle( 2, 3, e ) 
      end do
   end do

   allocate( BasisFunctionBC_Times_Difference_Position_Node_TE( 2, 3, Number_Element_ObjectiveFunction ) )

   !$omp parallel do default( none )   &
   !$omp private( e, i, j ) &
   !$omp shared( Number_Element_ObjectiveFunction, BasisFunctionBC_Times_Difference_Position_Node_TE ) &
   !$omp shared( BasisFunction_A, BasisFunction_B, BasisFunction_C, Difference_Position_Node_Element_Triangle )
   do e= 1, Number_Element_ObjectiveFunction
      do i= 1, 3
         do j= 1, 2
            BasisFunctionBC_Times_Difference_Position_Node_TE( j, i, e ) &
            = BasisFunction_B( i, e ) *Difference_Position_Node_Element_Triangle( 1, j, 3, e ) &
             +BasisFunction_C( i, e ) *Difference_Position_Node_Element_Triangle( 2, j, 3, e )
         end do
      end do
   end do
 
   deallocate( Difference_Position_Node_Element_Triangle )
 
   allocate( Coefficient_C( 3, 3, Number_Element_ObjectiveFunction ) )
   allocate( Coefficient_xi( 3, 3, Number_Element_ObjectiveFunction ) )
   allocate( Coefficient_eta( 3, 3, Number_Element_ObjectiveFunction ) )
   allocate( Coefficient_xi_eta( 3, 3, Number_Element_ObjectiveFunction ) )
   allocate( Coefficient_xi_xi( 3, 3, Number_Element_ObjectiveFunction ) )
   allocate( Coefficient_eta_eta( 3, 3, Number_Element_ObjectiveFunction ) )

 
   !$omp parallel do default( none )   &
   !$omp private( e, i, j ) &
   !$omp shared( Number_Element_ObjectiveFunction ) &
   !$omp shared( Basis_Function_Times_Position_Node_TE, Coefficient_C ) 
   do e= 1, Number_Element_ObjectiveFunction
      do i= 1, 3
         do j= 1, 3
            Coefficient_C( j, i, e ) & 
            = Basis_Function_Times_Position_Node_TE( i, e )*Basis_Function_Times_Position_Node_TE( j, e ) 
         end do
      end do
   end do

   !$omp parallel do default( none )   &
   !$omp private( e, i, j ) &
   !$omp shared( Number_Element_ObjectiveFunction, Coefficient_Xi ) & 
   !$omp shared( Basis_Function_Times_Position_Node_TE ) &
   !$omp shared( BasisFunctionBC_Times_Difference_Position_Node_TE ) 
   do e= 1, Number_Element_ObjectiveFunction
      do i= 1, 3
         do j= 1, 3
            Coefficient_Xi( j, i, e ) & 
            = Basis_Function_Times_Position_Node_TE( j, e ) &
             *BasisFunctionBC_Times_Difference_Position_Node_TE( 1, i, e ) &
             +Basis_Function_Times_Position_Node_TE( i, e ) &
             *BasisFunctionBC_Times_Difference_Position_Node_TE( 1, j, e ) 
         end do
      end do
   end do

   !$omp parallel do default( none )   &
   !$omp private( e, i, j ) &
   !$omp shared( Number_Element_ObjectiveFunction, Coefficient_Eta ) & 
   !$omp shared( Basis_Function_Times_Position_Node_TE ) &
   !$omp shared( BasisFunctionBC_Times_Difference_Position_Node_TE ) 
   do e= 1, Number_Element_ObjectiveFunction
      do i= 1, 3
         do j= 1, 3
            Coefficient_Eta( j, i, e ) &
            = Basis_Function_Times_Position_Node_TE( j, e ) &
             *BasisFunctionBC_Times_Difference_Position_Node_TE( 2, i, e ) &
             +Basis_Function_Times_Position_Node_TE( i, e ) &
             *BasisFunctionBC_Times_Difference_Position_Node_TE( 2, j, e ) 
         end do
      end do
   end do

   deallocate( Basis_Function_Times_Position_Node_TE ) 

   !$omp parallel do default( none )   &
   !$omp private( e, i, j ) &
   !$omp shared( Number_Element_ObjectiveFunction, Coefficient_Xi_Eta ) & 
   !$omp shared( BasisFunctionBC_Times_Difference_Position_Node_TE ) 
   do e= 1, Number_Element_ObjectiveFunction
      do i= 1, 3
         do j= 1, 3
            Coefficient_Xi_Eta( j, i, e ) & 
            = BasisFunctionBC_Times_Difference_Position_Node_TE( 2, i, e ) &
             *BasisFunctionBC_Times_Difference_Position_Node_TE( 1, j, e ) &
             +BasisFunctionBC_Times_Difference_Position_Node_TE( 1, i, e ) &
             *BasisFunctionBC_Times_Difference_Position_Node_TE( 2, j, e ) 
         end do
      end do
   end do

   !$omp parallel do default( none )   &
   !$omp private( e, i, j ) &
   !$omp shared( Number_Element_ObjectiveFunction, Coefficient_Xi_Xi ) & 
   !$omp shared( BasisFunctionBC_Times_Difference_Position_Node_TE ) 
   do e= 1, Number_Element_ObjectiveFunction
      do i= 1, 3
         do j= 1, 3
            Coefficient_Xi_Xi( j, i, e ) & 
            = BasisFunctionBC_Times_Difference_Position_Node_TE( 1, i, e ) &
             *BasisFunctionBC_Times_Difference_Position_Node_TE( 1, j, e ) 
         end do
      end do
   end do

   !$omp parallel do default( none )   &
   !$omp private( e, i, j ) &
   !$omp shared( Number_Element_ObjectiveFunction, Coefficient_Eta_Eta ) & 
   !$omp shared( BasisFunctionBC_Times_Difference_Position_Node_TE ) 
   do e= 1, Number_Element_ObjectiveFunction
      do i= 1, 3
         do j= 1, 3
            Coefficient_Eta_Eta( j, i, e ) & 
            = BasisFunctionBC_Times_Difference_Position_Node_TE( 2, i, e ) &
             *BasisFunctionBC_Times_Difference_Position_Node_TE( 2, j, e ) 
         end do
      end do
   end do
   
   deallocate( BasisFunctionBC_Times_Difference_Position_Node_TE )
 





   deallocate( BasisFunction_A )
   deallocate( BasisFunction_B )
   deallocate( BasisFunction_C )
  
   allocate( Integral_Element_Triangle( 3, 3, Number_Element_ObjectiveFunction ) )
 
   !$omp parallel do default( none )   &
   !$omp private( e, i, j ) &
   !$omp shared( Number_Element_ObjectiveFunction, Integral_Element_Triangle ) & 
   !$omp shared( Coefficient_C, Coefficient_xi, Coefficient_eta, Coefficient_xi_eta, Coefficient_xi_xi, Coefficient_eta_eta ) 
   do e= 1, Number_Element_ObjectiveFunction
      do i= 1, 3
         do j= 1, 3
            Integral_Element_Triangle( j, i, e ) & 
            = Coefficient_C( j, i, e ) +Coefficient_xi( j, i, e )/3.0d0 +Coefficient_eta( j, i, e )/3.0d0 &
             +Coefficient_xi_xi( j, i, e )/6.0d0 +Coefficient_eta_eta( j, i, e )/6.0d0 +Coefficient_xi_eta( j, i, e )/12.0d0 
         end do
      end do
   end do

   deallocate( Coefficient_C )
   deallocate( Coefficient_xi )
   deallocate( Coefficient_eta )
   deallocate( Coefficient_xi_eta )
   deallocate( Coefficient_xi_xi )
   deallocate( Coefficient_eta_eta ) 
  
   allocate( Objective_Function_Element_tmp( Number_Element_ObjectiveFunction ) )
   allocate( Intensity_Value_Objective_Function( 3, 3, Number_Element_ObjectiveFunction ) )

   !$omp parallel do default( none )   &
   !$omp private( e, i, j ) &
   !$omp shared( Number_Element_ObjectiveFunction, Intensity_Value_Objective_Function, Value_Objective_Function ) & 
   !$omp shared( Index_Element_2_Node_ObjectiveFunction, Core_Number ) 
   do e= 1, Number_Element_ObjectiveFunction
      do i= 1, 3  
         do j= 1, 3  
            Intensity_Value_Objective_Function( j, i, e ) &
            = Value_Objective_Function( Core_Number, Index_Element_2_Node_ObjectiveFunction( i, e ) ) & 
             *( Value_Objective_Function( Core_Number, Index_Element_2_Node_ObjectiveFunction( j, e ) ) ) 
         end do
      end do
   end do

   Objective_Function_Element_tmp= 0.0d0 
 
   !$omp parallel do default( none )   &
   !$omp private( e, i, j ) &
   !$omp shared( Number_Element_ObjectiveFunction, Objective_Function_Element_tmp ) &
   !$omp shared( Intensity_Value_Objective_Function, Integral_Element_Triangle ) 
   do e= 1, Number_Element_ObjectiveFunction
      do i= 1, 3
         do j= 1, 3
            Objective_Function_Element_tmp( e ) &
            = Objective_Function_Element_tmp( e ) &
              +Intensity_Value_Objective_Function( j, i, e ) *Integral_Element_Triangle( j, i, e ) 
         end do
      end do
   end do

   deallocate( Intensity_Value_Objective_Function )
   deallocate( Integral_Element_Triangle )
 
   allocate( Objective_Function_Element( Number_Element_ObjectiveFunction ) )

   !$omp parallel do default( none )   &
   !$omp private( e ) &
   !$omp shared( Number_Element_ObjectiveFunction, Objective_Function_Element, Area_Element, Objective_Function_Element_tmp ) 
   do e= 1, Number_Element_ObjectiveFunction
      Objective_Function_Element( e )= Area_Element( e ) *Objective_Function_Element_tmp( e )
   end do

   deallocate( Objective_Function_Element_tmp )

   Objective_Function_Value= 0.0d0

   do e= 1, Number_Element_ObjectiveFunction
      Objective_Function_Value= Objective_Function_Value +Objective_Function_Element( e )
   end do

   if( Flag_No_Structure==1 )then
      Objective_Function= Objective_Function_Value/Objective_Function_Value
   else
      Objective_Function= Objective_Function_Value/Value_ObjectiveFunction_Normalize
   end if

   deallocate( Position_Node_Element_Triangle )
   deallocate( Area_Element )
   deallocate( Objective_Function_Element )

   deallocate( Index_Element_2_Node_ObjectiveFunction ) 

   write(*,*)'==============================================================================' 
   write(*,*)'Objective Function/Value_ObjectiveFunction_Normalize=', & 
          Objective_Function
   write(*,*)'Objective Function=', Objective_Function_Value 
   write(*,*)'Value_ObjectiveFunction_Normalize=', Value_ObjectiveFunction_Normalize
   write(*,*)'==============================================================================' 

   return
end subroutine Compute_Objective_Function



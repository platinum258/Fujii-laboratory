
subroutine Compute_Poynting_Vector&
           ( Number_PV, Position_PV, Element_Number_PV, &
             Length_Maximum_Vector, Length_Minimum_Vector, &
             Number_Node, Position_Node, Temperature, &
             Index_Element_2_Node, Number_Element, &
             Thermal_Conductivity_PV, &
             !================================================
             Position_Arrow_PV, Level_PV )

   implicit none

   integer, intent(in) :: Number_PV
   double precision, intent(in) :: Position_PV( 2, Number_PV ) 
   integer, intent(in) :: Element_Number_PV( Number_PV ) 
   double precision, intent(in) :: Length_Maximum_Vector, Length_Minimum_Vector 

   integer, intent(in) :: Number_Node
   double precision, intent(in) :: Position_Node( 2, Number_Node ) 
   double precision, intent(in) :: Temperature( Number_Node ) 

   integer, intent(in) :: Number_Element 
   integer, intent(in) :: Index_Element_2_Node( 4, Number_Element )
   double precision, intent(in) :: Thermal_Conductivity_PV( Number_PV ) 

   double precision, intent(out) :: Position_Arrow_PV( 2, 2, Number_PV )
   double precision, intent(out) :: Level_PV( Number_PV )


   integer :: e, i, j, k
   integer, allocatable, dimension(:,:) :: Index_Element_2_Node_PV
   double precision, allocatable, dimension(:,:) :: BasisFunction_A, BasisFunction_B, BasisFunction_C

   double precision, allocatable, dimension(:,:) :: Temperature_PV 
   double precision, allocatable, dimension(:,:) :: Gradient_Temperature 
   double precision, allocatable, dimension(:,:) :: Heat_Flux, Heat_Flux_Plot 

   double precision, allocatable, dimension(:) :: Length_Poynting_Vector
   double precision :: Maximum_Length_Poynting_Vector, Minimum_Length_Poynting_Vector

   double precision, parameter :: Zero=0d0 !dcmplx( 0d0, 0d0 )
   double precision, parameter :: Complex_Unit=1d0 !dcmplx( 0d0, 1d0 )

   !====================================================================================================
   write(*,*)'   call Compute_Poynting_Vector '
   !====================================================================================================

   allocate( Index_Element_2_Node_PV( 3, Number_PV ) )

   do e= 1, Number_PV
      do i= 1, 3
         Index_Element_2_Node_PV( i, e )= Index_Element_2_Node( i, Element_Number_PV( e ) )
      end do
   end do

   allocate( Temperature_PV( 3, Number_PV ) )

   do e= 1, Number_PV
      do i= 1, 3
         Temperature_PV( i, e )= Temperature( Index_Element_2_Node_PV( i, e ) )
      end do
   end do

   !====================================================================================================
   write(*,*)'      Compute Basis Finctions'
   !====================================================================================================

   allocate( BasisFunction_A( 3, Number_PV ) )
   allocate( BasisFunction_B( 3, Number_PV ) )
   allocate( BasisFunction_C( 3, Number_PV ) )

   call Compute_Basis_Function_Triangle_Element &
        ( Number_Node, Position_Node, &
          Number_PV, Index_Element_2_Node_PV, &
          !====================================================================================================
          BasisFunction_A, BasisFunction_B, BasisFunction_C )

   allocate( Gradient_Temperature( 2, Number_PV ) )

   Gradient_Temperature = 0d0

   do e= 1, Number_PV
      do i= 1, 3
         Gradient_Temperature( 1, e ) &
         = Gradient_Temperature( 1, e ) +BasisFunction_B( i, e )*Temperature_PV( i, e ) 
         Gradient_Temperature( 2, e ) &
         = Gradient_Temperature( 2, e ) +BasisFunction_C( i, e )*Temperature_PV( i, e ) 
      end do
   end do

   allocate( Heat_Flux( 2, Number_PV ) )

   do e= 1, Number_PV
      do i= 1, 2
         Heat_Flux( i, e ) = -Thermal_Conductivity_PV( e ) *Gradient_Temperature( i, e )
      end do
   end do

   deallocate( Gradient_Temperature )
   allocate( Heat_Flux_Plot( 2, Number_PV ) )
   allocate( Length_Poynting_Vector( Number_PV ) )

   Heat_Flux_Plot = Heat_Flux

   do e= 1, Number_PV
      Length_Poynting_Vector( e )= sqrt( Heat_Flux_Plot( 1, e )**2 +Heat_Flux_Plot( 2, e )**2 )
   end do

   Maximum_Length_Poynting_Vector= -1d8
   Minimum_Length_Poynting_Vector= 1d8

   do e= 1, Number_PV
      if( Maximum_Length_Poynting_Vector < Length_Poynting_Vector( e ) )then
         Maximum_Length_Poynting_Vector = Length_Poynting_Vector( e )
      end if
      if( Minimum_Length_Poynting_Vector > Length_Poynting_Vector( e ) )then
         Minimum_Length_Poynting_Vector = Length_Poynting_Vector( e )
      end if
   end do

   if( Maximum_Length_Poynting_Vector -Minimum_Length_Poynting_Vector < Minimum_Length_Poynting_Vector*1d-8 )then
      Level_PV= 0.5d0
   else
      do e= 1, Number_PV
         Level_PV( e )&
         = ( Length_Poynting_Vector( e ) -Minimum_Length_Poynting_Vector ) &
          /( Maximum_Length_Poynting_Vector -Minimum_Length_Poynting_Vector )
      end do
   end if

   do e= 1, Number_PV
      if( Level_PV( e ) < -1d-8 .or. 1d0 +1d-8 < Level_PV( e ) )then
         write(*,*)'Compute_Poynting_Vector.f90  195'
         stop
      end if
   end do

   !if( Maximum_Length_Poynting_Vector -Minimum_Length_Poynting_Vector < Minimum_Length_Poynting_Vector*1d-8 )then
   !   do e= 1, Number_PV
   !      do i= 1, 2
   !         Heat_Flux_Plot( i, e ) &
   !         = Heat_Flux_Plot( i, e )/Minimum_Length_Poynting_Vector &
   !           *(Length_Maximum_Vector +Length_Minimum_Vector )/2d0/2d0
   !      end do
   !   end do
   !else
      do e= 1, Number_PV
         do i= 1, 2
            Heat_Flux_Plot( i, e )&
            = ( Level_PV( e )*( Length_Maximum_Vector -Length_Minimum_Vector ) +Length_Minimum_Vector ) & 
             *( Heat_Flux_Plot( i, e )/Length_Poynting_Vector( e ) )
         end do
      end do
   !end if

   deallocate( Length_Poynting_Vector )

   do e= 1, Number_PV
      do i= 1, 2
         Position_Arrow_PV( i, 1, e )= Position_PV( i, e ) -Heat_Flux_Plot( i, e )
         Position_Arrow_PV( i, 2, e )= Position_PV( i, e ) +Heat_Flux_Plot( i, e )
      end do
   end do
 
   !do e= 1, Number_PV
   !   Position_Arrow_PV( 1, 1, e )= Position_PV( 1, e ) -3d-2 !Heat_Flux_Plot( i, e )
   !   Position_Arrow_PV( 1, 2, e )= Position_PV( 1, e ) +3d-2 !Heat_Flux_Plot( i, e )
   !   Position_Arrow_PV( 2, 1, e )= Position_PV( 2, e ) !-Heat_Flux_Plot( i, e )
   !   Position_Arrow_PV( 2, 2, e )= Position_PV( 2, e ) !+Heat_Flux_Plot( i, e )
   !end do

   deallocate( Heat_Flux )
   deallocate( Heat_Flux_Plot )

   !====================================================================================================
   write(*,*)'   end Compute_Poynting_Vector '
   !====================================================================================================

   return
end subroutine Compute_Poynting_Vector 



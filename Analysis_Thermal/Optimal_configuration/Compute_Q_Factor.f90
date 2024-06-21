
subroutine Compute_Q_Factor &
       ( Number_Node, Number_Element_Triangle, ID_Element, Shape_Element, &
         Position_Node, Index_Element_2_Node_Triangle, Class_Element_Triangle, &
         Field_Electric, &
         Circular_Frequency_Normalized, IncidentAngle, Dielectric_Constant_Material )

   integer, intent(in) :: ID_Element
   integer, intent(in) :: Number_Node, Number_Element_Triangle, Shape_Element
   double precision, intent(in) :: Position_Node( 2, Number_Node )
   integer, intent(in) :: Index_Element_2_Node_Triangle( Shape_Element, Number_Element_Triangle ), Class_Element_Triangle( Number_Element_Triangle )
   double precision, intent(in) :: Circular_Frequency_Normalized, IncidentAngle
   complex( kind(0d0) ), intent(in) :: Dielectric_Constant_Material 
   complex( kind( 0d0 ) ) :: Field_Electric( Number_Node )

   integer :: e, i, n

   ! Edge
   integer :: Number_Element_1_on_Boundary_1, Number_Element_2_on_Boundary_1
   integer :: Number_Element_1_on_Boundary_2, Number_Element_2_on_Boundary_2
   integer, allocatable, dimension(:) :: Element_Number_1_Boundary_1, Element_Number_1_Boundary_1_tmp
   integer, allocatable, dimension(:) :: Element_Number_2_Boundary_1, Element_Number_2_Boundary_1_tmp
   integer, allocatable, dimension(:) :: Element_Number_1_Boundary_2, Element_Number_1_Boundary_2_tmp
   integer, allocatable, dimension(:) :: Element_Number_2_Boundary_2, Element_Number_2_Boundary_2_tmp
   integer, allocatable, dimension(:,:) :: Local_Node_Number_1_Boundary_1, Local_Node_Number_2_Boundary_1
   integer, allocatable, dimension(:,:) :: Local_Node_Number_1_Boundary_2, Local_Node_Number_2_Boundary_2
   integer, allocatable, dimension(:,:) :: Local_Node_Number_1_tmp, Local_Node_Number_2_tmp
   double precision, allocatable, dimension(:,:) :: Range_Position_Node ! ( x or y, min or max)

   complex( kind( 0d0 ) ), allocatable, dimension(:) :: Value_Compute_Line_Integral

   double precision :: Value_Integrated_Poynting_Vector_1, Value_Integrated_Poynting_Vector_2
   double precision :: Transmittivity

   !====================================================================
   write(*,*) ' Detect Elements on Boundary 1'
   !====================================================================

   allocate( Element_Number_1_Boundary_1_tmp( Number_Element_Triangle ) )
   allocate( Element_Number_2_Boundary_1_tmp( Number_Element_Triangle ) )
   allocate( Local_Node_Number_1_tmp( 2, Number_Element_Triangle ) )
   allocate( Local_Node_Number_2_tmp( 2, Number_Element_Triangle ) )

   allocate( Range_Position_Node( 2, 2 ) )

   Range_Position_Node( 1, 1 )= -1d0 -1d-8 
   Range_Position_Node( 1, 2 )= -1d0 +1d-8
   Range_Position_Node( 2, 1 )= -1d10
   Range_Position_Node( 2, 2 )= 1d10

   call Detect_Line_on_Waveguide &
      ( ID_Element, Range_Position_Node, &
        Number_Node, Position_Node, &
        Number_Element_Triangle, 3, Index_Element_2_Node_Triangle, Class_Element_Triangle, &
        2, & !Flag_Position_Detect_Boundary, &
      !====================================================================================================
        Number_Element_1_on_Boundary_1, Element_Number_1_Boundary_1_tmp, Local_Node_Number_1_tmp, &
        Number_Element_2_on_Boundary_1, Element_Number_2_Boundary_1_tmp, Local_Node_Number_2_tmp )

   allocate( Element_Number_1_Boundary_1( Number_Element_1_on_Boundary_1 ) )
   allocate( Element_Number_2_Boundary_1( Number_Element_2_on_Boundary_1 ) )
   allocate( Local_Node_Number_1_Boundary_1( 2, Number_Element_1_on_Boundary_1 ) )
   allocate( Local_Node_Number_2_Boundary_1( 2, Number_Element_2_on_Boundary_1 ) )

   do e= 1, Number_Element_1_on_Boundary_1
      Element_Number_1_Boundary_1( e )= Element_Number_1_Boundary_1_tmp( e )
      do i= 1, 2
         Local_Node_Number_1_Boundary_1( i, e )= Local_Node_Number_1_tmp( i, e )
      end do
   end do

   do e= 1, Number_Element_2_on_Boundary_1
      Element_Number_2_Boundary_1( e )= Element_Number_2_Boundary_1_tmp( e )
      do i= 1, 2
         Local_Node_Number_2_Boundary_1( i, e )= Local_Node_Number_2_tmp( i, e )
      end do
   end do

   deallocate( Element_Number_1_Boundary_1_tmp )
   deallocate( Element_Number_2_Boundary_1_tmp )
   deallocate( Local_Node_Number_1_tmp )
   deallocate( Local_Node_Number_2_tmp )

   !====================================================================
   write(*,*) ' Detect Elements on Boundary 2'
   !====================================================================

   allocate( Element_Number_1_Boundary_2_tmp( Number_Element_Triangle ) )
   allocate( Element_Number_2_Boundary_2_tmp( Number_Element_Triangle ) )
   allocate( Local_Node_Number_1_tmp( 2, Number_Element_Triangle ) )
   allocate( Local_Node_Number_2_tmp( 2, Number_Element_Triangle ) )

   Range_Position_Node( 1, 1 )= 1d0 -1d-8 
   Range_Position_Node( 1, 2 )= 1d0 +1d-8
   Range_Position_Node( 2, 1 )= -1d10
   Range_Position_Node( 2, 2 )= 1d10

   call Detect_Line_on_Waveguide &
      ( ID_Element, Range_Position_Node, & 
        Number_Node, Position_Node, &
        Number_Element_Triangle, 3, Index_Element_2_Node_Triangle, Class_Element_Triangle, &
        2, & !Flag_Position_Detect_Boundary, &
      !====================================================================================================
        Number_Element_1_on_Boundary_2, Element_Number_1_Boundary_2_tmp, Local_Node_Number_1_tmp, &
        Number_Element_2_on_Boundary_2, Element_Number_2_Boundary_2_tmp, Local_Node_Number_2_tmp )

   deallocate( Range_Position_Node )

   allocate( Element_Number_1_Boundary_2( Number_Element_1_on_Boundary_2 ) )
   allocate( Element_Number_2_Boundary_2( Number_Element_2_on_Boundary_2 ) )
   allocate( Local_Node_Number_1_Boundary_2( 2, Number_Element_1_on_Boundary_2 ) )
   allocate( Local_Node_Number_2_Boundary_2( 2, Number_Element_2_on_Boundary_2 ) )

   do e= 1, Number_Element_1_on_Boundary_2
      Element_Number_1_Boundary_2( e )= Element_Number_1_Boundary_2_tmp( e )
      do i= 1, 2
         Local_Node_Number_1_Boundary_2( i, e )= Local_Node_Number_1_tmp( i, e )
      end do
   end do

   !====================================================================================================
   write(*,*)' Compute Transmittivity'
   !====================================================================================================

   allocate( Value_Compute_Line_Integral( Number_Node ) )

   do n= 1, Number_Node
      Value_Compute_Line_Integral( n )= Field_Electric( n )
   end do

   call Compute_Line_Integral_Poynting_Vector &
      ( Number_Node, Position_Node, &
        Number_Element_Triangle, 3, Index_Element_2_Node_Triangle, &
        Number_Element_1_on_Boundary_1, Element_Number_1_Boundary_1, Local_Node_Number_1_Boundary_1, &
        Number_Element_2_on_Boundary_1, Element_Number_2_Boundary_1, Local_Node_Number_2_Boundary_1, &
        Value_Compute_Line_Integral, &
      !====================================================================================================
        Value_Integrated_Poynting_Vector_1 )

   deallocate( Value_Compute_Line_Integral )

   allocate( Value_Compute_Line_Integral( Number_Node ) )

   do n= 1, Number_Node
      Value_Compute_Line_Integral( n )= Field_Electric( n )
   end do

   call Compute_Line_Integral_Poynting_Vector &
      ( Number_Node, Position_Node, &
        Number_Element_Triangle, 3, Index_Element_2_Node_Triangle, &
        Number_Element_1_on_Boundary_2, Element_Number_1_Boundary_2, Local_Node_Number_1_Boundary_2, &
        Number_Element_2_on_Boundary_2, Element_Number_2_Boundary_2, Local_Node_Number_2_Boundary_2, &
        Value_Compute_Line_Integral, &
      !====================================================================================================
        Value_Integrated_Poynting_Vector_2 )

   deallocate( Value_Compute_Line_Integral )

   Transmittivity &
   = Value_Integrated_Poynting_Vector_2/Value_Integrated_Poynting_Vector_1 

   write(*,*) '=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-='
   write( *, * ) 'Transmittivity=', Transmittivity
   write(*,*) '=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-='
   open( 1000, file='./Transmission_Waveguide.dat', position='append' )
      write( 1000, '(2x,es15.8,2x,es15.8,2x,es15.8,2x,es15.8,2x,es15.8)' ) Circular_Frequency_Normalized, IncidentAngle, real( Dielectric_Constant_Material ), Transmittivity
   close( 1000 )

   return
end subroutine Compute_Q_Factor


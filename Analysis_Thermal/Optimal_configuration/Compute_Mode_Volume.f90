
subroutine Compute_Mode_Volume&
       ( Number_Node, Number_Element_Triangle, ID_Element_Target, &
         Position_Node, Index_Element_2_Node_Triangle, Class_Element_Triangle, Relative_Permittivity, &
         Electric_Field, &
         Circular_Frequency_Normalized, IncidentAngle, Dielectric_Constant_Material )

   implicit none

   integer, intent(in) :: Number_Node, Number_Element_Triangle, ID_Element_Target 
   double precision, intent(in) :: Position_Node( 2, Number_Node ), Relative_Permittivity( Number_Element_Triangle )
   integer, intent(in) :: Index_Element_2_Node_Triangle( 3, Number_Element_Triangle ), Class_Element_Triangle( Number_Element_Triangle )
   complex( kind( 0d0 ) ), intent(in) :: Electric_Field( Number_Node )
   double precision, intent(in) :: Circular_Frequency_Normalized, IncidentAngle 
   complex( kind( 0d0 ) ), intent(in) :: Dielectric_Constant_Material 

   integer :: e, i, j
  
   double precision, allocatable, dimension(:,:) :: Basis_Function_A, Basis_Function_B, Basis_Function_C 
   double precision, allocatable, dimension(:) :: Area_Element_Triangle
   !double precision, allocatable, dimension(:,:) :: XE, YE 
   double precision, allocatable, dimension(:,:,:) :: Position_Node_Triangle_Element
   complex( kind( 0d0 ) ), allocatable, dimension(:,:,:) :: Coefficient_C, Coefficient_xi, Coefficient_eta, Coefficient_xi_eta, Coefficient_xi_xi, Coefficient_eta_eta, Integration_tmp
  
   double precision, allocatable, dimension(:) :: Electric_Amplitude, Normalized_Electric_Amplitude
  
   double precision :: Mode_Volume, Denominator_Mode_Volume, Numerator_ModeVolume, Max_Value_Electric_Amplitude
  
   !=====================================================================================================================================
   write(*,*)'   call Compute_Mode_Volume'
   !=====================================================================================================================================
  
   allocate( Electric_Amplitude( Number_Node ) )
   allocate( Normalized_Electric_Amplitude( Number_Node ) )

   do i= 1, Number_Node 
      Electric_Amplitude( i )= abs( Electric_Field( i ) )
   end do

   allocate( Position_Node_Triangle_Element( 2, 3, Number_Element_Triangle) )
 
   do e= 1 , Number_Element_Triangle
      do i= 1, 3
         do j= 1, 2 
            Position_Node_Triangle_Element( j, i, e )= Position_Node( j, Index_Element_2_Node_Triangle( i, e ) )
         end do
      end do
   end do
  
   !================================
   write(*,*)'      Minimum & Maximum of Electric_Amplitude'
   !================================
   Denominator_Mode_Volume= 1d-100
   Max_Value_Electric_Amplitude= 1d-100
  
   do e= 1, Number_Element_Triangle
      do i= 1, 3
         if( Class_Element_Triangle( e )==ID_Element_Target )then
            if( Denominator_Mode_Volume < Relative_Permittivity( e ) *Electric_Amplitude( Index_Element_2_Node_Triangle( i, e ) )**2 )THEN
               Denominator_Mode_Volume= Relative_Permittivity( e ) *Electric_Amplitude( Index_Element_2_Node_Triangle( i, e ) )**2  
            end if
            if( Max_Value_Electric_Amplitude < Electric_Amplitude( Index_Element_2_Node_Triangle( i, e ) ) )then
               Max_Value_Electric_Amplitude= Electric_Amplitude( Index_Element_2_Node_Triangle( i, e ) )
            end if
         end if
      end do
   end do
 
   write(*,*)'      Denominator_Mode_Volume=', Denominator_Mode_Volume
 
   do i= 1, Number_Node 
      Normalized_Electric_Amplitude( i )= Electric_Amplitude( i )/Max_Value_Electric_Amplitude
   end do
 
   !================================
   write(*,*)'      Compute Element Area'
   !================================
   allocate( Area_Element_Triangle( Number_Element_Triangle ) )
    
   do e=1, Number_Element_Triangle
      Area_Element_Triangle( e ) &
      = abs( ( Position_Node_Triangle_Element( 1, 1, e )-Position_Node_Triangle_Element( 1, 3, e ) )*( Position_Node_Triangle_Element( 2, 2, e )-Position_Node_Triangle_Element( 2, 3, e ) ) & 
          -( Position_Node_Triangle_Element( 1, 2, e )-Position_Node_Triangle_Element( 1, 3, e ) )*( Position_Node_Triangle_Element( 2, 1, e )-Position_Node_Triangle_Element( 2, 3, e ) ) )/2d0
   end do

   allocate( Basis_Function_A( 3, Number_Element_Triangle ) )
   allocate( Basis_Function_B( 3, Number_Element_Triangle ) )
   allocate( Basis_Function_C( 3, Number_Element_Triangle ) )

   do e= 1, Number_Element_Triangle
      Basis_Function_A( 1, e )=( Position_Node_Triangle_Element( 1, 2, e )*Position_Node_Triangle_Element( 2, 3, e )-Position_Node_Triangle_Element( 1, 3, e )*Position_Node_Triangle_Element( 2, 2, e ) )/( 2d0*Area_Element_Triangle( e ) )
      Basis_Function_A( 2, e )=( Position_Node_Triangle_Element( 1, 3, e )*Position_Node_Triangle_Element( 2, 1, e )-Position_Node_Triangle_Element( 1, 1, e )*Position_Node_Triangle_Element( 2, 3, e ) )/( 2d0*Area_Element_Triangle( e ) )
      Basis_Function_A( 3, e )=( Position_Node_Triangle_Element( 1, 1, e )*Position_Node_Triangle_Element( 2, 2, e )-Position_Node_Triangle_Element( 1, 2, e )*Position_Node_Triangle_Element( 2, 1, e ) )/( 2d0*Area_Element_Triangle( e ) )
     
      Basis_Function_B( 1, e )=( Position_Node_Triangle_Element( 2, 2, e )-Position_Node_Triangle_Element( 2, 3, e ) )/( 2d0*Area_Element_Triangle( e ) )
      Basis_Function_B( 2, e )=( Position_Node_Triangle_Element( 2, 3, e )-Position_Node_Triangle_Element( 2, 1, e ) )/( 2d0*Area_Element_Triangle( e ) )
      Basis_Function_B( 3, e )=( Position_Node_Triangle_Element( 2, 1, e )-Position_Node_Triangle_Element( 2, 2, e ) )/( 2d0*Area_Element_Triangle( e ) )
     
      Basis_Function_C( 1, e )=( Position_Node_Triangle_Element( 1, 3, e )-Position_Node_Triangle_Element( 1, 2, e ) )/( 2d0*Area_Element_Triangle( e ) )
      Basis_Function_C( 2, e )=( Position_Node_Triangle_Element( 1, 1, e )-Position_Node_Triangle_Element( 1, 3, e ) )/( 2d0*Area_Element_Triangle( e ) )
      Basis_Function_C( 3, e )=( Position_Node_Triangle_Element( 1, 2, e )-Position_Node_Triangle_Element( 1, 1, e ) )/( 2d0*Area_Element_Triangle( e ) )
   end do

   allocate( Coefficient_C( 3, 3, Number_Element_Triangle ) )
   allocate( Coefficient_xi( 3, 3, Number_Element_Triangle ) )
   allocate( Coefficient_eta( 3, 3, Number_Element_Triangle ) )
   allocate( Coefficient_xi_eta( 3, 3, Number_Element_Triangle ) )
   allocate( Coefficient_xi_xi( 3, 3, Number_Element_Triangle ) )
   allocate( Coefficient_eta_eta( 3, 3, Number_Element_Triangle ) )
  
   !================================
   write(*,*)'      Compute Coefficient_C, Coefficient_xi, Coefficient_eta ...'
   !================================
   do e= 1, Number_Element_Triangle
      do i= 1, 3
         do j= 1, 3

            Coefficient_C( j, i, e ) & 
            =  ( Basis_Function_A( i, e )+ Basis_Function_B( i, e )*Position_Node_Triangle_Element( 1, 3, e )+ Basis_Function_C( i, e )*Position_Node_Triangle_Element( 2, 3, e ) ) &
              *( Basis_Function_A( j, e )+ Basis_Function_B( j, e )*Position_Node_Triangle_Element( 1, 3, e )+ Basis_Function_C( j, e )*Position_Node_Triangle_Element( 2, 3, e ) )
         
            Coefficient_xi( j, i, e ) &
            =  ( Basis_Function_A( j, e )+ Basis_Function_B( j, e )*Position_Node_Triangle_Element( 1, 3, e )+ Basis_Function_C( j, e )*Position_Node_Triangle_Element( 2, 3, e ) ) &
              *( Basis_Function_B( i, e )*( Position_Node_Triangle_Element( 1, 1, e )-Position_Node_Triangle_Element( 1, 3, e ) )+ Basis_Function_C( i, e )*( Position_Node_Triangle_Element( 2, 1, e )-Position_Node_Triangle_Element( 2, 3, e ) ) ) &
              +( Basis_Function_A( i, e )+ Basis_Function_B( i, e )*Position_Node_Triangle_Element( 1, 3, e )+ Basis_Function_C( i, e )*Position_Node_Triangle_Element( 2, 3, e ) ) &
              *( Basis_Function_B( j, e )*( Position_Node_Triangle_Element( 1, 1, e )-Position_Node_Triangle_Element( 1, 3, e ) )+ Basis_Function_C( j, e )*( Position_Node_Triangle_Element( 2, 1, e )-Position_Node_Triangle_Element( 2, 3, e ) ) )
         
            Coefficient_eta( j, i, e ) &
            =  ( Basis_Function_A( j, e )+ Basis_Function_B( j, e )*Position_Node_Triangle_Element( 1, 3, e )+ Basis_Function_C( j, e )*Position_Node_Triangle_Element( 2, 3, e ) ) &
              *( Basis_Function_B( i, e )*( Position_Node_Triangle_Element( 1, 2, e )-Position_Node_Triangle_Element( 1, 3, e ) )+ Basis_Function_C( i, e )*( Position_Node_Triangle_Element( 2, 2, e )-Position_Node_Triangle_Element( 2, 3, e ) ) ) &
              +( Basis_Function_A( i, e )+ Basis_Function_B( i, e )*Position_Node_Triangle_Element( 1, 3, e )+ Basis_Function_C( i, e )*Position_Node_Triangle_Element( 2, 3, e ) ) &
              *( Basis_Function_B( j, e )*( Position_Node_Triangle_Element( 1, 2, e )-Position_Node_Triangle_Element( 1, 3, e ) )+ Basis_Function_C( j, e )*( Position_Node_Triangle_Element( 2, 2, e )-Position_Node_Triangle_Element( 2, 3, e ) ) )
         
            Coefficient_xi_eta( j, i, e ) &
            =  ( Basis_Function_B( i, e )*( Position_Node_Triangle_Element( 1, 2, e )-Position_Node_Triangle_Element( 1, 3, e ) )+ Basis_Function_C( i, e )*( Position_Node_Triangle_Element( 2, 2, e )-Position_Node_Triangle_Element( 2, 3, e ) ) ) &
              *( Basis_Function_B( j, e )*( Position_Node_Triangle_Element( 1, 1, e )-Position_Node_Triangle_Element( 1, 3, e ) )+ Basis_Function_C( j, e )*( Position_Node_Triangle_Element( 2, 1, e )-Position_Node_Triangle_Element( 2, 3, e ) ) ) &
              +( Basis_Function_B( i, e )*( Position_Node_Triangle_Element( 1, 1, e )-Position_Node_Triangle_Element( 1, 3, e ) )+ Basis_Function_C( i, e )*( Position_Node_Triangle_Element( 2, 1, e )-Position_Node_Triangle_Element( 2, 3, e ) ) ) &
              *( Basis_Function_B( j, e )*( Position_Node_Triangle_Element( 1, 2, e )-Position_Node_Triangle_Element( 1, 3, e ) )+ Basis_Function_C( j, e )*( Position_Node_Triangle_Element( 2, 2, e )-Position_Node_Triangle_Element( 2, 3, e ) ) )
         
            Coefficient_xi_xi( j, i, e ) &
            =  ( Basis_Function_B( i, e )*( Position_Node_Triangle_Element( 1, 1, e )-Position_Node_Triangle_Element( 1, 3, e ) )+ Basis_Function_C( i, e )*( Position_Node_Triangle_Element( 2, 1, e )-Position_Node_Triangle_Element( 2, 3, e ) ) ) &
              *( Basis_Function_B( j, e )*( Position_Node_Triangle_Element( 1, 1, e )-Position_Node_Triangle_Element( 1, 3, e ) )+ Basis_Function_C( j, e )*( Position_Node_Triangle_Element( 2, 1, e )-Position_Node_Triangle_Element( 2, 3, e ) ) )
          
            Coefficient_eta_eta( j, i, e ) &
            =  ( Basis_Function_B( i, e )*( Position_Node_Triangle_Element( 1, 2, e )-Position_Node_Triangle_Element( 1, 3, e ) )+ Basis_Function_C( i, e )*( Position_Node_Triangle_Element( 2, 2, e )-Position_Node_Triangle_Element( 2, 3, e ) ) ) &
              *( Basis_Function_B( j, e )*( Position_Node_Triangle_Element( 1, 2, e )-Position_Node_Triangle_Element( 1, 3, e ) )+ Basis_Function_C( j, e )*( Position_Node_Triangle_Element( 2, 2, e )-Position_Node_Triangle_Element( 2, 3, e ) ) )

         end do
      end do
   end do
  
   deallocate( Basis_Function_A )
   deallocate( Basis_Function_B )
   deallocate( Basis_Function_C )
  
   !================================
   write(*,*)'      Compute Integration_tmp ...'
   !================================
   allocate( Integration_tmp(3, 3, Number_Element_Triangle) )
  
   do e= 1, Number_Element_Triangle
      do i= 1, 3
         do j= 1, 3
            Integration_tmp( j, i, e ) &
            = 1d0/2d0*Coefficient_C( j, i, e ) +1d0/6d0*Coefficient_xi( j, i, e )+1d0/6d0*Coefficient_eta( j, i, e ) &
             +1d0/12d0*Coefficient_xi_xi( j, i, e )+ 1d0/12d0*Coefficient_eta_eta( j, i, e )+1d0/24d0*Coefficient_xi_eta( j, i, e )
         end do
      end do
   end do
  
   deallocate( Coefficient_C )
   deallocate( Coefficient_xi )
   deallocate( Coefficient_eta )
   deallocate( Coefficient_xi_eta )
   deallocate( Coefficient_xi_xi )
   deallocate( Coefficient_eta_eta )
  
   !================================
   write(*,*)'      Compute Mode Volume'
   !================================
   Numerator_ModeVolume= 0d0

   do e= 1, Number_Element_Triangle
      if( Class_Element_Triangle( e )==ID_Element_Target )then
         do i= 1, 3
            do j= 1, 3
               Numerator_ModeVolume= Numerator_ModeVolume + Relative_Permittivity( e )*2d0*Area_Element_Triangle( e ) *real( Electric_Field( Index_Element_2_Node_Triangle( i, e ) )*CONJG( Electric_Field( Index_Element_2_Node_Triangle( j, e ) ) ) )*Integration_tmp( j, i, e ) 
            end do
         end do
      end if
   end do

   Mode_Volume= Numerator_ModeVolume/Denominator_Mode_Volume
 
   write(*,*)'      ModeVolume=', Mode_Volume, Numerator_ModeVolume, Denominator_Mode_Volume  !Numerator_ModeVolume/Denominator_Mode_Volume !( Relative_Permittivity( MaxClass_Element_Triangle_Triangle )*Denominator_Mode_Volume*Denominator_Mode_Volume ) *Omega*Omega

   open( 1000, file='./Mode_Volume.dat', position='append' ) 
      write( 1000, '(2x,es15.8,2x,es15.8,2x,es15.8,2x,es15.8,2x,es15.8)' ) Circular_Frequency_Normalized, IncidentAngle, real( Dielectric_Constant_Material ), Mode_Volume, Mode_Volume*Circular_Frequency_Normalized**2
   close( 1000 )
      write( 1001, '(2x,es15.8,2x,es15.8,2x,es15.8,2x,es15.8,2x,es15.8,2x,es15.8)' ) Circular_Frequency_Normalized, Mode_Volume, Numerator_ModeVolume, Denominator_Mode_Volume, Max_Value_Electric_Amplitude

   deallocate( Area_Element_Triangle )
  
   deallocate( Normalized_Electric_Amplitude )
   deallocate( Position_Node_Triangle_Element )
   deallocate( Integration_tmp )
  
end 


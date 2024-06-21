

subroutine Compute_Incident_Wave &
      ( Flag_Incident, Position_Node, Class_Node, Number_Node, & 
        Frequency_Normalized, Incident_Angle, &
        ID_Node_Material, ID_Node_Air, ID_Node_OpenRegion, ID_Node_Active, ID_Node_FixedDomain, &
        ID_Node_PML_X, ID_Node_PML_Y, ID_Node_PML_XY, &
        Dielectric_Constant_OpenRegion, &
        Position_Dipole_x, Position_Dipole_Y, &
       !==================================================================
        Field_Electric_Incident )

   !$ use omp_lib
   use Parameters
   implicit none

   integer, intent(in) :: Flag_Incident
   ! 'PW':plane wave, 'PS':point source, 'DR':dipole radiation 

   integer :: i 

   double precision, intent(in) :: Frequency_Normalized 
   double precision, intent(in) :: Incident_Angle 
   integer, intent(in) :: Number_Node, ID_Node_PML_X, ID_Node_PML_Y, ID_Node_PML_XY
   integer, intent(in) :: ID_Node_Material, ID_Node_Air, ID_Node_OpenRegion, ID_Node_Active, ID_Node_FixedDomain 
   double precision, intent(in) :: Position_Node( 2, Number_Node )
   double precision, intent(in) :: Position_Dipole_x, Position_Dipole_Y 
   integer, intent(in) :: Class_Node( Number_Node )
   double precision, intent(in) :: Dielectric_Constant_OpenRegion
   complex( kind( 0d0 ) ), intent(out) :: Field_Electric_Incident( Number_Node )

   double precision, allocatable, dimension(:) :: Distance_Dipole_2_Node, Position_Dipole
   double precision, allocatable, dimension(:) :: J0, Y0  
   complex( kind( 0d0 ) ), allocatable, dimension(:) :: H0 

   double precision :: Wave_Number_x, Wave_Number_y
   double precision :: Wave_Number_x_Lower, Wave_Number_y_Lower
   double precision :: Wave_Number_x_Higher, Wave_Number_y_Higher
   double precision :: dbesj0, dbesy0 

   double precision :: Distance_Dipole_2_Center_DesignDomain
   double precision :: J0_Center_DesignDomain, Y0_Center_DesignDomain, Amplitude_Dipole_Center_DesignDomain
   complex( kind( 0d0 ) ) :: H0_Center_DesignDomain 

   complex( kind( 0d0 ) ), allocatable, dimension(:) :: Complex_IncidentAngle_Tmp, Complex_Frequency_Tmp
   complex( kind( 0d0 ) ), allocatable, dimension(:,:) :: Exp_Frequency_Difference


   !========================================================================
   write(*,*)'         call Compute_Incident_Wave '
   !========================================================================

   !==================================
   if( Flag_Incident==0 )then
   !==================================
      write(*,*)'            Incident Wave : Plane Wave '
      write(*,*)'               Frequency_Normalized=', Frequency_Normalized

      Wave_Number_x= 2.0d0*Pi*Frequency_Normalized*dcos( Incident_Angle ) *sqrt( Dielectric_Constant_OpenRegion ) 
      Wave_Number_y= 2.0d0*Pi*Frequency_Normalized*dsin( Incident_Angle ) *sqrt( Dielectric_Constant_OpenRegion ) 

      do i= 1, Number_Node
         if( Class_Node( i ) == ID_Node_PML_X .or. &
             Class_Node( i ) == ID_Node_PML_Y .or. &
             Class_Node( i ) == ID_Node_PML_XY )then

            Field_Electric_Incident( i )= Zero
         else
            Field_Electric_Incident( i ) & 
            = Amplitude_Incident_Wave & 
              *exp( Imaginary_Unit &
                  *( Wave_Number_x*Position_Node( 1, i ) +Wave_Number_y*Position_Node( 2, i ) ) )
         end if
      end do

!   !==================================
!   else if( Flag_Incident==20 )then
!   !==================================
!      write(*,*)'            Incident Wave : Plane Wave for Wide Angle '
!      write(*,*)'               Incident_Angle_Radian_Lower=', Incident_Angle_Radian_Lower
!      write(*,*)'               Incident_Angle_Radian_Higher=', Incident_Angle_Radian_Higher
!
!      allocate( Complex_IncidentAngle_Tmp( Number_Node ) )
!      
!      !$omp parallel do default( none )   &
!      !$omp private( i ) &
!      !$omp shared( Number_Node, Field_Electric_Incident, Position_Node, Class_Node, Frequency_Normalized ) 
!      do i= 1, Number_Node
!         if( Class_Node( i ) == ID_Node_PML_X .or. &
!             Class_Node( i ) == ID_Node_PML_Y .or. &
!             Class_Node( i ) == ID_Node_PML_XY )then
!
!            Field_Electric_Incident( i )= Zero
!         else
!
!            call Trapezoidal_Rule_Complex &
!               ( Incident_Angle_Radian_Lower, Incident_Angle_Radian_Higher, 1000, & 
!                 Position_Node( 1, i ), Position_Node( 2, i ), & 
!                 Amplitude_Incident_Wave, 2.0d0*Pi*Frequency_Normalized, &
!               !==================================================
!                 Field_Electric_Incident( i ) )
!
!         end if
!      end do
!
!      deallocate( Complex_IncidentAngle_Tmp )

   !==================================
   else if( Flag_Incident==100 )then
   !==================================
      write(*,*)'            Incident Wave : Dipole Radiation '
      write(*,*)'               Position_Center_DesignDomain_X=', Position_Center_DesignDomain_X
      write(*,*)'               Position_Center_DesignDomain_Y=', Position_Center_DesignDomain_Y

      allocate( Position_Dipole( 2 ) )

      Position_Dipole( 1 )= Position_Dipole_x
      Position_Dipole( 2 )= Position_Dipole_y

      allocate( Distance_Dipole_2_Node( Number_Node ) )
      allocate( J0( Number_Node ) )
      allocate( Y0( Number_Node ) )
      allocate( H0( Number_Node ) )

      do i= 1, Number_Node
         Distance_Dipole_2_Node( i ) & 
         = sqrt( ( Position_Node( 1, i ) -Position_Dipole( 1 ) )**2 &
              +( Position_Node( 2, i ) -Position_Dipole( 2 ) )**2  )   
      
         if( Distance_Dipole_2_Node( i )==0.0d0 )then
            Distance_Dipole_2_Node( i )= Tolerance_Digit_Error 
         end if
      end do

      do i= 1, Number_Node
         if( Class_Node( i ) == ID_Node_PML_X .or. &
             Class_Node( i ) == ID_Node_PML_Y .or. &
             Class_Node( i ) == ID_Node_PML_XY )then

            Field_Electric_Incident( i )= Zero
         else
 
            J0( i )=dbesj0( 2.0D0*Pi*Frequency_Normalized*Distance_Dipole_2_Node( i ) )
            Y0( i )=dbesy0( 2.0D0*Pi*Frequency_Normalized*Distance_Dipole_2_Node( i ) ) 
            H0( i )= J0( i ) +Imaginary_Unit*Y0( i ) 
 
            Field_Electric_Incident( i )= -( 2.0D0*Pi*Frequency_Normalized )**2 *Imaginary_Unit*H0( i )/4.0D0
         end if
      end do

      Distance_Dipole_2_Center_DesignDomain & 
      = sqrt( ( Position_Center_DesignDomain_X -Position_Dipole( 1 ) )**2 &
           +( Position_Center_DesignDomain_Y -Position_Dipole( 2 ) )**2  )   

      J0_Center_DesignDomain=dbesj0( 2.0D0*Pi*Frequency_Normalized*Distance_Dipole_2_Center_DesignDomain )
      Y0_Center_DesignDomain=dbesy0( 2.0D0*Pi*Frequency_Normalized*Distance_Dipole_2_Center_DesignDomain ) 
      H0_Center_DesignDomain= J0_Center_DesignDomain +Imaginary_Unit*Y0_Center_DesignDomain 

      Amplitude_Dipole_Center_DesignDomain &
      = abs( -( 2d0*Pi*Frequency_Normalized )**2 *Imaginary_Unit*H0_Center_DesignDomain/4d0 )

      write(*,*)'   Amplitude_Dipole_Center_DesignDomain=', Amplitude_Dipole_Center_DesignDomain
      if( Amplitude_Dipole_Center_DesignDomain==0.0d0 )then
         call Output_Error( 'Compute_Incident_Wave', 156 )
      end if

      do i= 1, Number_Node
         Field_Electric_Incident( i )= Field_Electric_Incident( i )/Amplitude_Dipole_Center_DesignDomain 
      end do

      deallocate( Position_Dipole )
      deallocate( Distance_Dipole_2_Node )
      deallocate( J0 )
      deallocate( Y0 )
      deallocate( H0 )

   !==================================
   else if( Flag_Incident==200 )then
   !==================================
      write(*,*)'            Incident Wave : Posint Source '

      allocate( Position_Dipole( 2 ) )

      Position_Dipole( 1 )= Position_Dipole_x
      Position_Dipole( 2 )= Position_Dipole_y

      allocate( Distance_Dipole_2_Node( Number_Node ) )
      allocate( J0( Number_Node ) )
      allocate( Y0( Number_Node ) )
      allocate( H0( Number_Node ) )

      do i= 1, Number_Node
         Distance_Dipole_2_Node( i ) & 
         = sqrt( ( Position_Node( 1, i ) -Position_Dipole( 1 ) )**2 &
              +( Position_Node( 2, i ) -Position_Dipole( 2 ) )**2  )  
 
         if( Distance_Dipole_2_Node( i )==0.0d0 )then
            Distance_Dipole_2_Node( i )= Tolerance_Digit_Error 
         end if
      end do

      do i= 1, Number_Node
         if( Class_Node( i ) == ID_Node_Material .or. &
             Class_Node( i ) == ID_Node_Air      .or. &
             Class_Node( i ) == ID_Node_OpenRegion .or. & 
             Class_Node( i ) == ID_Node_Active   .or. &
             Class_Node( i ) == ID_Node_FixedDomain )then

            J0( i )=dbesj0( 2.0D0*Pi*Frequency_Normalized*Distance_Dipole_2_Node( i ) ) 
            Y0( i )=dbesy0( 2.0D0*Pi*Frequency_Normalized*Distance_Dipole_2_Node( i ) ) 
            H0( i )= J0( i ) +Imaginary_Unit*Y0( i ) 
 
            Field_Electric_Incident( i )= Imaginary_Unit*H0( i )/4.0D0

         else
            Field_Electric_Incident( i )= Zero
         end if
      end do

      deallocate( Position_Dipole )
      deallocate( Distance_Dipole_2_Node )
      deallocate( J0 )
      deallocate( Y0 )
      deallocate( H0 )


   !==================================
   else
   !==================================
      write(*,*)'ERROR'
      write(*,*)'  Flag_Incident', Flag_Incident
      write(*,*)'  Compute_Incident_Wave.f90 59 '
      stop
   end if

   return
end subroutine Compute_Incident_Wave


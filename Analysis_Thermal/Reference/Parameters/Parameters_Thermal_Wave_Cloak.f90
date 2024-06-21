
module Parameters  
   implicit none    

   double precision, parameter :: Pi=3.14159265358979d0
   double precision, parameter :: Light_Speed_Vacuum=299792458d0
   double precision, parameter :: Eps0=8.854187817D-12
   double precision, parameter :: Mu0=1.0d0/(Light_Speed_Vacuum*Light_Speed_Vacuum*Eps0)
   
   complex( kind(0d0) ), parameter :: Imaginary_Unit= dcmplx( 0.0d0, 1.0d0 )
   complex( kind(0d0) ), parameter :: One= dcmplx( 1.0d0, 0.0d0 )
   complex( kind(0d0) ), parameter :: Zero= dcmplx( 0.0d0, 0.0d0 ) 

   integer, parameter :: Flag_Thermal_Device= 0
   ! 0: Cloak
   ! 10: Carpet Cloak 
   ! 20: Cloak Concentrator
   integer, parameter :: Flag_Material_FixedDomain= 0
   integer, parameter :: Number_Objective_Function= 1 
   integer, parameter :: Flag_Refine_Element= 0 

   !integer, parameter ::  Flag_No_Structure= 1 ! 0: off, 1: On

   integer, parameter :: Flag_On= 1
   integer, parameter :: Flag_Off= 0

   !===============================================================================================
   ! Paramters for Plot 
   !===============================================================================================

   ! 0:Blue-Green_Red
   ! 1:Blue-Green_Red
   ! 2:darkBlue-darkGreen_darkRed
   ! 3:darkBlue-darkGreen_darkRed
   ! 5:DarkBlue-Blue-Green_Red-DarkRed
   integer, parameter :: Flag_RGB_Color= 0! 0
   integer, parameter :: Flag_RGB_Color_Diff= 1 
   integer, parameter :: Flag_RGB_Color_OC= 1 ! 1
   integer, parameter :: Flag_RGB_Color_Gradient= 5 

   integer, parameter :: Flag_Range_Value_Plot= 1
   integer, parameter :: Flag_Range_Value_Plot_Gradient= 1

      double precision, parameter :: Maximum_Value_Plotted_Fixed= 1d0 
      double precision, parameter :: Minimum_Value_Plotted_Fixed= 0d0

      double precision, parameter :: Maximum_Value_Plotted_Diff= 1d-2
      double precision, parameter :: Minimum_Value_Plotted_Diff= -Maximum_Value_Plotted_Diff 

      double precision, parameter :: Maximum_Value_Plotted_Gradient= 5d0 
      double precision, parameter :: Minimum_Value_Plotted_Gradient= 0d0

      !double precision, parameter :: RGB_Lower_R= 0.0d0
      !double precision, parameter :: RGB_Lower_G= 0.0d0
      !double precision, parameter :: RGB_Lower_B= 0.3d0

      !double precision, parameter :: RGB_Higher_R= 0.3d0
      !double precision, parameter :: RGB_Higher_G= 0.0d0
      !double precision, parameter :: RGB_Higher_B= 0.0d0 

   !===============================================================================================
   ! Paramters for Plot 
   !===============================================================================================
   integer, parameter :: Number_Circle_Initial_LSF= 3
   double precision, parameter :: Ratio_InitialConfiguration_2_DesignDomain= 1d0

   integer, parameter :: Number_Level_Plot_EM= 100
   integer, parameter :: Max_Number_Level_Element_Limit= Number_Level_Plot_EM +1

   integer, parameter :: Number_Grid_Scale_DesignDomain= 100

   double precision, parameter :: Ratio_DesignDomain_2_Scattering= 2d0

   integer, parameter :: Number_Grid_X_DesignDomain= Number_Grid_Scale_DesignDomain 
   integer, parameter :: Number_Grid_Y_DesignDomain= Number_Grid_Scale_DesignDomain

   integer, parameter :: Number_Grid_X_FixedDomain= Number_Grid_Scale_DesignDomain/3d0 
   integer, parameter :: Number_Grid_Y_FixedDomain= Number_Grid_Scale_DesignDomain/3d0
   !integer, parameter :: Number_Grid_X_FixedDomain= Number_Grid_Scale_DesignDomain*3d0/5d0 
   !integer, parameter :: Number_Grid_Y_FixedDomain= Number_Grid_Scale_DesignDomain*3d0/5d0

   integer, parameter :: Number_Grid_X_Scattering= Number_Grid_X_DesignDomain *Ratio_DesignDomain_2_Scattering
   !integer, parameter :: Number_Grid_Y_Scattering= Number_Grid_Y_DesignDomain *Ratio_DesignDomain_2_Scattering 
   integer, parameter :: Number_Grid_Y_Scattering= Number_Grid_Y_DesignDomain +Number_Grid_Y_FixedDomain

   double precision, parameter :: Length_Characteristic= dble( Number_Grid_Scale_DesignDomain )
   double precision, parameter :: Length_Characteristic_Normalized= 1d0 

   double precision, parameter :: Size_X_DesignDomain= dble( Number_Grid_X_DesignDomain ) 
   double precision, parameter :: Size_Y_DesignDomain= dble( Number_Grid_Y_DesignDomain ) 
   double precision, parameter :: Radius_DesignDomain= dble( Number_Grid_Scale_DesignDomain ) 
   
   double precision, parameter :: Size_X_Scattering= dble( Number_Grid_X_Scattering ) 
   double precision, parameter :: Size_Y_Scattering= dble( Number_Grid_Y_Scattering ) 
   
   double precision, parameter :: Size_X_FixedDomain= dble( Number_Grid_X_FixedDomain ) 
   double precision, parameter :: Size_Y_FixedDomain= dble( Number_Grid_Y_FixedDomain ) 
   double precision, parameter :: Radius_FixedDomain= Number_Grid_Scale_DesignDomain*3d0/5d0 
   !double precision, parameter :: Radius_FixedDomain= Size_X_DesignDomain/3d0 
   !double precision, parameter :: Radius_FixedDomain= 0d0 

   double precision, parameter :: Position_Center_FixedDomain_X= 0d0 
   double precision, parameter :: Position_Center_FixedDomain_Y= 0d0
   double precision, parameter :: Position_Center_DesignDomain_X= 0d0 
   double precision, parameter :: Position_Center_DesignDomain_Y= 0d0   

   !================================================================================================
   ! Fundamental Settings
   !================================================================================================
   !integer, parameter ::  Number_Loop_Frequency= 6 ! multiple of 6
      double precision, parameter :: Minimum_Temperature_Lower_Side_Boundary= 0d0 
      double precision, parameter :: Maximum_Temperature_Lower_Side_Boundary= 1d0
      double precision, parameter :: Normal_Temperature_Lower_Side_Boundary= 0d0

   !integer, parameter ::  Number_Loop_Temperature_Right_Side_Boundary= 0
      double precision, parameter :: Minimum_Temperature_Higher_Side_Boundary= 0d0 !-Pi/2.0d0 
      double precision, parameter :: Maximum_Temperature_Higher_Side_Boundary= 1d0 !Pi/2.0d0
      double precision, parameter :: Normal_Temperature_Higher_Side_Boundary= 1d0

   !integer, parameter :: Number_Loop_DielectricConstant= 0
      double precision, parameter :: Minimum_Thermal_Conductivity= 1d0 
      double precision, parameter :: Maximum_Thermal_Conductivity= 10d0
      double precision, parameter :: Normal_Thermal_Conductivity_Material= 204d0/67d0
 
      double precision, parameter :: Thermal_Conductivity_Base_Material= 1.0d0 
      double precision, parameter :: Thermal_Conductivity_OpenRegion= 1.0d0 
      double precision, parameter :: Thermal_Conductivity_FixedDomain= 0.0d0 
 
      double precision, parameter :: Minimum_Position_Source_x= -3d0
      double precision, parameter :: Maximum_Position_Source_x= -1d0
      double precision, parameter :: Normal_Position_Source_x= -2d0

      double precision, parameter :: Minimum_Position_Source_y= -3d0
      double precision, parameter :: Maximum_Position_Source_y=  3d0
      double precision, parameter :: Normal_Position_Source_y= 0d0 

   !character(len=2), parameter :: Type_Incident_Wave='PW' 

   double precision, parameter :: Amplitude_Incident_Wave= 1.0d0

   !===============================================================================================
   ! Postscript
   !===============================================================================================

   ! Check on Desktop
   integer, parameter :: WidthPS= 360
   integer, parameter :: WidthPS_X= WidthPS*Number_Grid_X_Scattering/Number_Grid_X_Scattering 
   !integer, parameter :: WidthPS_Y= WidthPS*Number_Grid_Y_Scattering/Number_Grid_X_Scattering
   integer, parameter :: WidthPS_Y= WidthPS_X*Number_Grid_Y_Scattering/Number_Grid_X_Scattering 
   integer, parameter :: Translation_X= 10
   integer, parameter :: Translation_Y= 20
   integer, parameter :: Width_Color_Var= WidthPS/20

   integer, parameter :: Flag_Color_Var= 2 
   ! 0:off, 1:Right Side, 2:Upper Side

   !===============================================================================================
   ! Paramters for Plot 
   !===============================================================================================
   double precision, parameter :: &
   Position_Center_FixedDomain_X_PS &
   = Position_Center_FixedDomain_X/( 2d0*Size_X_Scattering ) *dble( WidthPS_X ) &
    +dble( WidthPS_X )/2.0d0 +Translation_X

   double precision, parameter :: &
   Position_Center_FixedDomain_Y_PS &
   = Position_Center_FixedDomain_Y/( 2d0*Size_Y_Scattering ) *dble( WidthPS_Y ) &
    +dble( WidthPS_Y )/2.0d0 +Translation_Y

   double precision, parameter :: &
   Position_Center_DesignDomain_X_PS &
   = Position_Center_DesignDomain_X/( 2d0*Size_X_Scattering ) *dble( WidthPS_X ) &
    +dble( WidthPS_X )/2.0d0 +Translation_X

   double precision, parameter :: &
   Position_Center_DesignDomain_Y_PS &
   = Position_Center_DesignDomain_Y/( 2d0*Size_Y_Scattering ) *dble( WidthPS_Y ) &
    +dble( WidthPS_Y )/2.0d0 +Translation_Y

   double precision, parameter :: &
   Radius_FixedDomain_Plot_PS &
   = Radius_FixedDomain/( 2d0*Size_X_Scattering ) *dble( WidthPS_X )

   double precision, parameter :: &
   Radius_DesignDomain_Plot_PS &
   = Radius_DesignDomain/( 2d0*Size_X_Scattering ) *dble( WidthPS_X )


   double precision, parameter :: &
   Position_Center_WholeDomain_X_PS &
   = Position_Center_FixedDomain_X_PS 

   double precision, parameter :: &
   Position_Center_WholeDomain_Y_PS &
   = Position_Center_FixedDomain_Y_PS 

   double precision, parameter :: &
   Position_Corner_WholeDomain_X1_PS & 
   = Position_Center_WholeDomain_X_PS -Size_X_Scattering/( 2d0*Size_X_Scattering ) *dble( WidthPS_X ) 
 
   double precision, parameter :: &
   Position_Corner_WholeDomain_Y1_PS & 
   = Position_Center_WholeDomain_Y_PS -Size_Y_Scattering/( 2d0*Size_Y_Scattering ) *dble( WidthPS_Y ) 
 
   double precision, parameter :: &
   Position_Corner_WholeDomain_X2_PS & 
   = Position_Center_WholeDomain_X_PS +Size_X_Scattering/( 2d0*Size_X_Scattering ) *dble( WidthPS_X ) 
 
   double precision, parameter :: &
   Position_Corner_WholeDomain_Y2_PS & 
   = Position_Center_WholeDomain_Y_PS +Size_Y_Scattering/( 2d0*Size_Y_Scattering ) *dble( WidthPS_Y ) 





   double precision, parameter :: &
   Position_Corner_FixedDomain_X1_PS & 
   = Position_Center_FixedDomain_X_PS -Size_X_FixedDomain/( 2d0*Size_X_Scattering ) *dble( WidthPS_X ) 
 
   double precision, parameter :: &
   Position_Corner_FixedDomain_Y1_PS & 
   = Position_Center_FixedDomain_Y_PS -Size_Y_FixedDomain/( 2d0*Size_Y_Scattering ) *dble( WidthPS_Y ) 
 
   double precision, parameter :: &
   Position_Corner_FixedDomain_X2_PS & 
   = Position_Center_FixedDomain_X_PS +Size_X_FixedDomain/( 2d0*Size_X_Scattering ) *dble( WidthPS_X ) 
 
   double precision, parameter :: &
   Position_Corner_FixedDomain_Y2_PS & 
   = Position_Center_FixedDomain_Y_PS +Size_Y_FixedDomain/( 2d0*Size_Y_Scattering ) *dble( WidthPS_Y ) 

   double precision, parameter :: &
   Position_Corner_DesignDomain_X1_PS & 
   = Position_Center_DesignDomain_X_PS -Size_X_DesignDomain/( 2d0*Size_X_Scattering ) *dble( WidthPS_X ) 
 
   double precision, parameter :: &
   Position_Corner_DesignDomain_Y1_PS & 
   = Position_Center_DesignDomain_Y_PS -Size_Y_DesignDomain/( 2d0*Size_Y_Scattering ) *dble( WidthPS_Y ) 
 
   double precision, parameter :: &
   Position_Corner_DesignDomain_X2_PS & 
   = Position_Center_DesignDomain_X_PS +Size_X_DesignDomain/( 2d0*Size_X_Scattering ) *dble( WidthPS_X ) 
 
   double precision, parameter :: &
   Position_Corner_DesignDomain_Y2_PS & 
   = Position_Center_DesignDomain_Y_PS +Size_Y_DesignDomain/( 2d0*Size_Y_Scattering ) *dble( WidthPS_Y ) 

   integer, parameter :: Integer_Initialization= -10D8
   double precision, parameter :: Double_Precision_Initialization= -10D10

   integer, parameter :: Integer_Initialization_Plus= 10D8
   integer, parameter :: Integer_Initialization_Minus= -10D8

   double precision, parameter :: Double_Precision_Initialization_Plus= 10D10
   double precision, parameter :: Double_Precision_Initialization_Minus= -10D10
 
   integer, parameter :: Flag_DesignDomain= 1
   ! 1: Circle
   ! 2: Square
   
   integer, parameter :: Flag_FixedDomain= 1
   ! 0: Vacant
   ! 1: Circle
   ! 2: Square

   integer, parameter :: Number_Optimization_Step= 2000
   integer, parameter :: Length_Character_Optimization_Step= log10(dble(Number_Optimization_Step)) +1

   !===============================================================================================
   ! Tolerance
   !===============================================================================================

   double precision, parameter :: Tolerance_Position_Node= 1.0d-13 !10d-5 
   double precision, parameter :: Tolerance_Boundary_Group= 2.0d0*Tolerance_Position_Node 
   double precision, parameter :: Tolerance_LSF= 1.0d-8
   double precision, parameter :: Tolerance_Cross_Check= 1.0d-11

   double precision, parameter :: Tolerance_Digit_Error= 1.0d-8
   double precision, parameter :: Tolerance_Digit_Error_Plus= 1.0d0 +Tolerance_Digit_Error
   double precision, parameter :: Tolerance_Digit_Error_Minus= 1.0d0 -Tolerance_Digit_Error

   double precision, parameter :: Tolerance_Plot= 1.0d-2

end module Parameters 




module Parameters  
   implicit none    
 
   !================================================================================================
   ! Fundamental Settings
   !================================================================================================
   integer, parameter :: Number_Optimization_Step= 1d5 ! set
   character(len=5), parameter :: Flag_Read_LSF = 'false'
   character(len=5), parameter :: Flag_xmean = 'true' ! true

   integer, parameter :: Number_Grid_Scale_DesignDomain= 60 ! set ! 3N 
   double precision, parameter :: Radius_Search= 0.3d0 ! 0.3( x_max - x_min )=0.6

   integer, parameter :: Flag_Shape_Expression = 2
   ! 1 : nodal
   ! 2 : Radial Basis Function
   ! 3 : Spectral Level Set Method ( Fourier series expantion )

   character(len=9), parameter :: Type_CMA = 'Full'
   ! F : Normal CMA, D : Separable CMA ( Diagonal ), S : Sparse CMA, 
   ! ??? : Sepctral LSM based CMA ( Fourier series expantion) 

   integer, parameter :: Number_Neighboring_GP = 8

   integer, parameter :: Flag_Box_Constraint= 2 !set 
   ! 0:Off, 1:v1, 2:v2

   integer, parameter :: Number_Candidates = 200

   integer, parameter :: Number_Obj_Func= 1
   integer, parameter :: Loop_Cloak= 1 
   integer, parameter :: Loop_Flux_X= 2 
   integer, parameter :: Loop_Flux_Y= 3 

   integer, parameter :: Number_Obj_Func_Source= 1
   integer, parameter :: Number_Obj_Func_MP_T= 1 
   integer, parameter :: Number_Obj_Func_MP_V= 1 

   integer, parameter :: Number_Type_LSF= 1

   logical, parameter :: Flag_Ordering= .true.
   integer, parameter :: Seed_Number= 1
   integer, parameter :: Flag_CPU_Time= 1
   ! 0:Off, 1:On

   integer, parameter :: Flag_MultiPhysics_Device= 0 
   ! 0: Electric Current Cloak
   integer, parameter :: Flag_ObjectiveFunction= 11
   ! 0: intensity of Total Field
   ! 11: intensity of ( Scattered Field -Reference Field )
   integer, parameter :: Flag_Reference_OF= 1
   character(len=3), parameter :: Type_ObjectiveFunction= 'sum'
   ! sum: Normal
   ! max: max 
   ! div: Psi_1/Psi_2
   character(len=3), parameter :: Type_Fitness= 'sum'
   ! sum: summation 
   ! max: max

   integer, parameter :: Flag_Thermal_Insulation_BC= 1
   ! 0: off, 1: on
   integer, parameter :: Flag_Electrical_Insulation_BC= 1
   ! 0: off, 1: on

   integer, parameter :: Flag_Physics= 1 
   ! 0: Bifunctional, 1: Thermal flow, 2: DC electric
   integer, parameter :: Flag_Symmetry_LSF= 4 !2
   ! 0: off, 1: x-axis, 2: y-axis, 4: x & y axis, 8: x, y axes, y=x
   integer, parameter :: Flag_Symmetry_Sensitivity= 0
   ! 0: off, 1: x-axis, 2: Rotational, 3: Origin
   integer, parameter :: Flag_Incident_Wave= 0 
   ! 0: Plane Wave
   ! 1: Plane Wave - Multifrequency 
   ! 10: Plane Wave for Wide Range
   ! 11: Plane Wave for Wide Range Linearly Weighted Amplitude
   ! 12: Plane Wave for Wide Range Quadratically Weighted Amplitude
   ! 13: Plane Wave for Wide Range Linearly +1 Weighted Amplitude
   ! 14: Plane Wave for Wide Range Quadratically +1 Weighted Amplitude
   ! 20: Plane Wave for Wide Angle
   ! 100: Dipole Radiation
   ! 200: Point Source

   integer, parameter :: Flag_LSF_Piecewise_Constant= 1 
   ! 0: off, 1: on 

   double precision, parameter :: Ratio_Obj_Func_Thermal= 0.5d0
   double precision, parameter :: Ratio_Obj_Func_DC= 0.5d0

   double precision, parameter :: Sign_ObjectiveFunction= 1d0
   ! 1: minimum, -1: maximum

   !================================================================================================
   ! Update Level Set Function
   !================================================================================================
   double precision, parameter :: LSF_Maximum= 1d0
   double precision, parameter :: LSF_Minimum= -1d0

   !===============================================================================================
   ! Spectral Level Set Methodorogy (SLSM)
   !===============================================================================================
   integer, parameter :: Order_SLSM = 3
   double precision, parameter :: Extention_SLSM = 1.0d-1

   !===============================================================================================
   ! Radial Basis Function (RBF)
   !===============================================================================================
   integer, parameter :: Number_Grid_Scale_DesignDomain_RBF=10

   integer, parameter :: Number_Grid_X_DesignDomain_RBF = Number_Grid_Scale_DesignDomain_RBF 
   integer, parameter :: Number_Grid_Y_DesignDomain_RBF = Number_Grid_Scale_DesignDomain_RBF

   double precision, parameter :: Order_RBF = dble(Number_Grid_Scale_DesignDomain) / dble(Number_Grid_Scale_DesignDomain_RBF)

   integer, parameter :: Number_GridPoint_DesignDomain_RBF &
   = ( Number_Grid_X_DesignDomain_RBF*2 + 1 ) * ( Number_Grid_Y_DesignDomain_RBF*2 + 1 )

   !================================================================================================
   ! Constraints
   !================================================================================================
   double precision, parameter :: Coefficient_Complexity_Tau_Normal= 1d-2  
   double precision, parameter :: Coefficient_Additional_Obj_Func= 1.0d0
   double precision, parameter :: Coefficient_Additional_Obj_Func_2= 0.0d0
   integer, parameter :: Power_Obj_Func= 4

   integer, parameter :: Flag_Perimeter_Constraint= 1 
   ! 0:No, 1:Explicit, 2:Implicit

   integer, parameter :: Flag_Volume_Constraint= 0
   double precision, parameter :: Sign_VolumeConstraint= 1d0
   ! 1d0: V < Vc, -1d0: V > Vc
   double precision, parameter :: Ratio_Volume_Constraint= 0.4d0 
   double precision, parameter :: Coefficient_Decay_Volume_Constraint= 100d0

   integer, parameter :: Level_Constraint_Minimum_Length= 0 
   ! 0:Off, N(>0): Around N Grids
 
   !================================================================================================
   ! Output
   !================================================================================================

   double precision, parameter :: ObjectiveFunction_Threshold_Output= 1d-2 ! set 
   double precision, parameter :: Perimeter_Threshold_Output= 50d0 
   integer, parameter :: Flag_Update_Threshold_Plot= 1 

   integer, parameter :: Optimization_Step_Output= 0 
   !integer, parameter :: Optimization_Step_Output= Number_Optimization_Step +1 

   integer, parameter :: Flag_Output_LSF= 1 !Default: 1

   integer, parameter :: Flag_Output_Distribution_Result= 1 !Default: 2
      integer, parameter :: Number_Level_Plot_EM= 100
      integer, parameter :: Max_Number_Level_Element_Limit= Number_Level_Plot_EM +1
      integer, parameter :: Flag_RGB_Color= 4!1
 
   integer, parameter :: Flag_Output_PostScript= 1 ! Default:10, PML:6, Interval:7
   ! 0: Off, 1: On_Dual_Node, 2:On_No_Dual_Node, 3:On_No_Dual_Node_With_PML
      integer, parameter :: Number_Plot_1= 1 
      integer, parameter :: Number_Plot_2= 1 

   integer, parameter :: Flag_Range_Value_Plot= 1 
   ! 0: flexible, 1: fixed 
      ! Flag_Range_Value_Plot= 1
      double precision, parameter :: Maximum_Value_Plotted_Fixed= 1d0 
      double precision, parameter :: Minimum_Value_Plotted_Fixed= 0d0

      double precision, parameter :: RGB_Lower_R= 0d0
      double precision, parameter :: RGB_Lower_G= 0d0
      double precision, parameter :: RGB_Lower_B= 0.5d0

      double precision, parameter :: RGB_Higher_R= 0.5d0
      double precision, parameter :: RGB_Higher_G= 0d0
      double precision, parameter :: RGB_Higher_B= 0d0

   integer, parameter :: Flag_Output_Finite_Element_Data= 2
   ! Output Flag 0: off, 1: on, 2: OB(i-2) > OB(i-1) < OB(i)

      !Flag_Output_Finite_Element_Data= 1
      integer, parameter :: OptimiationStep_Output_Finite_Element_Data= Optimization_Step_Output 

   integer, parameter :: Flag_Plot_Node_Group= 0
   integer, parameter :: Flag_Plot_Boundary_Group= 0
   ! 0:Off, 1:On

   integer, parameter :: Flag_Find_Dual_Node= 1
   ! 0: Slow, 1:Fast   

   integer, parameter :: Flag_Output_Number_GirdPoint_LSF_P0M= 0
   ! Output Flag 0: off, 1: on
  
   integer, parameter :: Flag_Output_GridPoint= 0 
   ! 0: off, 1: on
  
   integer, parameter :: Length_Character_Optimization_Step= log10(dble(Number_Optimization_Step)) +1

   !===============================================================================================
   ! Physical Parameters 
   !===============================================================================================
   double precision, parameter :: Pi=3.14159265358979d0
   double precision, parameter :: Light_Speed_Vacuum=299792458d0
   double precision, parameter :: Eps0=8.854187817D-12
   double precision, parameter :: Mu0=1d0/(Light_Speed_Vacuum*Light_Speed_Vacuum*Eps0)

   complex( kind(0d0) ), parameter :: Imaginary_Unit= dcmplx( 0d0, 1d0 )
   complex( kind(0d0) ), parameter :: One= dcmplx( 1d0, 0d0 )
   complex( kind(0d0) ), parameter :: Zero= dcmplx( 0d0, 0d0 )

   !================================================================================================
   ! Meshing
   !================================================================================================
   integer, parameter :: Number_Grid_X_DesignDomain= Number_Grid_Scale_DesignDomain 
   integer, parameter :: Number_Grid_Y_DesignDomain= Number_Grid_Scale_DesignDomain
 
   double precision, parameter :: Length_Characteristic= dble( Number_Grid_Scale_DesignDomain )
   double precision, parameter :: Length_Characteristic_Normalized= 1d0 

   integer, parameter :: Width_PML= 10

   double precision, parameter :: Ratio_DesignDomain_2_Scattering= 2d0

   integer, parameter :: Number_Grid_X_Scattering= Number_Grid_X_DesignDomain *Ratio_DesignDomain_2_Scattering
   integer, parameter :: Number_Grid_Y_Scattering= Number_Grid_Y_DesignDomain *4d0/3d0
 
   integer, parameter :: Number_Grid_X_FixedDomain= Number_Grid_X_Scattering !Number_Grid_X_DesignDomain*1.5d0 
   integer, parameter :: Number_Grid_Y_FixedDomain= Number_Grid_Y_DesignDomain/5d0

   integer, parameter :: Number_Grid_Scattering= Number_Grid_X_Scattering*Number_Grid_Y_Scattering*4
   
   integer, parameter :: Number_Grid_X_InitialConfiguration= Number_Grid_X_DesignDomain*0.8d0
   integer, parameter :: Number_Grid_Y_InitialConfiguration= Number_Grid_Y_DesignDomain*0.8d0
   
   integer, parameter :: Number_GridPoint_All &
                 = ( ( Number_Grid_X_Scattering +Width_PML )*2 +1 ) &
                  *( ( Number_Grid_Y_Scattering +Width_PML )*2 +1 )
   
   integer, parameter :: Number_GridPoint_Scattering &
                 = ( Number_Grid_X_Scattering*2 +1  ) &
                  *( Number_Grid_Y_Scattering*2 +1  )
   
   integer, parameter :: Number_Node_tmp &
                 = ( ( Number_Grid_X_Scattering +Width_PML )*2 +1 ) &
                  *( ( Number_Grid_Y_Scattering +Width_PML )*2 +1 ) *8 
    
   integer, parameter :: Number_Element_tmp & 
                 = ( ( Number_Grid_X_Scattering +Width_PML )*2 +1 ) &
                  *( ( Number_Grid_Y_Scattering +Width_PML )*2 +1 ) *8 
   
   integer, parameter :: Number_Grid_PML&
                 = 2*Number_Grid_X_Scattering *Width_PML*2 &
                  +2*Number_Grid_Y_Scattering *Width_PML*2 &
                  +Width_PML*Width_PML*4
   
   integer, parameter :: Number_GridPoint_PML&
                 = ( ( Number_Grid_Y_Scattering+ Width_PML )*2 +1 )* ( Width_PML+1 ) *2 &
                   + ( ( Number_Grid_X_Scattering-1 )*2 +1 )* ( Width_PML+1 ) *2
   
   integer, parameter :: Max_Number_Node_OneGrid= 16  
   integer, parameter :: Max_Number_Element_OneGrid= 16 
   integer, parameter :: Max_Number_CrossPoint= 5
   
   integer, parameter :: Max_Number_Element_OneGrid_PML= 1
  
   double precision, parameter :: Size_X_DesignDomain= dble( Number_Grid_X_DesignDomain ) 
   double precision, parameter :: Size_Y_DesignDomain= dble( Number_Grid_Y_DesignDomain ) 
   double precision, parameter :: Radius_DesignDomain= dble( Number_Grid_Scale_DesignDomain ) 
 
   double precision, parameter :: Size_X_Scattering= dble( Number_Grid_X_Scattering ) 
   double precision, parameter :: Size_Y_Scattering= dble( Number_Grid_Y_Scattering )

   double precision, parameter :: Size_X_FixedDomain= dble( Number_Grid_X_FixedDomain ) 
   double precision, parameter :: Size_Y_FixedDomain= dble( Number_Grid_Y_FixedDomain ) 
   double precision, parameter :: Radius_FixedDomain= Radius_DesignDomain*3d0/5d0 
   !double precision, parameter :: Radius_FixedDomain= 0d0 ! set 

   double precision, parameter :: Position_Surface_Design_2_Fixed_X= Size_X_DesignDomain 
   double precision, parameter :: Position_Surface_Design_2_Fixed_X_PML= Size_X_Scattering *0.8 
   
   !===============================================================================================
   ! General settings
   !===============================================================================================

   integer, parameter :: Flag_DesignDomain= 1
   ! 1: Circle
   ! 2: Square
   
   integer, parameter :: Flag_FixedDomain= 1 !2 aho 
   ! 0: Vacant
   ! 1: Circle
   ! 2: Square

   integer, parameter :: Flag_Exterior_Domain= 0 
   ! 0: Vacant
   ! 1: Periodic Cylinder

   integer, parameter :: Flag_Lattice_Period_Exterior= 1
   ! 1: Triangular Lattice  
   ! 2: Square Lattice 

   double precision, parameter :: Position_Center_FixedDomain_X= 0d0 
   double precision, parameter :: Position_Center_FixedDomain_Y= 0d0! aho -Size_Y_Scattering -Size_Y_FixedDomain
   double precision, parameter :: Position_Center_FixedDomain_X_Cloaked_Region= 0d0
   double precision, parameter :: Position_Center_FixedDomain_Y_Cloaked_Region= 0d0 
 
   double precision, parameter :: Position_Center_DesignDomain_X= 0d0 
   double precision, parameter :: Position_Center_DesignDomain_Y= 0d0 

   double precision, parameter :: Position_Center_InitialConfiguration_X= Position_Center_FixedDomain_X_Cloaked_Region 
   double precision, parameter :: Position_Center_InitialConfiguration_Y= Position_Center_FixedDomain_Y_Cloaked_Region
   
   ! Periodic Cylinder in Exterior Domain 
   double precision, parameter :: Length_Edge_Hexagonal= Size_X_DesignDomain
   double precision, parameter :: Periodic_Length_Exterior_Domain = Length_Edge_Hexagonal/2d0 
   integer, parameter :: Number_Cylinder_Periodic_ExteriorDomain_X= Size_X_Scattering/Periodic_Length_Exterior_Domain 
   integer, parameter :: Number_Cylinder_Periodic_ExteriorDomain_Y= Number_Cylinder_Periodic_ExteriorDomain_X 
   integer, parameter :: Number_Cylinder_Periodic_ExteriorDomain_Maximum &
                = Number_Cylinder_Periodic_ExteriorDomain_X*Number_Cylinder_Periodic_ExteriorDomain_Y 

   double precision, parameter :: Radius_Cylinder_Periodic_ExteriorDomain= Periodic_Length_Exterior_Domain/3d0

   double precision, parameter :: Upper_Position_Waveguide_Y= Periodic_Length_Exterior_Domain *( -3d0 ) *dsin( 2d0*Pi*60d0/360d0 ) 
   double precision, parameter :: Lower_Position_Waveguide_Y= Periodic_Length_Exterior_Domain *( -4d0 ) *dsin( 2d0*Pi*60d0/360d0 )
 
   !===============================================================================================
   ! Grouping Nodes to Find Dual Nodes
   !===============================================================================================
   integer, parameter :: Number_Node_OneGrid_Assumed= 3 
   integer, parameter :: Number_Group_Node_Position_Axis & 
                 = SQRT( dble( Number_Node_OneGrid_Assumed*( Number_Grid_X_Scattering +Width_PML ) )/2d0 )
   
   integer, parameter :: Number_Group_Node_All= ( 2*Number_Group_Node_Position_Axis )**2
   
   integer, parameter :: Max_Number_Node_OneGroup= 3d0*( Number_Node_tmp +Number_GridPoint_PML ) /Number_Group_Node_All 
   
   !===============================================================================================
   ! Symmetry of Sensitivity 
   !===============================================================================================

   integer, parameter :: Number_Level_Distance_Symmetry_Sensitivity_Rotational= Number_Grid_X_Scattering*sqrt( 2d0 ) +1 
   integer, parameter :: Maximum_Number_GridPoint_OneLevel= Number_GridPoint_Scattering/Number_Level_Distance_Symmetry_Sensitivity_Rotational*10

   !===============================================================================================
   ! Configuration
   !===============================================================================================
   integer, parameter :: Flag_InitialConfiguration= 200 !set 
   ! 0: Vacant
   ! 1: Vacant Initial Configuration with Pediodic distribution of LSF
   ! 5: All Dielectrics
   ! 6: All Dielectrics Initial Configuration with Pediodic distribution of LSF
   ! 10: Circle
   ! 11: Three Circles
   ! 20: Square
   ! 30: Triangle
   ! 40: Random
   ! 50: Cylinders Random
   ! 60: Ellipse 
   ! 70: Double Circle
   ! 80: Double Elipse 
   ! 81: Double Elipse with y  
   ! 82: Outer Elipse & Inner Circle 
   ! 90: Rectangle Elipse  
   ! 100: Rectangle Elipse with point 
   ! 110: Cylinders Periodic
   ! 111: Cylinders Periodic in Circle
   ! 112: Cylinders Periodic in Ellipse
   ! 113: Cylinders Periodic in Whole Circular Design Domain
   ! 114: Cylinders Periodic in Whole Rectangular Design Domain
   ! 115: Cylinders Radiative 
   ! 116: Cylinders Pendulum 
   ! 117: Cylinders Gridpoint
   ! 120: Porous Periodic in Triangular Lattice
   ! 121: Porous Radiative 
   ! 122: Porous Radiative Waveguide 
   ! 123: DUal Periodic Porous in Triangular Lattice
   ! 124: DUal Periodic Porous in Triangular Lattice
   ! 140: Uneven
   ! 150: Multilayer
   ! 151: Angular Multilayer
   ! 152: Circular Multilayer
   ! 153: Radiative Multilayer
   ! 154: Splitted Multilayer
   ! 155: Splitted Multilayer shifted Half Perood
   ! 156: Circular Multilayer Splitted
   ! 157: Circular Multilayer Ellipse
   ! 160: Square fractal
   ! 170: Fishnet
   ! 200: CMAES

   ! 0 -- 10
   double precision, parameter :: Size_X_InitialConfiguration= dble( Number_Grid_X_InitialConfiguration ) 
   double precision, parameter :: Size_Y_InitialConfiguration= dble( Number_Grid_Y_InitialConfiguration ) 
   double precision, parameter :: Radius_InitialConfiguration= Radius_DesignDomain*0.95
 
   !===============================================================================================
   ! PML Material
   !===============================================================================================
   integer, parameter :: Number_Domain_Material_PML= 2

   double precision, parameter :: Boundary_PML_Material_X1_Low= -Size_X_Scattering -dble( Width_PML ) -1d0
   double precision, parameter :: Boundary_PML_Material_X1_High= Size_X_Scattering +dble( Width_PML ) +1d0
   double precision, parameter :: Boundary_PML_Material_Y1_Low= ( Position_Center_FixedDomain_Y -Size_Y_FixedDomain )/Length_Characteristic
   double precision, parameter :: Boundary_PML_Material_Y1_High= ( Position_Center_FixedDomain_Y +Size_Y_FixedDomain )/Length_Characteristic

   double precision, parameter :: Boundary_PML_Material_X2_Low= 1d10 
   double precision, parameter :: Boundary_PML_Material_X2_High= 1d10
   double precision, parameter :: Boundary_PML_Material_Y2_Low= 1d10
   double precision, parameter :: Boundary_PML_Material_Y2_High= 1d10
 
   integer, parameter :: ID_Element_PML_Base_Material= 0
   integer, parameter :: ID_Element_PML_Material= 1

   !===============================================================================================
   ! Postscript
   !===============================================================================================

   ! Check on Desktop
   integer, parameter :: WidthPS= 360
   integer, parameter :: WidthPS_X= WidthPS*Number_Grid_X_DesignDomain/Number_Grid_Scale_DesignDomain 
   integer, parameter :: WidthPS_Y= WidthPS_X*Number_Grid_Y_Scattering/Number_Grid_X_Scattering
   !integer, parameter :: WidthPS_Y= WidthPS*Number_Grid_Y_DesignDomain/Number_Grid_Scale_DesignDomain
   integer, parameter :: Translation_X= 10
   integer, parameter :: Translation_Y= 450
   integer, parameter :: Width_Color_Var= WidthPS/20

   integer, parameter :: Flag_Color_Var= 2 
   ! 0:off, 1:Right Side, 2:Upper Side

   !===============================================================================================
   ! Paramters for Plot 
   !===============================================================================================
   double precision, parameter :: &
   Position_Center_FixedDomain_X_PS &
   = Position_Center_FixedDomain_X/( 2d0*Size_X_Scattering ) *dble( WidthPS_X ) &
    +dble( WidthPS_X )/2d0 +Translation_X

   double precision, parameter :: &
   Position_Center_FixedDomain_Y_PS &
   = Position_Center_FixedDomain_Y/( 2d0*Size_Y_Scattering ) *dble( WidthPS_Y ) &
    +dble( WidthPS_Y )/2d0 +Translation_Y

   double precision, parameter :: &
   Position_Center_DesignDomain_X_PS &
   = Position_Center_DesignDomain_X/( 2d0*Size_X_Scattering ) *dble( WidthPS_X ) &
    +dble( WidthPS_X )/2d0 +Translation_X

   double precision, parameter :: &
   Position_Center_DesignDomain_Y_PS &
   = Position_Center_DesignDomain_Y/( 2d0*Size_Y_Scattering ) *dble( WidthPS_Y ) &
    +dble( WidthPS_Y )/2d0 +Translation_Y

   double precision, parameter :: &
   Radius_FixedDomain_Plot_PS &
   = Radius_FixedDomain/( 2d0*Size_X_Scattering ) *dble( WidthPS_X )

   double precision, parameter :: &
   Radius_DesignDomain_Plot_PS &
   = Radius_DesignDomain/( 2d0*Size_X_Scattering ) *dble( WidthPS_X )

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

   !===============================================================================================
   ! Create Input Data 
   !===============================================================================================
   
   integer, parameter :: Flag_Check_Dual_Node= 0 
   !0 : off, 1 : on

   integer, parameter :: Flag_Output_Number_Node_Boundary_PML= 0 !0
   integer, parameter :: Flag_Output_Node_Element_Renumbered= 0  
   ! Output Flag 0: off, 1: on

   !===============================================================================================
   ! Create Input Data 
   !===============================================================================================

   double precision, parameter :: Temperature_Boundary_Lower= 0d0 
      double precision, parameter :: Temperature_Boundary_Lower_Lower= 0d0
      double precision, parameter :: Temperature_Boundary_Lower_Higher= 5d0
      double precision, parameter :: Temperature_Boundary_Lower_Center &
                                     = ( Temperature_Boundary_Lower_Higher +Temperature_Boundary_Lower_Lower )/2d0
      double precision, parameter :: Width_Temperature_Boundary_Lower &
                                     = Temperature_Boundary_Lower_Higher -Temperature_Boundary_Lower_Lower 

   double precision, parameter :: Temperature_Boundary_Higher= 1d0 
      double precision, parameter :: Temperature_Boundary_Higher_Lower= 0d0 
      double precision, parameter :: Temperature_Boundary_Higher_Higher= 5d0 

   double precision, parameter :: Electric_Potential_Boundary_Lower= 0d0 
      double precision, parameter :: Electric_Potential_Boundary_Lower_Lower= 0d0
      double precision, parameter :: Electric_Potential_Boundary_Lower_Higher= 5d0
      double precision, parameter :: Electric_Potential_Boundary_Lower_Center &
                                     = ( Electric_Potential_Boundary_Lower_Higher +Electric_Potential_Boundary_Lower_Lower )/2d0
      double precision, parameter :: Width_Electric_Potential_Boundary_Lower &
                                     = Electric_Potential_Boundary_Lower_Higher -Electric_Potential_Boundary_Lower_Lower 

   double precision, parameter :: Electric_Potential_Boundary_Higher= 1d0 
      double precision, parameter :: Electric_Potential_Boundary_Higher_Lower= 0d0 
      double precision, parameter :: Electric_Potential_Boundary_Higher_Higher= 5d0 

   double precision, parameter :: Position_Source_x= -Size_X_Scattering/Length_Characteristic*0.9d0
      double precision, parameter :: Position_Source_x_1= -Size_X_Scattering/Length_Characteristic*0.9d0
      double precision, parameter :: Position_Source_x_2= 0.0d0

   double precision, parameter :: Position_Source_y= Size_Y_Scattering/Length_Characteristic -Size_X_Scattering*0.1d0/Length_Characteristic
      double precision, parameter :: Position_Source_y_1= Position_Source_y
      double precision, parameter :: Position_Source_y_2= -Position_Source_y

   !===============================================================================================
   ! Analyse_Light_Scattering
   !===============================================================================================
  
   ! http://www.hakko.co.jp/qa/qakit/html/h01020.htm 
   double precision, parameter :: Thermal_Conductivity_Material= 386d0/67d0 ! Pure Copper [W/m K]
   double precision, parameter :: Thermal_Conductivity_Active= 1d0
   double precision, parameter :: Thermal_Conductivity_Base_Material= 67d0/67d0 ! Pure Iron [W/m K]
   double precision, parameter :: Thermal_Conductivity_OuterDomain= Thermal_Conductivity_Base_Material 
   double precision, parameter :: Thermal_Conductivity_FixedDomain= Thermal_Conductivity_Material

   double precision, parameter :: Thermal_Conductivity_Material_Min= 4.0d0
   double precision, parameter :: Thermal_Conductivity_Material_Max= 8.0d0


   double precision, parameter :: Electric_Conductivity_Material= 59.0d0/9.9d0 ! Pure Copper [W/m K] 
   double precision, parameter :: Electric_Conductivity_Active= 1d0
   double precision, parameter :: Electric_Conductivity_Base_Material= 1d0 !9.9d0/9.9d0 ! Pure Iron [W/m K]
   double precision, parameter :: Electric_Conductivity_OuterDomain= Electric_Conductivity_Base_Material 
   double precision, parameter :: Electric_Conductivity_FixedDomain= Electric_Conductivity_Material

   double precision, parameter :: Electric_Conductivity_Material_Min= 4.0d0
   double precision, parameter :: Electric_Conductivity_Material_Max= 8.0d0


   !===============================================================================================
   ! ID
   !===============================================================================================
   
   integer, parameter :: ID_GridPoint_Initialization= -100000000
   integer, parameter :: ID_GridPoint_DesignDomain= 1 
   integer, parameter :: ID_GridPoint_ObjectiveFunction= 0 
   integer, parameter :: ID_GridPoint_PML= -1
   integer, parameter :: ID_GridPoint_FixedDomain= 2
   
   integer, parameter :: ID_GridPoint_Boundary_DesignDomain= 200 
   integer, parameter :: ID_GridPoint_Boundary_FixedDomain= 201 
   integer, parameter :: ID_GridPoint_Boundary_PML= 202
   
   integer, parameter :: ID_Element_Material= 1
   integer, parameter :: ID_Element_FixedDomain= 2
   integer, parameter :: ID_Element_Active= 3
   integer, parameter :: ID_Element_Base_Material= 4
   integer, parameter :: ID_Element_OuterDomain= 5
   
   integer, parameter :: ID_Element_Material_Exterior= 6
   integer, parameter :: ID_Element_Base_Material_Exterior= ID_Element_OuterDomain
   
   integer, parameter :: ID_Element_PML_Xp= 10
   integer, parameter :: ID_Element_PML_Xm= 11
   integer, parameter :: ID_Element_PML_Yp= 12
   integer, parameter :: ID_Element_PML_Ym= 13
   integer, parameter :: ID_Element_PML_XpYp= 14
   integer, parameter :: ID_Element_PML_XpYm= 15
   integer, parameter :: ID_Element_PML_XmYp= 16
   integer, parameter :: ID_Element_PML_XmYm= 17
  
   integer, parameter :: ID_Element_Obj_Func_1= ID_Element_OuterDomain 
   integer, parameter :: ID_Element_Obj_Func_2= ID_Element_FixedDomain 
 
   integer, parameter :: ID_LSF_Plus= 1
   integer, parameter :: ID_LSF_Zero= 0
   integer, parameter :: ID_LSF_Minus= -1

   integer, parameter :: ID_Node_Material= 1
   integer, parameter :: ID_Node_FixedDomain= 2
   integer, parameter :: ID_Node_Active= 3
   integer, parameter :: ID_Node_Base_Material= 4
   integer, parameter :: ID_Node_OuterDomain= 5
  
   integer, parameter :: ID_Node_Material_Exterior= 6
   integer, parameter :: ID_Node_Base_Material_Exterior= ID_Node_OuterDomain

   integer, parameter :: ID_Node_PML_X= 10
   integer, parameter :: ID_Node_PML_Y= 11
   integer, parameter :: ID_Node_PML_XY= 12
     
   !===============================================================================================
   ! MPI  
   !===============================================================================================
   integer, parameter :: Flag_Output_Data_MPI= 0 
   integer, parameter :: Flag_Check_Value= 0 

   !===============================================================================================
   ! CMA-ES Data 
   !===============================================================================================
   integer, parameter :: Flag_Output_xmean= 0 

   !===============================================================================================
   ! Initialization  
   !===============================================================================================

   integer, parameter :: Flag_On= 1
   integer, parameter :: Flag_Off= 0
   
   integer, parameter :: Integer_Initialization= -10D8
   double precision, parameter :: Double_Precision_Initialization= -10D10

   integer, parameter :: Integer_Initialization_Plus= 10D8
   integer, parameter :: Integer_Initialization_Minus= -10D8
   double precision, parameter :: Double_Precision_Initialization_Plus= 10D10
   double precision, parameter :: Double_Precision_Initialization_Minus= -10D10
   
   !===============================================================================================
   ! Tolerance
   !===============================================================================================

   double precision, parameter :: Tolerance_Position_Node= 1.0d-13 !10d-5 
   double precision, parameter :: Tolerance_Boundary_Group= 2d0*Tolerance_Position_Node 
   !double precision, parameter :: Tolerance_LSF= 1.0d-8
   double precision, parameter :: Tolerance_LSF= 1d-5
   double precision, parameter :: Tolerance_Cross_Check= 1.0d-11
   
   double precision, parameter :: Tolerance_Digit_Error= 1.0d-8
   double precision, parameter :: Tolerance_Digit_Error_Plus= 1d0 +Tolerance_Digit_Error
   double precision, parameter :: Tolerance_Digit_Error_Minus= 1d0 -Tolerance_Digit_Error

   double precision, parameter :: Tolerance_Plot= 1.0d-2

end module Parameters 



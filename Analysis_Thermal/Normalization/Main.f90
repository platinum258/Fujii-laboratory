
program Main

   !$ use omp_lib
   use Parameters
   implicit none

   character(len=100) :: File_Name_FEM_Data_Optimal
   character(len=100) :: File_Name_FEM_Data_Normalization
   character(len=100) :: File_Name_FEM_Data_Reference
   character(len=150) :: FilePass_Optimal
   character(len=150) :: FilePass_Normalization, FilePass_Reference 
   character(len=150) :: File_Name_Result_Normalization 
   character(len=150) :: File_Pass_Result_Normalization, File_Pass_Reference_Field 

   integer :: Number_Loop_Frequency, Number_Loop_IncidentAngle, Number_Loop_Material_Property
   integer :: Number_Loop_Position_Source_X, Number_Loop_Position_Source_Y
   integer :: Flag_Theoretical_Solution
   integer :: Flag_Structure, Flag_Output_PhysicalField, Flag_Plot_Configuration, Flag_Plot_Mesh_Configuration
   integer :: Flag_Plot_EFD_Configuration, Flag_Plot_EFD_OC_PoyntingVector 
   integer :: Flag_Reference_Device

   integer :: Number_Node_tmp
   integer :: Flag_Reference_Data

   integer :: Flag_Thermal_Insulation_BC  
   integer :: e, i, j, k, No, Counter
   integer :: NumThread, Loop_Core
 
   integer, allocatable, dimension(:) :: itg1
   integer, allocatable, dimension(:,:) :: itg2
   double precision, allocatable, dimension(:,:) :: dp2 

   integer :: LoopCounter_Frequency, LoopCounter_IncidentAngle, LoopCounter_Material_Property
   integer :: LoopCounter_Position_Source_X, LoopCounter_Position_Source_Y
   integer :: ID_Element_ObjectiveFunction
   integer :: ID_Element_Material, ID_Element_OpenRegion, ID_Element_FixedDomain, ID_Element_Base_Material 
   integer :: ID_Element_Material_Exterior, ID_Element_Base_Material_Exterior

   double precision, allocatable, dimension(:) :: Temperature_Left_Side_Boundary, Thermal_Conductivity_Material_Real
   double precision, allocatable, dimension(:) :: Temperature_Right_Side_Boundary 
   double precision, allocatable, dimension(:) :: Thermal_Conductivity_Material, Thermal_Conductivity_FixedDomain_Variable

   double precision, allocatable, dimension(:) :: Position_Source_X, Position_Source_Y

   integer :: Number_Node, Number_Element_Triangle
   integer :: Number_Node_on_PEC_Boundary, Number_Node_in_PEC 
   integer :: Max_Number_Element_Share_OneNode, Width_Matrix_LHS 
   double precision, allocatable, dimension(:,:) :: Position_Node 
   double precision, allocatable, dimension(:,:) :: Max_Position_Node
   integer, allocatable, dimension(:,:) :: Index_Element_2_Node_Triangle, Index_Element_2_Node_Plot
   integer, allocatable, dimension(:) :: Class_Element_Triangle
   integer, allocatable, dimension(:) :: Class_Element_Plot

   integer, allocatable, dimension(:,:) :: Element_and_LocalNode_Number_on_PEC_BC
   integer :: Number_Element_PEC_Boundary

   double precision, allocatable, dimension(:) :: Temperature_Solution, Temperature_Reference 
   double precision, allocatable, dimension(:) :: Thermal_Conductivity_Element, Thermal_Conductivity_PV 

   integer :: Number_Node_Reference, Number_Element_Plot, Counter_Element_Plot 
   !complex( kind( 0d0 ) ), allocatable, dimension(:,:) :: Value_Objective_Function
   double precision, allocatable, dimension(:,:) :: Value_Objective_Function
   double precision, allocatable, dimension(:) :: Value_Plot_Distribution

   double precision, allocatable, dimension(:,:,:,:,:,:) :: Normalized_Objective_Function, Original_OF_Value_without_Normalization, Value_OF_Normalize

   integer, allocatable, dimension(:,:) :: LoopCounter_Frequency_Range, LoopCounter_IncidentAngle_Range, LoopCounter_Material_Property_Range
   integer, allocatable, dimension(:,:) :: LoopCounter_Position_Source_X_Range, LoopCounter_Position_Source_Y_Range 

   ! Plot PV
   integer :: Number_Point_PV_X, Number_Point_PV_Y, Number_PV, Number_PV_tmp 
   double precision, allocatable, dimension(:,:) :: Position_PV, Position_PV_tmp 
   integer, allocatable, dimension(:) :: Element_Number_PV, Element_Number_PV_tmp 
   double precision :: Length_Maximum_Vector, Length_Minimum_Vector 
   double precision, allocatable, dimension(:,:,:) :: Position_Arrow_PV, Position_Arrow_tmp
   double precision, allocatable, dimension(:) :: Level_PV, Level_PV_tmp
        


   !***********************************************************************************
   !addition
   !***********************************************************************************

   double precision :: high,thickness
   integer :: penalty_density,g_min
   
   double precision, allocatable, dimension(:) :: Thermal_conductivity_by_density

   !================================================================================
   write(*,*)'====================================================================='
   write(*,*)' Heat Transfer Analysis Program for 2D Cloaking Structure'
   write(*,*)' Finite Element Method ( Node Base Element )'
   write(*,*)' Thermal Insulation Boundary Condition'
   write(*,*)' Multifrontal Method ( MA57 )'
   write(*,*)'  by Garuda Fujii'
   write(*,*)'   Ver.1: 16 Feb. 2017'
   write(*,*)'====================================================================='
   !================================================================================

   open( 100, file='./Input_Data' )
      read(100,*) File_Name_FEM_Data_Optimal
      read(100,*) File_Name_FEM_Data_Normalization
      read(100,*) File_Name_FEM_Data_Reference
      read(100,*) FilePass_Optimal 
      read(100,*) FilePass_Normalization 
      read(100,*) FilePass_Reference 
      read(100,*) Number_Loop_Frequency
      read(100,*) Number_Loop_IncidentAngle
      read(100,*) Number_Loop_Material_Property
      read(100,*) Number_Loop_Position_Source_X
      read(100,*) Number_Loop_Position_Source_Y
      read(100,*) Flag_Structure 
      read(100,*) Flag_Theoretical_Solution
      read(100,*) Flag_Plot_Configuration 
      read(100,*) Flag_Plot_Mesh_Configuration 
      read(100,*) Flag_Output_PhysicalField 
      read(100,*) Flag_Plot_EFD_Configuration 
      read(100,*) Flag_Plot_EFD_OC_PoyntingVector 
      read(100,*) File_Name_Result_Normalization 
      read(100,*) File_Pass_Result_Normalization 
      read(100,*) File_Pass_Reference_Field
      read(100,*) Flag_Reference_Device
   close( 100 )

   if( Flag_Structure==2 )then
      write(*,*)'open : ', trim(FilePass_Reference)//trim(File_Name_FEM_Data_Reference) 
      open( 10, file=trim(FilePass_Reference)//trim(File_Name_FEM_Data_Reference) ) 
   else if( Flag_Structure==1 )then
      write(*,*)'open : ', trim(FilePass_Normalization)//trim(File_Name_FEM_Data_Normalization) 
      open( 10, file=trim(FilePass_Normalization)//trim(File_Name_FEM_Data_Normalization) ) 
   else if( Flag_Structure==0  )then
      write(*,*)'open : ', trim(FilePass_Optimal)//trim(File_Name_FEM_Data_Optimal) 
      open( 10, file=trim(FilePass_Optimal)//trim(File_Name_FEM_Data_Optimal) ) 
   end if
      read(10,*) Flag_Thermal_Insulation_BC
      read(10,*) Number_Node
      read(10,*) Number_Element_Triangle
      read(10,*) Max_Number_Element_Share_OneNode
      read(10,*) Number_Node_on_PEC_Boundary  
      read(10,*) Number_Node_in_PEC 
      read(10,*) Number_Element_PEC_Boundary
      read(10,*) ID_Element_ObjectiveFunction
      read(10,*) ID_Element_Material, ID_Element_OpenRegion, ID_Element_FixedDomain, ID_Element_Base_Material
      read(10,*) ID_Element_Material_Exterior, ID_Element_Base_Material_Exterior
      read(10,*) high
      read(10,*) thickness
      read(10,*) penalty_density
      read(10,*) g_min


   allocate( Max_Position_Node( 2, 2 ) )

      read(10,*) Max_Position_Node( 1, 1 ), Max_Position_Node( 1, 2 ), Max_Position_Node( 2, 1 ), Max_Position_Node( 2, 2 )
 
   allocate( Position_Node( 2, Number_Node ) )

      do i= 1, Number_Node
         read(10,*) No, Position_Node( 1, i ), Position_Node( 2, i )
      end do

   allocate( Index_Element_2_Node_Triangle( 3, Number_Element_Triangle ) )
   allocate( Class_Element_Triangle( Number_Element_Triangle ) )
   allocate( Thermal_conductivity_by_density( Number_Element_Triangle ) )

      do i= 1, Number_Element_Triangle
         read(10,*) No, Class_Element_Triangle( i ), Thermal_conductivity_by_density( i ),&
                Index_Element_2_Node_Triangle( 1, i ), Index_Element_2_Node_Triangle( 2, i ), Index_Element_2_Node_Triangle( 3, i ) 
      end do

   allocate( Element_and_LocalNode_Number_on_PEC_BC( 3, Number_Element_PEC_Boundary ) )

      do i= 1, Number_Element_PEC_Boundary
         read(10,*) No, Element_and_LocalNode_Number_on_PEC_BC( 1, i ), &
                        Element_and_LocalNode_Number_on_PEC_BC( 2, i ), &
                        Element_and_LocalNode_Number_on_PEC_BC( 3, i )
      end do

   close( 10 )

                                                  open( 23,file='aho.dat')
                                                 do i =1,Number_Element_Triangle
                                                        write( 23,*)Class_Element_Triangle( i )
                                                 enddo
                                                  close(23)


   Width_Matrix_LHS= Max_Number_Element_Share_OneNode +1

   !====================================================================
   write(*,*) 'Modify Input Data for Reference Field'
   !====================================================================

   if( Flag_Refine_Element==1 )then
      allocate( itg2( 3, Number_Element_Triangle*3 ) )
      allocate( itg1( Number_Element_Triangle*3 ) )
      itg2= 99999999
      itg1= 88888888
   
      do e= 1, Number_Element_Triangle
         do i= 1, 3 
            itg2( i, e )= Index_Element_2_Node_Triangle( i, e )
         end do
         itg1( e )= Class_Element_Triangle( e )
      end do
   
      allocate( dp2( 2, Number_Node+Number_Element_Triangle ) )
   
      do i= 1, Number_Node
         do j= 1, 2
            dp2( j, i )= Position_Node( j, i )
         end do
      end do
   
      write(*,*)'Number_Node, Number_Element_Triangle=', Number_Node, Number_Element_Triangle
   
      call Refine_Finite_Elements&
           ( 999, & 
             !ID_Element_FixedDomain, &
             Number_Node, dp2,  &
             Number_Element_Triangle, itg2, itg1, &
             Width_Matrix_LHS )
   
      write(*,*)'Number_Node, Number_Element_Triangle=', Number_Node, Number_Element_Triangle
   
      deallocate( Index_Element_2_Node_Triangle )
      deallocate( Position_Node )
      deallocate( Class_Element_Triangle )
      deallocate( Thermal_conductivity_by_density )  
 
      allocate( Index_Element_2_Node_Triangle( 3, Number_Element_Triangle ) )
      allocate( Class_Element_Triangle( Number_Element_Triangle ) )
      allocate( Position_Node( 2, Number_Node ) )
   
      Position_Node( :, 1:Number_Node )= dp2( :, 1:Number_Node )
   
      do e= 1, Number_Element_Triangle
         Class_Element_Triangle( e )= itg1( e )
         do i= 1, 3
            Index_Element_2_Node_Triangle( i, e )= itg2( i, e )
         end do
      end do
   
      deallocate( itg2 )
      deallocate( itg1 )
      deallocate( dp2 )

   end if

   !====================================================================
   write(*,*) 'Modify Input Data for Reference Field'
   !====================================================================

   if( FilePass_Normalization/=FilePass_Optimal )then
      if( Flag_Structure==2 .or. Flag_Structure==0 )then
         Flag_Reference_Data=1
      else
         Flag_Reference_Data=0
      end if
   else
      if( File_Name_FEM_Data_Optimal==File_Name_FEM_Data_Reference .and. &
          FilePass_Optimal==FilePass_Reference )then
   
         if( Flag_Structure==2 )then 
            if( Flag_Thermal_Device==20 .or. Flag_Thermal_Device==21 .or. Flag_Thermal_Device==30 )then 
               do i= 1, Number_Element_Triangle
                  if( Class_Element_Triangle( i )/=ID_Element_FixedDomain )then
                     Class_Element_Triangle( i )=ID_Element_OpenRegion
                  end if
               end do
            else
               do i= 1, Number_Element_Triangle
                  if( Class_Element_Triangle( i )==ID_Element_Material .or. &
                      Class_Element_Triangle( i )==ID_Element_FixedDomain )then
   
                     Class_Element_Triangle( i )=ID_Element_Base_Material
                  end if
               end do
            end if
         end if
         Flag_Reference_Data=1
   
      else
         Flag_Reference_Data=0
      end if
   
      if( File_Name_FEM_Data_Optimal==File_Name_FEM_Data_Normalization .and. &
          FilePass_Optimal==FilePass_Normalization )then
   
         if( Flag_Structure==1 )then 
            do i= 1, Number_Element_Triangle
               if( Flag_Thermal_Device==20 .or. Flag_Thermal_Device==21 .or. Flag_Thermal_Device==30 )then 
                  if( Class_Element_Triangle( i )==ID_Element_Base_Material )then
                     Class_Element_Triangle( i )=ID_Element_Material
                  end if
               else
                  if( Class_Element_Triangle( i )==ID_Element_Material )then
                     Class_Element_Triangle( i )=ID_Element_Base_Material
                  end if
               end if
            end do
         end if
      end if
   end if

   !====================================================================
   write(*,*) 'Define Number of Cores for Parallel Computation '
   !====================================================================
   open(9999, file='Number_Thread')
   read(9999,*)NumThread
   !$ write(*,*)'      NumThread=', NumThread 
   !$ call omp_set_num_threads( NumThread )
   close(9999)

   allocate( LoopCounter_Frequency_Range( 2, NumThread ) )
   allocate( LoopCounter_IncidentAngle_Range( 2, NumThread ) )
   allocate( LoopCounter_Material_Property_Range( 2, NumThread ) )
   allocate( LoopCounter_Position_Source_X_Range( 2, NumThread ) )
   allocate( LoopCounter_Position_Source_Y_Range( 2, NumThread ) )

   if( Number_Loop_Frequency==0 )then
   
      do Loop_Core= 1, NumThread
         LoopCounter_Frequency_Range( 1, Loop_Core )= 0 
         LoopCounter_Frequency_Range( 2, Loop_Core )= 0 
   
         LoopCounter_IncidentAngle_Range( 1, Loop_Core )= ( Loop_Core -1 ) *( Number_Loop_IncidentAngle +1 )/NumThread 
         LoopCounter_IncidentAngle_Range( 2, Loop_Core )= Loop_Core*( Number_Loop_IncidentAngle +1 )/NumThread -1
   
         LoopCounter_Material_Property_Range( 1, Loop_Core )= 0 
         LoopCounter_Material_Property_Range( 2, Loop_Core )= Number_Loop_Material_Property 

         LoopCounter_Position_Source_X_Range( 1, Loop_Core )= 0
         LoopCounter_Position_Source_X_Range( 2, Loop_Core )= Number_Loop_Position_Source_X 

         LoopCounter_Position_Source_Y_Range( 1, Loop_Core )= 0
         LoopCounter_Position_Source_Y_Range( 2, Loop_Core )= Number_Loop_Position_Source_Y 

      end do

   else if( Number_Loop_IncidentAngle==0 )then
   
      do Loop_Core= 1, NumThread
         LoopCounter_Frequency_Range( 1, Loop_Core )= ( Loop_Core -1 )*( Number_Loop_Frequency +1 )/NumThread
         LoopCounter_Frequency_Range( 2, Loop_Core )= Loop_Core*( Number_Loop_Frequency +1 )/NumThread -1
   
         LoopCounter_IncidentAngle_Range( 1, Loop_Core )= 0 
         LoopCounter_IncidentAngle_Range( 2, Loop_Core )= 0 
   
         LoopCounter_Material_Property_Range( 1, Loop_Core )= 0 
         LoopCounter_Material_Property_Range( 2, Loop_Core )= Number_Loop_Material_Property 

         LoopCounter_Position_Source_X_Range( 1, Loop_Core )= 0
         LoopCounter_Position_Source_X_Range( 2, Loop_Core )= Number_Loop_Position_Source_X 

         LoopCounter_Position_Source_Y_Range( 1, Loop_Core )= 0
         LoopCounter_Position_Source_Y_Range( 2, Loop_Core )= Number_Loop_Position_Source_Y 

      end do

   else if( Number_Loop_Material_Property==0 )then 
   
      do Loop_Core= 1, NumThread
         LoopCounter_Frequency_Range( 1, Loop_Core )= ( Loop_Core -1 )*( Number_Loop_Frequency +1 )/NumThread
         LoopCounter_Frequency_Range( 2, Loop_Core )= Loop_Core*( Number_Loop_Frequency +1 )/NumThread -1
   
         LoopCounter_IncidentAngle_Range( 1, Loop_Core )= 0 
         LoopCounter_IncidentAngle_Range( 2, Loop_Core )= Number_Loop_IncidentAngle 
   
         LoopCounter_Material_Property_Range( 1, Loop_Core )= 0 
         LoopCounter_Material_Property_Range( 2, Loop_Core )= 0  

         LoopCounter_Position_Source_X_Range( 1, Loop_Core )= 0
         LoopCounter_Position_Source_X_Range( 2, Loop_Core )= Number_Loop_Position_Source_X 

         LoopCounter_Position_Source_Y_Range( 1, Loop_Core )= 0
         LoopCounter_Position_Source_Y_Range( 2, Loop_Core )= Number_Loop_Position_Source_Y 

      end do

   else if( Number_Loop_Position_Source_X==0 )then
   
      do Loop_Core= 1, NumThread
         LoopCounter_Frequency_Range( 1, Loop_Core )= ( Loop_Core -1 )*( Number_Loop_Frequency +1 )/NumThread
         LoopCounter_Frequency_Range( 2, Loop_Core )= Loop_Core*( Number_Loop_Frequency +1 )/NumThread -1
   
         LoopCounter_IncidentAngle_Range( 1, Loop_Core )= 0 
         LoopCounter_IncidentAngle_Range( 2, Loop_Core )= Number_Loop_IncidentAngle 
   
         LoopCounter_Material_Property_Range( 1, Loop_Core )= 0 
         LoopCounter_Material_Property_Range( 2, Loop_Core )= Number_Loop_Material_Property 

         LoopCounter_Position_Source_X_Range( 1, Loop_Core )= 0
         LoopCounter_Position_Source_X_Range( 2, Loop_Core )= 0 

         LoopCounter_Position_Source_Y_Range( 1, Loop_Core )= 0
         LoopCounter_Position_Source_Y_Range( 2, Loop_Core )= Number_Loop_Position_Source_Y 

      end do

   else if( Number_Loop_Position_Source_Y==0 )then
   
      do Loop_Core= 1, NumThread
         LoopCounter_Frequency_Range( 1, Loop_Core )= ( Loop_Core -1 )*( Number_Loop_Frequency +1 )/NumThread
         LoopCounter_Frequency_Range( 2, Loop_Core )= Loop_Core*( Number_Loop_Frequency +1 )/NumThread -1
   
         LoopCounter_IncidentAngle_Range( 1, Loop_Core )= 0 
         LoopCounter_IncidentAngle_Range( 2, Loop_Core )= Number_Loop_IncidentAngle 
   
         LoopCounter_Material_Property_Range( 1, Loop_Core )= 0 
         LoopCounter_Material_Property_Range( 2, Loop_Core )= Number_Loop_Material_Property 

         LoopCounter_Position_Source_X_Range( 1, Loop_Core )= 0
         LoopCounter_Position_Source_X_Range( 2, Loop_Core )= Number_Loop_Position_Source_X 

         LoopCounter_Position_Source_Y_Range( 1, Loop_Core )= 0
         LoopCounter_Position_Source_Y_Range( 2, Loop_Core )= 0 

      end do

   else
   
      do Loop_Core= 1, NumThread
         LoopCounter_Frequency_Range( 1, Loop_Core )= ( Loop_Core -1 )*( Number_Loop_Frequency +1 )/NumThread
         LoopCounter_Frequency_Range( 2, Loop_Core )= Loop_Core*( Number_Loop_Frequency +1 )/NumThread -1
   
         LoopCounter_IncidentAngle_Range( 1, Loop_Core )= 0 
         LoopCounter_IncidentAngle_Range( 2, Loop_Core )= Number_Loop_IncidentAngle 
   
         LoopCounter_Material_Property_Range( 1, Loop_Core )= 0 
         LoopCounter_Material_Property_Range( 2, Loop_Core )= Number_Loop_Material_Property 

         LoopCounter_Position_Source_X_Range( 1, Loop_Core )= 0
         LoopCounter_Position_Source_X_Range( 2, Loop_Core )= Number_Loop_Position_Source_X 

         LoopCounter_Position_Source_Y_Range( 1, Loop_Core )= 0
         LoopCounter_Position_Source_Y_Range( 2, Loop_Core )= Number_Loop_Position_Source_Y 
      end do

   end if

   allocate( Value_OF_Normalize &
             ( Number_Objective_Function, &
               0:Number_Loop_Position_Source_Y, 0:Number_Loop_Position_Source_X, & 
               0:Number_Loop_Material_Property, 0:Number_Loop_IncidentAngle, 0:Number_Loop_Frequency ) )

   if( Flag_Structure==0 )then
   !========================================================================================================================================
   write(*,*)'Read Normalizing Value'
   !========================================================================================================================================
 
      !open( 11, file='../No_Structure/Result_Normalization.dat' )
      open( 11, file=trim(File_Pass_Result_Normalization)//trim(File_Name_Result_Normalization) ) 

      !if( Number_Objective_Function==1 )then
         do LoopCounter_Frequency= 0, Number_Loop_Frequency
            do LoopCounter_IncidentAngle= 0, Number_Loop_IncidentAngle
               do LoopCounter_Material_Property= 0, Number_Loop_Material_Property 
                  do LoopCounter_Position_Source_X= 0, Number_Loop_Position_Source_X 
                     do LoopCounter_Position_Source_Y= 0, Number_Loop_Position_Source_Y
                        read(  11, * ) &
                        Value_OF_Normalize & 
                        ( 1, LoopCounter_Position_Source_Y, LoopCounter_Position_Source_X, &
                          LoopCounter_Material_Property, LoopCounter_IncidentAngle, LoopCounter_Frequency )
                     end do
                  end do
               end do
            end do
         end do
      !end if
      close( 11 )

      if( Number_Objective_Function==2 .and. &
          ( Flag_Thermal_Device==20 .or. Flag_Thermal_Device==21 .or. Flag_Thermal_Device==30 ) )then
         open( 11, file=trim(File_Pass_Reference_Field)//trim(File_Name_Result_Normalization) ) 
            do LoopCounter_Frequency= 0, Number_Loop_Frequency
               do LoopCounter_IncidentAngle= 0, Number_Loop_IncidentAngle
                  do LoopCounter_Material_Property= 0, Number_Loop_Material_Property 
                     do LoopCounter_Position_Source_X= 0, Number_Loop_Position_Source_X 
                        do LoopCounter_Position_Source_Y= 0, Number_Loop_Position_Source_Y
                           read(  11, * ) &
                           Value_OF_Normalize & 
                           ( 2, LoopCounter_Position_Source_Y, LoopCounter_Position_Source_X, &
                             LoopCounter_Material_Property, LoopCounter_IncidentAngle, LoopCounter_Frequency )
                        end do
                     end do
                  end do
               end do
            end do
         close( 11 )
      end if

   else
      Value_OF_Normalize= 1d0
   end if

   !========================================================================================================================================
   write(*,*)'Compute Range of Parameters'
   !========================================================================================================================================

   allocate( Temperature_Left_Side_Boundary( 0:Number_Loop_Frequency ) )
   allocate( Temperature_Right_Side_Boundary( 0:Number_Loop_IncidentAngle ) )
   allocate( Thermal_Conductivity_Material_Real( 0:Number_Loop_Material_Property ) )
   allocate( Thermal_Conductivity_Material( 0:Number_Loop_Material_Property ) )
   allocate( Thermal_Conductivity_FixedDomain_Variable( 0:Number_Loop_Material_Property ) )

   allocate( Position_Source_X( 0:Number_Loop_Position_Source_X ) )
   allocate( Position_Source_Y( 0:Number_Loop_Position_Source_Y ) )

   do LoopCounter_Frequency= 0, Number_Loop_Frequency 
   
      if( Number_Loop_Frequency==0 )then
         Temperature_Left_Side_Boundary( LoopCounter_Frequency )= Normal_Temperature_Lower_Side_Boundary
      else
         Temperature_Left_Side_Boundary( LoopCounter_Frequency ) &
         = Minimum_Temperature_Lower_Side_Boundary &
          +LoopCounter_Frequency* ( Maximum_Temperature_Lower_Side_Boundary -Minimum_Temperature_Lower_Side_Boundary )/Number_Loop_Frequency
      end if
   end do

   do LoopCounter_IncidentAngle= 0, Number_Loop_IncidentAngle  

      if( Number_Loop_IncidentAngle==0 )then
         Temperature_Right_Side_Boundary( LoopCounter_IncidentAngle )= Normal_Temperature_Higher_Side_Boundary
      else
         Temperature_Right_Side_Boundary( LoopCounter_IncidentAngle ) &
         = Minimum_Temperature_Higher_Side_Boundary &
          +LoopCounter_IncidentAngle*( Maximum_Temperature_Higher_Side_Boundary -Minimum_Temperature_Higher_Side_Boundary )/Number_Loop_IncidentAngle
      end if
   end do

   !====================================================================
   if( Flag_Thermal_Device==20 .or. Flag_Thermal_Device==21 .or. Flag_Thermal_Device==30 )then
   !====================================================================
      do LoopCounter_Material_Property= 0, Number_Loop_Material_Property 
   
         if( Number_Loop_Material_Property==0 )then
            Thermal_Conductivity_Material_Real( LoopCounter_Material_Property )= Thermal_Conductivity_FixedDomain !Normal_Thermal_Conductivity_Material
         else
            Thermal_Conductivity_Material_Real( LoopCounter_Material_Property ) &
            = Minimum_Thermal_Conductivity &
             +LoopCounter_Material_Property *( Maximum_Thermal_Conductivity -Minimum_Thermal_Conductivity )/Number_Loop_Material_Property 
         end if
      
         !Thermal_Conductivity_Material( LoopCounter_Material_Property )= Thermal_Conductivity_Material_Real( LoopCounter_Material_Property )
         Thermal_Conductivity_Material( LoopCounter_Material_Property )= Normal_Thermal_Conductivity_Material 
   
         if( Flag_Material_FixedDomain==0 )then
            Thermal_Conductivity_FixedDomain_Variable( LoopCounter_Material_Property )= Thermal_Conductivity_Material( LoopCounter_Material_Property )
   
         else if( Flag_Material_FixedDomain==1 )then
            Thermal_Conductivity_FixedDomain_Variable( LoopCounter_Material_Property )= Thermal_Conductivity_Base_Material
         
         else if( Flag_Material_FixedDomain==2 )then
            if( Flag_Structure==0 )then
               Thermal_Conductivity_FixedDomain_Variable( LoopCounter_Material_Property )= Thermal_Conductivity_Material_Real( LoopCounter_Material_Property ) 
            else 
               Thermal_Conductivity_FixedDomain_Variable( LoopCounter_Material_Property )= Thermal_Conductivity_FixedDomain 
            end if
         end if
      end do
   !====================================================================
   else
   !====================================================================
      do LoopCounter_Material_Property= 0, Number_Loop_Material_Property 
   
         if( Number_Loop_Material_Property==0 )then
            Thermal_Conductivity_Material_Real( LoopCounter_Material_Property )= Normal_Thermal_Conductivity_Material
         else
            Thermal_Conductivity_Material_Real( LoopCounter_Material_Property ) &
            = Minimum_Thermal_Conductivity &
             +LoopCounter_Material_Property *( Maximum_Thermal_Conductivity -Minimum_Thermal_Conductivity )/Number_Loop_Material_Property 
         end if
      
         Thermal_Conductivity_Material( LoopCounter_Material_Property )= Thermal_Conductivity_Material_Real( LoopCounter_Material_Property )
   
   
         if( Flag_Material_FixedDomain==0 )then
            Thermal_Conductivity_FixedDomain_Variable( LoopCounter_Material_Property )= Thermal_Conductivity_Material( LoopCounter_Material_Property )
   
         else if( Flag_Material_FixedDomain==1 )then
            Thermal_Conductivity_FixedDomain_Variable( LoopCounter_Material_Property )= Thermal_Conductivity_Base_Material
         
         else if( Flag_Material_FixedDomain==2 )then
            Thermal_Conductivity_FixedDomain_Variable( LoopCounter_Material_Property )= Thermal_Conductivity_FixedDomain 
         
         end if
      end do
   !====================================================================
   end if 
   !====================================================================

   do LoopCounter_Position_Source_X= 0, Number_Loop_Position_Source_X 

      if( Number_Loop_Position_Source_X==0 )then
         Position_Source_X( LoopCounter_Position_Source_X )= Normal_Position_Source_x
      else
         Position_Source_X( LoopCounter_Position_Source_X ) &
         = Minimum_Position_Source_x &
          +LoopCounter_Position_Source_X*( Maximum_Position_Source_x -Minimum_Position_Source_x )/Number_Loop_Position_Source_X
      end if
   end do

   do LoopCounter_Position_Source_Y= 0, Number_Loop_Position_Source_Y

      if( Number_Loop_Position_Source_Y==0 )then
         Position_Source_Y( LoopCounter_Position_Source_Y )= Normal_Position_Source_y
      else
         Position_Source_Y( LoopCounter_Position_Source_Y ) &
         = Minimum_Position_Source_y &
          +LoopCounter_Position_Source_Y*( Maximum_Position_Source_y -Minimum_Position_Source_y )/Number_Loop_Position_Source_Y
      end if
   end do

   allocate( Value_Objective_Function( NumThread, Number_Node ) )
   allocate( Value_Plot_Distribution( Number_Node ) )

   allocate( Normalized_Objective_Function&
             ( Number_Objective_Function, &
               0:Number_Loop_Position_Source_Y, 0:Number_Loop_Position_Source_X, & 
               0:Number_Loop_Material_Property, 0:Number_Loop_IncidentAngle, 0:Number_Loop_Frequency ) )

   allocate( Original_OF_Value_without_Normalization&
             ( Number_Objective_Function, &
               0:Number_Loop_Position_Source_Y, 0:Number_Loop_Position_Source_X, &
               0:Number_Loop_Material_Property, 0:Number_Loop_IncidentAngle, 0:Number_Loop_Frequency ) ) 


   !==========================
   do Loop_Core= 1, NumThread
   !==========================
 
      !==================================================================================================================
      do LoopCounter_Frequency= LoopCounter_Frequency_Range( 1, Loop_Core ), LoopCounter_Frequency_Range( 2, Loop_Core )
      !==================================================================================================================

         !==============================================================================================================================
         do LoopCounter_IncidentAngle= LoopCounter_IncidentAngle_Range( 1, Loop_Core ), LoopCounter_IncidentAngle_Range( 2, Loop_Core ) 
         !==============================================================================================================================
    
            !===============================================================================================================
            do LoopCounter_Material_Property &
               = LoopCounter_Material_Property_Range( 1, Loop_Core ), LoopCounter_Material_Property_Range( 2, Loop_Core ) 
            !===============================================================================================================

               !=============================================================================================================
               do LoopCounter_Position_Source_X &
                  = LoopCounter_Position_Source_X_Range( 1, Loop_Core ), LoopCounter_Position_Source_X_Range( 2, Loop_Core )
               !=============================================================================================================
   
                  !=============================================================================================================
                  do LoopCounter_Position_Source_Y &
                     = LoopCounter_Position_Source_Y_Range( 1, Loop_Core ), LoopCounter_Position_Source_Y_Range( 2, Loop_Core )
                  !=============================================================================================================
      
                     write(*,*)'LoopCounter_Frequency=', LoopCounter_Frequency
                     write(*,*)'Temperature_Left_Side_Boundary=', Temperature_Left_Side_Boundary( LoopCounter_Frequency ) 
                     write(*,*)'LoopCounter_IncidentAngle=', LoopCounter_IncidentAngle
                     write(*,*)'Temperature_Right_Side_Boundary=', Temperature_Right_Side_Boundary( LoopCounter_IncidentAngle )
                     write(*,*)'LoopCounter_Material_Property=', LoopCounter_Material_Property
                     write(*,*)'Thermal_Conductivity_Material=', Thermal_Conductivity_Material( LoopCounter_Material_Property )
                     write(*,*)'Position_Source_X=', Position_Source_X( LoopCounter_Position_Source_X ) 
                     write(*,*)'Position_Source_Y=', Position_Source_Y( LoopCounter_Position_Source_Y ) 

                     allocate( Temperature_Solution( Number_Node ) )
                     allocate( Thermal_Conductivity_Element( Number_Element_Triangle ) )

                     call Analyze_Steady_State_Heat_Conduction &
                        ( Flag_Structure,&
                          Temperature_Left_Side_Boundary( LoopCounter_Frequency ), &
                          Temperature_Right_Side_Boundary( LoopCounter_IncidentAngle ), &
                          Flag_Thermal_Insulation_BC, &
                          Thermal_Conductivity_Material( LoopCounter_Material_Property ), &
                          Thermal_Conductivity_Base_Material, & 
                          Thermal_Conductivity_OpenRegion, &
                          Thermal_Conductivity_FixedDomain_Variable( LoopCounter_Material_Property ), &
                          ID_Element_Material, ID_Element_OpenRegion, ID_Element_FixedDomain, ID_Element_Base_Material,  & 
                          ID_Element_Material_Exterior, ID_Element_Base_Material_Exterior, & 
                          Number_Node, Number_Element_Triangle, &
                          Position_Node, &
                          Width_Matrix_LHS, Max_Position_Node, &
                          Index_Element_2_Node_Triangle, Class_Element_Triangle, Thermal_conductivity_by_density, &
                          Position_Source_X( LoopCounter_Position_Source_X ), Position_Source_Y( LoopCounter_Position_Source_Y ), &
                          high, thickness, penalty_density, &
                          !============================================================================================
                          Temperature_Solution, Thermal_Conductivity_Element )

if ( Flag_Structure == 2 ) then
   open( 444, file='aho_T_ref.dat' )
        do i = 1, Number_Node
           write( 444,* ) Temperature_Solution( i )
        end do
   close( 444 )
end if

                     !========================================================================
                     write(*,*)'   Obtain Reference Data for Objective Function'
                     !========================================================================
 
                     if( Flag_Reference_Device==1 )then
   
                        call Counter_Number_Node_in_Element &
                           ( Number_Node, Number_Element_Triangle, &
                             Index_Element_2_Node_Triangle, Class_Element_Triangle, ID_Element_ObjectiveFunction, &
                             Number_Node_Reference )
  
                        Number_Node_tmp=Number_Node_Reference
                        if( Flag_Reference_Data==1 .or. Flag_Theoretical_Solution==1 ) Number_Node_tmp=Number_Node
                  
                        if( Flag_Theoretical_Solution==1 )then

                           allocate( Temperature_Reference( Number_Node_tmp ) )

                           do i= 1, Number_Node_tmp
                              Temperature_Reference( i ) &
                              = Position_Node( 1, i ) *( 1.0d0 -0.0d0 )/( Max_Position_Node( 1, 2 ) -Max_Position_Node( 1, 1 ) ) &
                               +( 1.0d0 +0.0d0 )/2d0
                           end do
                        else
 
                           allocate( Temperature_Reference( Number_Node_tmp ) )
write(*,*)'Flag_Structure= main.f90', Flag_Structure
                           call Operate_Data_Reference_Field &
                              ( Flag_Structure, File_Pass_Reference_Field, &
                                Number_Node, Number_Node_Reference, Number_Node_tmp, Temperature_Solution, &
                                dble(LoopCounter_Position_Source_X), dble(LoopCounter_Position_Source_Y), &
                                !Temperature_Left_Side_Boundary( LoopCounter_Frequency ), Temperature_Right_Side_Boundary( LoopCounter_IncidentAngle ), &
                                Temperature_Reference, Number_Node_Reference )
                        
                        end if 
                     end if
   
                     if( Flag_Thermal_Device==0 .or. Flag_Thermal_Device==10 .or. Flag_Thermal_Device==20 .or. &
                         Flag_Thermal_Device==21 .or. Flag_Thermal_Device==30 )then
                        !do i= 1, Number_Node
                        !   Value_Objective_Function( Loop_Core, i )= Temperature_Solution( i )/Amplitude_Incident_Wave
                        !end do

                        Value_Objective_Function( Loop_Core, : )= Temperature_Solution( : )
                        do i= 1, Number_Node_tmp 
                           Value_Objective_Function( Loop_Core, i )= ( Temperature_Solution( i ) -Temperature_Reference( i ) )
                        end do

                     else
                        call Output_Error( 'Main', 522 ) 
                     end if
   
                     call Compute_Objective_Function &
                        ( Flag_Structure, ID_Element_ObjectiveFunction, Loop_Core, &
                          Value_Objective_Function, & 
                          Position_Node, Number_Node, &
                          Index_Element_2_Node_Triangle, Class_Element_Triangle, Number_Element_Triangle, &
                          Value_OF_Normalize& 
                          ( 1, LoopCounter_Position_Source_Y, LoopCounter_Position_Source_X, &
                            LoopCounter_Material_Property, LoopCounter_IncidentAngle, LoopCounter_Frequency ),  &
                          !=================================================================
                          Normalized_Objective_Function&
                          ( 1, LoopCounter_Position_Source_Y, LoopCounter_Position_Source_X, &
                            LoopCounter_Material_Property, LoopCounter_IncidentAngle, LoopCounter_Frequency ), & 
                          Original_OF_Value_without_Normalization&
                          ( 1, LoopCounter_Position_Source_Y, LoopCounter_Position_Source_X, &
                            LoopCounter_Material_Property, LoopCounter_IncidentAngle, LoopCounter_Frequency ) )
          

                     if( ( Flag_Thermal_Device==20 .or. Flag_Thermal_Device==21  .or. Flag_Thermal_Device==30 ) .and. &
                           Number_Objective_Function==2 )then

                        call Compute_Objective_Function_Flux &
                             ( Flag_Structure, ID_Element_FixedDomain, &
                               Temperature_Solution, Number_Node, &
                               Position_Node, Number_Node, &
                               Index_Element_2_Node_Triangle, Class_Element_Triangle, Number_Element_Triangle, &
                               Thermal_Conductivity_FixedDomain_Variable( LoopCounter_Material_Property ), &
                               Value_OF_Normalize& 
                               ( 2, LoopCounter_Position_Source_Y, LoopCounter_Position_Source_X, &
                                 LoopCounter_Material_Property, LoopCounter_IncidentAngle, LoopCounter_Frequency ),  &
                               !=================================================================
                               Normalized_Objective_Function&
                               ( 2, LoopCounter_Position_Source_Y, LoopCounter_Position_Source_X, &
                                 LoopCounter_Material_Property, LoopCounter_IncidentAngle, LoopCounter_Frequency ), & 
                               Original_OF_Value_without_Normalization&
                               ( 2, LoopCounter_Position_Source_Y, LoopCounter_Position_Source_X, &
                                 LoopCounter_Material_Property, LoopCounter_IncidentAngle, LoopCounter_Frequency ) )

                     end if

                     !======================================================================================================== 
                     if( Flag_Plot_Configuration==1 )then
                     !======================================================================================================== 
          
                        allocate( Index_Element_2_Node_Plot( 4, Number_Element_Triangle ) )
         
                        do i= 1, Number_Element_Triangle
                           do j= 1, 3      
                              Index_Element_2_Node_Plot( j, i )= Index_Element_2_Node_Triangle( j, i )  
                           end do
                           Index_Element_2_Node_Plot( 4, i )= Index_Element_2_Node_Plot( 3, i ) 
                        end do
    
                        call Plot_Configuration &
                           ( Position_Node, Number_Node, &
                             Index_Element_2_Node_Plot, Number_Element_Triangle, &
                             Class_Element_Triangle, &
                             ID_Element_Material, ID_Element_OpenRegion, ID_Element_FixedDomain, ID_Element_Base_Material, &
                             ID_Element_Material_Exterior, ID_Element_Base_Material_Exterior )   
   
                        deallocate( Index_Element_2_Node_Plot )
                     end if
   
                     !======================================================================================================== 
                     if( Flag_Plot_Mesh_Configuration==1 )then
                     !======================================================================================================== 
          
                        allocate( Index_Element_2_Node_Plot( 4, Number_Element_Triangle ) )
         
                        do i= 1, Number_Element_Triangle
                           do j= 1, 3      
                              Index_Element_2_Node_Plot( j, i )= Index_Element_2_Node_Triangle( j, i )  
                           end do
                           Index_Element_2_Node_Plot( 4, i )= Index_Element_2_Node_Plot( 3, i ) 
                        end do
    
                        !call Plot_Mesh_Configuration &
                        !   ( Position_Node, Number_Node, &
                        !     Index_Element_2_Node_Plot, Number_Element_Triangle, &
                        !     Class_Element_Triangle, &
                        !     ID_Element_Material, ID_Element_OpenRegion, ID_Element_FixedDomain, ID_Element_Base_Material, &
                        !     ID_Element_Material_Exterior, ID_Element_Base_Material_Exterior, &
                        !     ID_Element_PML_Xp, ID_Element_PML_Xm, ID_Element_PML_Yp, ID_Element_PML_Ym, &
                        !     ID_Element_PML_XpYp,  ID_Element_PML_XmYp, ID_Element_PML_XpYm,  ID_Element_PML_XmYm )
   
                        deallocate( Index_Element_2_Node_Plot )
                     end if

                     if( Flag_Output_PhysicalField==1 )then
         
                        if( Flag_Thermal_Device==0 .or. Flag_Thermal_Device==10 )then 
                           do i= 1, Number_Node
                              Value_Plot_Distribution( i )= Temperature_Solution( i ) 
                           end do
                        else       
                           do i= 1, Number_Node
                              Value_Plot_Distribution( i )= Temperature_Solution( i )
                           end do
                        end if
   
                        !================================================
                        if( ( Flag_Thermal_Device==50 .or. Flag_Thermal_Device==51 ) .and. Flag_Thermal_Insulation_BC==1 )then
                        !================================================
   
                           Counter_Element_Plot= 0
                           do i= 1, Number_Element_Triangle
                              if( Class_Element_Triangle( i ) /= ID_Element_FixedDomain )then
                                 Counter_Element_Plot= Counter_Element_Plot +1
                              end if
                           end do
   
                           Number_Element_Plot= Counter_Element_Plot
   
                           allocate( Index_Element_2_Node_Plot( 4, Number_Element_Plot ) )
   
                           Counter_Element_Plot= 0
                           do i= 1, Number_Element_Triangle
                              if( Class_Element_Triangle( i ) /= ID_Element_FixedDomain )then
                                 Counter_Element_Plot= Counter_Element_Plot +1
                                 do j= 1, 3
                                    Index_Element_2_Node_Plot( j, Counter_Element_Plot )= Index_Element_2_Node_Triangle( j, i )
                                 end do
                              end if
                           end do
 
                        !================================================
                        else
                        !================================================
                           Number_Element_Plot= Number_Element_Triangle
   
                           allocate( Index_Element_2_Node_Plot( 4, Number_Element_Plot ) )
         
                           do i= 1, Number_Element_Plot
                              do j= 1, 3      
                                 Index_Element_2_Node_Plot( j, i )= Index_Element_2_Node_Triangle( j, i )  
                              end do
                           end do
                     
                        end if
   
                        do i= 1, Number_Element_Plot
                           Index_Element_2_Node_Plot( 4, i )= Index_Element_2_Node_Plot( 3, i )
                        end do
         
                        call Plot_Result_Element_Triangle &
                           ( 100, &!200 +LoopCounter_Frequency +LoopCounter_IncidentAngle +LoopCounter_Material_Property &
                                  !+LoopCounter_Position_Source_X + LoopCounter_Position_Source_Y, &
                             Value_Plot_Distribution, Number_Level_Plot_EM, & 
                             LoopCounter_Frequency +LoopCounter_IncidentAngle +LoopCounter_Material_Property & 
                             +LoopCounter_Position_Source_X +LoopCounter_Position_Source_Y , &
                             Position_Node, Number_Node,  &
                             Index_Element_2_Node_Plot, Number_Element_Plot, &
                             Normalized_Objective_Function&
                             ( 1, LoopCounter_Position_Source_Y, LoopCounter_Position_Source_X, & 
                               LoopCounter_Material_Property, LoopCounter_IncidentAngle, LoopCounter_Frequency ) ) 
         
                        deallocate( Index_Element_2_Node_Plot )
         
                     end if
   
                     !======================================================================================================== 
                     if( Flag_Plot_EFD_Configuration==1 )then
                     !======================================================================================================== 
           
                        if( Flag_Thermal_Device==0 .or. Flag_Thermal_Device==10 )then 
                           do i= 1, Number_Node
                              Value_Plot_Distribution( i )= Temperature_Solution( i )
                           end do
                        else       
                           do i= 1, Number_Node
                              Value_Plot_Distribution( i )= Temperature_Solution( i ) 
                           end do
                        end if
   
                        !================================================
                        if( ( Flag_Thermal_Device==50 .or. Flag_Thermal_Device==51 ) .and. Flag_Thermal_Insulation_BC==1 )then
                        !================================================
   
                           Counter_Element_Plot= 0
                           do i= 1, Number_Element_Triangle
                              if( Class_Element_Triangle( i ) /= ID_Element_FixedDomain )then
                                 Counter_Element_Plot= Counter_Element_Plot +1
                              end if
                           end do
   
                           Number_Element_Plot= Counter_Element_Plot
   
                           allocate( Index_Element_2_Node_Plot( 4, Number_Element_Plot ) )
                           allocate( Class_Element_Plot( Number_Element_Plot ) )
   
                           Counter_Element_Plot= 0
                           do i= 1, Number_Element_Triangle
                              if( Class_Element_Triangle( i ) /= ID_Element_FixedDomain )then
                                 Counter_Element_Plot= Counter_Element_Plot +1
                                 do j= 1, 3
                                    Index_Element_2_Node_Plot( j, Counter_Element_Plot )= Index_Element_2_Node_Triangle( j, i )
                                 end do
                                 Class_Element_Plot( Counter_Element_Plot )= Class_Element_Triangle( i )
                              end if
                           end do
                        !================================================
                        else
                        !================================================
                           Number_Element_Plot= Number_Element_Triangle
   
                           allocate( Index_Element_2_Node_Plot( 4, Number_Element_Plot ) )
                           allocate( Class_Element_Plot( Number_Element_Plot ) )
         
                           do i= 1, Number_Element_Plot
                              do j= 1, 3      
                                 Index_Element_2_Node_Plot( j, i )= Index_Element_2_Node_Triangle( j, i )  
                              end do
                              Class_Element_Plot( i )= Class_Element_Triangle( i )
                           end do
                     
                        end if
   
                        do i= 1, Number_Element_Plot
                           Index_Element_2_Node_Plot( 4, i )= Index_Element_2_Node_Plot( 3, i )
                        end do
         
                        call Plot_Result_and_Configuration_Element_Triangle &
                           ( 200 +LoopCounter_Frequency +LoopCounter_IncidentAngle +LoopCounter_Material_Property &
                                 +LoopCounter_Position_Source_X + LoopCounter_Position_Source_Y, &
                             Value_Plot_Distribution, Number_Level_Plot_EM, & 
                             LoopCounter_Frequency +LoopCounter_IncidentAngle +LoopCounter_Material_Property & 
                             +LoopCounter_Position_Source_X +LoopCounter_Position_Source_Y , &
                             Position_Node, Number_Node,  &
                             Index_Element_2_Node_Plot, Number_Element_Plot, &
                             Class_Element_Plot, ID_Element_Material, &
                             !Class_Element_Plot, ID_Element_Base_Material, &
                             !Class_Element_Plot, ID_Element_OpenRegion, &
                             Normalized_Objective_Function&
                             ( 1, LoopCounter_Position_Source_Y, LoopCounter_Position_Source_X, & 
                               LoopCounter_Material_Property, LoopCounter_IncidentAngle, LoopCounter_Frequency ) ) 
         
                        deallocate( Index_Element_2_Node_Plot )
                        deallocate( Class_Element_Plot )
         
                     end if
   
                     !======================================================================================================== 
                     if( Flag_Plot_EFD_OC_PoyntingVector==1 )then
                     !======================================================================================================== 
    
                        if( Flag_Thermal_Device==0 .or. Flag_Thermal_Device==10 )then 
                           do i= 1, Number_Node
                              Value_Plot_Distribution( i )= Temperature_Solution( i ) 
                           end do
                        else       
                           do i= 1, Number_Node
                              Value_Plot_Distribution( i )= Temperature_Solution( i )
                           end do
                        end if
   
                        !================================================
                        if( ( Flag_Thermal_Device==50 .or. Flag_Thermal_Device==51 ) .and. Flag_Thermal_Insulation_BC==1 )then
                        !================================================
   
                           Counter_Element_Plot= 0
                           do i= 1, Number_Element_Triangle
                              if( Class_Element_Triangle( i ) /= ID_Element_FixedDomain )then
                                 Counter_Element_Plot= Counter_Element_Plot +1
                              end if
                           end do
   
                           Number_Element_Plot= Counter_Element_Plot
   
                           allocate( Index_Element_2_Node_Plot( 4, Number_Element_Plot ) )
                           allocate( Class_Element_Plot( Number_Element_Plot ) )
   
                           Counter_Element_Plot= 0
                           do i= 1, Number_Element_Triangle
                              if( Class_Element_Triangle( i ) /= ID_Element_FixedDomain )then
                                 Counter_Element_Plot= Counter_Element_Plot +1
                                 do j= 1, 3
                                    Index_Element_2_Node_Plot( j, Counter_Element_Plot )= Index_Element_2_Node_Triangle( j, i )
                                 end do
                                 Class_Element_Plot( Counter_Element_Plot )= Class_Element_Triangle( i )
                              end if
                           end do
                        !================================================
                        else
                        !================================================
                           Number_Element_Plot= Number_Element_Triangle
   
                           allocate( Index_Element_2_Node_Plot( 4, Number_Element_Plot ) )
                           allocate( Class_Element_Plot( Number_Element_Plot ) )
         
                           do i= 1, Number_Element_Plot
                              do j= 1, 3      
                                 Index_Element_2_Node_Plot( j, i )= Index_Element_2_Node_Triangle( j, i )  
                              end do
                              Class_Element_Plot( i )= Class_Element_Triangle( i )
                           end do
                     
                        end if
   
                        do i= 1, Number_Element_Plot
                           Index_Element_2_Node_Plot( 4, i )= Index_Element_2_Node_Plot( 3, i )
                        end do
   
                        !================================================
                        write(*,*)'   Compute PV Data' 
                        !================================================
   
                        Number_Point_PV_Y= 20 ! >= 20 
                        Number_Point_PV_X= Number_Point_PV_Y*2
   
                        allocate( Position_PV_tmp( 2, Number_Point_PV_X*Number_Point_PV_Y ) ) 
                        allocate( Element_Number_PV_tmp( Number_Point_PV_X*Number_Point_PV_Y ) ) 
   
                        write(*,*)'Number_Element_Plot=', Number_Element_Plot

                        call Compute_Position_and_Element_PV&
                           ( Number_Point_PV_X, Number_Point_PV_Y, &
                             Position_Node, Number_Node,  &
                             Index_Element_2_Node_Plot, Number_Element_Plot, &
                             Position_PV_tmp, Element_Number_PV_tmp, Number_PV )
   
                        allocate( Position_PV( 2, Number_PV ) ) 
   
                        do i= 1, Number_PV 
                           do j= 1, 2
                              Position_PV( j, i )= Position_PV_tmp( j, i )
                           end do
                        end do
   
                        deallocate( Position_PV_tmp ) 
                        allocate( Element_Number_PV( Number_PV ) ) 
   
                        do i= 1, Number_PV 
                           Element_Number_PV( i )= Element_Number_PV_tmp( i )
                        end do
                        deallocate( Element_Number_PV_tmp ) 
   
                        !do i= 1, Number_PV 
                        !   write(888,*) Position_PV( 1, i ), Position_PV( 2, i )
                        !end do
   
                        !do i= 1, Number_PV 
                        !   do k= 1, 3 
                        !      write(889,*) Position_Node( 1, Index_Element_2_Node_Plot( k, Element_Number_PV( i ) ) ), & 
                        !               Position_Node( 2, Index_Element_2_Node_Plot( k, Element_Number_PV( i ) ) )  
                        !   end do
                        !end do
   
                        write(*,*)'   end Compute PV Data' 
   
                        !================================================
                        write(*,*)'   Compute PV' 
                        !================================================
    
                        !Length_Maximum_Vector= 2d0/Number_Point_PV_X*1.5d0
                        !Length_Minimum_Vector= 2d0/Number_Point_PV_X*0.5d0 
                        Length_Maximum_Vector= 2d0/Number_Point_PV_X*1.8d0
                        Length_Minimum_Vector= 2d0/Number_Point_PV_X*0.1d0 
                        !Length_Maximum_Vector= 2d0/Number_Point_PV_Y
                        !Length_Maximum_Vector= 2d0/Number_Point_PV_Y*5d0
                        !Length_Minimum_Vector= 2d0/Number_Point_PV_Y/2d0 
   
                        allocate( Position_Arrow_PV( 2, 2, Number_PV ) ) 
                        allocate( Level_PV( Number_PV ) ) 
   
                        if( Flag_Thermal_Device==0 .or. Flag_Thermal_Device==10 )then 
  
                           allocate( Thermal_Conductivity_PV( Number_PV ) )

                           do i= 1, Number_PV 
                              Thermal_Conductivity_PV( i ) = Thermal_Conductivity_Element( Element_Number_PV( i ) )
                           end do
 
                           call Compute_Poynting_Vector &
                                ( Number_PV, Position_PV, Element_Number_PV, &
                                  Length_Maximum_Vector, Length_Minimum_Vector, &
                                  Number_Node, Position_Node, Temperature_Solution, &
                                  Index_Element_2_Node_Plot, Number_Element_Plot, &
                                  Thermal_Conductivity_PV, &
                                  !================================================
                                  Position_Arrow_PV, Level_PV )

                           allocate( Position_Arrow_tmp( 2, 2, Number_PV ) ) 
                           allocate( Level_PV_tmp( Number_PV ) ) 

                           Counter= 0
                           if( Flag_Thermal_Insulation_BC==1 )then 
                              do e= 1, Number_PV
                                 if( Class_Element_Triangle( Element_Number_PV( e ) )/=ID_Element_FixedDomain )then 
                                    Counter= Counter +1
                                    Position_Arrow_tmp(:,:,Counter) = Position_Arrow_PV(:,:,e)
                                    Level_PV_tmp( Counter ) = Level_PV( e )
                                 end if
                              end do

                              deallocate( Position_Arrow_PV ) 
                              deallocate( Level_PV ) 
                              Number_PV= Counter
                              allocate( Position_Arrow_PV( 2, 2, Number_PV ) ) 
                              allocate( Level_PV( Number_PV ) ) 
                              Position_Arrow_PV( :, :, 1:Number_PV )=Position_Arrow_tmp( :,:,1:Number_PV )
                              Level_PV( 1:Number_PV )=Level_PV_tmp( 1:Number_PV )

                           end if 

                           deallocate( Position_Arrow_tmp ) 
                           deallocate( Level_PV_tmp ) 
                           deallocate( Thermal_Conductivity_PV )

                        end if
   
                        deallocate( Element_Number_PV ) 
                        deallocate( Position_PV ) 
   
                        !do e= 1, Number_PV
                        !   write(1469,1475) Position_Arrow_PV( 1, 1, e ), Position_Arrow_PV( 2, 1, e ), &
                        !                    Position_Arrow_PV( 1, 2, e ), Position_Arrow_PV( 2, 2, e )
                        !   !write(1469,1475) ( Position_Arrow_PV( 1, 1, e ) +Position_Arrow_PV( 1, 2, e ) )/2d0, &
                        !   !         ( Position_Arrow_PV( 2, 1, e ) +Position_Arrow_PV( 2, 2, e ) )/2d0, &
                        !   !         Position_Arrow_PV( 1, 2, e ) -Position_Arrow_PV( 1, 1, e ), &
                        !   !         Position_Arrow_PV( 2, 2, e ) -Position_Arrow_PV( 2, 1, e )
                        !   1475 format(es15.8,1x,es15.8,1x,es15.8,1x,es15.8) 
                        !end do
   
                        !================================================
                        write(*,*)'   Subroutine Plot_HF_T_OC' 
                        !================================================
         
                        !call Plot_HF_T_OC &
                        !   ( 200 +LoopCounter_Frequency +LoopCounter_IncidentAngle +LoopCounter_Material_Property &
                        !          +LoopCounter_Position_Source_X + LoopCounter_Position_Source_Y, &
                        !     Value_Plot_Distribution, Number_Level_Plot_EM, & 
                        !     LoopCounter_Frequency +LoopCounter_IncidentAngle +LoopCounter_Material_Property & 
                        !     +LoopCounter_Position_Source_X +LoopCounter_Position_Source_Y , &
                        !     Position_Node, Number_Node,  &
                        !     Index_Element_2_Node_Plot, Number_Element_Plot, &
                        !     Class_Element_Plot, ID_Element_Material, ID_Element_FixedDomain, &
                        !     Number_PV, Position_Arrow_PV, Length_Maximum_Vector, Level_PV, &
                        !     Normalized_Objective_Function&
                        !     ( 1, LoopCounter_Position_Source_Y, LoopCounter_Position_Source_X, & 
                        !       LoopCounter_Material_Property, LoopCounter_IncidentAngle, LoopCounter_Frequency ) ) 

                        call Plot_Contour_T_OC &
                           ( 200 +LoopCounter_Frequency +LoopCounter_IncidentAngle +LoopCounter_Material_Property &
                                  +LoopCounter_Position_Source_X + LoopCounter_Position_Source_Y, &
                             Flag_Thermal_Insulation_BC, Value_Plot_Distribution, Number_Level_Plot_EM, & 
                             LoopCounter_Frequency +LoopCounter_IncidentAngle +LoopCounter_Material_Property & 
                             +LoopCounter_Position_Source_X +LoopCounter_Position_Source_Y , &
                             Position_Node, Number_Node,  &
                             Index_Element_2_Node_Plot, Number_Element_Plot, &
                             Class_Element_Plot, ID_Element_Material, ID_Element_FixedDomain, &
                             Number_PV, Position_Arrow_PV, Length_Maximum_Vector, Level_PV, &
                             Normalized_Objective_Function&
                             ( 1, LoopCounter_Position_Source_Y, LoopCounter_Position_Source_X, & 
                               LoopCounter_Material_Property, LoopCounter_IncidentAngle, LoopCounter_Frequency ) ) 
 
    !                    call Plot_Gradient_T_OC &
     !                      ( 1064 +LoopCounter_Frequency +LoopCounter_IncidentAngle +LoopCounter_Material_Property &
      !                            +LoopCounter_Position_Source_X + LoopCounter_Position_Source_Y, &
       !                      Flag_Thermal_Insulation_BC, Value_Plot_Distribution, Thermal_Conductivity_Element, &
        !                     Number_Level_Plot_EM, & 
         !                    LoopCounter_Frequency +LoopCounter_IncidentAngle +LoopCounter_Material_Property & 
          !                   +LoopCounter_Position_Source_X +LoopCounter_Position_Source_Y , &
    !                         Position_Node, Number_Node,  &
     !                        Index_Element_2_Node_Plot, Number_Element_Plot, &
      !                       Class_Element_Plot, ID_Element_Material, ID_Element_FixedDomain, &
       !                      Normalized_Objective_Function&
        !                     ( 1, LoopCounter_Position_Source_Y, LoopCounter_Position_Source_X, & 
         !                      LoopCounter_Material_Property, LoopCounter_IncidentAngle, LoopCounter_Frequency ) ) 

 
!                        call Plot_Gradient_Arrow_T_OC &
 !                          ( 1177 +LoopCounter_Frequency +LoopCounter_IncidentAngle +LoopCounter_Material_Property &
  !                                +LoopCounter_Position_Source_X + LoopCounter_Position_Source_Y, &
   !!                          Flag_Thermal_Insulation_BC, Value_Plot_Distribution, Thermal_Conductivity_Element, &
    !                         Number_Level_Plot_EM, & 
     !                        LoopCounter_Frequency +LoopCounter_IncidentAngle +LoopCounter_Material_Property & 
      !                       +LoopCounter_Position_Source_X +LoopCounter_Position_Source_Y , &
       !                      Position_Node, Number_Node,  &
        !                     Index_Element_2_Node_Plot, Number_Element_Plot, &
        !                     Class_Element_Plot, ID_Element_Material, ID_Element_FixedDomain, &
         !                    Normalized_Objective_Function&
        !                     ( 1, LoopCounter_Position_Source_Y, LoopCounter_Position_Source_X, & 
         !                      LoopCounter_Material_Property, LoopCounter_IncidentAngle, LoopCounter_Frequency ) ) 
         

         
                        if( Flag_Thermal_Device==0 .or. Flag_Thermal_Device==10 )then 
                           do i= 1, Number_Node
                              Value_Plot_Distribution( i )= Temperature_Solution( i ) -Temperature_Reference( i ) 
                           end do
                        else
                           do i= 1, Number_Node
                              Value_Plot_Distribution( i )= Temperature_Solution( i ) -Temperature_Reference( i )
                           end do
                        end if
   
                        call Plot_DTFD_OC &
                           ( 2, 200 +LoopCounter_Frequency +LoopCounter_IncidentAngle +LoopCounter_Material_Property &
                                  +LoopCounter_Position_Source_X + LoopCounter_Position_Source_Y, &
                             Value_Plot_Distribution, Number_Level_Plot_EM, & 
                             LoopCounter_Frequency +LoopCounter_IncidentAngle +LoopCounter_Material_Property & 
                             +LoopCounter_Position_Source_X +LoopCounter_Position_Source_Y , &
                             Position_Node, Number_Node,  &
                             Index_Element_2_Node_Plot, Number_Element_Plot, &
                             Class_Element_Plot, ID_Element_Material, &
                             ID_Element_Material, ID_Element_FixedDomain, ID_Element_Base_Material, &
                             Number_PV, Position_Arrow_PV, Length_Maximum_Vector, &
                             Normalized_Objective_Function&
                             ( 1, LoopCounter_Position_Source_Y, LoopCounter_Position_Source_X, & 
                              LoopCounter_Material_Property, LoopCounter_IncidentAngle, LoopCounter_Frequency ) ) 
                
           
                        !call Plot_OC_PV &
                        !   ( 1, 200 +LoopCounter_Frequency +LoopCounter_IncidentAngle +LoopCounter_Material_Property &
                        !          +LoopCounter_Position_Source_X + LoopCounter_Position_Source_Y, &
                        !     Value_Plot_Distribution, Number_Level_Plot_EM, & 
                        !     LoopCounter_Frequency +LoopCounter_IncidentAngle +LoopCounter_Material_Property & 
                        !     +LoopCounter_Position_Source_X +LoopCounter_Position_Source_Y , &
                        !     Position_Node, Number_Node,  &
                        !     Index_Element_2_Node_Plot, Number_Element_Plot, &
                        !     Class_Element_Plot, ID_Element_Material, &
                        !     Number_PV, Position_Arrow_PV, Level_PV, Length_Maximum_Vector, &
                        !     Normalized_Objective_Function&
                        !     ( 1, LoopCounter_Position_Source_Y, LoopCounter_Position_Source_X, & 
                        !       LoopCounter_Material_Property, LoopCounter_IncidentAngle, LoopCounter_Frequency ) ) 
    
                        deallocate( Position_Arrow_PV )
                        deallocate( Level_PV ) 
                        deallocate( Index_Element_2_Node_Plot )
                        deallocate( Class_Element_Plot )
   
                     end if

                     deallocate( Thermal_Conductivity_Element ) 
                     deallocate( Temperature_Solution )
      
                     if( Flag_Reference_Device==1 )then
                        deallocate( Temperature_Reference )
                     end if
   
                  end do ! Loop_Position_Source_Y
               end do ! Loop_Position_Source_X
            end do ! Loop_Material_Property
         end do ! Loop_IncidentAngle
      end do ! Loop_Frequency
   end do ! Loop_Core

!   if( Number_Loop_Frequency/=0 .and. Number_Loop_IncidentAngle==0 .and. Number_Loop_Material_Property==0 )then
!      open( 100, file='Result_Frequency.dat', position='append' ) 
!   else if( Number_Loop_Frequency==0 .and. Number_Loop_IncidentAngle/=0 .and. Number_Loop_Material_Property==0 )then
!      open( 100, file='Result_IncidentAngle.dat', position='append' ) 
!   else if( Number_Loop_Frequency==0 .and. Number_Loop_IncidentAngle==0 .and. Number_Loop_Material_Property/=0 )then
!      open( 100, file='Result_Material_Property.dat', position='append' ) 
!   else
      open( 100, file='Result.dat', position='append' ) 
!   end if


   if( Number_Objective_Function==1 )then 
      write(100,*)'# T_left, T_right, Thermal_Conductivity_Material, Objective Function'
      do LoopCounter_Frequency= 0, Number_Loop_Frequency 
         do LoopCounter_IncidentAngle= 0, Number_Loop_IncidentAngle
            do LoopCounter_Material_Property= 0, Number_Loop_Material_Property
               do LoopCounter_Position_Source_X= 0, Number_Loop_Position_Source_X 
                  do LoopCounter_Position_Source_Y= 0, Number_Loop_Position_Source_Y
      
                     write(100, '(2x,es15.6,2x,es15.6,2x,es15.6,2x,es15.6,2x,es15.6,2x,es15.6,2x,es15.6)') & 
                     Position_Source_X( LoopCounter_Position_Source_X ), &
                     Position_Source_Y( LoopCounter_Position_Source_Y ), &
                     Thermal_Conductivity_Material_Real( LoopCounter_Material_Property ), &
                     Normalized_Objective_Function & 
                     ( 1, LoopCounter_Position_Source_Y, LoopCounter_Position_Source_X, & 
                       LoopCounter_Material_Property, LoopCounter_IncidentAngle, LoopCounter_Frequency )
   
                  end do ! LoopCounter_Position_Source_Y
               end do ! LoopCounter_Position_Source_X
            end do ! Loop_Material_Property
         end do ! Loop_IncidentAngle
      end do ! Loop_Frequency

   else if( Number_Objective_Function==2 )then
      write(100,'(a115)')'#       T_left,         T_right, Thermal_Conductivity_Material, Objective Function 1, Objective Function 2 '
      do LoopCounter_Frequency= 0, Number_Loop_Frequency 
         do LoopCounter_IncidentAngle= 0, Number_Loop_IncidentAngle
            do LoopCounter_Material_Property= 0, Number_Loop_Material_Property
               do LoopCounter_Position_Source_X= 0, Number_Loop_Position_Source_X 
                  do LoopCounter_Position_Source_Y= 0, Number_Loop_Position_Source_Y
      
                     write(100, '(2x,es15.6,2x,es15.6,2x,es15.6,2x,es15.6,2x,es15.6,2x,es15.6,2x,es15.6)') & 
                     Position_Source_X( LoopCounter_Position_Source_X ), &
                     Position_Source_Y( LoopCounter_Position_Source_Y ), &
                     Thermal_Conductivity_Material_Real( LoopCounter_Material_Property ), &
                     Normalized_Objective_Function & 
                     ( 1, LoopCounter_Position_Source_Y, LoopCounter_Position_Source_X, & 
                       LoopCounter_Material_Property, LoopCounter_IncidentAngle, LoopCounter_Frequency ), & 
                     Normalized_Objective_Function & 
                     ( 2, LoopCounter_Position_Source_Y, LoopCounter_Position_Source_X, & 
                       LoopCounter_Material_Property, LoopCounter_IncidentAngle, LoopCounter_Frequency ) 
   
                  end do ! LoopCounter_Position_Source_Y
               end do ! LoopCounter_Position_Source_X
            end do ! Loop_Material_Property
         end do ! Loop_IncidentAngle
      end do ! Loop_Frequency
   end if

   close( 100 )

   if( Flag_Structure==1 )then
      !open( 101, file='Result_Normalization.dat', position='append' )
      open( 101, file=trim(File_Name_Result_Normalization), position='append' ) 

      !if( Number_Objective_Function==1 )then
         do LoopCounter_Frequency= 0, Number_Loop_Frequency 
            do LoopCounter_IncidentAngle= 0, Number_Loop_IncidentAngle
               do LoopCounter_Material_Property= 0, Number_Loop_Material_Property
                  do LoopCounter_Position_Source_X= 0, Number_Loop_Position_Source_X 
                     do LoopCounter_Position_Source_Y= 0, Number_Loop_Position_Source_Y
         
                        write( 101, * ) &
                        Original_OF_Value_without_Normalization& 
                        ( 1, LoopCounter_Position_Source_Y, LoopCounter_Position_Source_X, &
                          LoopCounter_Material_Property, LoopCounter_IncidentAngle, LoopCounter_Frequency ) 
      
                     end do ! LoopCounter_Position_Source_Y
                  end do ! LoopCounter_Position_Source_X
               end do ! Loop_Material_Property
            end do ! Loop_IncidentAngle
         end do ! Loop_Frequency
      !end if

         write( 101, * ) 'Number_Loop_Frequency=', Number_Loop_Frequency
         write( 101, * ) 'Number_Loop_IncidentAngle=', Number_Loop_IncidentAngle
         write( 101, * ) 'Number_Loop_Material_Property=', Number_Loop_Material_Property
      close( 101 )

   else if( Flag_Structure==2 .and. Number_Objective_Function==2 .and. &
            ( Flag_Thermal_Device==20 .or. Flag_Thermal_Device==21 .or. Flag_Thermal_Device==30 ) )then

      open( 101, file=trim(File_Name_Result_Normalization), position='append' ) 

         do LoopCounter_Frequency= 0, Number_Loop_Frequency 
            do LoopCounter_IncidentAngle= 0, Number_Loop_IncidentAngle
               do LoopCounter_Material_Property= 0, Number_Loop_Material_Property
                  do LoopCounter_Position_Source_X= 0, Number_Loop_Position_Source_X 
                     do LoopCounter_Position_Source_Y= 0, Number_Loop_Position_Source_Y
         
                        write( 101, * ) &
                        Original_OF_Value_without_Normalization& 
                        ( 2, LoopCounter_Position_Source_Y, LoopCounter_Position_Source_X, &
                          LoopCounter_Material_Property, LoopCounter_IncidentAngle, LoopCounter_Frequency ) 
      
                     end do ! LoopCounter_Position_Source_Y
                  end do ! LoopCounter_Position_Source_X
               end do ! Loop_Material_Property
            end do ! Loop_IncidentAngle
         end do ! Loop_Frequency

         write( 101, * ) 'Number_Loop_Frequency=', Number_Loop_Frequency
         write( 101, * ) 'Number_Loop_IncidentAngle=', Number_Loop_IncidentAngle
         write( 101, * ) 'Number_Loop_Material_Property=', Number_Loop_Material_Property
      close( 101 )

   end if

   deallocate( Value_Plot_Distribution )
   deallocate( Value_Objective_Function )

   deallocate( LoopCounter_Frequency_Range )
   deallocate( LoopCounter_IncidentAngle_Range )
   deallocate( LoopCounter_Material_Property_Range )
   deallocate( LoopCounter_Position_Source_X_Range )
   deallocate( LoopCounter_Position_Source_Y_Range )

   deallocate( Temperature_Left_Side_Boundary )
   deallocate( Temperature_Right_Side_Boundary )
   deallocate( Thermal_Conductivity_Material_Real )
   deallocate( Thermal_Conductivity_Material )
   deallocate( Thermal_Conductivity_FixedDomain_Variable )
   deallocate( Position_Source_X )
   deallocate( Position_Source_Y )

   deallocate( Normalized_Objective_Function )
   deallocate( Original_OF_Value_without_Normalization )
   deallocate( Value_OF_Normalize )

   deallocate( Position_Node )
   deallocate( Max_Position_Node )
   deallocate( Index_Element_2_Node_Triangle )

   deallocate( Class_Element_Triangle )

   deallocate( Element_and_LocalNode_Number_on_PEC_BC )

   write(*,*)'============================================================='
   write(*,*)'============================================================='
   write(*,*)'=====      ===   ===   ===   ===   ====     ===   ===   ====='
   write(*,*)'=====   ======   ===    ==   ===   ===    =====   ===   ====='
   write(*,*)'=====      ===   ===         ===   ====    ====         ====='
   write(*,*)'=====   ======   ===   ==    ===   =====    ===   ===   ====='
   write(*,*)'=====   ======   ===   ===   ===   ===     ====   ===   ====='
   write(*,*)'============================================================='
   write(*,*)'============================================================='
 
end program Main
 


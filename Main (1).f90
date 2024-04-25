
program Main

   !$use omp_lib
   use Parameters
   !use Cmaes, only : CmaParam !, Cmaes_Manager, Cmaes_Prestep
   use structure

   implicit none
   include 'mpif.h'

   !type(CmaParam) :: cma                   !CHANGE!
   type(FEM_Data) :: femdata                   !CHANGE!
   type(Temporary_Data) :: tmp                  !CHANGE!
   integer :: e, i, j, k, Counter, No
   integer :: Device_Number 
   integer :: NumThread_MKL, NumThread_OMP, NumThread_MPI 
   integer :: Loop_Solution, Loop_Start 
   integer :: Loop_Source, Loop_MatParam_All, Loop_MatParam_T, Loop_MatParam_V, Loop_Obj_Func 
   integer :: Optimization_Step 
   integer :: info
   character(len=256) :: Filename_Log
   character(len=256) :: Filename_ObjectiveFunction
   character(len=256) :: Filename_VolumeConstraint
   character(len=256) :: Filename_FEM_Data
   character(len=256) :: Filename_Time_Step
   character(len=256) :: Filename_Coefficient_Complexity
   character(len=256) :: Filename_Structure
   character(len=256) :: FileName_LSF 

   integer :: Number_Sampling
   integer :: Number_Obj_Func_MaterialParameter
   double precision, allocatable, dimension(:,:,:) :: Value_Obj_Func_Normalize_Thermal, Value_Obj_Func_Normalize_DC 

   double precision, allocatable, dimension(:) :: Position_X_Source, Position_Y_Source
   double precision, allocatable, dimension(:) :: Relative_Thermal_Conductivity, Relative_Electrical_Conductivity 
   double precision :: Relative_Thermal_Conductivity_Evaluated, Relative_Electrical_Conductivity_Evaluated 
  
   double precision, allocatable, dimension(:) :: Fitness_CMAES, Fitness_Convergence
   double precision, allocatable, dimension(:) :: Convergence_Ratio, Convergence_Ratio_Average, Convergence_Ratio_xmean
   double precision, allocatable, dimension(:,:,:,:,:) :: Obj_Func_Multi
   integer :: Number_Physics, Loop_Physics
   double precision, allocatable, dimension(:) :: Obj_Func_Real
   double precision Obj_Func_tmp

   double precision :: Convergence_Ratio_xmean_dowhile, Convergence_Error
   character(len=1) :: Termination_dowhile
   character(len=1), allocatable, dimension(:) :: Termination

   integer, allocatable, dimension(:) :: Individual_Optimal

   ! Grid Point
   integer :: Number_Quadrant
   integer :: Number_GridPoint_Optimization, Number_Design_Variable, n_dv
   integer, allocatable, dimension(:,:) :: Index_DV_Symmetric_2_GP_All 
   double precision, allocatable, dimension(:,:) :: LSF_GridPoint
   double precision, allocatable, dimension(:) :: LSF_GridPoint_tmp 
   integer, allocatable, dimension(:) :: Class_GridPoint_Scattering
   integer, allocatable, dimension(:) :: Correspondence_Number 
   double precision, allocatable, dimension(:,:) :: Position_GridPoint_Scattering
   double precision, allocatable, dimension(:) :: LSF_FixedDomain_GridPoint
   double precision, allocatable, dimension(:) :: LSF_DesignDomain_GridPoint
   double precision, allocatable, dimension(:) :: LSF_ExteriorDomain_GridPoint
   double precision, allocatable, dimension(:,:) :: Design_Variable
   double precision, allocatable, dimension(:,:) :: Design_Variable_Blackbox
   double precision, allocatable, dimension(:) :: Design_Variable_Output
   double precision, allocatable, dimension(:,:) :: LSF_Piecewise_Constant
   double precision, allocatable, dimension(:,:) :: Diff_LSF_PC
   double precision, allocatable, dimension(:) :: Average_Diff_LSF_PC
 
   integer, allocatable, dimension(:,:) :: Index_Grid_2_GridPoint
   integer, allocatable, dimension(:,:) :: Index_Grid_total
   integer, allocatable, dimension(:,:) :: Index_Grid_2_GridPoint_twice
   integer, allocatable, dimension(:,:) :: Index_Grid_2_GridPoint_3
   integer, allocatable, dimension(:,:) :: Index_Grid_2_GridPoint_4
   integer, allocatable, dimension(:,:) :: Index_Grid_2_GridPoint_10
   integer, allocatable, dimension(:,:) :: Index_Grid_2_GridPoint_15

   !! Grid Point in PML 
   !integer, allocatable, dimension(:) :: Class_GridPoint_PML
   !double precision, allocatable, dimension(:,:) :: Position_GridPoint_PML
   !integer, allocatable, dimension(:,:) :: Index_Grid_PML_2_GridPoint_PML  

   !! Mesh in PML
   !double precision, allocatable, dimension(:,:,:)  :: Position_Node_OneGrid_PML 
   !integer, allocatable, dimension(:) :: Number_Node_OneGrid_PML
   !integer, allocatable, dimension(:) :: Class_Grid_PML 

   !integer, allocatable, dimension(:) :: Number_Element_OneGrid_PML
   !integer, allocatable, dimension(:,:) :: Class_Element_OneGrid_PML
   !integer, allocatable, dimension(:,:,:) :: Index_Element_PML_2_Node_PML

   ! Node, Element
   double precision, allocatable, dimension(:,:) :: Position_Node_Preserve
 
   integer, allocatable, dimension(:) :: Class_Node_Preserve
   integer, allocatable, dimension(:,:) :: Index_Element_2_Node_Triangle_Preserve 
   integer, allocatable, dimension(:) :: Class_Element_Triangle_Preserve 

   integer, allocatable, dimension(:) :: Index_GridPointNumber_2_NodeNumber, Index_GridPointNumber_2_NodeNumber_Renumbered

   double precision, allocatable, dimension(:) :: Perimeter_Structure_All, Perimeter_Implicit_All
   double precision, allocatable, dimension(:) :: Volume_Structure_All
   double precision :: Perimeter_Structure
   double precision :: Perimeter_Structure_Minimum, Obj_Func_Real_Minimum, Volume_Structure_Minimum

   double precision :: Fractal_Structure
   double precision, allocatable, dimension(:) :: Fractal_MPI, Fractal_Structure_All
   double precision,parameter :: Perimeter_i=0.0d0 ,Perimeter_a=7.0d0
   double precision,parameter :: Fractal_i=0.0d0 ,Fractal_a=1.0d0
   double precision,parameter :: Obj_Func_i=0.0d0 ,Obj_Func_a=1.0d-3
   double precision ::Fractal_trade,Obj_Func_trade,Perimeter_trade

   integer :: Number_Node_on_Electrical_Insulation_Boundary, Number_Node_in_Electrical_Insulation

   !Create Input Data
   integer :: Max_Number_Element_Share_OneNode
   integer :: Width_Matrix_LHS 

   double precision, allocatable, dimension(:,:,:) :: Temperature_Solution 
   double precision, allocatable, dimension(:) :: Temperature_tmp 
   double precision, allocatable, dimension(:,:,:) :: Electric_Potential_Solution 
   double precision, allocatable, dimension(:) :: Electric_Potential_tmp 
   double precision, allocatable, dimension(:,:,:) :: Field_Plot_1, Field_Plot_2 
   integer :: Flag_Physics_Plot 

   integer :: Number_Node_Reference 
   double precision, allocatable, dimension(:,:,:) :: Temperature_Reference
   double precision, allocatable, dimension(:,:,:) :: Electric_Potential_Reference

   double precision, allocatable, dimension(:,:) :: Max_Position_Node, Max_Position_Node_Preserve

   integer :: Number_Element_Electrical_Insulation_Boundary
   integer, allocatable, dimension(:,:) :: Element_and_LocalNode_Number_on_Electrical_Insulation_BC
   integer, allocatable, dimension(:,:) :: Element_and_LocalNode_Num_on_Electrical_Insulation_BC_Preserve

   ! Plot Result
   double precision, allocatable, dimension(:) :: Value_Plot
   double precision, allocatable, dimension(:,:,:) :: Value_Plot_Preserve_Thermal, Value_Plot_Preserve_DC
   integer :: Interval_Plot_Optimization_Step
   integer :: Generation_Convergence_xmean 
   double precision :: Threshold_Sigma_xmean
   double precision, allocatable, dimension(:) :: Threshold_Convergence, Threshold_Convergence_Average, Threshold_Convergence_xmean
   double precision :: Threshold_Convergence_tmp, Threshold_Convergence_xmean_read
   integer :: Flag_Symmetry_LSF_read
   integer :: Number_Threshold_Convergence, Number_Threshold_Convergence_Average, Number_Threshold_Convergence_xmean

   ! Objective Function
   double precision, allocatable, dimension(:) :: Obj_Func_All 
   double precision, allocatable, dimension(:,:,:,:) :: Obj_Func_Single_Thermal, Obj_Func_Single_DC
   !complex( kind( 0d0 ) ), allocatable, dimension(:) :: Value_Obj_Func
   double precision, allocatable, dimension(:) :: Value_Obj_Func
   integer :: ID_Element_Obj_Func
   ! Output FEM Data
   double precision :: Fitness_Output 
   double precision :: Obj_Func_Output  
   double precision :: Perimeter_Output
   integer :: Generation_Output, Generation_Output_Read 
   character(len=Length_Character_Optimization_Step) :: Generation_Character
   character(len=8) :: Format_Filenumber
   integer :: Number_Node_Preserve 
   integer :: Number_Element_Triangle_Preserve  
   integer :: Number_Node_on_Electrical_Insulation_Boundary_Preserve, Number_Node_in_Electrical_Insulation_Preserve
   integer :: Number_Element_Electrical_Insulation_Boundary_Preserve 
   integer :: Max_Number_Element_Share_OneNode_Preserve

   integer :: Number_Node_Output
   integer :: Number_Element_Triangle_Output 
   integer :: Number_Node_on_Electrical_Insulation_Boundary_Output, Number_Node_in_Electrical_Insulation_Output
   integer :: Number_Element_Electrical_Insulation_Boundary_Output 
   integer :: Max_Number_Element_Share_OneNode_Output

   double precision, allocatable, dimension(:,:) :: Position_Node_Output
   integer, allocatable, dimension(:) :: Class_Node_Output
   integer, allocatable, dimension(:,:) :: Index_Element_2_Node_Triangle_Output 
   integer, allocatable, dimension(:) :: Class_Element_Triangle_Output 
   double precision, allocatable, dimension(:,:) :: Max_Position_Node_Output
   integer, allocatable, dimension(:,:) :: Element_and_LocalNode_Number_on_Electrical_Insulation_BC_Output
   double precision, allocatable, dimension(:,:,:) :: Value_Plot_Output_Thermal, Value_Plot_Output_DC

   integer :: Flag_Output_PSFile, Flag_Output_EFD, Flag_Output_Interval_LSF 

   double precision, allocatable, dimension(:,:) :: Position_Node_Plot 
   integer, allocatable, dimension(:,:) :: Index_Element_2_Node_Plot 
   integer, allocatable, dimension(:) :: Class_Element_Plot 
   !Optimization

   ! Volume Constraint
   double precision :: Volume_Material, Volume_DesignDomain

   !CPU Time
   double precision, allocatable, dimension(:) :: Time_CPU
   !double precision :: Obj_Func_CPU_Threshold
   double precision, allocatable, dimension(:) :: Fitness_CPU_Threshold, Obj_Func_CPU_Threshold
   integer :: Number_Fitness_CPU_Threshold, Number_Obj_Func_CPU_Threshold
   double precision :: Time_CPU_Start
   integer, allocatable, dimension(:) :: Generation_CPU
 
   !MPI
   integer :: ierr, MyRank, Nproc, Number_Loop_MPI, MyRank_Output_tmp, MyRank_Output, Individual_MPI
   integer, allocatable, dimension(:) :: Loop_Num_MPI
   double precision, allocatable, dimension(:) :: Perimeter_MPI, Volume_MPI, Obj_Func_MPI, Fitness_MPI 
   double precision, allocatable, dimension(:) :: Perimeter_Implicit_MPI, Average_Diff_LSF_PC_MPI
   double precision, allocatable, dimension(:) :: Penalty_Volume_Constraint_MPI, Penalty_Volume_Constraint 
   double precision, allocatable, dimension(:) :: Vector_Data_MPI, Vector_Data
   !complex( kind( 0d0 ) ), allocatable, dimension(:) :: Vector_Data_Complex_MPI
   double precision, allocatable, dimension(:) :: Vector_Data_Real_MPI

   ! Sparse CMA-ES
   integer, allocatable, dimension(:,:) :: Neighboring_LSF
   integer, allocatable, dimension(:,:) :: dep_dv 
   integer, allocatable, dimension(:) :: Index_GP_All_2_DV_Symmetric

   double precision :: sigma, flag_eigen

   integer :: Number_Point_ATD( 2 )
   character(len=1) :: match
   integer :: Number_Grid_quarter,Number_Grid_16,Number_Grid_9
   integer :: Number_Grid_100,Number_Grid_225,Grid_total
   integer,parameter :: Box_Number=6
   integer,allocatable,dimension(:)::plot_Number
   !================================================================
   ! MPI Setting
   !================================================================
   CALL MPI_Init(ierr)
   CALL MPI_Comm_size(MPI_COMM_WORLD, Nproc, ierr)
   CALL MPI_Comm_rank(MPI_COMM_WORLD, MyRank, ierr)
   call MPI_Barrier( MPI_COMM_WORLD, ierr )

   IF(MyRank.EQ.0) WRITE(*,*) ' Nproc = ',Nproc
   WRITE(*,*) 'Hello world ! from processor ', MyRank

   open(302, file='../src/Device_Number')
      read(302,*)Device_Number
      write(*,*)'Device_Number=', Device_Number
   close(302)

  if( Device_Number==999 )then
      Number_Sampling= 1
   else
      Number_Sampling= Number_Candidates
      if( mod( Number_Sampling, Nproc )/=0 )then
         write(*,*)'Number_Sampling=', Number_Sampling
         write(*,*)'Nproc=', Nproc
         call Output_Error( 'Main', 179 )
      end if
   end if

   !if( mod( Number_Grid_Scale_DesignDomain, 3 )/=0 )then
   !   write(*,*)'Number_Grid_Scale_DesignDomain=', Number_Grid_Scale_DesignDomain
   !   write(*,*)'mod( Number_Grid_Scale_DesignDomain, 3 )=', mod( Number_Grid_Scale_DesignDomain, 3 )
   !   call Output_Error( 'Main', 196 ) 
   !end if

   if( MyRank==0 )then
      open(178, file='MPI_Info.dat', position='append')
         write(178,*)'Nproc=', Nproc 
         write(178,*)'MyRank=', MyRank 
      close(178)
   end if

   !================================================================
   ! OpenMP Setting
   !================================================================
   open(9999, file='Number_Thread')
   read(9999,*)NumThread_OMP, NumThread_MKL, NumThread_MPI
   if( MyRank==0 ) write(*,*)'NumThread_OMP=', NumThread_OMP
   if( MyRank==0 ) write(*,*)'NumThread_MKL=', NumThread_MKL
   if( MyRank==0 ) write(*,*)'NumThread_MPI=', NumThread_MPI
   !$ call omp_set_num_threads( NumThread_OMP )
   !!$ call mkl_set_num_threads( NumThread_MKL )
   close(9999)

   open(229, file='Convergence_Error')
      read(229,*) Convergence_Error 
      if( MyRank==0 ) write(*,*)'Convergence_Error=', Convergence_Error 
   close(229)

   !================================================================
   if( MyRank==0 ) write(*,*)'====================================================================='
   if( MyRank==0 ) write(*,*)'Topology Optimization Program based on Level Set Method' 
   if( MyRank==0 ) write(*,*)'for 2D Dielectric Structure'
   if( MyRank==0 ) write(*,*)' Finite Element Method ( Node Base Element: TM mode )'
   if( MyRank==0 ) write(*,*)' PML Boundary Condition'
   if( MyRank==0 ) write(*,*)' Multifrontal Method (ZMUMPS,DMUMPS)'
   if( MyRank==0 ) write(*,*)'  by Garuda Fujii'
   if( MyRank==0 ) write(*,*)'   Ver.1: 13 June   2012'
   if( MyRank==0 ) write(*,*)'   Ver.2:  6 February 2013'
   if( MyRank==0 ) write(*,*)'====================================================================='
   !================================================================
   
   if( MyRank==0 ) write(*,*)'                       ' 
   if( MyRank==0 ) write(*,*)'Number_Optimization_Step=', Number_Optimization_Step 
   if( MyRank==0 ) write(*,*)'Size_X_Scattering=', Size_X_Scattering
   if( MyRank==0 ) write(*,*)'Size_Y_Scattering=', Size_Y_Scattering 
   if( MyRank==0 ) write(*,*)'Radius_FixedDomain', Radius_FixedDomain
   if( MyRank==0 ) write(*,*)'Radius_DesignDomain',Radius_DesignDomain
   if( MyRank==0 ) write(*,*)'====================================='

   !========================================================================
   if( MyRank==0 ) write(*,*) ' '
   if( MyRank==0 ) write(*,*) '<<<<<<<<<<<<<<< Pre-Process Starts >>>>>>>>>>>>>>>'
   if( MyRank==0 ) write(*,*) ' '
   !========================================================================
  
   call Check_Parameters( info )
  
   if( info/=0 )then
      write(*,*)'ERROR: info/=0 at Check_Parameters'
      write(*,*)'info=', info
      stop
   end if
   
   call Define_Filename &
      ( Filename_Log, &
        Filename_ObjectiveFunction, &
        Filename_VolumeConstraint, Filename_FEM_Data, Filename_Time_Step, Filename_Coefficient_Complexity, &
        Filename_Structure )

   call Output_Parameters( Filename_Log )
 
   !======================================================
   if( MyRank==0 ) write(*,*)'   Parameter Setting'
   !======================================================

   if( Flag_Physics==1 .or. Flag_Physics==2  )then
      Number_Physics= 1
   else if( Flag_Physics==0 )then
      Number_Physics= 2
   else
      write(*,*)'Flag_Physics=', Flag_Physics
      call Output_Error( 'Main', 266 ) 
   end if

   !=========================================================================
   if( MyRank==0 ) write(*,*)'   Parameter Thresholds for Convergence'
   !=========================================================================

   Number_Threshold_Convergence= 4
   Number_Threshold_Convergence_Average= Number_Threshold_Convergence
   Number_Fitness_CPU_Threshold= 4
   Number_Obj_Func_CPU_Threshold= 4
   Number_Threshold_Convergence_xmean= 4
 
   !=========================================================================
   if( MyRank==0 ) write(*,*)'   Structural Symmetry'
   !=========================================================================

   if( Flag_Symmetry_LSF==0 )then
      Number_Quadrant= 1 
   else if( Flag_Symmetry_LSF==1 .or. Flag_Symmetry_LSF==2 .or. Flag_Symmetry_LSF==12 )then
      Number_Quadrant= 2 
   else if( Flag_Symmetry_LSF==4 .or. Flag_Symmetry_LSF==14 )then
      Number_Quadrant= 4 
   else if( Flag_Symmetry_LSF==8 )then
      Number_Quadrant= 8
   else
      write(*,*)'Flag_Symmetry_LSF=', Flag_Symmetry_LSF
      call Output_Error( 'Main', 221 ) 
   end if

   !=========================================================================
   if( MyRank==0 ) write(*,*)'   Relative Relations of Constraints'
   !=========================================================================

   !if( Flag_Volume_Constraint==1 .and. Coefficient_Complexity_Tau_Normal >= 1d-2 )then
   !   write(*,*)'Coefficient_Complexity_Tau_Normal=', Coefficient_Complexity_Tau_Normal
   !   write(*,*)'Flag_Volume_Constraint=', Flag_Volume_Constraint
   !   call Output_Error( 'Main', 313 ) 
   !end if

   !=========================================================================
   if( MyRank==0 ) write(*,*)'   Optimization for Min. or Max'
   !=========================================================================

   if( Flag_Ordering )then
      if( Sign_ObjectiveFunction /= 1d0 )then 
         write(*,*)'Flag_Ordering=', Flag_Ordering
         write(*,*)'Sign_ObjectiveFunction=', Sign_ObjectiveFunction
         call Output_Error( 'Main', 214 ) 
      end if
   else
      if( Sign_ObjectiveFunction /=- 1d0 )then 
         write(*,*)'Flag_Ordering=', Flag_Ordering
         write(*,*)'Sign_ObjectiveFunction=', Sign_ObjectiveFunction
         call Output_Error( 'Main', 222 ) 
      end if
   end if

   !======================================================
   if( MyRank==0 ) write(*,*)'   Setting -- Optical Device'
   !======================================================
  
   allocate( Class_GridPoint_Scattering( Number_GridPoint_Scattering ) )
   allocate( Position_GridPoint_Scattering( 2, Number_GridPoint_Scattering ) )
   allocate( LSF_FixedDomain_GridPoint( Number_GridPoint_Scattering ) )
   allocate( LSF_DesignDomain_GridPoint( Number_GridPoint_Scattering ) )
   allocate( LSF_ExteriorDomain_GridPoint( Number_GridPoint_Scattering ) )
   allocate( tmp%itg1( Number_GridPoint_Scattering ) )
   allocate( tmp%itg2( Number_Quadrant, Number_GridPoint_Scattering ) )

   call Select_Laplacian_Device & 
      ( Number_GridPoint_Scattering, Number_Grid_X_Scattering, Number_Grid_Y_Scattering, &
       !================================================================================
        Position_GridPoint_Scattering, Class_GridPoint_Scattering, &
        Number_Quadrant, Number_GridPoint_Optimization, tmp%itg2, &
        LSF_FixedDomain_GridPoint, LSF_DesignDomain_GridPoint, LSF_ExteriorDomain_GridPoint, tmp%itg1 )

   allocate( Index_DV_Symmetric_2_GP_All( Number_Quadrant, Number_GridPoint_Optimization ) )
   allocate( Correspondence_Number( Number_GridPoint_Optimization ) )

   if( MyRank==0 ) write(*,*)'   Class_GridPoint_Scattering( 100 )=', Class_GridPoint_Scattering( 100 ) 
   if( MyRank==0 ) write(*,*)'   Position_GridPoint_Scattering( 2, 100 )=', Position_GridPoint_Scattering( 2, 100 ) 
   if( MyRank==0 ) write(*,*)'   Position_Center_DesignDomain_Y=', Position_Center_DesignDomain_Y 
   if( MyRank==0 ) write(*,*)'   Number_GridPoint_Optimization=', Number_GridPoint_Optimization
   Index_DV_Symmetric_2_GP_All( :, 1:Number_GridPoint_Optimization )= tmp%itg2( :, 1:Number_GridPoint_Optimization )
   Correspondence_Number( 1:Number_GridPoint_Optimization )= tmp%itg1( 1:Number_GridPoint_Optimization )

   do i = 1, Number_GridPoint_Optimization
      do j= 1, Number_Quadrant 
         if( 1 <= Index_DV_Symmetric_2_GP_All( j, i ) .and. &
           Index_DV_Symmetric_2_GP_All( j, i ) <= Number_GridPoint_Scattering  )then
           284 continue
         else
         write( *, * )'===================================================='
         write( *, * )'Index_DV_Symmetric_2_GP_All( j, i )=', Index_DV_Symmetric_2_GP_All( j, i )
         write( *, * )'ERROR'
         write( *, * )'Main', '.f90', 289
         write( *, * )'===================================================='
         stop
         end if
      end do
   end do

   allocate( Index_GP_All_2_DV_Symmetric( Number_GridPoint_Scattering ) )

   Index_GP_All_2_DV_Symmetric = 0
   do i= 1, Number_GridPoint_Optimization
      do j= 1, Number_Quadrant
         Index_GP_All_2_DV_Symmetric( Index_DV_Symmetric_2_GP_All( j, i ) ) = i
      end do
   end do

   deallocate( tmp%itg2 )
   deallocate( tmp%itg1 )

   !======================================================
   !if( MyRank==0 ) write(*,*)'   Classify Grid Points in PML'
   !======================================================
  
   !allocate( Position_GridPoint_PML( 2, Number_GridPoint_PML+(Width_PML+1)*4 ) )
   !allocate( Class_GridPoint_PML( Number_GridPoint_PML+(Width_PML+1)*4 ) )
   !allocate( Class_Grid_PML( Number_Grid_PML ) )
   !allocate( Index_Grid_PML_2_GridPoint_PML( 4, Number_Grid_PML ) )
  
   !call Classify_GridPoint_PML( Position_GridPoint_PML, Class_GridPoint_PML, Class_Grid_PML, Index_Grid_PML_2_GridPoint_PML )
  
   !======================================================
   !if( MyRank==0 ) write(*,*)'   Generate Mesh in PML'
   !======================================================
  
   !allocate( Position_Node_OneGrid_PML( 2, 4, Number_Grid_PML ) )
   !allocate( Number_Node_OneGrid_PML( Number_Grid_PML ) )
   !allocate( Number_Element_OneGrid_PML( Number_Grid_PML ) )
   !allocate( Class_Element_OneGrid_PML( Max_Number_Element_OneGrid_PML, Number_Grid_PML ) )
   !allocate( Index_Element_PML_2_Node_PML( 4, Max_Number_Element_OneGrid_PML, Number_Grid_PML ) )
  
   !call Meshing_Grid_PML &
   !   ( Position_GridPoint_PML, Class_Grid_PML, Index_Grid_PML_2_GridPoint_PML, &
   !    !=========================================================================
   !     Position_Node_OneGrid_PML, Number_Node_OneGrid_PML, Number_Element_OneGrid_PML, & 
   !     Class_Element_OneGrid_PML, Index_Element_PML_2_Node_PML )

   !======================================================
   if( MyRank==0 ) write(*,*)'   Initialize Leve Set Function'
   !======================================================

   Number_Design_Variable= Number_GridPoint_Optimization*Number_Type_LSF
   n_dv = Number_Design_Variable

   allocate( Design_Variable( Number_Design_Variable, Number_Sampling ) )
   allocate( Design_Variable_Blackbox( Number_Design_Variable, Number_Sampling) ) 
   allocate( Design_Variable_Output( Number_Design_Variable ) )
   allocate( LSF_Piecewise_Constant( Number_Design_Variable, Number_Sampling ) ) 
   allocate( Diff_LSF_PC( Number_Design_Variable, Number_Sampling) ) 

   if( MyRank==0 ) write(*,*)'============================================' 
   if( MyRank==0 ) write(*,*)'LSF_Minimum          =',    LSF_Minimum
   if( MyRank==0 ) write(*,*)'LSF_Maximum          =',    LSF_Maximum
   if( MyRank==0 ) write(*,*)'Number_GridPoint_Scattering  =', Number_GridPoint_Scattering
   if( MyRank==0 ) write(*,*)'Number_Design_Variable=', Number_Design_Variable
   if( MyRank==0 ) write(*,*)'Number_GridPoint_Optimization=', Number_GridPoint_Optimization
   if( MyRank==0 ) write(*,*)'Number_Type_LSF=', Number_Type_LSF 
   if( MyRank==0 ) write(*,*)'Number_GridPoint_PML       =', Number_GridPoint_PML
   if( MyRank==0 ) write(*,*)'============================================' 

   open(302, file='../src/Device_Number')
      read(302,*)Device_Number
      write(*,*)'Device_Number=', Device_Number
   close(302)

   if( Device_Number==999 )then
      Number_Sampling= 1
   else
      Number_Sampling= Number_Candidates
      if( mod( Number_Sampling, Nproc )/=0 )then
         write(*,*)'Number_Sampling=', Number_Sampling
         write(*,*)'Nproc=', Nproc
         call Output_Error( 'Main', 179 )
      end if
   end if

   if( Number_Sampling < 4 +3*int( log( dble( Number_Design_Variable ) ) ) )then
      write(*,*)'4 +3*log( Number_Design_Variable )=', int( 4d0 +3d0*log( dble( Number_Design_Variable ) ) ) +1 
      write(*,*)'Number_Sampling', Number_Sampling

      if( Device_Number/=999 .and. Flag_InitialConfiguration/=0 )then
         call Output_Error( 'Main', 382 ) ! aho 
      end if
   end if

   !========================================================================
   if( MyRank==0 ) write(*,*) ' '
   if( MyRank==0 ) write(*,*) '<<<<<<<<<<<<<<< Optimization Process Starts >>>>>>>>>>>>>>>'
   if( MyRank==0 ) write(*,*) ' '
   !========================================================================

   allocate( Obj_Func_All( 0:Number_Optimization_Step ) )
   allocate( Obj_Func_Single_Thermal( 0:Number_Optimization_Step, Number_Obj_Func_Source, Number_Obj_Func_MP_T, Number_Obj_Func ) )
   allocate( Obj_Func_Single_DC( 0:Number_Optimization_Step, Number_Obj_Func_Source, Number_Obj_Func_MP_V, Number_Obj_Func ) )

   Obj_Func_All= 0d0
   Obj_Func_Single_Thermal= 0d0
   Obj_Func_Single_DC= 0d0

   !Number_Grid_quarter=Number_Grid_Scattering/4
   !Number_Grid_16=Number_Grid_Scattering/16
   !Number_Grid_9=Number_Grid_Scattering/9
   !Number_Grid_100=Number_Grid_Scattering/100
   !Number_Grid_225=Number_Grid_Scattering/225
   Grid_total=0
   do i=1,30
        Grid_total=Grid_total+Number_Grid_Scattering/(i**2)
   enddo
   
   allocate( Index_Grid_2_GridPoint( 4, Number_Grid_Scattering ) )
   allocate( Index_Grid_total( 961, Grid_total))
   allocate( plot_Number(Box_Number))
   !allocate( Index_Grid_2_GridPoint_twice( 9, Number_Grid_quarter ) )
   !allocate( Index_Grid_2_GridPoint_3( 16, Number_Grid_9 ) )
   !allocate( Index_Grid_2_GridPoint_4( 25, Number_Grid_16 ) )
   !allocate( Index_Grid_2_GridPoint_10( 121, Number_Grid_100 ) )
   !allocate( Index_Grid_2_GridPoint_15( 256, Number_Grid_225 ) )
   !allocate( Index_rid_2_GridPoint_15times( 4, Number_Grid_Scattering ) )

   call Create_Index_Grid_2_GridPoint( Index_Grid_2_GridPoint )
   !call Grid_twice(Number_Grid_quarter,Number_Grid_9,Number_Grid_16,Number_Grid_100,Number_Grid_225,Size_X_Scattering,&
   !        Index_Grid_2_GridPoint_twice,Index_Grid_2_GridPoint_3,Index_Grid_2_GridPoint_4,Index_Grid_2_GridPoint_10,Index_Grid_2_GridPoint_15)
   call Create_Box(Grid_total,Size_X_Scattering,Size_Y_Scattering,Box_Number,Number_Grid_Scattering,Number_GridPoint_Scattering,&
                        plot_Number,Index_Grid_total)
   !call Grid_4times(Number_Grid_16,Size_X_Scattering,Index_Grid2_GridPoint_4times)
   allocate( Position_Node_Preserve( 2, 3*Number_GridPoint_Scattering ) )
   allocate( Index_Element_2_Node_Triangle_Preserve( 3, 5*Number_GridPoint_Scattering ) )
   allocate( Class_Element_Triangle_Preserve( 5*Number_GridPoint_Scattering ) )
   allocate( Value_Plot_Preserve_Thermal( 3*Number_GridPoint_Scattering, Number_Obj_Func_Source, Number_Obj_Func_MP_T ) )
   allocate( Value_Plot_Preserve_DC( 3*Number_GridPoint_Scattering, Number_Obj_Func_Source, Number_Obj_Func_MP_V ) )

   allocate( Class_Node_Preserve( 3*Number_GridPoint_Scattering ) )
   allocate( Max_Position_Node_Preserve( 2, 2 ) )
   allocate( Element_and_LocalNode_Num_on_Electrical_Insulation_BC_Preserve( 3, 5*Number_GridPoint_Scattering ) )
 
   allocate( Position_Node_Output( 2, 3*Number_GridPoint_Scattering ) )
   allocate( Index_Element_2_Node_Triangle_Output( 3, 5*Number_GridPoint_Scattering ) )
   allocate( Class_Element_Triangle_Output( 5*Number_GridPoint_Scattering ) )
   allocate( Value_Plot_Output_Thermal( 3*Number_GridPoint_Scattering, Number_Obj_Func_Source, Number_Obj_Func_MP_T ) )
   allocate( Value_Plot_Output_DC( 3*Number_GridPoint_Scattering, Number_Obj_Func_Source, Number_Obj_Func_MP_V ) )
   allocate( Class_Node_Output( 3*Number_GridPoint_Scattering ) )
   allocate( Max_Position_Node_Output( 2, 2 ) )
   allocate( Element_and_LocalNode_Number_on_Electrical_Insulation_BC_Output( 3, 5*Number_GridPoint_Scattering ) )

   !==========================================================================================================================
   ! Position_X_Source, Position_Y_Source
   !==========================================================================================================================
   allocate( Position_X_Source( Number_Obj_Func_Source ) )
   allocate( Position_Y_Source( Number_Obj_Func_Source ) )

   if( Number_Obj_Func_Source==1 )then
      Position_X_Source( 1 )= Position_Source_x
      Position_Y_Source( 1 )= Position_Source_y
   else if( Number_Obj_Func_Source==2 )then
      Position_X_Source( 1 ) = Position_Source_x
      Position_X_Source( 2 ) = Position_Source_x
      Position_Y_Source( 1 ) = -Position_Source_y
      Position_Y_Source( 2 ) = Position_Source_y
   else if( Number_Obj_Func_Source==3 )then
      Position_X_Source( 1 ) = Position_Source_x
      Position_X_Source( 2 ) = Position_Source_x
      Position_X_Source( 3 ) = 0.0d0
      Position_Y_Source( 1 ) = -Position_Source_y
      Position_Y_Source( 2 ) = Position_Source_y
      Position_Y_Source( 3 ) = Position_Source_y
   else
      do Loop_Source= 1, Number_Obj_Func_Source
         Position_X_Source( Loop_Source ) &
         = Position_Source_x_1 &
          +( Position_Source_x_2 - Position_Source_x_1 )/( Number_Obj_Func_Source-1 )*( Loop_Source -1 )
      end do ! Loop_Source
      do Loop_Source= 1, Number_Obj_Func_Source
         Position_Y_Source( Loop_Source) &
         = Position_Source_y_1 &
          +( Position_Source_y_2 -Position_Source_y_1 )/( Number_Obj_Func_Source -1 )*( Loop_Source -1 )
      end do ! Loop_Source 
   end if

   !==========================================================================================================================
   ! Relative_Thermal_Conductivity, Relative_Electrical_Conductivity
   !==========================================================================================================================
   allocate( Relative_Thermal_Conductivity( Number_Obj_Func_MP_T ) )
   allocate( Relative_Electrical_Conductivity( Number_Obj_Func_MP_V ) )

   if( Device_Number==999 .and. ( Flag_MultiPhysics_Device==20 .or. Flag_MultiPhysics_Device==30 .or. Flag_MultiPhysics_Device==40 ) )then
      if( Flag_Physics==0 .or. Flag_Physics==1 )then
         Relative_Thermal_Conductivity( : )= Thermal_Conductivity_FixedDomain 
      end if

   else if( Flag_MultiPhysics_Device==20 .or. Flag_MultiPhysics_Device==30 .or. Flag_MultiPhysics_Device==40 )then
           !Flag_InitialConfiguration==200 .or. Flag_InitialConfiguration==5

      if( Flag_Physics==0 .or. Flag_Physics==1 )then
         if( Number_Obj_Func_MP_T==1 )then
            Relative_Thermal_Conductivity( 1 )= Thermal_Conductivity_FixedDomain
         else if( Number_Obj_Func_MP_T==2 )then
            Relative_Thermal_Conductivity( 1 )= Thermal_Conductivity_Material_Min
            Relative_Thermal_Conductivity( 2 )= Thermal_Conductivity_Material_Max
         else
            do Loop_MatParam_T= 1, Number_Obj_Func_MP_T
               Relative_Thermal_Conductivity( Loop_MatParam_T ) &
               = Thermal_Conductivity_Material_Min &
                +( Thermal_Conductivity_Material_Max - Thermal_Conductivity_Material_Min )/( Number_Obj_Func_MP_T -1 )*( Loop_MatParam_T -1 )
            end do ! Loop_MatParam_T
         end if
      end if
   else
      if( Number_Obj_Func_MP_T==1 )then
         Relative_Thermal_Conductivity( 1 )= Thermal_Conductivity_Material
      else if( Number_Obj_Func_MP_T==2 )then
         Relative_Thermal_Conductivity( 1 )= Thermal_Conductivity_Material_Min
         Relative_Thermal_Conductivity( 2 )= Thermal_Conductivity_Material_Max
      else
         do Loop_MatParam_T= 1, Number_Obj_Func_MP_T
            Relative_Thermal_Conductivity( Loop_MatParam_T ) &
            = Thermal_Conductivity_Material_Min &
             +( Thermal_Conductivity_Material_Max - Thermal_Conductivity_Material_Min )/( Number_Obj_Func_MP_T -1 )*( Loop_MatParam_T -1 )
         end do ! Loop_MatParam_T
      end if
   end if

   if( Device_Number==999 .and. ( Flag_MultiPhysics_Device==20 .or. Flag_MultiPhysics_Device==30 .or. Flag_MultiPhysics_Device==40 ) )then
      if( Flag_Physics==0 .or. Flag_Physics==2 )then
         Relative_Electrical_Conductivity( : )= Electric_Conductivity_FixedDomain 
      end if

   else if( Flag_MultiPhysics_Device==20 .or. Flag_MultiPhysics_Device==30 .or. Flag_MultiPhysics_Device==40 )then !Flag_InitialConfiguration==200 .or. Flag_InitialConfiguration==5

      if( Flag_Physics==0 .or. Flag_Physics==2 )then

         if( Number_Obj_Func_MP_V==1 )then
            Relative_Electrical_Conductivity( 1 )= Electric_Conductivity_FixedDomain 
         else if( Number_Obj_Func_MP_V==2 )then
            Relative_Electrical_Conductivity( 1 )= Electric_Conductivity_Material_Min 
            Relative_Electrical_Conductivity( 2 )= Electric_Conductivity_Material_Max
         else
            do Loop_MatParam_V= 1, Number_Obj_Func_MP_V 
               Relative_Electrical_Conductivity( Loop_MatParam_V ) &
               = Electric_Conductivity_Material_Min &
                +( Electric_Conductivity_Material_Max -Electric_Conductivity_Material_Min )/( Number_Obj_Func_MP_V -1 )*( Loop_MatParam_V -1 )
            end do ! Loop_MatParam_V 
         end if
      end if
   else

      if( Number_Obj_Func_MP_V==1 )then
         Relative_Electrical_Conductivity( 1 )= Electric_Conductivity_Material
      else if( Number_Obj_Func_MP_V==2 )then
         Relative_Electrical_Conductivity( 1 )= Electric_Conductivity_Material_Min 
         Relative_Electrical_Conductivity( 2 )= Electric_Conductivity_Material_Max
      else
         do Loop_MatParam_V= 1, Number_Obj_Func_MP_V 
            Relative_Electrical_Conductivity( Loop_MatParam_V ) &
            = Electric_Conductivity_Material_Min &
             +( Electric_Conductivity_Material_Max -Electric_Conductivity_Material_Min )/( Number_Obj_Func_MP_V -1 )*( Loop_MatParam_V -1 )
         end do ! Loop_MatParam_V 
      end if
   end if

   !==========================================================================================================================
   if( Flag_Physics==0 .or. Flag_Physics==1 )then
   !==========================================================================================================================

      allocate( Value_Obj_Func_Normalize_Thermal( Number_Obj_Func_Source, Number_Obj_Func_MP_T, Number_Obj_Func ) )
    
      if( Device_Number==999 .and. ( Flag_InitialConfiguration==0 .or. Flag_InitialConfiguration==5 ) )then
         Value_Obj_Func_Normalize_Thermal= 1.0d0
      else
         do Loop_Source= 1, Number_Obj_Func_Source
            do Loop_MatParam_T= 1, Number_Obj_Func_MP_T
               do Loop_Obj_Func= 1, Number_Obj_Func
                  call Define_Value_ObjectiveFunction_Normalize_HT &
                     ( Position_X_Source( Loop_Source ), Position_Y_Source( Loop_Source ), Relative_Thermal_Conductivity( Loop_MatParam_T ), Loop_Obj_Func, &
                       Value_Obj_Func_Normalize_Thermal( Loop_Source, Loop_MatParam_T, Loop_Obj_Func ) )
               end do ! Loop_Obj_Func
            end do ! Loop_Source
         end do ! Loop_MatParam_T
      end if

      call Check_Value_Matrix_3_Dimension &
           ( Number_Obj_Func_Source, Number_Obj_Func_MP_T, Number_Obj_Func, &
             Value_Obj_Func_Normalize_Thermal, 0.0d0+1d-8, 'Main', 686, 'out' )
   
      call Check_Value_Matrix_3_Dimension &
           ( Number_Obj_Func_Source, Number_Obj_Func_MP_T, Number_Obj_Func, &
             Value_Obj_Func_Normalize_Thermal, 1.0d8, 'Main', 650, 'in' )

   end if

   !==========================================================================================================================
   if( Flag_Physics==0 .or. Flag_Physics==2 )then
   !==========================================================================================================================

      allocate( Value_Obj_Func_Normalize_DC( Number_Obj_Func_Source, Number_Obj_Func_MP_V, Number_Obj_Func ) )
    
      if( Device_Number==999 .and. ( Flag_InitialConfiguration==0 .or. Flag_InitialConfiguration==5 ) )then
         Value_Obj_Func_Normalize_DC= 1.0d0
      else
         do Loop_Source= 1, Number_Obj_Func_Source
            do Loop_MatParam_V= 1, Number_Obj_Func_MP_V
               do Loop_Obj_Func= 1, Number_Obj_Func
                  call Define_Value_ObjectiveFunction_Normalize_DC &
                     ( Position_X_Source( Loop_Source ), Position_Y_Source( Loop_Source ), Relative_Electrical_Conductivity( Loop_MatParam_V ), Loop_Obj_Func,  &
                       Value_Obj_Func_Normalize_DC( Loop_Source, Loop_MatParam_V, Loop_Obj_Func ) )
               end do ! Loop_Obj_Func
            end do ! Loop_Source
         end do ! Loop_MatParam_V
      end if

      call Check_Value_Matrix_3_Dimension &
           ( Number_Obj_Func_Source, Number_Obj_Func_MP_V, Number_Obj_Func, &
             Value_Obj_Func_Normalize_DC, 0.0d0+1d-8, 'Main', 675, 'out' )

      call Check_Value_Matrix_3_Dimension &
           ( Number_Obj_Func_Source, Number_Obj_Func_MP_V, Number_Obj_Func, &
             Value_Obj_Func_Normalize_DC, 1.0d8, 'Main', 679, 'in' )

   end if

   !==========================================================================================================================
   ! allocate area
   !==========================================================================================================================
   allocate( Fitness_CMAES( Number_Sampling ) )
   allocate( Obj_Func_Real( Number_Sampling ) )
   allocate( Fractal_Structure_All( Number_Sampling ) )
   allocate( Perimeter_Structure_All( Number_Sampling ) )
   if( Flag_Perimeter_Constraint==2 ) allocate( Perimeter_Implicit_All( Number_Sampling ) )
   allocate( Volume_Structure_All( Number_Sampling ) )
   allocate( Penalty_Volume_Constraint( Number_Sampling ) )
   allocate( Individual_Optimal( 0:Number_Optimization_Step ) )

   if( Flag_Physics==0 )then
      Number_Obj_Func_MaterialParameter = max( Number_Obj_Func_MP_T, Number_Obj_Func_MP_V )
   else if( Flag_Physics==1 )then
      Number_Obj_Func_MaterialParameter = Number_Obj_Func_MP_T 
   else if( Flag_Physics==2 )then
      Number_Obj_Func_MaterialParameter = Number_Obj_Func_MP_V
   end if
   allocate( Obj_Func_Multi( Number_Obj_Func, Number_Obj_Func_MaterialParameter, Number_Obj_Func_Source, Number_Physics, Number_Sampling ) )

   allocate( Fitness_Convergence( 0:Number_Optimization_Step ) )
   allocate( Convergence_Ratio( Number_Optimization_Step ) )
   allocate( Convergence_Ratio_Average( Number_Optimization_Step ) )
   allocate( Convergence_Ratio_xmean( 0:Number_Optimization_Step ) )

   allocate( Termination( 0:Number_Optimization_Step ) )

   allocate( Threshold_Convergence( Number_Threshold_Convergence ) )
   allocate( Threshold_Convergence_Average( Number_Threshold_Convergence_Average ) )
   allocate( Threshold_Convergence_xmean( Number_Threshold_Convergence_xmean ) )
   allocate( Average_Diff_LSF_PC( Number_Sampling) ) 

   if( Flag_CPU_Time==1 )then
      allocate( Time_CPU(0:Number_Optimization_Step) )
      allocate( Generation_CPU( 2 ) )
      allocate( Fitness_CPU_Threshold( Number_Fitness_CPU_Threshold ) )
      allocate( Obj_Func_CPU_Threshold( Number_Obj_Func_CPU_Threshold ) )

      Fitness_CPU_Threshold( 1 )= 1d0
      Fitness_CPU_Threshold( 2 )= 0.75d0
      Fitness_CPU_Threshold( 3 )= 0.5d0
      Fitness_CPU_Threshold( 4 )= 0.25d0

      Obj_Func_CPU_Threshold( 1 )= 1d0
      Obj_Func_CPU_Threshold( 2 )= 0.75d0
      Obj_Func_CPU_Threshold( 3 )= 0.5d0
      Obj_Func_CPU_Threshold( 4 )= 0.25d0

      if( MyRank==0 )then
         call cpu_time( Time_CPU_Start )
      end if
   end if

   call MPI_Barrier( MPI_COMM_WORLD, ierr )

   Optimization_Step=0
   Termination_dowhile='n'
   Termination(:)='n'
   Convergence_Ratio_xmean_dowhile=1.0d0
   Convergence_Ratio_xmean= 1.0d8

   if( Type_CMA=='Sparse' .and. Flag_Symmetry_LSF/=0 .and. Flag_Symmetry_LSF/=360 )then 

      allocate( Neighboring_LSF( 1 +Number_Neighboring_GP, Number_Design_Variable ) ) !( 4, )

      call Detect_Neighboring_LSF&
         ( Number_Grid_Scattering, Number_GridPoint_Scattering, Index_Grid_2_GridPoint, &
           Number_Neighboring_GP, Number_Design_Variable, Number_GridPoint_Optimization, Number_Type_LSF, Number_Quadrant, &
           Index_GP_All_2_DV_Symmetric, Index_DV_Symmetric_2_GP_All, &
         !==============================================================================
         Neighboring_LSF )

      !do i= 1, Number_Design_Variable
      !   write(809,'(i8,a3)', advance='no') i, ' : ' 
      !   do j = 1, 1 +Number_Neighboring_GP -1
      !      write(809,'(i8)', advance='no') Neighboring_LSF( j, i )
      !   enddo
      !   write(809,'(i8)') Neighboring_LSF( 1 +Number_Neighboring_GP, i )
      !enddo
      !stop

      allocate( dep_dv( 1+Number_Neighboring_GP, Number_Design_Variable ) ) !( 4, )
      dep_dv = Neighboring_LSF

      open( 68, file='dep_dv.dat', status='replace')
         write(68,*) 1+Number_Neighboring_GP
         !write(68,*) dep_dv
         do j= 1, Number_Design_Variable
            do i= 1, Number_Neighboring_GP
               write(68,'(I8)', advance='no') dep_dv( i, j )
            end do
            write(68,'(I8)')  dep_dv( Number_Neighboring_GP+1, j )
         end do
      close( 68 )

      deallocate( Neighboring_LSF )
   end if


   !==========================================================================================================================
   !do Optimization_Step= 0, Number_Optimization_Step
   !do while( Convergence_Ratio_xmean_dowhile > Convergence_Error ) 
   do while( Termination_dowhile/='y' )
   !==========================================================================================================================
    
      call MPI_Barrier( MPI_COMM_WORLD, ierr )

      if( MyRank==0 )then
         if( Flag_Read_LSF=='true' .and. Flag_xmean=='true' .and. Optimization_Step==0 )then
             FileName_LSF= trim( '../../../Optimization/tau_1e-4/result/LSF_Convergence_001000' )
            write(*,*)'     Read Level Set Function and set it as initial mean vector'

            open( 794, file=FileName_LSF, action='read' ) 
            read( 794, * ) Generation_Output_Read 
            read( 794, * ) Threshold_Convergence_xmean_read
            read( 794, * ) Flag_Symmetry_LSF_read
            read( 794, * ) Number_GridPoint_Optimization
            do j= 1, Number_GridPoint_Optimization
               read( 794, * ) No, Design_Variable( j, 1 ) 
            end do
            close( 794 )
 
            if( Flag_Symmetry_LSF_read/=Flag_Symmetry_LSF )then
               write(*,*)'========================================'
               write(*,*)'Flag_Symmetry_LSF_read/=Flag_Symmetry_LSF'
               write(*,*)'Flag_Symmetry_LSF_read=', Flag_Symmetry_LSF_read
               write(*,*)'Flag_Symmetry_LSF=', Flag_Symmetry_LSF
               write(*,*)'========================================'
               call Output_Error( 'Main', 810 ) 
            end if

            do j= 1, Number_GridPoint_Optimization
               if( Design_Variable( j, 1 ) < -0.3d0 )then
                  Design_Variable( j, 1 )= -0.3d0
               else if( Design_Variable( j, 1 ) > 0.3d0 )then
                  Design_Variable( j, 1 )= 0.3d0
               end if
            end do

            call cmaes( Design_Variable, Termination( Optimization_Step ), sigma, flag_eigen, &
                        'min', Type_CMA, Optimization_Step, Number_Design_Variable, Number_Sampling, -1d0, 1d0, Fitness_CMAES )
            if( Optimization_Step >= 1 ) Termination_dowhile= Termination( Optimization_Step )
            !call cmaes( Design_Variable, Convergence_Ratio_xmean( Optimization_Step ), sigma, flag_eigen, &
            !            'min', Type_CMA, Optimization_Step, Number_Design_Variable, Number_Sampling, -1d0, 1d0, Fitness_CMAES )
            !if( Optimization_Step >= 1 ) Convergence_Ratio_xmean_dowhile= Convergence_Ratio_xmean( Optimization_Step )

         else if( Flag_Read_LSF=='true' .and. Number_Optimization_Step==1 .and. Number_Sampling==1 )then
             FileName_LSF= trim( '../../../Optimization/tau_1e-4/result/LSF_Convergence_001000' )
            open( 732, file=FileName_LSF, action='read' ) 
            read( 732, * ) Generation_Output_Read 
            read( 732, * ) Threshold_Convergence_xmean_read
            read( 732, * ) Flag_Symmetry_LSF_read
            read( 732, * ) Number_GridPoint_Optimization
            do j= 1, Number_GridPoint_Optimization
               read( 732, * ) No, Design_Variable( j, Number_Sampling ) 
            end do
            close( 732 )
 
           if( Flag_Symmetry_LSF_read/=Flag_Symmetry_LSF )then
              write(*,*)'========================================'
              write(*,*)'Flag_Symmetry_LSF_read/=Flag_Symmetry_LSF'
              write(*,*)'Flag_Symmetry_LSF_read=', Flag_Symmetry_LSF_read
              write(*,*)'Flag_Symmetry_LSF=', Flag_Symmetry_LSF
              write(*,*)'========================================'
              call Output_Error( 'Main', 752 ) 
           end if
         else if( Device_Number==999 .and. Number_Optimization_Step==1 .and. Number_Sampling==1 )then
            850 continue
         else
              !===========================================================================
              write(*,*)'     Update Level Set Function, CMA :: ', Type_CMA 
              !===========================================================================

            call cmaes( Design_Variable, Termination( Optimization_Step ), sigma, flag_eigen, &
                        'min', Type_CMA, Optimization_Step, Number_Design_Variable, Number_Sampling, -1d0, 1d0, Fitness_CMAES )
            if( Optimization_Step >= 1 ) Termination_dowhile= Termination( Optimization_Step )
            !call cmaes( Design_Variable, Convergence_Ratio_xmean( Optimization_Step ), sigma, flag_eigen, &
            !            'min', Type_CMA, Optimization_Step, Number_Design_Variable, Number_Sampling, -1d0, 1d0, Fitness_CMAES )
            !if( Optimization_Step >= 1 ) Convergence_Ratio_xmean_dowhile= Convergence_Ratio_xmean( Optimization_Step )
         end if
      end if

      call MPI_Barrier( MPI_COMM_WORLD, ierr )
      call MPI_Bcast( Convergence_Ratio_xmean(0), 1+Optimization_Step, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
      call MPI_Bcast( Termination(0), 1+Optimization_Step, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)

      allocate( Vector_Data_MPI( Number_Design_Variable*Number_Sampling ) )

      call Matrix_2_Vector&
           ( Design_Variable, Number_Design_Variable, Number_Sampling, &
             Vector_Data_MPI, Number_Design_Variable*Number_Sampling )
      call MPI_Barrier( MPI_COMM_WORLD, ierr )
      call MPI_Bcast( Vector_Data_MPI( 1 ), Number_Design_Variable*Number_Sampling, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
      call Vector_2_Matrix&
           ( Vector_Data_MPI, Number_Design_Variable*Number_Sampling, & 
             Design_Variable, Number_Design_Variable, Number_Sampling )

      deallocate( Vector_Data_MPI )

      if( MyRank==0 ) write(*,*) '====================================================================='
      if( MyRank==0 ) write(*,*) 'Optimization Step=', Optimization_Step, '/', Number_Optimization_Step
      if( MyRank==0 ) write(*,*) '====================================================================='

      if( Optimization_Step==0 )then
         Loop_Start= 1
      else
         !Loop_Start= Number_Solution_Survived +1
         Loop_Start= 1
      end if

      allocate( Loop_Num_MPI( 2 ) )

      call MPI_Loop_Number( Loop_Start, Number_Sampling, Nproc, MyRank, Loop_Num_MPI, Number_Loop_MPI )

      if( MyRank==0 )then
         open(560, file='MPI_Info.dat', position='append')
            write(560,*)'Number_Loop_MPI=', Number_Loop_MPI, MyRank
         close(560)
      end if

      if( Optimization_Step <= 10 )then
         Interval_Plot_Optimization_Step= 1
      else if( Optimization_Step <= 100 )then
         Interval_Plot_Optimization_Step= 10
      else if( Optimization_Step <= 1000 )then
         Interval_Plot_Optimization_Step= 100
      else if( Optimization_Step <= 3000 )then
         Interval_Plot_Optimization_Step= 200
      else if( Optimization_Step <= 10000 )then
         Interval_Plot_Optimization_Step= 500
      else if( Optimization_Step > 10000 )then
         Interval_Plot_Optimization_Step= 1000
      end if


      allocate( Fractal_MPI( Loop_Num_MPI( 1 ):Loop_Num_MPI( 2 ) ) )
      allocate( Perimeter_MPI( Loop_Num_MPI( 1 ):Loop_Num_MPI( 2 ) ) )
      allocate( Perimeter_Implicit_MPI( Loop_Num_MPI( 1 ):Loop_Num_MPI( 2 ) ) )
      allocate( Volume_MPI( Loop_Num_MPI( 1 ):Loop_Num_MPI( 2 ) ) )
      allocate( Obj_Func_MPI( Loop_Num_MPI( 1 ):Loop_Num_MPI( 2 ) ) )
      allocate( Fitness_MPI( Loop_Num_MPI( 1 ):Loop_Num_MPI( 2 ) ) )
      allocate( Average_Diff_LSF_PC_MPI( Loop_Num_MPI( 1 ):Loop_Num_MPI( 2 ) ) )
      allocate( Penalty_Volume_Constraint_MPI( Loop_Num_MPI( 1 ):Loop_Num_MPI( 2 ) ) )

      !do Loop_Solution= Loop_Start, Number_Sampling 
      do Loop_Solution= Loop_Num_MPI( 1 ), Loop_Num_MPI( 2 ) 

         write(*,*) '   ====================================================================='
         write(*,*) '   Loop_Solution=', Loop_Solution -Loop_Start +1, '/',  Number_Sampling -Loop_Start +1
         write(*,*) '   ====================================================================='

         if( MyRank==0 )then
            open(497, file='Process.dat', position='append')
               write(497,*) Loop_Solution, '/', Loop_Num_MPI( 2 ), ':', Optimization_Step
            close(497)
         end if

         allocate( LSF_GridPoint( Number_GridPoint_Scattering, Number_Type_LSF ) )

         if( Flag_InitialConfiguration==200 .or. Flag_InitialConfiguration==10 )then
            LSF_GridPoint= 0d0
         else if( Flag_InitialConfiguration==0 )then
            LSF_GridPoint= -1d0
         else if( Flag_InitialConfiguration==5 )then
            LSF_GridPoint= 0d0
            Design_Variable= 1d0
         else
            write(*,*)'Flag_InitialConfiguration=', Flag_InitialConfiguration
            call Output_Error( 'Main', 406 ) 
         end if

         if( Flag_Check_Value==1 ) &
         call Check_Value_Matrix( Number_Design_Variable, Number_Sampling, Design_Variable, 1d0 +1d-3, 'Main', 686 )

         if( Flag_InitialConfiguration==10 )then
            do i = 1, Number_Type_LSF
               do j= 1, Number_GridPoint_Scattering
                  LSF_GridPoint( j, i )&
                  = Radius_InitialConfiguration &
                    -sqrt( ( Position_GridPoint_Scattering( 1, j )-Position_Center_FixedDomain_X_Cloaked_Region )**2 &
                          +( Position_GridPoint_Scattering( 2, j )-Position_Center_FixedDomain_Y_Cloaked_Region )**2 )
               end do
            end do

            do i = 1, Number_Type_LSF
               do j= 1, Number_GridPoint_Scattering
                  if( LSF_GridPoint( j, i ) < -1.0d0 )then
                     LSF_GridPoint( j, i ) = -1.0d0
                  else if( LSF_GridPoint( j, i ) > 1.0d0 )then
                     LSF_GridPoint( j, i ) = 1.0d0
                  end if
               end do
            end do
         else if( Flag_InitialConfiguration/=0 )then
            if( Flag_Symmetry_LSF==0 )then
               do i = 1, Number_Type_LSF 
                  do j = 1, Number_GridPoint_Optimization
                     LSF_GridPoint( Correspondence_Number( j ), i )=  Design_Variable( j*i, Loop_Solution)
                  end do    ! set
               end do    ! set
            else
               do i = 1, Number_Type_LSF 
                  do j = 1, Number_GridPoint_Optimization
                     do k= 1, Number_Quadrant 
                        LSF_GridPoint( Index_DV_Symmetric_2_GP_All( k, j ), i )= Design_Variable( j*i, Loop_Solution )
                     end do
                  end do
               end do
            end if
         end if

         if( Level_Constraint_Minimum_Length >= 1 )then 
            allocate( LSF_GridPoint_tmp( Number_GridPoint_Scattering ) )
            do i = 1, Number_Type_LSF 
               LSF_GridPoint_tmp( : )= LSF_GridPoint( :, i )

               call Constraint_Minimum_Length &
                  ( Number_Grid_Scattering, Number_GridPoint_Scattering, Index_Grid_2_GridPoint, &
                    LSF_FixedDomain_GridPoint, LSF_DesignDomain_GridPoint, & 
                    Level_Constraint_Minimum_Length, Position_GridPoint_Scattering, &
                  !==============================================================================
                    LSF_GridPoint_tmp )

               LSF_GridPoint( :, i )= LSF_GridPoint_tmp( : )
            end do
            deallocate( LSF_GridPoint_tmp )
         end if

         ! aho
         if( Flag_InitialConfiguration==200 )then
            call Resize_LSF &
               ( Index_Grid_2_GridPoint, Position_GridPoint_Scattering,  &
               !==============================================================================
                 LSF_GridPoint, LSF_FixedDomain_GridPoint, LSF_DesignDomain_GridPoint ) 
         end if

         if( Flag_Check_Value==1 ) &
         call Check_Value_Vector( Number_GridPoint_Scattering, LSF_GridPoint, 1d0 +1d-3, 'Main', 576 )

         do i = 1, Number_Type_LSF 
            do j = 1, Number_GridPoint_Optimization
               LSF_Piecewise_Constant( j*i, Loop_Solution )= LSF_GridPoint( Index_DV_Symmetric_2_GP_All( 1, j ), i )
            end do
         end do

         do i = 1, Number_Type_LSF 
            do j = 1, Number_GridPoint_Optimization
               Diff_LSF_PC( j*i, Loop_Solution )= Design_Variable( j*i, Loop_Solution ) -LSF_Piecewise_Constant( j*i, Loop_Solution )
            end do
         end do

         Average_Diff_LSF_PC_MPI( Loop_Solution )= 0d0
         do i = 1, Number_Type_LSF 
            do j = 1, Number_GridPoint_Optimization
               Average_Diff_LSF_PC_MPI( Loop_Solution )&
               = Average_Diff_LSF_PC_MPI( Loop_Solution ) &
                +abs( Diff_LSF_PC( j*i, Loop_Solution ) )/dble( Number_Design_Variable ) 
            end do
         end do

         if( Flag_Check_Value==1 ) &
         call Check_Value_Matrix( Number_Design_Variable, Number_Sampling, Design_Variable, 1d0 +1d-3, 'Main', 577 )

         !========================================================================
         if( MyRank==0 ) write(*,*) '   Generate Mesh'
         !========================================================================
 
         allocate( tmp%dp2( 2, Number_Node_tmp ) )
         allocate( tmp%itg2( 4, Number_Element_tmp ) )
         allocate( tmp%itg1( Number_Element_tmp ) )

         allocate( Index_GridPointNumber_2_NodeNumber( Number_GridPoint_Scattering ) )

         call Generate_Mesh( Optimization_Step, &
              LSF_GridPoint, Position_GridPoint_Scattering, &
              !Position_Node_OneGrid_PML, Number_Node_OneGrid_PML, Number_Element_OneGrid_PML, &
              !Class_Element_OneGrid_PML, Index_Element_PML_2_Node_PML, & 
              LSF_FixedDomain_GridPoint, LSF_DesignDomain_GridPoint, LSF_ExteriorDomain_GridPoint, &
              Index_Grid_2_GridPoint, &
              !================================================================================
              tmp%dp2, tmp%itg2, tmp%itg1, & 
              femdata%number_node, femdata%number_element_all, &
              Index_GridPointNumber_2_NodeNumber )
 
         if( MyRank==0 ) write(*,*)'   =================================================== '
         if( MyRank==0 ) write(*,*)'   Number_Node=', femdata%number_node 
         if( MyRank==0 ) write(*,*)'   Number_Element=', femdata%number_element_all 
         if( MyRank==0 ) write(*,*)'   =================================================== '
    
         if( MyRank==0 )then
            open( 934, file= Filename_FEM_Data, position='append')
               write(934,*)'Optimization_Step=', Optimization_Step 
               write(934,*)'Number_Node=', femdata%number_node 
               write(934,*)'Number_Element=', femdata%number_element_all 
               write(934,*)' '
            close( 934 )
         end if
    
         !========================================================================
         if( MyRank==0 ) write(*,*) '   Temporaly Data --> Input Data '
         !========================================================================
    
         allocate( femdata%position_node( 2, femdata%number_node ) )
         allocate( femdata%index_element_2_node( 4, femdata%number_element_all ) )
         allocate( femdata%type_element( femdata%number_element_all ) )
    
         femdata%position_node( :, 1:femdata%number_node )= tmp%dp2( :, 1:femdata%number_node )
         femdata%index_element_2_node( :, 1:femdata%number_element_all )= tmp%itg2( :, 1:femdata%number_element_all )
         femdata%type_element( 1:femdata%number_element_all )= tmp%itg1( 1:femdata%number_element_all ) 

         deallocate( tmp%dp2 )
         deallocate( tmp%itg2 )
         deallocate( tmp%itg1 )

         !========================================================================
         if( MyRank==0 ) write(*,*) '   Create Input Data '
         !========================================================================
   
         allocate( femdata%position_node_renumber( 2, femdata%number_node ) )
         allocate( femdata%type_node( femdata%number_node ) )
         allocate( tmp%itg2( 3, femdata%number_element_all ) )
         allocate( tmp%itg1( femdata%number_element_all ) )
    
         allocate( Max_Position_Node( 2, 2 ) )
    
         allocate( Index_GridPointNumber_2_NodeNumber_Renumbered( Number_GridPoint_Scattering ) )

         call Create_Input_Data & 
            ( femdata%number_node, femdata%number_element_all, &
              femdata%position_node, femdata%index_element_2_node, femdata%type_element, &
              Length_Characteristic, Index_GridPointNumber_2_NodeNumber, & 
              ID_Element_FixedDomain, &
             !==============================================================================================
              femdata%position_node_renumber, femdata%type_node,  &
              tmp%itg2, tmp%itg1, femdata%number_element_triangle, & 
              Number_Node_Reference, & 
              Number_Node_on_Electrical_Insulation_Boundary, Number_Node_in_Electrical_Insulation, & 
              Number_Element_Electrical_Insulation_Boundary, &
              Max_Number_Element_Share_OneNode, Max_Position_Node, &
              Index_GridPointNumber_2_NodeNumber_Renumbered )

         allocate( Element_and_LocalNode_Number_on_Electrical_Insulation_BC( 3, Number_Element_Electrical_Insulation_Boundary ) )
    
         call ReEdit_InputData_Electrical_Insulation_Boundary_Condition &
            ( tmp%itg2, femdata%number_element_triangle, &
              Number_Node_on_Electrical_Insulation_Boundary, Number_Node_in_Electrical_Insulation, & 
              Element_and_LocalNode_Number_on_Electrical_Insulation_BC )
    
         femdata%number_element_all= femdata%number_element_triangle 
   
         deallocate( Index_GridPointNumber_2_NodeNumber )
         deallocate( femdata%position_node )
         deallocate( femdata%type_element )
         deallocate( femdata%index_element_2_node )

         !========================================================================
         if( MyRank==0 ) write(*,*) '   Temporaly Element Data --> Element Data '
         !========================================================================
    
         allocate( femdata%index_element_triangle_2_node( 3, femdata%number_element_triangle ) )
         allocate( femdata%type_element_triangle( femdata%number_element_triangle ) )
   
         femdata%index_element_triangle_2_node( :, 1:femdata%number_element_triangle ) &
         = tmp%itg2( :, 1:femdata%number_element_triangle )
         femdata%type_element_triangle( 1:femdata%number_element_triangle ) &
         = tmp%itg1( 1:femdata%number_element_triangle )

         deallocate( tmp%itg2 )
         deallocate( tmp%itg1 )
 
         if( Flag_Output_Finite_Element_Data/=0 .and. Loop_Solution==1 .and. &
             Optimization_Step==Optimization_Step_Output .and. Flag_InitialConfiguration==0 .and. &
             MyRank==0 )then

            !========================================================================
            write(*,*) '   Output FEM Data for Normalization '
            !========================================================================
    
            call Output_Finite_Element_Data &
               ( femdata%number_node, femdata%position_node_renumber, &
                 Max_Number_Element_Share_OneNode, Max_Position_Node, &
                 femdata%number_element_triangle, femdata%index_element_triangle_2_node, &
                 femdata%type_element_triangle, &
                 Element_and_LocalNode_Number_on_Electrical_Insulation_BC, &
                 Number_Node_on_Electrical_Insulation_Boundary, Number_Node_in_Electrical_Insulation, &
                 Number_Element_Electrical_Insulation_Boundary, &
                 0, Optimization_Step )
         end if

         !========================================================================
         if( MyRank==0 ) write(*,*) '   Modify Local Node Number'
         !========================================================================
    
         call Modify_Local_Node_Number &
            ( femdata%number_element_triangle, 3, femdata%index_element_triangle_2_node, & 
              femdata%number_node, femdata%position_node_renumber )
    
         !call Modify_Local_Node_Number &
         !   ( femdata%number_element_square, 4, femdata%index_element_square_2_node, & 
         !     femdata%number_node, femdata%position_node_renumber )
   
         !========================================================================
         if( MyRank==0 ) write(*,*) '   Compute Perimeter of Dielectric Strucuture '
         !========================================================================
         if( Device_Number==999 )then
            Fractal_Structure = 1.0d0 
         else
          if(Type_fractal=='area')then  
          call fractal_area &
                 ( LSF_GridPoint, Position_GridPoint_Scattering, Index_Grid_2_GridPoint,&
                   Number_GridPoint_Scattering, Number_Type_LSF, Number_Grid_Scattering,&
                   Size_X_Scattering, Box_Number, Index_Grid_total,Grid_total,plot_Number,&
                   Fractal_Structure ) 
         
          else if (Type_fractal=='line')then
           call fractal_line&
                 ( LSF_GridPoint, Position_GridPoint_Scattering,Index_Grid_2_GridPoint,&
                   Number_GridPoint_Scattering, Number_Type_LSF,Number_Grid_Scattering,&
                   Size_X_Scattering, Box_Number,Index_Grid_total,Grid_total,plot_Number,&
                   Fractal_Structure )
          end if
         end if
 
         
         call Compute_Perimeter &
              ( femdata%position_node_renumber, femdata%number_node, & 
                femdata%index_element_triangle_2_node, femdata%type_element_triangle, femdata%number_element_triangle, &
                ID_Element_Material, ID_Element_Base_Material, &
                !======================================================
                Perimeter_Structure )
   
         if( Flag_Perimeter_Constraint==2 )then 
            !========================================================================
            if( MyRank==0 ) write(*,*) '   Compute Implicit Perimeter of Dielectric Strucuture, |Nabla Phi|'
            !========================================================================

            call Compute_Implicit_Perimeter &
               ( Index_Grid_2_GridPoint, Position_GridPoint_Scattering, &
                 LSF_GridPoint, LSF_FixedDomain_GridPoint, LSF_DesignDomain_GridPoint, &
               !==============================================================================
                 Perimeter_Implicit_MPI( Loop_Solution ) ) 
         end if

         !========================================================================
         if( MyRank==0 ) write(*,*) '   Compute Volume of Dielectric Strucuture '
         !========================================================================

         call Compute_Volume_Constraint &
            ( femdata%position_node_renumber, femdata%number_node, &
              femdata%index_element_triangle_2_node, &
              femdata%type_element_triangle, femdata%number_element_triangle, &
              Volume_Material, Volume_DesignDomain )

         if( Flag_Check_Value==1 ) &
         call Check_Value_Matrix( Number_Design_Variable, Number_Sampling, Design_Variable, 1d0 +1d-3, 'Main', 962 )

         !========================================================================
         if( MyRank==0 ) write(*,*) '   Compute Averaged Thermal Diffusivity '
         !========================================================================

         !Number_Point_ATD( 2 )= 100
         !if( mod( Number_Point_ATD( 2 ), 2 )==0 )then
         !   Number_Point_ATD( 1 )= 2*Number_Point_ATD( 2 )
         !else if( mod( Number_Point_ATD( 2 ), 2 )==1 )then
         !   Number_Point_ATD( 1 )= 2*( Number_Point_ATD( 2 ) -1 ) +1
         !else
         !   write(*,*)'mod( Number_Point_ATD( 2 ), 2 )=', mod( Number_Point_ATD( 2 ), 2 )
         !   call Output_Error( 'Main', 1190 )
         !end if
       
         !allocate( tmp%dp3( 2, 3, femdata%number_element_triangle ) )

         !do e= 1, femdata%number_element_triangle
         !   do i= 1, 3 
         !      do j= 1, 2 
         !         tmp%dp3( j, i, e )= femdata%position_node_renumber( j, femdata%index_element_triangle_2_node( i, e ) )
         !      end do
         !   end do
         !end do

         !call Compute_Averaged_Material_Properties &
         !   ( femdata%number_element_triangle, Number_Point_ATD, &
         !     Max_Position_Node, & !Position_Minimum_Plot, Position_Maximum_Plot, &
         !     tmp%dp3, femdata%type_element_triangle, & !Position_Node_Plot_Element, Class_Element, &
         !     !ID_Element_Material, ID_Element_OuterDomain, ID_Element_FixedDomain, ID_Element_Base_Material, &
         !     !ID_Element_Material_Exterior, ID_Element_Base_Material_Exterior, &
         !     Optimization_Step )

         !deallocate( tmp%dp3 )

         !========================================================================
         if( MyRank==0 ) write(*,*) '   Analysis for Steady State Electric Potential '
         !========================================================================
   
         allocate( Temperature_Solution( femdata%number_node, Number_Obj_Func_Source, Number_Obj_Func_MP_T ) )
         allocate( Electric_Potential_Solution( femdata%number_node, Number_Obj_Func_Source, Number_Obj_Func_MP_V ) )
         Width_Matrix_LHS= Max_Number_Element_Share_OneNode +1

         do Loop_Source= 1, Number_Obj_Func_Source
            if( Flag_Physics==0 .or. Flag_Physics==1 )then
               do Loop_MatParam_T= 1, Number_Obj_Func_MP_T
      
                  allocate( Temperature_tmp( femdata%number_node ) )
      
                  call Analyze_Steady_State_Heat_Conduction &
                     ( 0.0d0, 1.0d0, & 
                       Position_X_Source( Loop_Source ), Position_Y_Source( Loop_Source ), &
                       Relative_Thermal_Conductivity( Loop_MatParam_T ), &
                       femdata%number_node, femdata%number_element_triangle, &
                       femdata%position_node_renumber, & 
                       Width_Matrix_LHS, Max_Position_Node, &
                       femdata%index_element_triangle_2_node, femdata%type_element_triangle, &
                       !============================================================================================
                       Temperature_tmp )
      
                  Temperature_Solution( :, Loop_Source, Loop_MatParam_T )= Temperature_tmp( : )
      
                  deallocate( Temperature_tmp )
      
               end do ! Loop_MatParam_T
            end if
 
            if( Flag_Physics==0 .or. Flag_Physics==2 )then
               do Loop_MatParam_V= 1, Number_Obj_Func_MP_V
      
                  allocate( Electric_Potential_tmp( femdata%number_node ) )
      
                  call Analyze_Steady_State_Electric_Potential &
                     ( 0.0d0, 1.0d0, & 
                       Position_X_Source( Loop_Source ), Position_Y_Source( Loop_Source ), &
                       Relative_Electrical_Conductivity( Loop_MatParam_V ), &
                       femdata%number_node, femdata%number_element_triangle, &
                       femdata%position_node_renumber, & 
                       Width_Matrix_LHS, Max_Position_Node, &
                       femdata%index_element_triangle_2_node, femdata%type_element_triangle, &
                       !============================================================================================
                       Electric_Potential_tmp )
      
                  Electric_Potential_Solution( :, Loop_Source, Loop_MatParam_V )= Electric_Potential_tmp( : )
      
                  deallocate( Electric_Potential_tmp )
      
               end do ! Loop_MatParam_V
            end if
         end do ! Loop_Source
   
         !========================================================================
         if( MyRank==0 ) write(*,*)'   Obtain Reference Data for Objective Function'
         !========================================================================

         if( Optimization_Step==0 .and. ( Flag_Physics==0 .or. Flag_Physics==1 ) )then
            if( Flag_ObjectiveFunction==11 .and. Loop_Solution==Loop_Num_MPI( 1 ) )then
               allocate( Temperature_Reference( Number_Node_Reference, Number_Obj_Func_Source, Number_Obj_Func_MP_T ) )
      
               do Loop_Source= 1, Number_Obj_Func_Source
                  do Loop_MatParam_T= 1, Number_Obj_Func_MP_T
      
                     allocate( Temperature_tmp( femdata%number_node ) )
      
                     Temperature_tmp( : )= Temperature_Solution( :, Loop_Source, Loop_MatParam_T )
      
                     allocate( tmp%dp1( Number_Node_Reference ) )
       
                     call Operate_Data_Reference_Field &
                        ( Optimization_Step, 1, femdata%number_node, Number_Node_Reference, Temperature_tmp, &
                          dble( Loop_Source ), dble( Loop_MatParam_T ), &
                          tmp%dp1 )
   
                     deallocate( Temperature_tmp )
       
                     Temperature_Reference( :, Loop_Source, Loop_MatParam_T )= tmp%dp1( : )
       
                     deallocate( tmp%dp1 )
                  end do
               end do
            end if
   
            call MPI_Barrier( MPI_COMM_WORLD, ierr )
   
            if( Flag_ObjectiveFunction==11 .and. Loop_Solution==Loop_Num_MPI( 1 ) )then
   
               allocate( Vector_Data_Real_MPI( Number_Node_Reference*Number_Obj_Func_Source*Number_Obj_Func_MP_T ) )
               call Matrix3_2_Vector_Real &
                  ( Temperature_Reference, Number_Node_Reference, Number_Obj_Func_Source, Number_Obj_Func_MP_T, &
                    Vector_Data_Real_MPI, Number_Node_Reference*Number_Obj_Func_Source*Number_Obj_Func_MP_T )
               call MPI_Barrier( MPI_COMM_WORLD, ierr )
               call MPI_Bcast( Vector_Data_Real_MPI( 1 ), Number_Node_Reference*Number_Obj_Func_Source*Number_Obj_Func_MP_T, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )
               call Vector_2_Matrix3_Real&
                  ( Vector_Data_Real_MPI, Number_Node_Reference*Number_Obj_Func_Source*Number_Obj_Func_MP_T, &
                    Temperature_Reference, Number_Node_Reference, Number_Obj_Func_Source, Number_Obj_Func_MP_T )
               deallocate( Vector_Data_Real_MPI )
            end if
    
            call MPI_Barrier( MPI_COMM_WORLD, ierr )

         end if

         if( Optimization_Step==0 .and. ( Flag_Physics==0 .or. Flag_Physics==2 ) )then
            if( Flag_ObjectiveFunction==11 .and. Loop_Solution==Loop_Num_MPI( 1 ) )then
               allocate( Electric_Potential_Reference( Number_Node_Reference, Number_Obj_Func_Source, Number_Obj_Func_MP_V ) )
   
               do Loop_Source= 1, Number_Obj_Func_Source
                  do Loop_MatParam_V= 1, Number_Obj_Func_MP_V
   
                     allocate( Electric_Potential_tmp( femdata%number_node ) )
   
                     Electric_Potential_tmp( : )= Electric_Potential_Solution( :, Loop_Source, Loop_MatParam_V )
   
                     allocate( tmp%dp1( Number_Node_Reference ) )
    
                     call Operate_Data_Reference_Field &
                        ( Optimization_Step, 2, femdata%number_node, Number_Node_Reference, Electric_Potential_tmp, &
                          dble( Loop_Source ), dble( Loop_MatParam_V ), &
                          tmp%dp1 )

                     deallocate( Electric_Potential_tmp )
    
                     Electric_Potential_Reference( :, Loop_Source, Loop_MatParam_V )= tmp%dp1( : )
    
                     deallocate( tmp%dp1 )
                  end do
               end do
            end if

            call MPI_Barrier( MPI_COMM_WORLD, ierr )

            if( Flag_ObjectiveFunction==11 .and. Loop_Solution==Loop_Num_MPI( 1 ) )then

               allocate( Vector_Data_Real_MPI( Number_Node_Reference*Number_Obj_Func_Source*Number_Obj_Func_MP_V ) )
               call Matrix3_2_Vector_Real &
                  ( Electric_Potential_Reference, Number_Node_Reference, Number_Obj_Func_Source, Number_Obj_Func_MP_V, &
                    Vector_Data_Real_MPI, Number_Node_Reference*Number_Obj_Func_Source*Number_Obj_Func_MP_V )
               call MPI_Barrier( MPI_COMM_WORLD, ierr )
               call MPI_Bcast( Vector_Data_Real_MPI( 1 ), Number_Node_Reference*Number_Obj_Func_Source*Number_Obj_Func_MP_V, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )
               call Vector_2_Matrix3_Real&
                  ( Vector_Data_Real_MPI, Number_Node_Reference*Number_Obj_Func_Source*Number_Obj_Func_MP_V, &
                    Electric_Potential_Reference, Number_Node_Reference, Number_Obj_Func_Source, Number_Obj_Func_MP_V )
               deallocate( Vector_Data_Real_MPI )
            end if
 
            call MPI_Barrier( MPI_COMM_WORLD, ierr )
         end if
         !========================================================================
         if( MyRank==0 ) write(*,*) '   Compute Objective Function'
         !========================================================================
   
         Obj_Func_All( Optimization_Step )= 0d0
         Counter= 0
   
         do Loop_Source= 1, Number_Obj_Func_Source
            !========================================================================
            if( Flag_Physics==0 .or. Flag_Physics==1 )then
            !========================================================================
                
               do Loop_MatParam_T= 1, Number_Obj_Func_MP_T
                  do Loop_Obj_Func= 1, Number_Obj_Func
 
                     allocate( Value_Obj_Func( femdata%number_node ) )
             
                     !=========================================
                     if( Loop_Obj_Func==Loop_Cloak )then
                     !=========================================
                        if( Flag_ObjectiveFunction==0 )then
                           Value_Obj_Func( : )&
                           = Temperature_Solution( :, Loop_Source, Loop_MatParam_T )
                        else if( Flag_ObjectiveFunction==11 )then
                           Value_Obj_Func( 1:Number_Node_Reference )&
                           = ( Temperature_Solution( 1:Number_Node_Reference, Loop_Source, Loop_MatParam_T ) &
                              -Temperature_Reference( 1:Number_Node_Reference, Loop_Source, Loop_MatParam_T ) )
      
                           Value_Obj_Func( Number_Node_Reference+1:femdata%number_node )&
                           = Temperature_Solution( Number_Node_Reference+1:femdata%number_node, Loop_Source, Loop_MatParam_T ) 
                        end if

                        if( Flag_Check_Value==1 ) &
                        !call Check_Value_Vector( femdata%number_node, abs( Value_Obj_Func ), 1.7d308, 'Main', 1160 )
                        call Check_Value_Vector( femdata%number_node, Value_Obj_Func, 1.7d308, 'Main', 1160 )
   
                        ID_Element_Obj_Func = ID_Element_Obj_Func_1

                        call Compute_Objective_Function &
                           ( ID_Element_Obj_Func, &
                             Value_Obj_Func, femdata%number_node, & 
                             femdata%position_node_renumber, femdata%number_node, &
                             femdata%index_element_triangle_2_node, femdata%type_element_triangle, femdata%number_element_triangle, &
                             Optimization_Step, Loop_Solution, &
                             Perimeter_Structure, &
                             !=================================================================
                             Obj_Func_Single_Thermal( Optimization_Step, Loop_Source, Loop_MatParam_T, Loop_Obj_Func ) )

                     !=========================================
                     else if( Loop_Obj_Func==Loop_Flux_X )then
                     !=========================================
                        Value_Obj_Func( : )&
                        = Temperature_Solution( :, Loop_Source, Loop_MatParam_T )

                        if( Flag_Check_Value==1 ) &
                        !call Check_Value_Vector( femdata%number_node, abs( Value_Obj_Func ), 1.7d308, 'Main', 1160 )
                        call Check_Value_Vector( femdata%number_node, Value_Obj_Func, 1.7d308, 'Main', 1160 )
   
                        ID_Element_Obj_Func = ID_Element_Obj_Func_2
                        Relative_Thermal_Conductivity_Evaluated = Relative_Thermal_Conductivity( Loop_MatParam_T )
 
                        call Compute_Objective_Function_Flux &
                           ( ID_Element_Obj_Func, 'nx', 'inner', 'sig', &
                             Value_Obj_Func, femdata%number_node, & 
                             femdata%position_node_renumber, femdata%number_node, &
                             femdata%index_element_triangle_2_node, femdata%type_element_triangle, femdata%number_element_triangle, &
                             Optimization_Step, Loop_Solution, &
                             Perimeter_Structure, &
                             Relative_Thermal_Conductivity_Evaluated, &
                             !=================================================================
                             Obj_Func_Single_Thermal( Optimization_Step, Loop_Source, Loop_MatParam_T, Loop_Obj_Func ) )

                     !=========================================
                     else if( Loop_Obj_Func==Loop_Flux_Y )then
                     !=========================================
                        Value_Obj_Func( : )&
                        = Temperature_Solution( :, Loop_Source, Loop_MatParam_T )

                        if( Flag_Check_Value==1 ) &
                        !call Check_Value_Vector( femdata%number_node, abs( Value_Obj_Func ), 1.7d308, 'Main', 1160 )
                        call Check_Value_Vector( femdata%number_node, Value_Obj_Func, 1.7d308, 'Main', 1160 )
   
                        ID_Element_Obj_Func = ID_Element_Obj_Func_2
                        Relative_Thermal_Conductivity_Evaluated = Relative_Thermal_Conductivity( Loop_MatParam_T )
 
                        call Compute_Objective_Function_Flux &
                           ( ID_Element_Obj_Func, 'nx', 'outer', 'abs', &
                           !( ID_Element_Obj_Func, 'nx', 'outer', 'sqn', &
                             Value_Obj_Func, femdata%number_node, & 
                             femdata%position_node_renumber, femdata%number_node, &
                             femdata%index_element_triangle_2_node, femdata%type_element_triangle, femdata%number_element_triangle, &
                             Optimization_Step, Loop_Solution, &
                             Perimeter_Structure, &
                             Relative_Thermal_Conductivity_Evaluated, &
                             !=================================================================
                             Obj_Func_Single_Thermal( Optimization_Step, Loop_Source, Loop_MatParam_T, Loop_Obj_Func ) )

                     !============================================ 
                     else 
                     !============================================ 
                        write(*,*)'Loop_Obj_Func=', Loop_Obj_Func
                        call Output_Error( 'Main', 1257 ) 
                     !=========================================
                     end if
                     !=========================================

                     deallocate( Value_Obj_Func )
             
                  end do ! Loop_Obj_Func
               end do ! Loop_MatParam_T
            end if
    
            !========================================================================
            if( Flag_Physics==0 .or. Flag_Physics==2 )then
            !========================================================================
                
               do Loop_MatParam_V= 1, Number_Obj_Func_MP_V
                  do Loop_Obj_Func= 1, Number_Obj_Func
 
                     allocate( Value_Obj_Func( femdata%number_node ) )
             
                     !============================================ 
                     if( Loop_Obj_Func==Loop_Cloak )then
                     !============================================ 
                        if( Flag_ObjectiveFunction==1 )then
                           Value_Obj_Func( : )&
                           = Electric_Potential_Solution( :, Loop_Source, Loop_MatParam_V )
                        else if( Flag_ObjectiveFunction==11 )then
                           Value_Obj_Func( 1:Number_Node_Reference )&
                           = ( Electric_Potential_Solution( 1:Number_Node_Reference, Loop_Source, Loop_MatParam_V ) &
                              -Electric_Potential_Reference( 1:Number_Node_Reference, Loop_Source, Loop_MatParam_V ) )
      
                           Value_Obj_Func( Number_Node_Reference+1:femdata%number_node )&
                           = Electric_Potential_Solution( Number_Node_Reference+1:femdata%number_node, Loop_Source, Loop_MatParam_V ) 
                        end if

                        if( Flag_Check_Value==1 ) &
                        !call Check_Value_Vector( femdata%number_node, abs( Value_Obj_Func ), 1.7d308, 'Main', 1160 )
                        call Check_Value_Vector( femdata%number_node, Value_Obj_Func, 1.7d308, 'Main', 1160 )
   
                        ID_Element_Obj_Func = ID_Element_Obj_Func_1
 
                        call Compute_Objective_Function &
                           ( ID_Element_Obj_Func, &
                             Value_Obj_Func, femdata%number_node, & 
                             femdata%position_node_renumber, femdata%number_node, &
                             femdata%index_element_triangle_2_node, femdata%type_element_triangle, femdata%number_element_triangle, &
                             Optimization_Step, Loop_Solution, &
                             Perimeter_Structure, &
                             !=================================================================
                             Obj_Func_Single_DC( Optimization_Step, Loop_Source, Loop_MatParam_V, Loop_Obj_Func ) )
    
                     !============================================ 
                     else if( Loop_Obj_Func==Loop_Flux_X )then
                     !============================================ 
                        Value_Obj_Func( : )&
                        = Electric_Potential_Solution( :, Loop_Source, Loop_MatParam_V )

                        if( Flag_Check_Value==1 ) &
                        !call Check_Value_Vector( femdata%number_node, abs( Value_Obj_Func ), 1.7d308, 'Main', 1160 )
                        call Check_Value_Vector( femdata%number_node, Value_Obj_Func, 1.7d308, 'Main', 1160 )
   
                        ID_Element_Obj_Func = ID_Element_Obj_Func_2
                        Relative_Electrical_Conductivity_Evaluated = Relative_Electrical_Conductivity( Loop_MatParam_V ) 
 
                        call Compute_Objective_Function_Flux &
                           ( ID_Element_Obj_Func, 'nx', 'inner', 'sig', &
                             Value_Obj_Func, femdata%number_node, & 
                             femdata%position_node_renumber, femdata%number_node, &
                             femdata%index_element_triangle_2_node, femdata%type_element_triangle, femdata%number_element_triangle, &
                             Optimization_Step, Loop_Solution, &
                             Perimeter_Structure, &
                             Relative_Electrical_Conductivity_Evaluated, &
                             !=================================================================
                             Obj_Func_Single_DC( Optimization_Step, Loop_Source, Loop_MatParam_V, Loop_Obj_Func ) )

                     !============================================ 
                     else if( Loop_Obj_Func==Loop_Flux_Y )then
                     !============================================ 
                        Value_Obj_Func( : )&
                        = Electric_Potential_Solution( :, Loop_Source, Loop_MatParam_V )

                        if( Flag_Check_Value==1 ) &
                        !call Check_Value_Vector( femdata%number_node, abs( Value_Obj_Func ), 1.7d308, 'Main', 1160 )
                        call Check_Value_Vector( femdata%number_node, Value_Obj_Func, 1.7d308, 'Main', 1160 )
   
                        ID_Element_Obj_Func = ID_Element_Obj_Func_2
                        Relative_Electrical_Conductivity_Evaluated = Relative_Electrical_Conductivity( Loop_MatParam_V ) 
 
                        call Compute_Objective_Function_Flux &
                           ( ID_Element_Obj_Func, 'nx', 'outer', 'abs', &
                             Value_Obj_Func, femdata%number_node, & 
                             femdata%position_node_renumber, femdata%number_node, &
                             femdata%index_element_triangle_2_node, femdata%type_element_triangle, femdata%number_element_triangle, &
                             Optimization_Step, Loop_Solution, &
                             Perimeter_Structure, &
                             Relative_Electrical_Conductivity_Evaluated, &
                             !=================================================================
                             Obj_Func_Single_DC( Optimization_Step, Loop_Source, Loop_MatParam_V, Loop_Obj_Func ) )

                     !============================================ 
                     else 
                     !============================================ 
                        write(*,*)'Loop_Obj_Func=', Loop_Obj_Func
                        call Output_Error( 'Main', 1332 ) 
                     !============================================ 
                     end if
                     !============================================ 

                     deallocate( Value_Obj_Func )
             
                  end do ! Loop_Obj_Func
               end do ! Loop_MatParam_V
            end if
 
!write(2138,*) Obj_Func_Single_Thermal
!write(2139,*) Obj_Func_Single_DC 

            !========================================================================
            if( Flag_Physics==0 )then
            !========================================================================
               do Loop_MatParam_T= 1, Number_Obj_Func_MP_T
                  !if( Type_ObjectiveFunction/='div' .and. Type_ObjectiveFunction/='smx' )then
                  if( Flag_MultiPhysics_Device==0 .or. Flag_MultiPhysics_Device==10 .or. &
                      Flag_MultiPhysics_Device==20 .or. Flag_MultiPhysics_Device==40 )then
                     do Loop_Obj_Func= 1, Number_Obj_Func
   
                        if( Loop_Obj_Func==Loop_Cloak )then
                           if( Type_ObjectiveFunction=='sum' )then
                              Obj_Func_tmp &
                              = Obj_Func_Single_Thermal( Optimization_Step, Loop_Source, Loop_MatParam_T, Loop_Obj_Func )&
                               /Value_Obj_Func_Normalize_Thermal( Loop_Source, Loop_MatParam_T, Loop_Obj_Func )&
                               /dble( Number_Obj_Func_Source )/dble( Number_Obj_Func_MP_T + Number_Obj_Func_MP_V )
                          
                              Obj_Func_All( Optimization_Step ) &
                              = Obj_Func_All( Optimization_Step ) +Obj_Func_tmp            
                           else if( Type_ObjectiveFunction=='max' )then
                              Obj_Func_tmp &
                              = Obj_Func_Single_Thermal( Optimization_Step, Loop_Source, Loop_MatParam_T, Loop_Obj_Func ) &
                               /Value_Obj_Func_Normalize_Thermal( Loop_Source, Loop_MatParam_T, Loop_Obj_Func )
           
                              Obj_Func_All( Optimization_Step ) &
                              = max( Obj_Func_All( Optimization_Step ), Obj_Func_tmp )             

                           else if( Type_ObjectiveFunction=='mxd' )then
                              Obj_Func_tmp &
                              = Obj_Func_Single_Thermal( Optimization_Step, Loop_Source, Loop_MatParam_T, Loop_Obj_Func ) &
                               /Value_Obj_Func_Normalize_Thermal( Loop_Source, Loop_MatParam_T, Loop_Obj_Func )
            
                              Obj_Func_All( Optimization_Step ) &
                              = max( Obj_Func_All( Optimization_Step ), Obj_Func_tmp ) 
            
                           end if
                        else if( Loop_Obj_Func==Loop_Flux_X )then
                           if( Type_ObjectiveFunction=='sum' )then
                              ! Maximize F = Minimize  -F
                              Obj_Func_tmp &
                              = -Obj_Func_Single_Thermal( Optimization_Step, Loop_Source, Loop_MatParam_T, Loop_Obj_Func )&
                                /Value_Obj_Func_Normalize_Thermal( Loop_Source, Loop_MatParam_T, Loop_Obj_Func )&
                                /dble( Number_Obj_Func_Source )/dble( Number_Obj_Func_MP_T + Number_Obj_Func_MP_V )
 
                              Obj_Func_All( Optimization_Step ) &
                              = Obj_Func_All( Optimization_Step ) +Obj_Func_tmp
            
                           else if( Type_ObjectiveFunction=='max' )then
                              Obj_Func_tmp &
                              = 1.0d0/( Obj_Func_Single_Thermal( Optimization_Step, Loop_Source, Loop_MatParam_T, Loop_Obj_Func ) &
                                        /Value_Obj_Func_Normalize_Thermal( Loop_Source, Loop_MatParam_T, Loop_Obj_Func ) )
            
                              Obj_Func_All( Optimization_Step ) &
                              = max( Obj_Func_All( Optimization_Step ), Obj_Func_tmp ) 
             
                           else if( Type_ObjectiveFunction=='mxd' )then
                              Obj_Func_tmp &
                              = 1.0d0/( Obj_Func_Single_Thermal( Optimization_Step, Loop_Source, Loop_MatParam_T, Loop_Obj_Func ) &
                               /Value_Obj_Func_Normalize_Thermal( Loop_Source, Loop_MatParam_T, Loop_Obj_Func ) )**Power_Obj_Func
            
                              Obj_Func_All( Optimization_Step ) &
                              = max( Obj_Func_All( Optimization_Step ), Obj_Func_tmp )
            
                           end if
   
                        end if
                     end do ! Loop_Obj_Func

                  else if( Type_ObjectiveFunction=='smx' .and. Flag_MultiPhysics_Device==20 .and. Loop_Cloak==1 .and. Loop_Flux_X==2 )then

                     do Loop_MatParam_V= 1, Number_Obj_Func_MP_V ! for Type_ObjectiveFunction=='smx' ONLY
                        do Loop_Obj_Func= 1, Number_Obj_Func
                           if( Loop_Obj_Func==Loop_Cloak )then

                              Obj_Func_All( Optimization_Step ) &
                              = Obj_Func_All( Optimization_Step ) & 
                                +max( Obj_Func_Single_Thermal( Optimization_Step, Loop_Source, Loop_MatParam_T, Loop_Cloak )&
                                      /Value_Obj_Func_Normalize_Thermal( Loop_Source, Loop_MatParam_T, Loop_Cloak ), &
                                      Obj_Func_Single_DC( Optimization_Step, Loop_Source, Loop_MatParam_V, Loop_Cloak )&
                                      /Value_Obj_Func_Normalize_DC( Loop_Source, Loop_MatParam_V, Loop_Cloak ) &
                                    ) 

                           else if( Loop_Obj_Func==Loop_Flux_X )then

                              Obj_Func_All( Optimization_Step ) &
                              = Obj_Func_All( Optimization_Step ) & 
                                +max( 1.0d0 &
                                      !----------------------------------------------------------------------------------------------
                                      /( Obj_Func_Single_Thermal( Optimization_Step, Loop_Source, Loop_MatParam_T, Loop_Flux_X ) &
                                         /Value_Obj_Func_Normalize_Thermal( Loop_Source, Loop_MatParam_T, Loop_Flux_X ) ), &
                                      1.0d0 &
                                      !----------------------------------------------------------------------------------------------
                                      /( Obj_Func_Single_DC( Optimization_Step, Loop_Source, Loop_MatParam_V, Loop_Flux_X ) &
                                         /Value_Obj_Func_Normalize_DC( Loop_Source, Loop_MatParam_V, Loop_Flux_X ) ) &
                                    ) 
                           end if
                        end do ! Loop_Obj_Func
                     end do ! Loop_MatParam_V

                  else if( Type_ObjectiveFunction=='div' .and. Flag_MultiPhysics_Device==20 .and. Loop_Cloak==1 .and. Loop_Flux_X==2 )then

                     do Loop_Obj_Func= 1, Number_Obj_Func
                        if( Loop_Obj_Func==Loop_Cloak )then

                           Obj_Func_All( Optimization_Step ) &
                           = Obj_Func_All( Optimization_Step ) & 
                             +Obj_Func_Single_Thermal( Optimization_Step, Loop_Source, Loop_MatParam_T, Loop_Cloak )&
                               /Value_Obj_Func_Normalize_Thermal( Loop_Source, Loop_MatParam_T, Loop_Cloak ) &
                               /( dble( Number_Obj_Func_Source )*dble( Number_Obj_Func_MP_T + Number_Obj_Func_MP_V ) ) 

                        else if( Loop_Obj_Func==Loop_Flux_X )then

                           Obj_Func_All( Optimization_Step ) &
                           = Obj_Func_All( Optimization_Step ) & 
                             +1.0d0 &
                             !----------------------------------------------------------------------------------------------
                             /( Obj_Func_Single_Thermal( Optimization_Step, Loop_Source, Loop_MatParam_T, Loop_Flux_X ) &
                                /Value_Obj_Func_Normalize_Thermal( Loop_Source, Loop_MatParam_T, Loop_Flux_X ) )**Power_Obj_Func &
                             /( dble( Number_Obj_Func_Source )*dble( Number_Obj_Func_MP_T + Number_Obj_Func_MP_V ) ) 

                        end if
                     end do ! Loop_Obj_Func

                  else if( Type_ObjectiveFunction=='div' .and. Flag_MultiPhysics_Device==20 .and. Loop_Flux_X==1 )then

                     Obj_Func_All( Optimization_Step ) &
                     = Obj_Func_All( Optimization_Step ) & 
                       +1.0d0 &
                       !----------------------------------------------------------------------------------------------
                       /( Obj_Func_Single_Thermal( Optimization_Step, Loop_Source, Loop_MatParam_T, Loop_Flux_X ) &
                          /Value_Obj_Func_Normalize_Thermal( Loop_Source, Loop_MatParam_T, Loop_Flux_X ) )**Power_Obj_Func &
                       /( dble( Number_Obj_Func_Source )*dble( Number_Obj_Func_MP_T + Number_Obj_Func_MP_V ) ) 

                  ! aho hoge
                  else if( Flag_MultiPhysics_Device==30 .and. Number_Obj_Func==3 .and. Flag_Read_LSF=='true' .and. Flag_xmean=='true' .and. &
                           Loop_Cloak==1 .and. Loop_Flux_X==2 .and. Loop_Flux_Y==3 .and. Type_ObjectiveFunction=='sum' )then

                     Obj_Func_All( Optimization_Step ) &
                     = Obj_Func_All( Optimization_Step ) &
                      +Coefficient_Additional_Obj_Func*( & 
                         Obj_Func_Single_Thermal( Optimization_Step, Loop_Source, Loop_MatParam_T, Loop_Cloak )&
                         /Value_Obj_Func_Normalize_Thermal( Loop_Source, Loop_MatParam_T, Loop_Cloak )&
                         /dble( Number_Obj_Func_Source )/dble( Number_Obj_Func_MP_T + Number_Obj_Func_MP_V ) & 
                        +Obj_Func_Single_Thermal( Optimization_Step, Loop_Source, Loop_MatParam_T, Loop_Flux_X )&
                         /Value_Obj_Func_Normalize_Thermal( Loop_Source, Loop_MatParam_T, Loop_Flux_X )&
                         /dble( Number_Obj_Func_Source )/dble( Number_Obj_Func_MP_T + Number_Obj_Func_MP_V ) ) & 
                      +Obj_Func_Single_Thermal( Optimization_Step, Loop_Source, Loop_MatParam_T, Loop_Flux_Y )&
                       /Value_Obj_Func_Normalize_Thermal( Loop_Source, Loop_MatParam_T, Loop_Flux_Y )&
                       /dble( Number_Obj_Func_Source )/dble( Number_Obj_Func_MP_T + Number_Obj_Func_MP_V ) 

                  else if( Flag_MultiPhysics_Device==30 .and. Loop_Cloak==1 .and. Loop_Flux_X==2 .and. Type_ObjectiveFunction=='sum' )then

                     Obj_Func_All( Optimization_Step ) &
                     = Obj_Func_All( Optimization_Step ) &
                      +Obj_Func_Single_Thermal( Optimization_Step, Loop_Source, Loop_MatParam_T, Loop_Cloak )&
                       /Value_Obj_Func_Normalize_Thermal( Loop_Source, Loop_MatParam_T, Loop_Cloak )&
                       /dble( Number_Obj_Func_Source )/dble( Number_Obj_Func_MP_T + Number_Obj_Func_MP_V ) & 
                      +Obj_Func_Single_Thermal( Optimization_Step, Loop_Source, Loop_MatParam_T, Loop_Flux_X )&
                       /Value_Obj_Func_Normalize_Thermal( Loop_Source, Loop_MatParam_T, Loop_Flux_X )&
                       /dble( Number_Obj_Func_Source )/dble( Number_Obj_Func_MP_T + Number_Obj_Func_MP_V ) 

                  else if( Flag_MultiPhysics_Device==30 .and. Type_ObjectiveFunction=='sum' )then

                     do Loop_Obj_Func= 1, Number_Obj_Func
                        if( Loop_Obj_Func==Loop_Flux_Y )then
                           Obj_Func_All( Optimization_Step ) &
                           = Obj_Func_All( Optimization_Step ) &
                            +Coefficient_Additional_Obj_Func &
                             *Obj_Func_Single_Thermal( Optimization_Step, Loop_Source, Loop_MatParam_T, Loop_Obj_Func )&
                             /Value_Obj_Func_Normalize_Thermal( Loop_Source, Loop_MatParam_T, Loop_Obj_Func )&
                             /dble( Number_Obj_Func_Source )/dble( Number_Obj_Func_MP_T + Number_Obj_Func_MP_V ) 
                        else
                           Obj_Func_All( Optimization_Step ) &
                           = Obj_Func_All( Optimization_Step ) &
                            +Obj_Func_Single_Thermal( Optimization_Step, Loop_Source, Loop_MatParam_T, Loop_Obj_Func )&
                             /Value_Obj_Func_Normalize_Thermal( Loop_Source, Loop_MatParam_T, Loop_Obj_Func )&
                             /dble( Number_Obj_Func_Source )/dble( Number_Obj_Func_MP_T + Number_Obj_Func_MP_V ) 
                        end if
                     end do
                 
                  else if( Flag_MultiPhysics_Device==30 .and. Type_ObjectiveFunction=='sdv' .and. &
                           Loop_Cloak==1 .and. Loop_Flux_X==2 .and. Loop_Flux_Y==3 )then

                     Obj_Func_All( Optimization_Step ) &
                     = Obj_Func_Single_Thermal( Optimization_Step, Loop_Source, Loop_MatParam_T, Loop_Cloak )&
                        /Value_Obj_Func_Normalize_Thermal( Loop_Source, Loop_MatParam_T, Loop_Cloak )&
                        /dble( Number_Obj_Func_Source )/dble( Number_Obj_Func_MP_T +Number_Obj_Func_MP_V ) &
                       ! 
                       +( Obj_Func_Single_Thermal( Optimization_Step, Loop_Source, Loop_MatParam_T, Loop_Flux_X )&
                          /Value_Obj_Func_Normalize_Thermal( Loop_Source, Loop_MatParam_T, Loop_Flux_X ) )&
                       /( Obj_Func_Single_Thermal( Optimization_Step, Loop_Source, Loop_MatParam_T, Loop_Flux_Y )&
                          /Value_Obj_Func_Normalize_Thermal( Loop_Source, Loop_MatParam_T, Loop_Flux_Y ) )&
                          /dble( Number_Obj_Func_Source )/dble( Number_Obj_Func_MP_T +Number_Obj_Func_MP_V )  
                 
                  else if( Flag_MultiPhysics_Device==30 .and. Type_ObjectiveFunction=='sdv' .and. &
                           Loop_Flux_X==1 .and. Loop_Flux_Y==2 )then

                     Obj_Func_All( Optimization_Step ) &
                     = ( Obj_Func_Single_Thermal( Optimization_Step, Loop_Source, Loop_MatParam_T, Loop_Flux_X )&
                          /Value_Obj_Func_Normalize_Thermal( Loop_Source, Loop_MatParam_T, Loop_Flux_X ) )&
                       /( Obj_Func_Single_Thermal( Optimization_Step, Loop_Source, Loop_MatParam_T, Loop_Flux_Y )&
                          /Value_Obj_Func_Normalize_Thermal( Loop_Source, Loop_MatParam_T, Loop_Flux_Y ) )&
                          /dble( Number_Obj_Func_Source )/dble( Number_Obj_Func_MP_T +Number_Obj_Func_MP_V )  

                  else if( Flag_MultiPhysics_Device==30 .and. Loop_Flux_X==1 .and. Type_ObjectiveFunction=='exp' )then

                        Obj_Func_All( Optimization_Step ) &
                        = Obj_Func_All( Optimization_Step ) &
                         +exp( &
                          Obj_Func_Single_Thermal( Optimization_Step, Loop_Source, Loop_MatParam_T, Loop_Flux_X )&
                          /Value_Obj_Func_Normalize_Thermal( Loop_Source, Loop_MatParam_T, Loop_Flux_X ) &
                          ) &
                          /dble( Number_Obj_Func_Source )/dble( Number_Obj_Func_MP_T + Number_Obj_Func_MP_V ) 
                  else
                     write(*,*)'Type_ObjectiveFunction=', Type_ObjectiveFunction
                     write(*,*)'Flag_MultiPhysics_Device=', Flag_MultiPhysics_Device
                     call Output_Error( 'Main', 1653 ) 
                  end if
               end do ! Loop_MatParam_T

               do Loop_MatParam_V= 1, Number_Obj_Func_MP_V
                  !if( Type_ObjectiveFunction/='div' .and. Type_ObjectiveFunction/='smx' )then
                  if( Flag_MultiPhysics_Device==0 .or. Flag_MultiPhysics_Device==10 .or. &
                      Flag_MultiPhysics_Device==20 .or. Flag_MultiPhysics_Device==40 )then
                     do Loop_Obj_Func= 1, Number_Obj_Func
   
                        if( Loop_Obj_Func==Loop_Cloak )then
                           if( Type_ObjectiveFunction=='sum' )then
                          
                              Obj_Func_All( Optimization_Step ) &
                              = Obj_Func_All( Optimization_Step ) &
                               +Obj_Func_Single_DC( Optimization_Step, Loop_Source, Loop_MatParam_V, Loop_Obj_Func )&
                               /Value_Obj_Func_Normalize_DC( Loop_Source, Loop_MatParam_V, Loop_Obj_Func )&
                               /dble( Number_Obj_Func_Source )/dble( Number_Obj_Func_MP_T + Number_Obj_Func_MP_V ) 
            
                           else if( Type_ObjectiveFunction=='max' )then
            
                              Obj_Func_All( Optimization_Step ) &
                              = max( &
                                Obj_Func_All( Optimization_Step ), &
                                Obj_Func_Single_DC( Optimization_Step, Loop_Source, Loop_MatParam_V, Loop_Obj_Func ) &
                               /Value_Obj_Func_Normalize_DC( Loop_Source, Loop_MatParam_V, Loop_Obj_Func ) &
                                )
             
                           else if( Type_ObjectiveFunction=='mxd' )then
            
                              Obj_Func_All( Optimization_Step ) &
                              = max( &
                                Obj_Func_All( Optimization_Step ), &
                                Obj_Func_Single_DC( Optimization_Step, Loop_Source, Loop_MatParam_V, Loop_Obj_Func ) &
                               /Value_Obj_Func_Normalize_DC( Loop_Source, Loop_MatParam_V, Loop_Obj_Func ) &
                                )

                           end if
                        else if( Loop_Obj_Func==Loop_Flux_X )then
                           if( Type_ObjectiveFunction=='sum' )then
                       
                              Obj_Func_All( Optimization_Step ) &
                              = Obj_Func_All( Optimization_Step ) &
                               +Obj_Func_Single_DC( Optimization_Step, Loop_Source, Loop_MatParam_V, Loop_Obj_Func )&
                               /Value_Obj_Func_Normalize_DC( Loop_Source, Loop_MatParam_V, Loop_Obj_Func )&
                               /dble( Number_Obj_Func_Source )/dble( Number_Obj_Func_MP_T + Number_Obj_Func_MP_V ) 
            
                           else if( Type_ObjectiveFunction=='max' )then
            
                              Obj_Func_All( Optimization_Step ) &
                              = max( &
                                Obj_Func_All( Optimization_Step ), &
                                1.0d0/( Obj_Func_Single_DC( Optimization_Step, Loop_Source, Loop_MatParam_V, Loop_Obj_Func ) &
                               /Value_Obj_Func_Normalize_DC( Loop_Source, Loop_MatParam_V, Loop_Obj_Func ) ) &
                                )
             
                           else if( Type_ObjectiveFunction=='mxd' )then
            
                              Obj_Func_All( Optimization_Step ) &
                              = max( &
                                Obj_Func_All( Optimization_Step ), &
                                abs( ( Obj_Func_Single_DC( Optimization_Step, Loop_Source, Loop_MatParam_V, Loop_Obj_Func ) &
                               /Value_Obj_Func_Normalize_DC( Loop_Source, Loop_MatParam_V, Loop_Obj_Func ) ) ) &
                                )
                              != max( &
                              !  Obj_Func_All( Optimization_Step ), &
                              !  1.0d0/( Obj_Func_Single_DC( Optimization_Step, Loop_Source, Loop_MatParam_V, Loop_Obj_Func ) &
                              ! /Value_Obj_Func_Normalize_DC( Loop_Source, Loop_MatParam_V, Loop_Obj_Func ) )**Power_Obj_Func &
                              !  )

                           end if
   
                        end if
                     end do ! Loop_Obj_Func

                  else if( Type_ObjectiveFunction=='div' .and. Flag_MultiPhysics_Device==20  .and. Loop_Cloak==1 .and. Loop_Flux_X==2 )then
                     do Loop_Obj_Func= 1, Number_Obj_Func
                        if( Loop_Obj_Func==Loop_Cloak )then
                           Obj_Func_All( Optimization_Step ) &
                           = Obj_Func_All( Optimization_Step ) & 
                             + Obj_Func_Single_DC( Optimization_Step, Loop_Source, Loop_MatParam_V, Loop_Cloak )&
                               /Value_Obj_Func_Normalize_DC( Loop_Source, Loop_MatParam_V, Loop_Cloak ) &
                               /( dble( Number_Obj_Func_Source )*dble( Number_Obj_Func_MP_T + Number_Obj_Func_MP_V ) ) 

                        else if( Loop_Obj_Func==Loop_Flux_X )then
                           Obj_Func_All( Optimization_Step ) &
                           = Obj_Func_All( Optimization_Step ) & 
                             +1.0d0 &
                             !----------------------------------------------------------------------------------------------
                             /( Obj_Func_Single_DC( Optimization_Step, Loop_Source, Loop_MatParam_V, Loop_Flux_X ) &
                                /Value_Obj_Func_Normalize_DC( Loop_Source, Loop_MatParam_V, Loop_Flux_X ) )**Power_Obj_Func &
                             /( dble( Number_Obj_Func_Source )*dble( Number_Obj_Func_MP_T + Number_Obj_Func_MP_V ) ) 

                        end if
                     end do ! Loop_Obj_Func

                  !else if( Type_ObjectiveFunction=='div' .and. Flag_MultiPhysics_Device==20  .and. Loop_Cloak==1 .and. Loop_Flux_X==2 )then
                  !   Obj_Func_All( Optimization_Step ) &
                  !   = Obj_Func_All( Optimization_Step ) & 
                  !     +( Obj_Func_Single_DC( Optimization_Step, Loop_Source, Loop_MatParam_V, Loop_Cloak )&
                  !       /Value_Obj_Func_Normalize_DC( Loop_Source, Loop_MatParam_V, Loop_Cloak ) ) &
                  !     !----------------------------------------------------------------------------------------------
                  !     /( Obj_Func_Single_DC( Optimization_Step, Loop_Source, Loop_MatParam_V, Loop_Flux_X ) &
                  !        /Value_Obj_Func_Normalize_DC( Loop_Source, Loop_MatParam_V, Loop_Flux_X ) ) &
                  !     /( dble( Number_Obj_Func_Source )*dble( Number_Obj_Func_MP_T + Number_Obj_Func_MP_V ) ) 

                  else if( Type_ObjectiveFunction=='div' .and. Flag_MultiPhysics_Device==20  .and. Loop_Flux_X==1 )then
                     Obj_Func_All( Optimization_Step ) &
                     = Obj_Func_All( Optimization_Step ) & 
                       +1.0d0 &
                       !----------------------------------------------------------------------------------------------
                       /( Obj_Func_Single_DC( Optimization_Step, Loop_Source, Loop_MatParam_V, Loop_Flux_X ) &
                          /Value_Obj_Func_Normalize_DC( Loop_Source, Loop_MatParam_V, Loop_Flux_X ) )**Power_Obj_Func &
                       /( dble( Number_Obj_Func_Source )*dble( Number_Obj_Func_MP_T + Number_Obj_Func_MP_V ) ) 

                  ! aho hoge
                  else if( Flag_MultiPhysics_Device==30 .and. Number_Obj_Func==3 .and. Flag_Read_LSF=='true' .and. Flag_xmean=='true' .and. &
                           Loop_Cloak==1 .and. Loop_Flux_X==2 .and. Loop_Flux_Y==3 .and. Type_ObjectiveFunction=='sum' )then

                     Obj_Func_All( Optimization_Step ) &
                     = Obj_Func_All( Optimization_Step ) &
                      +Coefficient_Additional_Obj_Func*( & 
                         Obj_Func_Single_DC( Optimization_Step, Loop_Source, Loop_MatParam_V, Loop_Cloak )&
                         /Value_Obj_Func_Normalize_DC( Loop_Source, Loop_MatParam_V, Loop_Cloak )&
                         /dble( Number_Obj_Func_Source )/dble( Number_Obj_Func_MP_T + Number_Obj_Func_MP_V ) & 
                        +Obj_Func_Single_DC( Optimization_Step, Loop_Source, Loop_MatParam_V, Loop_Flux_X )&
                         /Value_Obj_Func_Normalize_DC( Loop_Source, Loop_MatParam_V, Loop_Flux_X )&
                         /dble( Number_Obj_Func_Source )/dble( Number_Obj_Func_MP_T + Number_Obj_Func_MP_V ) ) & 
                      +Obj_Func_Single_DC( Optimization_Step, Loop_Source, Loop_MatParam_V, Loop_Flux_Y )&
                       /Value_Obj_Func_Normalize_DC( Loop_Source, Loop_MatParam_V, Loop_Flux_Y )&
                       /dble( Number_Obj_Func_Source )/dble( Number_Obj_Func_MP_T + Number_Obj_Func_MP_V ) 

                  else if( Flag_MultiPhysics_Device==30 .and. Type_ObjectiveFunction=='sum' )then
                     
                     do Loop_Obj_Func= 1, Number_Obj_Func
                        if( Loop_Obj_Func==Loop_Flux_Y )then
                           Obj_Func_All( Optimization_Step ) &
                           = Obj_Func_All( Optimization_Step ) &
                            +Coefficient_Additional_Obj_Func &
                             *Obj_Func_Single_DC( Optimization_Step, Loop_Source, Loop_MatParam_V, Loop_Obj_Func )&
                             /Value_Obj_Func_Normalize_DC( Loop_Source, Loop_MatParam_V, Loop_Obj_Func )&
                             /dble( Number_Obj_Func_Source )/dble( Number_Obj_Func_MP_T + Number_Obj_Func_MP_V ) 
                        else
                           Obj_Func_All( Optimization_Step ) &
                           = Obj_Func_All( Optimization_Step ) &
                            +Obj_Func_Single_DC( Optimization_Step, Loop_Source, Loop_MatParam_V, Loop_Obj_Func )&
                             /Value_Obj_Func_Normalize_DC( Loop_Source, Loop_MatParam_V, Loop_Obj_Func )&
                             /dble( Number_Obj_Func_Source )/dble( Number_Obj_Func_MP_T + Number_Obj_Func_MP_V ) 
                        end if
                     end do

                  else if( Flag_MultiPhysics_Device==30 .and. Type_ObjectiveFunction=='sdv' .and. &
                           Loop_Cloak==1 .and. Loop_Flux_X==2 .and. Loop_Flux_Y==3 )then

                     Obj_Func_All( Optimization_Step ) &
                     = Obj_Func_Single_DC( Optimization_Step, Loop_Source, Loop_MatParam_V, Loop_Cloak )&
                        /Value_Obj_Func_Normalize_DC( Loop_Source, Loop_MatParam_V, Loop_Cloak )&
                        /dble( Number_Obj_Func_Source )/dble( Number_Obj_Func_MP_T +Number_Obj_Func_MP_V ) &
                       ! 
                       +( Obj_Func_Single_DC( Optimization_Step, Loop_Source, Loop_MatParam_V, Loop_Flux_X )&
                          /Value_Obj_Func_Normalize_DC( Loop_Source, Loop_MatParam_V, Loop_Flux_X ) )&
                       /( Obj_Func_Single_DC( Optimization_Step, Loop_Source, Loop_MatParam_V, Loop_Flux_Y )&
                          /Value_Obj_Func_Normalize_DC( Loop_Source, Loop_MatParam_V, Loop_Flux_Y ) )&
                          /dble( Number_Obj_Func_Source )/dble( Number_Obj_Func_MP_T +Number_Obj_Func_MP_V )  

                  else if( Flag_MultiPhysics_Device==30 .and. Type_ObjectiveFunction=='sdv' .and. &
                           Loop_Flux_X==1 .and. Loop_Flux_Y==2 )then

                     Obj_Func_All( Optimization_Step ) &
                     = ( Obj_Func_Single_DC( Optimization_Step, Loop_Source, Loop_MatParam_V, Loop_Flux_X )&
                          /Value_Obj_Func_Normalize_DC( Loop_Source, Loop_MatParam_V, Loop_Flux_X ) )&
                       /( Obj_Func_Single_DC( Optimization_Step, Loop_Source, Loop_MatParam_V, Loop_Flux_Y )&
                          /Value_Obj_Func_Normalize_DC( Loop_Source, Loop_MatParam_V, Loop_Flux_Y ) )&
                          /dble( Number_Obj_Func_Source )/dble( Number_Obj_Func_MP_T +Number_Obj_Func_MP_V )  

                  else if( Flag_MultiPhysics_Device==30 .and. Loop_Flux_X==1 .and. Type_ObjectiveFunction=='exp' )then

                        Obj_Func_All( Optimization_Step ) &
                        = Obj_Func_All( Optimization_Step ) &
                         +exp( &
                          Obj_Func_Single_DC( Optimization_Step, Loop_Source, Loop_MatParam_V, Loop_Flux_X )&
                          /Value_Obj_Func_Normalize_DC( Loop_Source, Loop_MatParam_V, Loop_Flux_X )&
                          ) &
                          /dble( Number_Obj_Func_Source )/dble( Number_Obj_Func_MP_T + Number_Obj_Func_MP_V ) 

                  else if( Type_ObjectiveFunction/='smx' )then
                     write(*,*)'Type_ObjectiveFunction=', Type_ObjectiveFunction
                     write(*,*)'Flag_MultiPhysics_Device=', Flag_MultiPhysics_Device
                     call Output_Error( 'Main', 1800 ) 
                  end if
               end do ! Loop_MatParam_V

            !========================================================================
            else if( Flag_Physics==1 )then
            !========================================================================

               do Loop_MatParam_T= 1, Number_Obj_Func_MP_T
                  !if( Type_ObjectiveFunction/='div' .and. Type_ObjectiveFunction/='smx' )then
                  if( Flag_MultiPhysics_Device==0 .or. Flag_MultiPhysics_Device==10 )then
                     do Loop_Obj_Func= 1, Number_Obj_Func
    
                        if( Loop_Obj_Func==Loop_Cloak )then
                           if( Type_ObjectiveFunction=='sum' )then
                          
                              Obj_Func_All( Optimization_Step ) &
                              = Obj_Func_All( Optimization_Step ) &
                               +Obj_Func_Single_Thermal( Optimization_Step, Loop_Source, Loop_MatParam_T, Loop_Obj_Func )&
                               /Value_Obj_Func_Normalize_Thermal( Loop_Source, Loop_MatParam_T, Loop_Obj_Func )&
                               /Number_Obj_Func_Source/Number_Obj_Func_MP_T
 
                           else if( Type_ObjectiveFunction=='max' )then
            
                              Obj_Func_All( Optimization_Step ) &
                              = max( &
                                Obj_Func_All( Optimization_Step ), &
                                Obj_Func_Single_Thermal( Optimization_Step, Loop_Source, Loop_MatParam_T, Loop_Obj_Func ) &
                               /Value_Obj_Func_Normalize_Thermal( Loop_Source, Loop_MatParam_T, Loop_Obj_Func ) &
                                )

                           else                   
                              write(*,*)'Type_ObjectiveFunction=', Type_ObjectiveFunction
                              call Output_Error( 'Main', 1361 ) 
                           end if
                        else if( Loop_Obj_Func==Loop_Flux_X )then
                           if( Type_ObjectiveFunction=='sum' )then
                          
                              Obj_Func_All( Optimization_Step ) &
                              = Obj_Func_All( Optimization_Step ) &
                               +Obj_Func_Single_Thermal( Optimization_Step, Loop_Source, Loop_MatParam_T, Loop_Obj_Func )&
                               /Value_Obj_Func_Normalize_Thermal( Loop_Source, Loop_MatParam_T, Loop_Obj_Func )&
                               /Number_Obj_Func_Source/Number_Obj_Func_MP_T
            
                           else if( Type_ObjectiveFunction=='mxd' )then
            
                              Obj_Func_All( Optimization_Step ) &
                              = max( &
                                Obj_Func_All( Optimization_Step ), &
                                1.0d0/( Obj_Func_Single_Thermal( Optimization_Step, Loop_Source, Loop_MatParam_T, Loop_Obj_Func ) &
                               /Value_Obj_Func_Normalize_Thermal( Loop_Source, Loop_MatParam_T, Loop_Obj_Func ) )**Power_Obj_Func &
                                )
 
                           else                   
                              write(*,*)'Type_ObjectiveFunction=', Type_ObjectiveFunction
                              call Output_Error( 'Main', 1374 ) 
                           end if
                        end if
                     end do ! Loop_Obj_Func

                  else if( Type_ObjectiveFunction=='amd' .and. Flag_MultiPhysics_Device==20 .and. Loop_Cloak==1 .and. Loop_Flux_X==2 )then

                     do Loop_MatParam_V= 1, Number_Obj_Func_MP_V ! for Type_ObjectiveFunction=='smx' ONLY
                        do Loop_Obj_Func= 1, Number_Obj_Func
                           if( Loop_Obj_Func==Loop_Cloak )then

                              Obj_Func_All( Optimization_Step ) &
                              = max( Obj_Func_All( Optimization_Step ), &
                                     Obj_Func_Single_Thermal( Optimization_Step, Loop_Source, Loop_MatParam_T, Loop_Cloak )&
                                     /Value_Obj_Func_Normalize_Thermal( Loop_Source, Loop_MatParam_T, Loop_Cloak ) &
                                    ) 

                           else if( Loop_Obj_Func==Loop_Flux_X )then

                              Obj_Func_All( Optimization_Step ) &
                              = max( Obj_Func_All( Optimization_Step ), &
                                      1.0d0 &
                                      !----------------------------------------------------------------------------------------------
                                      /( Obj_Func_Single_Thermal( Optimization_Step, Loop_Source, Loop_MatParam_T, Loop_Flux_X ) &
                                         /Value_Obj_Func_Normalize_Thermal( Loop_Source, Loop_MatParam_T, Loop_Flux_X ) )**Power_Obj_Func &
                                    ) 
                           end if
                        end do ! Loop_Obj_Func
                     end do ! Loop_MatParam_V

                  else if( Type_ObjectiveFunction=='div' .and. Flag_MultiPhysics_Device==20 .and. Loop_Cloak==1 .and. Loop_Flux_X==2 )then

                     do Loop_Obj_Func= 1, Number_Obj_Func
                        if( Loop_Obj_Func==Loop_Cloak )then

                           Obj_Func_All( Optimization_Step ) &
                           = Obj_Func_All( Optimization_Step ) & 
                             +Obj_Func_Single_Thermal( Optimization_Step, Loop_Source, Loop_MatParam_T, Loop_Cloak )&
                               /Value_Obj_Func_Normalize_Thermal( Loop_Source, Loop_MatParam_T, Loop_Cloak ) &
                               /( dble( Number_Obj_Func_Source )*dble( Number_Obj_Func_MP_T ) ) 

                        else if( Loop_Obj_Func==Loop_Flux_X )then

                           Obj_Func_All( Optimization_Step ) &
                           = Obj_Func_All( Optimization_Step ) & 
                             +1.0d0 &
                             !----------------------------------------------------------------------------------------------
                             /( Obj_Func_Single_Thermal( Optimization_Step, Loop_Source, Loop_MatParam_T, Loop_Flux_X ) &
                                /Value_Obj_Func_Normalize_Thermal( Loop_Source, Loop_MatParam_T, Loop_Flux_X ) )**Power_Obj_Func & !hoge
                             /( dble( Number_Obj_Func_Source )*dble( Number_Obj_Func_MP_T ) ) 

                        end if
                     end do ! Loop_Obj_Func

                  !else if( Type_ObjectiveFunction=='div' .and. Flag_MultiPhysics_Device==20 .and. Loop_Cloak==1 .and. Loop_Flux_X==2 )then
                  !   Obj_Func_All( Optimization_Step ) &
                  !   = Obj_Func_All( Optimization_Step ) & 
                  !     +( Obj_Func_Single_Thermal( Optimization_Step, Loop_Source, Loop_MatParam_T, Loop_Cloak )&
                  !       /Value_Obj_Func_Normalize_Thermal( Loop_Source, Loop_MatParam_T, Loop_Cloak ) ) &
                  !     !----------------------------------------------------------------------------------------------
                  !     /( Obj_Func_Single_Thermal( Optimization_Step, Loop_Source, Loop_MatParam_T, Loop_Flux_X ) &
                  !        /Value_Obj_Func_Normalize_Thermal( Loop_Source, Loop_MatParam_T, Loop_Flux_X ) ) &
                  !     /( dble( Number_Obj_Func_Source )*dble( Number_Obj_Func_MP_T ) ) 

                  else if( Type_ObjectiveFunction=='div' .and. Flag_MultiPhysics_Device==20 .and. Loop_Flux_X==1 )then
                     Obj_Func_All( Optimization_Step ) &
                     = Obj_Func_All( Optimization_Step ) & 
                       +1.0d0 &
                       !----------------------------------------------------------------------------------------------
                       /( Obj_Func_Single_Thermal( Optimization_Step, Loop_Source, Loop_MatParam_T, Loop_Flux_X ) &
                          /Value_Obj_Func_Normalize_Thermal( Loop_Source, Loop_MatParam_T, Loop_Flux_X ) )**Power_Obj_Func &
                       /( dble( Number_Obj_Func_Source )*dble( Number_Obj_Func_MP_T ) ) 

                  ! aho hoge
                  else if( Flag_MultiPhysics_Device==30 .and. Number_Obj_Func==3 .and. Flag_Read_LSF=='true' .and. Flag_xmean=='true' .and. &
                           Loop_Cloak==1 .and. Loop_Flux_X==2 .and. Loop_Flux_Y==3 .and. Type_ObjectiveFunction=='sum' )then

                     Obj_Func_All( Optimization_Step ) &
                     = Obj_Func_All( Optimization_Step ) &
                      +Coefficient_Additional_Obj_Func*( & 
                         Obj_Func_Single_Thermal( Optimization_Step, Loop_Source, Loop_MatParam_T, Loop_Cloak )&
                         /Value_Obj_Func_Normalize_Thermal( Loop_Source, Loop_MatParam_T, Loop_Cloak )&
                         /dble( Number_Obj_Func_Source )/dble( Number_Obj_Func_MP_T ) & 
                        +Obj_Func_Single_Thermal( Optimization_Step, Loop_Source, Loop_MatParam_T, Loop_Flux_X )&
                         /Value_Obj_Func_Normalize_Thermal( Loop_Source, Loop_MatParam_T, Loop_Flux_X )&
                         /dble( Number_Obj_Func_Source )/dble( Number_Obj_Func_MP_T ) ) & 
                      ! diferent term
                      +Obj_Func_Single_Thermal( Optimization_Step, Loop_Source, Loop_MatParam_T, Loop_Flux_Y )&
                       /Value_Obj_Func_Normalize_Thermal( Loop_Source, Loop_MatParam_T, Loop_Flux_Y )&
                       /dble( Number_Obj_Func_Source )/dble( Number_Obj_Func_MP_T ) 

                  ! aho hoge
                  else if( Flag_MultiPhysics_Device==30 .and. Number_Obj_Func==3 .and. Flag_Read_LSF=='false' .and. &
                           Loop_Cloak==1 .and. Loop_Flux_X==2 .and. Loop_Flux_Y==3 .and. Type_ObjectiveFunction=='sum' )then

                     Obj_Func_All( Optimization_Step ) &
                     = Obj_Func_All( Optimization_Step ) &
                      ! diferent term
                      +Obj_Func_Single_Thermal( Optimization_Step, Loop_Source, Loop_MatParam_T, Loop_Cloak )&
                       /Value_Obj_Func_Normalize_Thermal( Loop_Source, Loop_MatParam_T, Loop_Cloak )&
                       /dble( Number_Obj_Func_Source )/dble( Number_Obj_Func_MP_T ) & 
                      ! diferent term
                      +Coefficient_Additional_Obj_Func* & 
                       Obj_Func_Single_Thermal( Optimization_Step, Loop_Source, Loop_MatParam_T, Loop_Flux_X )&
                       /Value_Obj_Func_Normalize_Thermal( Loop_Source, Loop_MatParam_T, Loop_Flux_X )&
                       /dble( Number_Obj_Func_Source )/dble( Number_Obj_Func_MP_T ) & 
                      ! diferent term
                      +Coefficient_Additional_Obj_Func_2* & 
                      +Obj_Func_Single_Thermal( Optimization_Step, Loop_Source, Loop_MatParam_T, Loop_Flux_Y )&
                       /Value_Obj_Func_Normalize_Thermal( Loop_Source, Loop_MatParam_T, Loop_Flux_Y )&
                       /dble( Number_Obj_Func_Source )/dble( Number_Obj_Func_MP_T ) 

                  else if( Flag_MultiPhysics_Device==30 .and. Type_ObjectiveFunction=='sum' )then

                     do Loop_Obj_Func= 1, Number_Obj_Func
                        if( Loop_Obj_Func==Loop_Flux_Y )then
                           Obj_Func_All( Optimization_Step ) &
                           = Obj_Func_All( Optimization_Step ) &
                            +Coefficient_Additional_Obj_Func &
                             *Obj_Func_Single_Thermal( Optimization_Step, Loop_Source, Loop_MatParam_T, Loop_Obj_Func )&
                             /Value_Obj_Func_Normalize_Thermal( Loop_Source, Loop_MatParam_T, Loop_Obj_Func )&
                             /dble( Number_Obj_Func_Source )/dble( Number_Obj_Func_MP_T ) 
                        else
                           Obj_Func_All( Optimization_Step ) &
                           = Obj_Func_All( Optimization_Step ) &
                            +Obj_Func_Single_Thermal( Optimization_Step, Loop_Source, Loop_MatParam_T, Loop_Obj_Func )&
                             /Value_Obj_Func_Normalize_Thermal( Loop_Source, Loop_MatParam_T, Loop_Obj_Func )&
                             /dble( Number_Obj_Func_Source )/dble( Number_Obj_Func_MP_T ) 
                        end if
                     end do

                  else if( Flag_MultiPhysics_Device==30 .and. Type_ObjectiveFunction=='sdv' .and. &
                           Loop_Cloak==1 .and. Loop_Flux_X==2 .and. Loop_Flux_Y==3 )then

                     Obj_Func_All( Optimization_Step ) &
                     = Obj_Func_Single_Thermal( Optimization_Step, Loop_Source, Loop_MatParam_T, Loop_Cloak )&
                        /Value_Obj_Func_Normalize_Thermal( Loop_Source, Loop_MatParam_T, Loop_Cloak )&
                        /dble( Number_Obj_Func_Source )/dble( Number_Obj_Func_MP_T ) &
                       ! 
                       +( Obj_Func_Single_Thermal( Optimization_Step, Loop_Source, Loop_MatParam_T, Loop_Flux_X )&
                          /Value_Obj_Func_Normalize_Thermal( Loop_Source, Loop_MatParam_T, Loop_Flux_X ) )&
                       /( Obj_Func_Single_Thermal( Optimization_Step, Loop_Source, Loop_MatParam_T, Loop_Flux_Y )&
                          /Value_Obj_Func_Normalize_Thermal( Loop_Source, Loop_MatParam_T, Loop_Flux_Y ) )&
                          /dble( Number_Obj_Func_Source )/dble( Number_Obj_Func_MP_T )  

                  else if( Flag_MultiPhysics_Device==30 .and. Type_ObjectiveFunction=='sdv' .and. &
                           Loop_Flux_X==1 .and. Loop_Flux_Y==2 )then

                     Obj_Func_All( Optimization_Step ) &
                     = ( Obj_Func_Single_Thermal( Optimization_Step, Loop_Source, Loop_MatParam_T, Loop_Flux_X )&
                          /Value_Obj_Func_Normalize_Thermal( Loop_Source, Loop_MatParam_T, Loop_Flux_X ) )&
                       /( Obj_Func_Single_Thermal( Optimization_Step, Loop_Source, Loop_MatParam_T, Loop_Flux_Y )&
                          /Value_Obj_Func_Normalize_Thermal( Loop_Source, Loop_MatParam_T, Loop_Flux_Y ) )&
                          /dble( Number_Obj_Func_Source )/dble( Number_Obj_Func_MP_T )  

                  else if( Flag_MultiPhysics_Device==30 .and. Loop_Flux_X==1 .and. Type_ObjectiveFunction=='exp' )then

                        Obj_Func_All( Optimization_Step ) &
                        = Obj_Func_All( Optimization_Step ) &
                         +exp( &
                          Obj_Func_Single_Thermal( Optimization_Step, Loop_Source, Loop_MatParam_T, Loop_Flux_X )&
                          /Value_Obj_Func_Normalize_Thermal( Loop_Source, Loop_MatParam_T, Loop_Flux_X )&
                          ) &
                          /dble( Number_Obj_Func_Source )/dble( Number_Obj_Func_MP_T ) 

                  else
                     write(*,*)'Type_ObjectiveFunction=', Type_ObjectiveFunction
                     write(*,*)'Flag_MultiPhysics_Device=', Flag_MultiPhysics_Device
                     call Output_Error( 'Main', 1959 ) 
                  end if

               end do ! Loop_MatParam_T

            !========================================================================
            else if( Flag_Physics==2 )then
            !========================================================================
               do Loop_MatParam_V= 1, Number_Obj_Func_MP_V
                  !if( Type_ObjectiveFunction/='div' )then
                  if( Flag_MultiPhysics_Device==0 .or. Flag_MultiPhysics_Device==10 )then
                     do Loop_Obj_Func= 1, Number_Obj_Func
    
                        if( Loop_Obj_Func==Loop_Cloak )then
                           if( Type_ObjectiveFunction=='sum' )then
                          
                              Obj_Func_All( Optimization_Step ) &
                              = Obj_Func_All( Optimization_Step ) &
                               +Obj_Func_Single_DC( Optimization_Step, Loop_Source, Loop_MatParam_V, Loop_Obj_Func )&
                               /Value_Obj_Func_Normalize_DC( Loop_Source, Loop_MatParam_V, Loop_Obj_Func )&
                               /Number_Obj_Func_Source/Number_Obj_Func_MP_V
                           else                   
                              write(*,*)'Type_ObjectiveFunction=', Type_ObjectiveFunction
                              call Output_Error( 'Main', 1395 ) 
                           end if
                        else if( Loop_Obj_Func==Loop_Flux_X )then
                           if( Type_ObjectiveFunction=='sum' )then
                          
                              Obj_Func_All( Optimization_Step ) &
                              = Obj_Func_All( Optimization_Step ) &
                               +Obj_Func_Single_DC( Optimization_Step, Loop_Source, Loop_MatParam_V, Loop_Obj_Func )&
                               /Value_Obj_Func_Normalize_DC( Loop_Source, Loop_MatParam_V, Loop_Obj_Func )&
                               /Number_Obj_Func_Source/Number_Obj_Func_MP_V
             
                           else if( Type_ObjectiveFunction=='mxd' )then
            
                              Obj_Func_All( Optimization_Step ) &
                              = max( &
                                Obj_Func_All( Optimization_Step ), &
                                1.0d0/( Obj_Func_Single_DC( Optimization_Step, Loop_Source, Loop_MatParam_V, Loop_Obj_Func ) &
                               /Value_Obj_Func_Normalize_DC( Loop_Source, Loop_MatParam_V, Loop_Obj_Func ) )**Power_Obj_Func &
                                )
 
                           else                   
                              write(*,*)'Type_ObjectiveFunction=', Type_ObjectiveFunction
                              call Output_Error( 'Main', 1407 ) 
                           end if
   
                        end if
                     end do ! Loop_Obj_Func

                  else if( Type_ObjectiveFunction=='div' .and. Flag_MultiPhysics_Device==20  .and. Loop_Cloak==1 .and. Loop_Flux_X==2 )then

                     do Loop_Obj_Func= 1, Number_Obj_Func
                        if( Loop_Obj_Func==Loop_Cloak )then
                           Obj_Func_All( Optimization_Step ) &
                           = Obj_Func_All( Optimization_Step ) & 
                             + Obj_Func_Single_DC( Optimization_Step, Loop_Source, Loop_MatParam_V, Loop_Cloak )&
                               /Value_Obj_Func_Normalize_DC( Loop_Source, Loop_MatParam_V, Loop_Cloak ) &
                               /( dble( Number_Obj_Func_Source )*dble( Number_Obj_Func_MP_V ) ) 

                        else if( Loop_Obj_Func==Loop_Flux_X )then
                           Obj_Func_All( Optimization_Step ) &
                           = Obj_Func_All( Optimization_Step ) & 
                             +1.0d0 &
                             !----------------------------------------------------------------------------------------------
                             /( Obj_Func_Single_DC( Optimization_Step, Loop_Source, Loop_MatParam_V, Loop_Flux_X ) &
                                /Value_Obj_Func_Normalize_DC( Loop_Source, Loop_MatParam_V, Loop_Flux_X ) )**Power_Obj_Func &
                             /( dble( Number_Obj_Func_Source )*dble( Number_Obj_Func_MP_V ) ) 

                        end if
                     end do ! Loop_Obj_Func

                  else if( Type_ObjectiveFunction=='div' .and. Flag_MultiPhysics_Device==20 .and. Loop_Flux_X==1 )then

                     Obj_Func_All( Optimization_Step ) &
                     = Obj_Func_All( Optimization_Step ) & 
                       +1.0d0 &
                       !----------------------------------------------------------------------------------------------
                       /( Obj_Func_Single_DC( Optimization_Step, Loop_Source, Loop_MatParam_V, Loop_Flux_X ) &
                          /Value_Obj_Func_Normalize_DC( Loop_Source, Loop_MatParam_V, Loop_Flux_X ) )**Power_Obj_Func &
                       /( dble( Number_Obj_Func_Source )*dble( Number_Obj_Func_MP_V ) ) 

                  ! aho hoge
                  else if( Flag_MultiPhysics_Device==30 .and. Number_Obj_Func==3 .and. Flag_Read_LSF=='true' .and. Flag_xmean=='true' .and. &
                           Loop_Cloak==1 .and. Loop_Flux_X==2 .and. Loop_Flux_Y==3 .and. Type_ObjectiveFunction=='sum' )then

                     Obj_Func_All( Optimization_Step ) &
                     = Obj_Func_All( Optimization_Step ) &
                      +Coefficient_Additional_Obj_Func*( & 
                         Obj_Func_Single_DC( Optimization_Step, Loop_Source, Loop_MatParam_V, Loop_Cloak )&
                         /Value_Obj_Func_Normalize_DC( Loop_Source, Loop_MatParam_V, Loop_Cloak )&
                         /dble( Number_Obj_Func_Source )/dble( Number_Obj_Func_MP_V ) & 
                        +Obj_Func_Single_DC( Optimization_Step, Loop_Source, Loop_MatParam_V, Loop_Flux_X )&
                         /Value_Obj_Func_Normalize_DC( Loop_Source, Loop_MatParam_V, Loop_Flux_X )&
                         /dble( Number_Obj_Func_Source )/dble( Number_Obj_Func_MP_V ) ) & 
                      +Obj_Func_Single_DC( Optimization_Step, Loop_Source, Loop_MatParam_V, Loop_Flux_Y )&
                       /Value_Obj_Func_Normalize_DC( Loop_Source, Loop_MatParam_V, Loop_Flux_Y )&
                       /dble( Number_Obj_Func_Source )/dble( Number_Obj_Func_MP_V ) 

                  else if( Flag_MultiPhysics_Device==30 .and. Type_ObjectiveFunction=='sum' )then

                     do Loop_Obj_Func= 1, Number_Obj_Func
                        if( Loop_Obj_Func==Loop_Flux_Y )then
                           Obj_Func_All( Optimization_Step ) &
                           = Obj_Func_All( Optimization_Step ) &
                            +Coefficient_Additional_Obj_Func &
                             *Obj_Func_Single_DC( Optimization_Step, Loop_Source, Loop_MatParam_V, Loop_Obj_Func )&
                             /Value_Obj_Func_Normalize_DC( Loop_Source, Loop_MatParam_V, Loop_Obj_Func )&
                             /dble( Number_Obj_Func_Source )/dble( Number_Obj_Func_MP_V ) 
                        else
                           Obj_Func_All( Optimization_Step ) &
                           = Obj_Func_All( Optimization_Step ) &
                            +Obj_Func_Single_DC( Optimization_Step, Loop_Source, Loop_MatParam_V, Loop_Obj_Func )&
                             /Value_Obj_Func_Normalize_DC( Loop_Source, Loop_MatParam_V, Loop_Obj_Func )&
                             /dble( Number_Obj_Func_Source )/dble( Number_Obj_Func_MP_V ) 
                        end if
                     end do

                  else if( Flag_MultiPhysics_Device==30 .and. Type_ObjectiveFunction=='sdv' .and. &
                           Loop_Cloak==1 .and. Loop_Flux_X==2 .and. Loop_Flux_Y==3 )then

                     Obj_Func_All( Optimization_Step ) &
                     = Obj_Func_Single_DC( Optimization_Step, Loop_Source, Loop_MatParam_V, Loop_Cloak )&
                        /Value_Obj_Func_Normalize_DC( Loop_Source, Loop_MatParam_V, Loop_Cloak )&
                        /dble( Number_Obj_Func_Source )/dble( Number_Obj_Func_MP_V ) &
                       ! 
                       +( Obj_Func_Single_DC( Optimization_Step, Loop_Source, Loop_MatParam_V, Loop_Flux_X )&
                          /Value_Obj_Func_Normalize_DC( Loop_Source, Loop_MatParam_V, Loop_Flux_X ) )&
                       /( Obj_Func_Single_DC( Optimization_Step, Loop_Source, Loop_MatParam_V, Loop_Flux_Y )&
                          /Value_Obj_Func_Normalize_DC( Loop_Source, Loop_MatParam_V, Loop_Flux_Y ) )&
                          /dble( Number_Obj_Func_Source )/dble( Number_Obj_Func_MP_V )  

                  else if( Flag_MultiPhysics_Device==30 .and. Type_ObjectiveFunction=='sdv' .and. &
                           Loop_Flux_X==1 .and. Loop_Flux_Y==2 )then

                     Obj_Func_All( Optimization_Step ) &
                     = ( Obj_Func_Single_DC( Optimization_Step, Loop_Source, Loop_MatParam_V, Loop_Flux_X )&
                          /Value_Obj_Func_Normalize_DC( Loop_Source, Loop_MatParam_V, Loop_Flux_X ) )&
                       /( Obj_Func_Single_DC( Optimization_Step, Loop_Source, Loop_MatParam_V, Loop_Flux_Y )&
                          /Value_Obj_Func_Normalize_DC( Loop_Source, Loop_MatParam_V, Loop_Flux_Y ) )&
                          /dble( Number_Obj_Func_Source )/dble( Number_Obj_Func_MP_V )  


                  else if( Flag_MultiPhysics_Device==30 .and. Loop_Flux_X==1 .and. Type_ObjectiveFunction=='exp' )then

                        Obj_Func_All( Optimization_Step ) &
                        = Obj_Func_All( Optimization_Step ) &
                         +exp( &
                          Obj_Func_Single_DC( Optimization_Step, Loop_Source, Loop_MatParam_V, Loop_Flux_X )&
                          /Value_Obj_Func_Normalize_DC( Loop_Source, Loop_MatParam_V, Loop_Flux_X )&
                          ) &
                          /dble( Number_Obj_Func_Source )/dble( Number_Obj_Func_MP_V ) 

                  else
                     write(*,*)'Type_ObjectiveFunction=', Type_ObjectiveFunction
                     write(*,*)'Flag_MultiPhysics_Device=', Flag_MultiPhysics_Device
                     call Output_Error( 'Main', 2085 ) 
                  end if

               end do ! Loop_MatParam_V
            else                   
               write(*,*)'Flag_Physics=', Flag_Physics
               call Output_Error( 'Main', 1415 ) 
            end if

            write(*,*)'================================================='
            write(*,*)'Obj_Func_All=', Obj_Func_All( Optimization_Step )
            write(*,*)'================================================='
   
            if( Flag_Physics==0 .or. Flag_Physics==1 )then
               do Loop_MatParam_T= 1, Number_Obj_Func_MP_T
                  do Loop_Obj_Func= 1, Number_Obj_Func
 
                     write(*,*)'Normalized Obj_Func_Single_Thermal=', &
                     Obj_Func_Single_Thermal( Optimization_Step, Loop_Source, Loop_MatParam_T, Loop_Obj_Func )&
                     /Value_Obj_Func_Normalize_Thermal( Loop_Source, Loop_MatParam_T, Loop_Obj_Func )
                     write(*,*)'Obj_Func_Single_Thermal=', &
                     Obj_Func_Single_Thermal( Optimization_Step, Loop_Source, Loop_MatParam_T, Loop_Obj_Func )
                     write(*,*)'Value_Obj_Func_Normalize_Thermal=', &
                     Value_Obj_Func_Normalize_Thermal( Loop_Source, Loop_MatParam_T, Loop_Obj_Func ) 

                  end do ! Loop_Obj_Func
               end do ! Loop_MatParam_V
            end if
   
            if( Flag_Physics==0 .or. Flag_Physics==2 )then
               do Loop_MatParam_V= 1, Number_Obj_Func_MP_V
                  do Loop_Obj_Func= 1, Number_Obj_Func
 
                     write(*,*)'Normalized Obj_Func_Single_DC=', &
                     Obj_Func_Single_DC( Optimization_Step, Loop_Source, Loop_MatParam_V, Loop_Obj_Func )&
                     /Value_Obj_Func_Normalize_DC( Loop_Source, Loop_MatParam_V, Loop_Obj_Func )
                     write(*,*)'Obj_Func_Single_DC=', &
                     Obj_Func_Single_DC( Optimization_Step, Loop_Source, Loop_MatParam_V, Loop_Obj_Func )
                     write(*,*)'Value_Obj_Func_Normalize_DC=', &
                     Value_Obj_Func_Normalize_DC( Loop_Source, Loop_MatParam_V, Loop_Obj_Func )

                  end do ! Loop_Obj_Func
               end do ! Loop_MatParam_V
            end if
            write(*,*)'================================================='
 
         end do ! Loop_Source

         if( Flag_InitialConfiguration==0 .and. & 
             ( Flag_MultiPhysics_Device==0 .or. Flag_MultiPhysics_Device==10 ) .and. &
             Radius_FixedDomain > 1d-8 .and. &
             Optimization_Step==0 .and. Loop_Solution==1 )then

            if( Flag_Physics==0 .or. Flag_Physics==1 )then
               call Output_Value_Normalization_Thermal&
                    ( Obj_Func_Single_Thermal, Position_X_Source, Position_Y_Source, Relative_Thermal_Conductivity, Loop_Cloak ) 
            end if
   
            if( Flag_Physics==0 .or. Flag_Physics==2 )then
               call Output_Value_Normalization_DC&
                    ( Obj_Func_Single_DC, Position_X_Source, Position_Y_Source, Relative_Electrical_Conductivity, Loop_Cloak ) 
            end if

            !Convergence_Ratio_xmean_dowhile=0.0d0
            Termination_dowhile='y'
         end if

         if( Flag_InitialConfiguration/=200 .and. Flag_MultiPhysics_Device==20 .and. & ! 40 for grep
             Optimization_Step==0 .and. Loop_Solution==1 )then

            if( Flag_InitialConfiguration==0 )then
               if( Flag_Physics==0 .or. Flag_Physics==1 )then
                  call Output_Value_Normalization_Thermal&
                       ( Obj_Func_Single_Thermal, Position_X_Source, Position_Y_Source, Relative_Thermal_Conductivity, Loop_Flux_X ) 
               end if
      
               if( Flag_Physics==0 .or. Flag_Physics==2 )then
                  call Output_Value_Normalization_DC&
                       ( Obj_Func_Single_DC, Position_X_Source, Position_Y_Source, Relative_Electrical_Conductivity, Loop_Flux_X ) 
               end if
            else  if( Flag_InitialConfiguration==5 )then
               if( Flag_Physics==0 .or. Flag_Physics==1 )then
                  call Output_Value_Normalization_Thermal&
                       ( Obj_Func_Single_Thermal, Position_X_Source, Position_Y_Source, Relative_Thermal_Conductivity, Loop_Cloak ) 
               end if
      
               if( Flag_Physics==0 .or. Flag_Physics==2 )then
                  call Output_Value_Normalization_DC&
                       ( Obj_Func_Single_DC, Position_X_Source, Position_Y_Source, Relative_Electrical_Conductivity, Loop_Cloak ) 
               end if
            end if

            !Convergence_Ratio_xmean_dowhile=0.0d0
            Termination_dowhile='y'
         end if

         if( Flag_InitialConfiguration/=200 .and. Flag_MultiPhysics_Device==30 .and. &
             Optimization_Step==0 .and. Loop_Solution==1 )then

            if( Flag_InitialConfiguration==0 )then
               if( Flag_Physics==0 .or. Flag_Physics==1 )then
                  call Output_Value_Normalization_Thermal&
                       ( Obj_Func_Single_Thermal, Position_X_Source, Position_Y_Source, Relative_Thermal_Conductivity, Loop_Flux_X ) 

                  if( Loop_Flux_Y<=Number_Obj_Func ) &
                  call Output_Value_Normalization_Thermal&
                       ( Obj_Func_Single_Thermal, Position_X_Source, Position_Y_Source, Relative_Thermal_Conductivity, Loop_Flux_Y ) 
               end if
      
               if( Flag_Physics==0 .or. Flag_Physics==2 )then
                  call Output_Value_Normalization_DC&
                       ( Obj_Func_Single_DC, Position_X_Source, Position_Y_Source, Relative_Electrical_Conductivity, Loop_Flux_X ) 

                  if( Loop_Flux_Y<=Number_Obj_Func ) &
                  call Output_Value_Normalization_DC&
                       ( Obj_Func_Single_DC, Position_X_Source, Position_Y_Source, Relative_Electrical_Conductivity, Loop_Flux_Y ) 
               end if
            else  if( Flag_InitialConfiguration==5 )then
               if( Flag_Physics==0 .or. Flag_Physics==1 )then
                  call Output_Value_Normalization_Thermal&
                       ( Obj_Func_Single_Thermal, Position_X_Source, Position_Y_Source, Relative_Thermal_Conductivity, Loop_Cloak ) 
               end if
      
               if( Flag_Physics==0 .or. Flag_Physics==2 )then
                  call Output_Value_Normalization_DC&
                       ( Obj_Func_Single_DC, Position_X_Source, Position_Y_Source, Relative_Electrical_Conductivity, Loop_Cloak ) 
               end if
            end if

            !Convergence_Ratio_xmean_dowhile=0.0d0
            Termination_dowhile='y'
         end if

         if( Flag_InitialConfiguration/=200 .and. Flag_MultiPhysics_Device==40 .and. & 
             Optimization_Step==0 .and. Loop_Solution==1 )then

            if( Flag_InitialConfiguration==0 )then
               if( match( Thermal_Conductivity_Base_Material, 1.0d0 )=='y' )then 
                  if( Flag_Physics==0 .or. Flag_Physics==1 )then
                     call Output_Value_Normalization_Thermal&
                          ( Obj_Func_Single_Thermal, Position_X_Source, Position_Y_Source, Relative_Thermal_Conductivity, Loop_Flux_X ) 
                  end if
               else
                  if( Flag_Physics==0 .or. Flag_Physics==1 )then
                     call Output_Value_Normalization_Thermal&
                          ( Obj_Func_Single_Thermal, Position_X_Source, Position_Y_Source, Relative_Thermal_Conductivity, Loop_Cloak ) 
                  end if
               end if
      
               if( match( Electric_Conductivity_Base_Material, 1.0d0 )=='y' )then 
                  if( Flag_Physics==0 .or. Flag_Physics==2 )then
                     call Output_Value_Normalization_DC&
                          ( Obj_Func_Single_DC, Position_X_Source, Position_Y_Source, Relative_Electrical_Conductivity, Loop_Flux_X ) 
                  end if
               else
                  if( Flag_Physics==0 .or. Flag_Physics==2 )then
                     call Output_Value_Normalization_DC&
                          ( Obj_Func_Single_DC, Position_X_Source, Position_Y_Source, Relative_Electrical_Conductivity, Loop_Cloak ) 
                  end if
               end if
            end if

            !Convergence_Ratio_xmean_dowhile=0.0d0
            Termination_dowhile='y'
         end if



         !========================================================================
         if( MyRank==0 ) write(*,*) '   Preserve Individual Data'
         !========================================================================

         Fractal_MPI( Loop_Solution )= Fractal_Structure 
         Perimeter_MPI( Loop_Solution )= Perimeter_Structure

         Volume_MPI( Loop_Solution )= Volume_Material/Volume_DesignDomain

         Obj_Func_MPI( Loop_Solution )= Obj_Func_All( Optimization_Step )
         
         if( Volume_MPI( Loop_Solution )*Sign_VolumeConstraint > Ratio_Volume_Constraint*Sign_VolumeConstraint )then
            Penalty_Volume_Constraint_MPI( Loop_Solution ) &
            = exp( Sign_VolumeConstraint*Coefficient_Decay_Volume_Constraint &
                *( Volume_MPI( Loop_Solution ) -Ratio_Volume_Constraint ) ) &
             -1d0 ! exp( 0d0 ) 
         else
            Penalty_Volume_Constraint_MPI( Loop_Solution )= 0d0
         end if

         !======================================
         if( Type_Fitness=='sum' )then
         !======================================
            if( Flag_Volume_Constraint==1 .and. Flag_Perimeter_Constraint==1 )then
               Fitness_MPI( Loop_Solution )&
               = Obj_Func_MPI( Loop_Solution ) &
               + Fractal_MPI( Loop_Solution )*Coefficient_Complexity_Tau_Normal & 
               + Penalty_Volume_Constraint_MPI( Loop_Solution ) 
   
            else if( Flag_Volume_Constraint==0 .and. Flag_Perimeter_Constraint==1 )then
               Fitness_MPI( Loop_Solution )&
               = Obj_Func_MPI( Loop_Solution ) &
               + Fractal_MPI( Loop_Solution )*Coefficient_Complexity_Tau_Normal

            else if( Flag_Volume_Constraint==1 .and. Flag_Perimeter_Constraint==2 )then
               Fitness_MPI( Loop_Solution )&
               = Obj_Func_MPI( Loop_Solution ) &
               + Perimeter_Implicit_MPI( Loop_Solution )*Coefficient_Complexity_Tau_Normal & 
               + Penalty_Volume_Constraint_MPI( Loop_Solution ) 
   
            else if( Flag_Volume_Constraint==0 .and. Flag_Perimeter_Constraint==2 )then
               Fitness_MPI( Loop_Solution )&
               = Obj_Func_MPI( Loop_Solution ) &
               + Perimeter_Implicit_MPI( Loop_Solution )*Coefficient_Complexity_Tau_Normal 
            else
               call Output_Error( 'Main', 1219 ) 
            end if
         !======================================
         else if( Type_Fitness=='max' )then
         !======================================
            if( Flag_Volume_Constraint==1 .and. Flag_Perimeter_Constraint==1 )then
               Fitness_MPI( Loop_Solution )&
               = max( &
                 Obj_Func_MPI( Loop_Solution ), &
                 Fractal_MPI( Loop_Solution )*Coefficient_Complexity_Tau_Normal, & 
                 Penalty_Volume_Constraint_MPI( Loop_Solution ) &
                 ) 
   
            else if( Flag_Volume_Constraint==0 .and. Flag_Perimeter_Constraint==1 )then
               Fitness_MPI( Loop_Solution )&
               = max( &
                 Obj_Func_MPI( Loop_Solution ),  &
                 Fractal_MPI( Loop_Solution )*Coefficient_Complexity_Tau_Normal & 
                 ) 
   
            else if( Flag_Volume_Constraint==1 .and. Flag_Perimeter_Constraint==2 )then
               Fitness_MPI( Loop_Solution )&
               = max( &
                 Obj_Func_MPI( Loop_Solution ), &
                 Perimeter_Implicit_MPI( Loop_Solution )*Coefficient_Complexity_Tau_Normal, & 
                 Penalty_Volume_Constraint_MPI( Loop_Solution ) &
                 ) 
   
            else if( Flag_Volume_Constraint==0 .and. Flag_Perimeter_Constraint==2 )then
               Fitness_MPI( Loop_Solution )&
               = max( &
                 Obj_Func_MPI( Loop_Solution ), &
                 Perimeter_Implicit_MPI( Loop_Solution )*Coefficient_Complexity_Tau_Normal &
                 )
            else
              call Output_Error( 'Main', 1737 ) 
            end if
        !=============================== 
        else if(Type_Fitness=='tra' )then
        !================================
               Fractal_trade=( Fractal_MPI( Loop_Solution ) - Fractal_i )/( Fractal_a - Fractal_i )

               Obj_Func_trade=(Obj_Func_MPI(Loop_Solution)-Obj_Func_i)/(Obj_Func_a-Obj_Func_i)

               !Perimeter_trade=(Perimeter_MPI(Loop_Solution)-Perimeter_i)/(Perimeter_a-Perimeter_i)               

               Fitness_MPI(Loop_Solution)&
               =max( Obj_Func_trade, Fractal_trade )

         !======================================
         else
         !======================================
            call Output_Error( 'Main', 1719 ) 
         end if

         if( MyRank==0 ) write(*,*)'==============================================================================' 
         if( MyRank==0 ) write(*,*)'Obj_Func_MPI=', Obj_Func_MPI( Loop_Solution ) 
         if( MyRank==0 ) write(*,*)'Perimeter_MPI=', Perimeter_MPI( Loop_Solution )
         if( MyRank==0 .and. Flag_Perimeter_Constraint==2 ) write(*,*)'Perimeter_Implicit_MPI=', Perimeter_Implicit_MPI( Loop_Solution )
         if( MyRank==0 ) write(*,*)'Volume_MPI=', Volume_MPI( Loop_Solution ) 
         if( MyRank==0 ) write(*,*)'Fitness_MPI=', Fitness_MPI( Loop_Solution ) 
         if( MyRank==0 ) write(*,*)'Optimization_Step=', Optimization_Step 
         if( MyRank==0 ) write(*,*) 'Loop_Solution=', Loop_Solution -Loop_Start +1, '/', Number_Sampling -Loop_Start +1
         if( MyRank==0 ) write(*,*)'==============================================================================' 

         if( Flag_Physics==1 )then
            do Loop_Source= 1, Number_Obj_Func_Source
               do Loop_MatParam_T= 1, Number_Obj_Func_MP_T
                  do Loop_Obj_Func= 1, Number_Obj_Func
                     Obj_Func_Multi( Loop_Obj_Func, Loop_MatParam_T, Loop_Source, 1, Loop_Solution )&
                     = Obj_Func_Single_Thermal( Optimization_Step, Loop_Source, Loop_MatParam_T, Loop_Obj_Func )&
                       /Value_Obj_Func_Normalize_Thermal( Loop_Source, Loop_MatParam_T, Loop_Obj_Func ) 
                  end do
               end do
            end do
         else if( Flag_Physics==2 )then
            do Loop_Source= 1, Number_Obj_Func_Source
               do Loop_MatParam_V= 1, Number_Obj_Func_MP_V
                  do Loop_Obj_Func= 1, Number_Obj_Func
                     Obj_Func_Multi( Loop_Obj_Func, Loop_MatParam_V, Loop_Source, 1, Loop_Solution )&
                     = Obj_Func_Single_DC( Optimization_Step, Loop_Source, Loop_MatParam_V, Loop_Obj_Func )&
                       /Value_Obj_Func_Normalize_DC( Loop_Source, Loop_MatParam_V, Loop_Obj_Func ) 
                  end do
               end do
            end do
         else if( Flag_Physics==0 )then
            do Loop_Source= 1, Number_Obj_Func_Source
               do Loop_MatParam_T= 1, Number_Obj_Func_MP_T
                  do Loop_Obj_Func= 1, Number_Obj_Func
                     Obj_Func_Multi( Loop_Obj_Func, Loop_MatParam_T, Loop_Source, 1, Loop_Solution )&
                     = Obj_Func_Single_Thermal( Optimization_Step, Loop_Source, Loop_MatParam_T, Loop_Obj_Func )&
                       /Value_Obj_Func_Normalize_Thermal( Loop_Source, Loop_MatParam_T, Loop_Obj_Func ) 
                  end do
               end do
            end do 
            do Loop_Source= 1, Number_Obj_Func_Source
               do Loop_MatParam_V= 1,  Number_Obj_Func_MP_V
                  do Loop_Obj_Func= 1, Number_Obj_Func
                     Obj_Func_Multi( Loop_Obj_Func, Loop_MatParam_V, Loop_Source, 2, Loop_Solution )&
                     = Obj_Func_Single_DC( Optimization_Step, Loop_Source, Loop_MatParam_V, Loop_Obj_Func )&
                       /Value_Obj_Func_Normalize_DC( Loop_Source, Loop_MatParam_V, Loop_Obj_Func ) 
                  end do
               end do
            end do
         else
             write(*,*)'Flag_Physics=', Flag_Physics
             call Output_Error( 'Main', 1430 )
         end if

!write(*,*)'stop at Main.f90 2173'
!stop

         !========================================================================
         if( MyRank==0 ) write(*,*) '   Preserve FEM Data'
         !========================================================================
   
         if( Loop_Solution==Loop_Num_MPI( 1 ) )then
            Perimeter_Structure_Minimum = Perimeter_Structure
            Volume_Structure_Minimum= Volume_Structure_All( Loop_Solution )
            Obj_Func_Real_Minimum= Obj_Func_Real( Loop_Solution ) 
            Individual_MPI= Loop_Solution

            allocate( Field_Plot_1( femdata%number_node, Number_Obj_Func_Source, Number_Obj_Func_MP_T ) )
            allocate( Field_Plot_2( femdata%number_node, Number_Obj_Func_Source, Number_Obj_Func_MP_V ) )

            if( Flag_MultiPhysics_Device==0 .or. Flag_MultiPhysics_Device==10 .or. &
                Flag_MultiPhysics_Device==20 .or. Flag_MultiPhysics_Device==30 .or. Flag_MultiPhysics_Device==40 )then
               if( Flag_Physics==0 .or. Flag_Physics==1 ) Field_Plot_1= Temperature_Solution
               if( Flag_Physics==0 .or. Flag_Physics==2 ) Field_Plot_2= Electric_Potential_Solution
            else
               write(*,*)'Flag_MultiPhysics_Device=', Flag_MultiPhysics_Device
               call Output_Error( 'Main', 1300 ) 
            end if

            call Current_Data_2_Preserve_Data &
               ( femdata%number_node, &
                 femdata%number_element_triangle, & 
                 Number_Node_on_Electrical_Insulation_Boundary, Number_Node_in_Electrical_Insulation, Number_Element_Electrical_Insulation_Boundary, &
                 Max_Number_Element_Share_OneNode, &
                 femdata%position_node_renumber, femdata%type_node, &
                 femdata%index_element_triangle_2_node, &
                 femdata%type_element_triangle, &
                 Element_and_LocalNode_Number_on_Electrical_Insulation_BC, Max_Position_Node, &
                 Number_Obj_Func_Source, Number_Obj_Func_MP_T, Number_Obj_Func_MP_V, &
                 Field_Plot_1, Field_Plot_2,  & 
               !==========================================================================
                 Number_Node_Preserve, &
                 Number_Element_Triangle_Preserve, &
                 Number_Node_on_Electrical_Insulation_Boundary_Preserve, Number_Node_in_Electrical_Insulation_Preserve, &
                 Number_Element_Electrical_Insulation_Boundary_Preserve, &
                 Max_Number_Element_Share_OneNode_Preserve, &
                 Position_Node_Preserve, Class_Node_Preserve, &  
                 Index_Element_2_Node_Triangle_Preserve, &
                 Class_Element_Triangle_Preserve, &
                 Element_and_LocalNode_Num_on_Electrical_Insulation_BC_Preserve, Max_Position_Node_Preserve, &
                 Value_Plot_Preserve_Thermal, Value_Plot_Preserve_DC )

            deallocate( Field_Plot_1 )
            deallocate( Field_Plot_2 )

         else
            do i= Loop_Num_MPI( 1 ), Loop_Solution-1
               if( Fitness_MPI( i )*Sign_ObjectiveFunction < Fitness_MPI( Loop_Solution )*Sign_ObjectiveFunction )then 
                  goto 1141 
               end if
            end do

            Perimeter_Structure_Minimum= Perimeter_Structure
            Volume_Structure_Minimum= Volume_Structure_All( Loop_Solution )
            Obj_Func_Real_Minimum= Obj_Func_Real( Loop_Solution ) 
            Individual_MPI= Loop_Solution

            allocate( Field_Plot_1( femdata%number_node, Number_Obj_Func_Source, Number_Obj_Func_MP_T ) )
            allocate( Field_Plot_2( femdata%number_node, Number_Obj_Func_Source, Number_Obj_Func_MP_V ) )

            if( Flag_MultiPhysics_Device==0 .or. Flag_MultiPhysics_Device==10 .or. &
                Flag_MultiPhysics_Device==20 .or. Flag_MultiPhysics_Device==30 .or. Flag_MultiPhysics_Device==40 )then
               if( Flag_Physics==0 .or. Flag_Physics==1 ) Field_Plot_1= Temperature_Solution
               if( Flag_Physics==0 .or. Flag_Physics==2 ) Field_Plot_2= Electric_Potential_Solution
            else
               write(*,*)'Flag_MultiPhysics_Device=', Flag_MultiPhysics_Device
               call Output_Error( 'Main', 1300 ) 
            end if

            call Current_Data_2_Preserve_Data &
               ( femdata%number_node, &
                 femdata%number_element_triangle, & 
                 Number_Node_on_Electrical_Insulation_Boundary, Number_Node_in_Electrical_Insulation, Number_Element_Electrical_Insulation_Boundary, &
                 Max_Number_Element_Share_OneNode, &
                 femdata%position_node_renumber, femdata%type_node, &
                 femdata%index_element_triangle_2_node, &
                 femdata%type_element_triangle, &
                 Element_and_LocalNode_Number_on_Electrical_Insulation_BC, Max_Position_Node, &
                 Number_Obj_Func_Source, Number_Obj_Func_MP_T, Number_Obj_Func_MP_V, &
                 Field_Plot_1, Field_Plot_2, & 
               !==========================================================================
                 Number_Node_Preserve, &
                 Number_Element_Triangle_Preserve, &
                 Number_Node_on_Electrical_Insulation_Boundary_Preserve, Number_Node_in_Electrical_Insulation_Preserve, &
                 Number_Element_Electrical_Insulation_Boundary_Preserve, &
                 Max_Number_Element_Share_OneNode_Preserve, &
                 Position_Node_Preserve, Class_Node_Preserve, &  
                 Index_Element_2_Node_Triangle_Preserve, &
                 Class_Element_Triangle_Preserve, &
                 Element_and_LocalNode_Num_on_Electrical_Insulation_BC_Preserve, Max_Position_Node_Preserve, &
                 Value_Plot_Preserve_Thermal, Value_Plot_Preserve_DC )

            deallocate( Field_Plot_1 )
            deallocate( Field_Plot_2 )

            1141 continue
         end if

         deallocate( LSF_GridPoint )
         deallocate( Temperature_Solution )
         deallocate( Electric_Potential_Solution )
         deallocate( femdata%position_node_renumber )
         deallocate( Index_GridPointNumber_2_NodeNumber_Renumbered )
 
         deallocate( femdata%type_node )
         deallocate( femdata%index_element_triangle_2_node )
         deallocate( femdata%type_element_triangle )
         deallocate( Max_Position_Node )
 
         deallocate( Element_and_LocalNode_Number_on_Electrical_Insulation_BC )

         if( Flag_ObjectiveFunction==11 .and. abs( Radius_FixedDomain)<=1d-8 .and. Flag_InitialConfiguration==0 )then
            !====================================================================================
            write(*,*)'   Output FEM_Data_0 for Flat plane, Main.f90 1050'
            write(*,*)'   See ./data/'
            !====================================================================================
            if( Flag_Output_Distribution_Result/=999 .and. Flag_Output_PostScript/=999 ) stop
         end if

      end do ! Loop_Solution

      call MPI_Barrier( MPI_COMM_WORLD, ierr )

      !===========================================================================
      if( MyRank==0 ) write(*,*)'   MPI Data Communication : Line 1640' 
      !===========================================================================

write(*,*)'MPI_Allgather: Fractal_MPI --> Fractal_Structure_All'
      call MPI_Allgather&
         ( Fractal_MPI( Loop_Num_MPI( 1 ) ), Number_Loop_MPI, MPI_DOUBLE_PRECISION, &
           Fractal_Structure_All( 1 ), Number_Loop_MPI, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierr)

write(*,*)'MPI_Allgather: Perimeter_MPI --> Perimeter_Structure_All'
      call MPI_Allgather&
         ( Perimeter_MPI( Loop_Num_MPI( 1 ) ), Number_Loop_MPI, MPI_DOUBLE_PRECISION, &
           Perimeter_Structure_All( 1 ), Number_Loop_MPI, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierr)

      if( Flag_Perimeter_Constraint==2 )then
         call MPI_Allgather&
            ( Perimeter_Implicit_MPI( Loop_Num_MPI( 1 ) ), Number_Loop_MPI, MPI_DOUBLE_PRECISION, &
              Perimeter_Implicit_All( 1 ), Number_Loop_MPI, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierr)
      end if

      call MPI_Allgather&
         ( Volume_MPI( Loop_Num_MPI( 1 ) ), Number_Loop_MPI, MPI_DOUBLE_PRECISION, &
           Volume_Structure_All( 1 ), Number_Loop_MPI, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierr)

      call MPI_Allgather&
         ( Penalty_Volume_Constraint_MPI( Loop_Num_MPI( 1 ) ), Number_Loop_MPI, MPI_DOUBLE_PRECISION, &
           Penalty_Volume_Constraint( 1 ), Number_Loop_MPI, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierr)

      call MPI_Allgather&
         ( Obj_Func_MPI( Loop_Num_MPI( 1 ) ), Number_Loop_MPI, MPI_DOUBLE_PRECISION, &
           Obj_Func_Real( 1 ), Number_Loop_MPI, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierr)

      call MPI_Allgather&
         ( Fitness_MPI( Loop_Num_MPI( 1 ) ), Number_Loop_MPI, MPI_DOUBLE_PRECISION, &
           Fitness_CMAES( 1 ), Number_Loop_MPI, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierr)

      call MPI_Allgather&
         ( Average_Diff_LSF_PC_MPI( Loop_Num_MPI( 1 ) ), Number_Loop_MPI, MPI_DOUBLE_PRECISION, &
           Average_Diff_LSF_PC( 1 ), Number_Loop_MPI, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierr)

      !========================================================================================================
      !allocate( tmp%itg1( 3 ) )

      !tmp%itg1( 1 )= Number_Design_Variable*Number_Sampling
      !tmp%itg1( 2 )= Number_Design_Variable*( Loop_Num_MPI( 2 )-Loop_Num_MPI( 1 )+1 )
      !tmp%itg1( 3 )= Number_Design_Variable*( Loop_Num_MPI( 1 ) -1 ) +1

      !allocate( Vector_Data_MPI( tmp%itg1( 1 ) ) )

      !call Matrix_2_Vector&
      !   ( Design_Variable, Number_Design_Variable, Number_Sampling, &
      !     Vector_Data_MPI, tmp%itg1( 1 ) )

      !allocate( Vector_Data( tmp%itg1( 1 ) ) )

      !call MPI_Allgather&
      !   ( Vector_Data_MPI( tmp%itg1( 3 ) ), tmp%itg1( 2 ), MPI_DOUBLE_PRECISION, &
      !     Vector_Data( 1 ), tmp%itg1( 2 ), MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierr)

      !call Vector_2_Matrix&
      !   ( Vector_Data, tmp%itg1( 1 ), & 
      !     Design_Variable, Number_Design_Variable, Number_Sampling )

      !deallocate( Vector_Data_MPI )
      !deallocate( Vector_Data )
      !deallocate( tmp%itg1 )

      !========================================================================================================
      allocate( tmp%itg1( 3 ) )
      tmp%itg1( 1 )= Number_Obj_Func*Number_Obj_Func_MaterialParameter*Number_Obj_Func_Source*Number_Physics*Number_Sampling
      tmp%itg1( 2 )= Number_Obj_Func*Number_Obj_Func_MaterialParameter*Number_Obj_Func_Source*Number_Physics*( Loop_Num_MPI( 2 )-Loop_Num_MPI( 1 )+1 )
      tmp%itg1( 3 )= Number_Obj_Func*Number_Obj_Func_MaterialParameter*Number_Obj_Func_Source*Number_Physics*( Loop_Num_MPI( 1 ) -1 ) +1

      allocate( Vector_Data_MPI( tmp%itg1( 1 ) ) )
       
      !call Matrix4_2_Vector_Real&
      call Matrix5_2_Vector_Real&
         ( Obj_Func_Multi, Number_Obj_Func, Number_Obj_Func_MaterialParameter, Number_Obj_Func_Source, Number_Physics, Number_Sampling, &
           Vector_Data_MPI, tmp%itg1( 1 ) )

      allocate( Vector_Data( tmp%itg1( 1 ) ) )

      call MPI_Allgather&
         ( Vector_Data_MPI( tmp%itg1( 3 ) ), tmp%itg1( 2 ), MPI_DOUBLE_PRECISION, &
           Vector_Data( 1 ), tmp%itg1( 2 ), MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierr)

      !call Vector_2_Matrix4_Real&
      call Vector_2_Matrix5_Real&
         ( Vector_Data, tmp%itg1( 1 ), &
           Obj_Func_Multi, Number_Obj_Func, Number_Obj_Func_MaterialParameter, Number_Obj_Func_Source, Number_Physics, Number_Sampling )

      deallocate( Vector_Data_MPI )
      deallocate( Vector_Data )
      deallocate( tmp%itg1 )

      !========================================================================================================
      allocate( tmp%itg1( 3 ) )

      tmp%itg1( 1 )= Number_Design_Variable*Number_Sampling
      tmp%itg1( 2 )= Number_Design_Variable*( Loop_Num_MPI( 2 )-Loop_Num_MPI( 1 )+1 )
      tmp%itg1( 3 )= Number_Design_Variable*( Loop_Num_MPI( 1 ) -1 ) +1

      allocate( Vector_Data_MPI( tmp%itg1( 1 ) ) )

      call Matrix_2_Vector&
         ( Diff_LSF_PC, Number_Design_Variable, Number_Sampling, &
           Vector_Data_MPI, tmp%itg1( 1 ) )

      allocate( Vector_Data( tmp%itg1( 1 ) ) )

      call MPI_Allgather&
         ( Vector_Data_MPI( tmp%itg1( 3 ) ), tmp%itg1( 2 ), MPI_DOUBLE_PRECISION, &
           Vector_Data( 1 ), tmp%itg1( 2 ), MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, ierr)

      call Vector_2_Matrix&
         ( Vector_Data, tmp%itg1( 1 ), & 
           Diff_LSF_PC, Number_Design_Variable, Number_Sampling )

      deallocate( Vector_Data_MPI )
      deallocate( Vector_Data )
      deallocate( tmp%itg1 )

      if( Flag_Check_Value==1 )then
         call Check_Value_Vector( Number_Sampling, Obj_Func_Real, 1.7d308, 'Main', 1201 )
         call Check_Value_Vector( Number_Sampling, Fitness_CMAES, 1.7d308, 'Main', 1202 )
         call Check_Value_Vector( Number_Sampling, Fractal_Structure_All, 1.7d308, 'Main', 1203 )
         call Check_Value_Vector( Number_Sampling, Perimeter_Structure_All, 1.7d308, 'Main', 1203 )
         if( Flag_Perimeter_Constraint==2 ) call Check_Value_Vector( Number_Sampling, Perimeter_Implicit_All, 1.7d308, 'Main', 1204 )
         call Check_Value_Vector( Number_Sampling, Volume_Structure_All, 1.7d308, 'Main', 1205 )
      end if

      deallocate( Loop_Num_MPI )
      deallocate( Fractal_MPI )
      deallocate( Perimeter_MPI )
      deallocate( Perimeter_Implicit_MPI )
      deallocate( Volume_MPI )
      deallocate( Penalty_Volume_Constraint_MPI )
      deallocate( Obj_Func_MPI )
      deallocate( Fitness_MPI )
      deallocate( Average_Diff_LSF_PC_MPI )

      call MPI_Barrier( MPI_COMM_WORLD, ierr )


      if( MyRank==0 .and. Flag_Output_Data_MPI==1 )then
         open( 1382, file='./Obj_Func_ALL.dat', position='append' ) 
            write(1382,*)'Optimization_Step=', Optimization_Step 
            write(1382,*) Obj_Func_Real 
         close( 1382 )

         open( 1383, file='./Fitness_ALL.dat', position='append' ) 
            write(1383,*)'Optimization_Step=', Optimization_Step 
            write(1383,*) Fitness_CMAES 
         close( 1383 )

         open( 1384, file='./Fitness_ALL.dat', position='append' ) 
            write(1384,*)'Optimization_Step=', Optimization_Step 
            write(1384,*) Design_Variable( 1, : ) 
         close( 1384 )
      end if

      !===========================================================================
      if( MyRank==0 ) write(*,*)'Output the Results of Optimization' 
      !===========================================================================

      call Find_Optimal_Individual&
         ( Fitness_CMAES, size( Fitness_CMAES ), &
           Flag_Ordering, Individual_Optimal( Optimization_Step ) )

      Fitness_Convergence( Optimization_Step )= Fitness_CMAES( Individual_Optimal( Optimization_Step ) )

      if( MyRank==0 )then
         if( Device_Number==999 .and. Number_Optimization_Step==1 .and. Number_Sampling==1 )then
            !======================================================================================
            open( 2211, file='cmaes_ndv.f90', status='replace' )
            !======================================================================================
               write(2211,'(a29,1x,i8)') '   integer, parameter :: ndv=', Number_Design_Variable
               write(2211,'(a30,1x,i8)') '   integer, parameter :: nsmp=', Number_Candidates
            close( 2211 )
         end if
         !======================================================================================
         open( 58, file='Objective_Function.dat', position='append' ) 
         !======================================================================================
            if( Optimization_Step==0 )then
               write( 58, * )'# ========== Topology Optimization Parameters =========='
               write( 58, * )'# tau=', Coefficient_Complexity_Tau_Normal, ' Number_Design_Variable=', n_dv 
               write( 58, * )'# Flag_Symmetry_LSF=', Flag_Symmetry_LSF 
               write( 58, * )'# ========== Physical Parameters =========='
               write( 58, * )'# Flag_Physics=', Flag_Physics
               write( 58, * )'# Number_Obj_Func_Source=', Number_Obj_Func_Source 
               write( 58, * )'# Number_Obj_Func_MP_T=', Number_Obj_Func_MP_T ,' Number_Obj_Func_MP_V=', Number_Obj_Func_MP_V
               write( 58, * )'# ========== CMA-ES Parameters =========='
               write( 58, * )'# Flag_Box_Constraint=', Flag_Box_Constraint, ' Radius_Search=', Radius_Search 
               write( 58, * )'# n_cs=', Number_Sampling,' Nproc=', Nproc 
               write( 58, * )'# Flag_Volume_Constraint=', Flag_Volume_Constraint 
                         if( Flag_Volume_Constraint==1 ) write( 58, '(1x,a26,f15.8,2x,a22,f15.8)' )&
                                                         '# Ratio_Volume_Constraint=', Ratio_Volume_Constraint, &
                                                         'Sign_VolumeConstraint=', Sign_VolumeConstraint
               write( 58, * )'# ===================='
               write( 58, '(a111)' )' # Generation,  Obj. Func., Fractal,           Perimeter,       Fitness,        Volume,      cp_sigma, Opt.Sol.Num., Param. ED' 
            end if

            write( 58, 525 ) Optimization_Step, &
                       Obj_Func_Real( Individual_Optimal( Optimization_Step ) ), &
                       Fractal_Structure_All( Individual_Optimal( Optimization_Step ) ), &
                       Perimeter_Structure_All( Individual_Optimal( Optimization_Step ) ),&
                       Fitness_CMAES( Individual_Optimal( Optimization_Step ) ), &
                       Volume_Structure_All( Individual_Optimal( Optimization_Step ) ), &
                       Individual_Optimal( Optimization_Step ), &
                       sigma, flag_eigen
            525 format( 1x, I10, es15.7, es15.7, es15.7, es15.7, f15.7, I10, es15.7, f15.7 )
         close( 58 )
         
        !close(290)

         !if( Number_Obj_Func_Source >= 2 .or. Number_Obj_Func_MP_T >= 2 )then
            !======================================================================================
            open( 1128, file='Multi_OF.dat', position='append' )
            !======================================================================================
               if( Optimization_Step==0 )then
                  write( 1128, 1153 )'# Generation,      OF All,  '
                  1153 format( 1x, a30, $ )
                  do Loop_Physics= 1, Number_Physics
                     if( ( Loop_Physics==1 .and. Flag_Physics==0 ) .or. ( Loop_Physics==1 .and. Flag_Physics==1 ) )then
                        do Loop_MatParam_T= 1, Number_Obj_Func_MP_T
                           do Loop_Source= 1, Number_Obj_Func_Source
                              write( 1128, 1152 )'Src:', Loop_Source, 'Conductivity:', Loop_MatParam_T, 'Phys:', Loop_Physics, ' | '
                              1152 format( a4, i1, 1x, a13, i1, 1x, a5, i1, a3, $ )
                           end do
                        end do
                     else if( ( Loop_Physics==2 .and. Flag_Physics==0 ) .or. ( Loop_Physics==1 .and. Flag_Physics==2 ) )then
                        do Loop_MatParam_V= 1, Number_Obj_Func_MP_V
                           do Loop_Source= 1, Number_Obj_Func_Source
                              write( 1128, 1986 )'Src:', Loop_Source, 'Conductivity:', Loop_MatParam_V, 'Phys:', Loop_Physics, ' | '
                              1986 format( a4, i1, 1x, a13, i1, 1x, a5, i1, a3, $ )
                           end do
                        end do
                     end if
                  end do
                  write( 1128, 1154 )'Perimeter'
                  1154 format( 6x, a9 )
               end if
               write( 1128, 1148 ) Optimization_Step
               1148 format( 1x, I10, $  )
               write( 1128, 1134 ) Obj_Func_Real( Individual_Optimal( Optimization_Step ) )
               do Loop_Physics= 1, Number_Physics
                  if( ( Loop_Physics==1 .and. Flag_Physics==0 ) .or. ( Loop_Physics==1 .and. Flag_Physics==1 ) )then
                     do Loop_MatParam_T= 1, Number_Obj_Func_MP_T
                        do Loop_Source= 1, Number_Obj_Func_Source
                           do Loop_Obj_Func= 1, Number_Obj_Func
                              write( 1128, 1135 ) &
                              Obj_Func_Multi( Loop_Obj_Func, Loop_MatParam_T, Loop_Source, &
                                        Loop_Physics, Individual_Optimal( Optimization_Step ) )
                           end do
                        end do
                     end do
                  else if( ( Loop_Physics==2 .and. Flag_Physics==0 ) .or. ( Loop_Physics==1 .and. Flag_Physics==2 ) )then
                     do Loop_MatParam_V= 1, Number_Obj_Func_MP_V
                        do Loop_Source= 1, Number_Obj_Func_Source
                           do Loop_Obj_Func= 1, Number_Obj_Func
                              write( 1128, 1135 ) &
                              Obj_Func_Multi( Loop_Obj_Func, Loop_MatParam_V, Loop_Source, &
                                        Loop_Physics, Individual_Optimal( Optimization_Step ) )
                           end do
                        end do
                     end do
                  end if 
               end do
               write( 1128, 1134 ) Perimeter_Structure_All( Individual_Optimal( Optimization_Step ) )
               1134 format( 2x, es15.5, $ )
               1135 format( 6x, es15.5, $ )
               write( 1128, * ) ' '
            close( 1128 )
         !end if
         

         open( 458, file=Filename_Structure, position='append' ) 
            write( 458, 459 ) Optimization_Step, Perimeter_Structure_All( Individual_Optimal( Optimization_Step ) ), &
                        Volume_Structure_All( Individual_Optimal( Optimization_Step ) ), &
                        Penalty_Volume_Constraint( Individual_Optimal( Optimization_Step ) ) 
            459 format( I10, 2x, es15.8, 2x, es15.8, 2x, es15.8 )
         close( 458 )

         open( 576, file='./Penalty_function_Uncertainty.dat', position='append')
            if( Optimization_Step==0 )then
               write( 576,'(a118)') &
               'Generation, minval( Average_Diff_LSF_PC ), maxval( Average_Diff_LSF_PC ), &
                minval( Diff_LSF_PC ), maxval( Diff_LSF_PC )'
            end if
            write( 576,582) Optimization_Step, minval( Average_Diff_LSF_PC ), maxval( Average_Diff_LSF_PC ), &
                      minval( abs( Diff_LSF_PC ) ), maxval( abs( Diff_LSF_PC ) )
            582 format( I10, 1x, es15.8, 1x, es15.8, 1x, es15.8, 1x, es15.8 )
         close( 576 )

      end if
  
      call MPI_Barrier( MPI_COMM_WORLD, ierr )
      !===========================================================================
      if( MyRank==0 ) write(*,*)'Output Configuration' 
      !===========================================================================

      call Judge_Output_PSFile &
         ( Flag_Output_PostScript, Optimization_Step, Interval_Plot_Optimization_Step, &
           Individual_MPI, Individual_Optimal( Optimization_Step ), &
           Flag_Output_PSFile )

      if( Flag_Output_PSFile==1 )then

         allocate( Position_Node_Plot( 2, Number_Node_Preserve ) )

         Position_Node_Plot( :, 1:Number_Node_Preserve )= Position_Node_Preserve( :, 1:Number_Node_Preserve )

         allocate( Index_Element_2_Node_Plot( 4, Number_Element_Triangle_Preserve ) )
         allocate( Class_Element_Plot( Number_Element_Triangle_Preserve ) )

         Class_Element_Plot( 1:Number_Element_Triangle_Preserve )&
         = Class_Element_Triangle_Preserve( 1:Number_Element_Triangle_Preserve )
         Index_Element_2_Node_Plot( 1:3, 1:Number_Element_Triangle_Preserve )&
         = Index_Element_2_Node_Triangle_Preserve( 1:3, 1:Number_Element_Triangle_Preserve )
         Index_Element_2_Node_Plot( 4, : )= Index_Element_2_Node_Plot( 3, : )

         call Plot_FEM_Model &
            ( Number_Node_Preserve, Number_Element_Triangle_Preserve, &
              Optimization_Step, &
              Position_Node_Plot, Index_Element_2_Node_Plot, Class_Element_Plot )

         deallocate( Position_Node_Plot )
         deallocate( Index_Element_2_Node_Plot )
         deallocate( Class_Element_Plot )
      end if

      call Judge_Output_PSFile &
         ( Flag_Output_LSF, Optimization_Step, Interval_Plot_Optimization_Step, &
           Individual_MPI, Individual_Optimal( Optimization_Step ), &
           Flag_Output_Interval_LSF )

      if( Flag_Output_Interval_LSF==1 )then
         write( Format_Filenumber, '(A2,I1,A1,I1,A1)' ) & 
         '(i', Length_Character_Optimization_Step,'.',Length_Character_Optimization_Step,')'
        
         Generation_Output= Optimization_Step
         write( Generation_Character, Format_Filenumber ) Generation_Output
      
         FileName_LSF= trim( "LSF_Convergence_"//trim(Generation_Character) )
         open( 3030, file=FileName_LSF, position='append' ) 
            write( 3030, * ) Generation_Output, '# Generation_Output'
            !write( 3030, * ) Threshold_Convergence_xmean( i ), '# Threshold_Convergence_xmean'
            write( 3030, * ) Interval_Plot_Optimization_Step, '# Interval_Plot_Optimization_Step'
            write( 3030, * ) Flag_Symmetry_LSF, '# Flag_Symmetry_LSF'
            write( 3030, * ) Number_GridPoint_Optimization, '# Number_GridPoint_Optimization'
            do j= 1, Number_GridPoint_Optimization
               write( 3030, * ) j, Design_Variable_Output( j ) 
            end do
         close( 3030 )

         !FileName_LSF= trim( "LSF_Plot_"//trim(Generation_Character) )
         !open( 3287, file=FileName_LSF, position='append' ) 
         !   do j= 1, Number_GridPoint_Optimization
         !      write( 3287, * ) j, Design_Variable_Output( j ) 
         !   end do
         !close( 3287 )

      end if

      call Judge_Output_PSFile &
         ( Flag_Output_Distribution_Result, Optimization_Step, Interval_Plot_Optimization_Step, &
           Individual_MPI, Individual_Optimal( Optimization_Step ), &
           Flag_Output_EFD )

      !if( Flag_Output_EFD==1 .and. MyRank==0 )then
      if( Flag_Output_EFD==1 )then
 
         allocate( Position_Node_Plot( 2, Number_Node_Preserve ) )

         Position_Node_Plot( :, 1:Number_Node_Preserve )= Position_Node_Preserve( :, 1:Number_Node_Preserve )
         allocate( Index_Element_2_Node_Plot( 4, Number_Element_Triangle_Preserve ) )

         Index_Element_2_Node_Plot( 1:3, 1:Number_Element_Triangle_Preserve )&
         = Index_Element_2_Node_Triangle_Preserve( 1:3, 1:Number_Element_Triangle_Preserve )
         Index_Element_2_Node_Plot( 4, : )= Index_Element_2_Node_Plot( 3, : )

         !=============================================
         if( Flag_Physics==0 .or. Flag_Physics==1 )then
         !=============================================
            allocate( Value_Plot( Number_Node_Preserve ) ) 

            Value_Plot( 1:Number_Node_Preserve )= Value_Plot_Preserve_Thermal( 1:Number_Node_Preserve, 1, 1 )

            call Plot_Result_Element_Triangle &
               ( Optimization_Step, 1, &
                 Value_Plot, Number_Level_Plot_EM, & 
                 Optimization_Step, 1, 1, &
                 Position_Node_Plot, Number_Node_Preserve,  &
                 Index_Element_2_Node_Plot, Number_Element_Triangle_Preserve, &
                 Obj_Func_Real( Individual_Optimal( Optimization_Step ) ) ) 

            deallocate( Value_Plot ) 
         end if

         !=============================================
         if( Flag_Physics==0 .or. Flag_Physics==2 )then
         !=============================================
            allocate( Value_Plot( Number_Node_Preserve ) ) 

            Value_Plot( 1:Number_Node_Preserve )= Value_Plot_Preserve_DC( 1:Number_Node_Preserve, 1, 1 )

            call Plot_Result_Element_Triangle &
               ( Optimization_Step, 2, &
                 Value_Plot, Number_Level_Plot_EM, & 
                 Optimization_Step, 1, 1, &
                 Position_Node_Plot, Number_Node_Preserve,  &
                 Index_Element_2_Node_Plot, Number_Element_Triangle_Preserve, &
                 Obj_Func_Real( Individual_Optimal( Optimization_Step ) ) ) 

            deallocate( Value_Plot ) 
         end if

         deallocate( Position_Node_Plot )
         deallocate( Index_Element_2_Node_Plot )

      end if

      call MPI_Barrier( MPI_COMM_WORLD, ierr )
      !===========================================================================
      if( MyRank==0 ) write(*,*)'Update & Output FEM Data'
      !===========================================================================

      if( Optimization_Step==0 )then
         Fitness_Output= Fitness_Convergence( Optimization_Step ) 

         Obj_Func_Output= Obj_Func_Real( Individual_Optimal( Optimization_Step ) )
         Perimeter_Output= Perimeter_Structure_All( Individual_Optimal( Optimization_Step ) )

         Generation_Output= Optimization_Step

         MyRank_Output_tmp= 0
         if( Individual_MPI==Individual_Optimal( Optimization_Step ) )then
            MyRank_Output_tmp= MyRank
         end if 
         call MPI_Allreduce( MyRank_Output_tmp, MyRank_Output, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
         call MPI_Barrier( MPI_COMM_WORLD, ierr )

         Design_Variable_Output( : )=Design_Variable( :, Individual_Optimal( Optimization_Step ) ) 

         call Preserve_Data_2_Ouput_Data &
            ( Number_Node_Preserve, &
              Number_Element_Triangle_Preserve, &
              Number_Node_on_Electrical_Insulation_Boundary_Preserve, Number_Node_in_Electrical_Insulation_Preserve, Number_Element_Electrical_Insulation_Boundary_Preserve, &
              Max_Number_Element_Share_OneNode_Preserve, &
              Position_Node_Preserve, Class_Node_Preserve, &  
              Index_Element_2_Node_Triangle_Preserve, &
              Class_Element_Triangle_Preserve, &
              Element_and_LocalNode_Num_on_Electrical_Insulation_BC_Preserve, Max_Position_Node_Preserve, &
              Number_Obj_Func_Source, Number_Obj_Func_MP_T, Number_Obj_Func_MP_V, &
              Value_Plot_Preserve_Thermal, Value_Plot_Preserve_DC, &
              !==========================================================================
              Number_Node_Output, &
              Number_Element_Triangle_Output, &
              Number_Node_on_Electrical_Insulation_Boundary_Output, Number_Node_in_Electrical_Insulation_Output, Number_Element_Electrical_Insulation_Boundary_Output, &
              Max_Number_Element_Share_OneNode_Output, &
              Position_Node_Output, Class_Node_Output, &  
              Index_Element_2_Node_Triangle_Output, &
              Class_Element_Triangle_Output, & 
              Element_and_LocalNode_Number_on_Electrical_Insulation_BC_Output, Max_Position_Node_Output, &
              Value_Plot_Output_Thermal, Value_Plot_Output_DC )

      else if( Fitness_Convergence( Optimization_Step )*Sign_ObjectiveFunction < Fitness_Output*Sign_ObjectiveFunction )then
             Fitness_Output= Fitness_Convergence( Optimization_Step )

         Obj_Func_Output= Obj_Func_Real( Individual_Optimal( Optimization_Step ) )
         Perimeter_Output= Perimeter_Structure_All( Individual_Optimal( Optimization_Step ) )

         Generation_Output= Optimization_Step

         MyRank_Output_tmp= 0
         if( Individual_MPI==Individual_Optimal( Optimization_Step ) )then
            MyRank_Output_tmp= MyRank
         end if 
         call MPI_Allreduce( MyRank_Output_tmp, MyRank_Output, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
         call MPI_Barrier( MPI_COMM_WORLD, ierr )

         call Preserve_Data_2_Ouput_Data &
            ( Number_Node_Preserve, &
              Number_Element_Triangle_Preserve, &
              Number_Node_on_Electrical_Insulation_Boundary_Preserve, Number_Node_in_Electrical_Insulation_Preserve, Number_Element_Electrical_Insulation_Boundary_Preserve, &
              Max_Number_Element_Share_OneNode_Preserve, &
              Position_Node_Preserve, Class_Node_Preserve, &  
              Index_Element_2_Node_Triangle_Preserve, &
              Class_Element_Triangle_Preserve, &
              Element_and_LocalNode_Num_on_Electrical_Insulation_BC_Preserve, Max_Position_Node_Preserve, &
              Number_Obj_Func_Source, Number_Obj_Func_MP_T, Number_Obj_Func_MP_V, &
              Value_Plot_Preserve_Thermal, Value_Plot_Preserve_DC, &
              !==========================================================================
              Number_Node_Output, &
              Number_Element_Triangle_Output, &
              Number_Node_on_Electrical_Insulation_Boundary_Output, Number_Node_in_Electrical_Insulation_Output, Number_Element_Electrical_Insulation_Boundary_Output, &
              Max_Number_Element_Share_OneNode_Output, &
              Position_Node_Output, Class_Node_Output, &  
              Index_Element_2_Node_Triangle_Output, &
              Class_Element_Triangle_Output, &
              Element_and_LocalNode_Number_on_Electrical_Insulation_BC_Output, Max_Position_Node_Output, &
              Value_Plot_Output_Thermal, Value_Plot_Output_DC )

         Design_Variable_Output( : )=Design_Variable( :, Individual_Optimal( Optimization_Step ) ) 
      end if

      call MPI_Barrier( MPI_COMM_WORLD, ierr )
      !===========================================================================
      if( MyRank==0 ) write(*,*)'Output FEM Data : Generation Interval'
      !===========================================================================

      if( Flag_Output_Finite_Element_Data/=0 .and.  &
          ( mod( Optimization_Step, Interval_Plot_Optimization_Step ) == 0 .or. &
          Optimization_Step==Optimization_Step_Output ) .and. &
          MyRank_Output==MyRank )then

         call Output_Finite_Element_Data &
            ( Number_Node_Output, Position_Node_Output, &
              Max_Number_Element_Share_OneNode_Output, Max_Position_Node_Output, &
              Number_Element_Triangle_Output, Index_Element_2_Node_Triangle_Output, &
              Class_Element_Triangle_Output, &
              Element_and_LocalNode_Number_on_Electrical_Insulation_BC_Output, &
              Number_Node_on_Electrical_Insulation_Boundary_Output, Number_Node_in_Electrical_Insulation_Output, &
              Number_Element_Electrical_Insulation_Boundary_Output, &
              2, Optimization_Step )
      end if

      call MPI_Barrier( MPI_COMM_WORLD, ierr )

      !===========================================================================
      if( MyRank==0 ) write(*,*)'Output CPU Time' 
      !===========================================================================

      if( MyRank==0 )then
         call cpu_time( Time_CPU( Optimization_Step ) )
         open( 1609, file='./CPU_Time.dat', position='append' ) 
         if( Optimization_Step==0 )then
            write(1609,*) Optimization_Step, Time_CPU( Optimization_Step )-Time_CPU_Start, Time_CPU( Optimization_Step )-Time_CPU_Start
         else
            write(1609,*) Optimization_Step, Time_CPU( Optimization_Step )-Time_CPU_Start, Time_CPU( Optimization_Step )-Time_CPU( Optimization_Step-1 )
         end if
         close( 1609 )

         do i= 1, Number_Fitness_CPU_Threshold
            if( Flag_CPU_Time==1 .and. minval( Fitness_CMAES ) <= Fitness_CPU_Threshold( i ) )then
               !===========================================================================
               write(*,*)'Output CPU Time for Fitness' 
               !===========================================================================
   
               Generation_CPU( 1 )= Optimization_Step
               open( 1617, file='./CPU_Time_Fitness.dat', position='append' ) 
               write(1617,1630) Optimization_Step, &
                          Time_CPU( Optimization_Step )-Time_CPU_Start, Time_CPU( Optimization_Step )-Time_CPU( Generation_CPU( 1 ) ), &
                          minval( Fitness_CMAES ), Fitness_CPU_Threshold( i )
               1630 format( I10, 1x, es15.8, 1x, es15.8, 1x, es15.8, 1x, es15.8 )
               close( 1617 )
   
               Fitness_CPU_Threshold( i )= Fitness_CPU_Threshold( i ) *1d-1
               goto 1514
            end if
         end do

         1514 continue

         do i= 1, Number_Obj_Func_CPU_Threshold
            if( Flag_CPU_Time==1 .and. Obj_Func_Real_Minimum <= Obj_Func_CPU_Threshold( i ) )then
               !===========================================================================
               write(*,*)'Output CPU Time for OF' 
               !===========================================================================

               Generation_CPU( 2 )= Optimization_Step
               open( 1632, file='./CPU_Time_OF.dat', position='append' ) 
               write( 1632,1645) Optimization_Step, &
                           Time_CPU( Optimization_Step )-Time_CPU_Start, Time_CPU( Optimization_Step )-Time_CPU( Generation_CPU( 2 ) ), &
                           Obj_Func_Real_Minimum, Obj_Func_CPU_Threshold( i )
               1645 format( I10, 1x, es15.8, 1x, es15.8, 1x, es15.8, 1x, es15.8 )
               close( 1632 )

               Obj_Func_CPU_Threshold( i )= Obj_Func_CPU_Threshold( i ) *1d-1
               goto 1542
            end if
         end do

         1542 continue
      end if

      call MPI_Barrier( MPI_COMM_WORLD, ierr )

      !===========================================================================
      if( MyRank==0 ) write(*,*)'Output FEM Data : the Convergence of xmean'
      !===========================================================================

      if( Optimization_Step==0 )then
         Threshold_Convergence_xmean= 0d0
         Threshold_Convergence_xmean( 1 )= 1d-2
         Threshold_Convergence_xmean( 2 )= 2d-2
         Threshold_Convergence_xmean( 3 )= 3d-2
         Threshold_Convergence_xmean( 4 )= 5d-2

         Threshold_Sigma_xmean= Radius_Search *1d0
      else if( MyRank_Output==MyRank )then 
         do i= 1, Number_Threshold_Convergence_xmean 
            !if( abs( Convergence_Ratio_xmean( Optimization_Step ) ) <= Threshold_Convergence_xmean( i ) .and. &
            !    cma%sigma <= Threshold_Sigma_xmean )then
            !if( abs( Convergence_Ratio_xmean( Optimization_Step ) ) <= Threshold_Convergence_xmean( i ) )then
            if( Termination( Optimization_Step )=='y' )then

               Threshold_Convergence_tmp= Threshold_Convergence_xmean( i ) 

               Generation_Convergence_xmean= Optimization_Step

               !======================================================================================
               open( 2157, file='./Convergence_xmean_Generation.dat', position='append' ) 
               !======================================================================================
               write( 2157, 2161 ) Generation_Convergence_xmean, Generation_Output, &
                           Threshold_Convergence_xmean( i ), &
                           abs( Convergence_Ratio_xmean( Optimization_Step ) ), Convergence_Ratio_xmean( Optimization_Step )
               write( 2157, 2162 ) 'Before : ', Threshold_Convergence_xmean( 1 ), Threshold_Convergence_xmean( 2 ), &
                           Threshold_Convergence_xmean( 3 ), Threshold_Convergence_xmean( 4 ) 
               2161 format( i10, 1x, i10, 1x, es15.7, 1x, es15.7, 1x, es15.7  )
               2162 format( a9, 1x, es15.7, es15.7, es15.7, es15.7 )

               call Output_Finite_Element_Data &
                  ( Number_Node_Output, Position_Node_Output, &
                    Max_Number_Element_Share_OneNode_Output, Max_Position_Node_Output, &
                    Number_Element_Triangle_Output, Index_Element_2_Node_Triangle_Output, &
                    Class_Element_Triangle_Output, &
                    Element_and_LocalNode_Number_on_Electrical_Insulation_BC_Output, &
                    Number_Node_on_Electrical_Insulation_Boundary_Output, Number_Node_in_Electrical_Insulation_Output, &
                    Number_Element_Electrical_Insulation_Boundary_Output, &
                    4, Generation_Output )

               !==================================================
               if( Flag_Output_PostScript/=0 )then
               !==================================================
      
                  allocate( Position_Node_Plot( 2, Number_Node_Output ) )
                  allocate( Index_Element_2_Node_Plot( 4, Number_Element_Triangle_Output ) )
                  allocate( Class_Element_Plot( Number_Element_Triangle_Output ) )
      
                  Position_Node_Plot( :, 1:Number_Node_Output )= Position_Node_Output( :, 1:Number_Node_Output )
      
                  Class_Element_Plot( 1:Number_Element_Triangle_Output )&
                  = Class_Element_Triangle_Output( 1:Number_Element_Triangle_Output )
      
                  Index_Element_2_Node_Plot( 1:3, 1:Number_Element_Triangle_Output )&
                  = Index_Element_2_Node_Triangle_Output( 1:3, 1:Number_Element_Triangle_Output )
                  Index_Element_2_Node_Plot( 4, : )= Index_Element_2_Node_Plot( 3, : )
      
                  call Plot_FEM_Model &
                    ( Number_Node_Output, Number_Element_Triangle_Output, &
                      Generation_Output, &
                      Position_Node_Plot, Index_Element_2_Node_Plot, Class_Element_Plot )
      
                  deallocate( Position_Node_Plot )
                  deallocate( Index_Element_2_Node_Plot )
                  deallocate( Class_Element_Plot )
                !==================================================
                end if !Flag_Output_PostScript/=0
                !==================================================
  
                !==================================================
                if( Flag_Output_Distribution_Result/=0 )then
                !==================================================
      
                   allocate( Position_Node_Plot( 2, Number_Node_Output ) )
                   allocate( Index_Element_2_Node_Plot( 4, Number_Element_Triangle_Output ) )
                   allocate( Value_Plot( Number_Node_Output ) )
       
                   Position_Node_Plot( :, 1:Number_Node_Output )= Position_Node_Output( :, 1:Number_Node_Output )
       
                   Index_Element_2_Node_Plot( 1:3, 1:Number_Element_Triangle_Output )&
                   = Index_Element_2_Node_Triangle_Output( 1:3, 1:Number_Element_Triangle_Output )
                   Index_Element_2_Node_Plot( 4, : )= Index_Element_2_Node_Plot( 3, : )
       
                   do Loop_Source= 1, Number_Obj_Func_Source
       
                     if( Flag_Physics==0 .or. Flag_Physics==1 )then
                       do Loop_MatParam_T= 1, Number_Obj_Func_MP_T
       
                          Value_Plot( 1:Number_Node_Output )= Value_Plot_Output_Thermal( 1:Number_Node_Output, Loop_Source, Loop_MatParam_T )
                          Flag_Physics_Plot = 1
       
                          call Plot_Result_Element_Triangle &
                              ( Optimization_Step, Flag_Physics_Plot, &
                                Value_Plot, Number_Level_Plot_EM, &
                                Generation_Output, Loop_Source, Loop_MatParam_T, &
                                Position_Node_Plot, Number_Node_Output,  &
                                Index_Element_2_Node_Plot, Number_Element_Triangle_Output, &
                                Obj_Func_Output )
       
                       end do
                     end if
       
                     if( Flag_Physics==0 .or. Flag_Physics==2 )then
                       do Loop_MatParam_V= 1, Number_Obj_Func_MP_V
                          Value_Plot( 1:Number_Node_Output )= Value_Plot_Output_DC( 1:Number_Node_Output, Loop_Source, Loop_MatParam_V )
                          Flag_Physics_Plot = 2
       
                       call Plot_Result_Element_Triangle &
                           ( Optimization_Step, Flag_Physics_Plot, &
                             Value_Plot, Number_Level_Plot_EM, &
                             Generation_Output, Loop_Source, Loop_MatParam_V, &
                             Position_Node_Plot, Number_Node_Output,  &
                             Index_Element_2_Node_Plot, Number_Element_Triangle_Output, &
                             Obj_Func_Output )
                       end do
                     end if
       
                   end do
       
                   deallocate( Value_Plot )
                   deallocate( Position_Node_Plot )
                   deallocate( Index_Element_2_Node_Plot )
                !====================================================
                end if !Flag_Output_Distribution_Result/=0
                !====================================================

                !====================================================
                if( Flag_Output_LSF/=0 )then
                !====================================================
                   write( Format_Filenumber, '(A2,I1,A1,I1,A1)' ) & 
                   '(i', Length_Character_Optimization_Step,'.',Length_Character_Optimization_Step,')'
        
                   write( Generation_Character, Format_Filenumber ) Generation_Output
      
                   FileName_LSF= trim( "LSF_Convergence_"//trim(Generation_Character) )
                   open( 2666, file=FileName_LSF, position='append' ) 
                      write( 2666, * ) Generation_Output, '# Generation_Output'
                      write( 2666, * ) Threshold_Convergence_xmean( i ), '# Threshold_Convergence_xmean'
                      write( 2666, * ) Flag_Symmetry_LSF, '# Flag_Symmetry_LSF'
                      write( 2666, * ) Number_GridPoint_Optimization, '# Number_GridPoint_Optimization'
                      do j= 1, Number_GridPoint_Optimization
                         write( 2666, * ) j, Design_Variable_Output( j ) 
                      end do
                   close( 2666 )
                !====================================================
                end if !Flag_Output_LSF==1
                !====================================================

                !======================================================================================
                write(*,*)'Update the value for Convergence Criterion'
                !======================================================================================

                do j= 1, Number_Threshold_Convergence_xmean 
                   if( Threshold_Convergence_xmean( j )*( 1d0+1d-8 ) >= Threshold_Convergence_tmp )then
                      Threshold_Convergence_xmean( j )= Threshold_Convergence_xmean( j )*0.1d0
                   end if
                end do

                write( 2157, 2163 ) 'After  : ', Threshold_Convergence_xmean( 1 ), Threshold_Convergence_xmean( 2 ), &
                            Threshold_Convergence_xmean( 3 ), Threshold_Convergence_xmean( 4 ) 
                2163 format( a9, 1x, es15.7, es15.7, es15.7, es15.7 )
                close( 2157 )

                go to 2193 
            end if
         end do 
         2193 continue
      end if

      call MPI_Bcast( Threshold_Convergence_xmean(1), Number_Threshold_Convergence_xmean, MPI_DOUBLE_PRECISION, MyRank_Output, MPI_COMM_WORLD, ierr)
      call MPI_Barrier( MPI_COMM_WORLD, ierr )

      if( MyRank==0 )then
         !======================================================================================
         open( 1430, file='Convergence_Ratio.dat', position='append' ) 
         !======================================================================================
            if( Optimization_Step==0 )then
               write( 1430, '(a101)' )' # Generation, C. R. xmean,    abs(C. R.),    Ave(C. R.),    cp_sigma,     Fitness,   Param. ED' 
            else
               write( 1430, 1441 ) Optimization_Step, &
                           Convergence_Ratio_xmean( Optimization_Step ), &
                           Fitness_CMAES( Individual_Optimal( Optimization_Step ) )
               1441 format( 1x, I10, es15.5, es15.5 )
            end if
         close( 1430 )
      end if
  
      call MPI_Barrier( MPI_COMM_WORLD, ierr )

      !Update Generation
      Optimization_Step=Optimization_Step+1

      open(2816, file='Convergence_Error')
         read(2816,*) Convergence_Error 
         if( MyRank==0 ) write(*,*)'Convergence_Error=', Convergence_Error 
      close(2816)

   !========================================================================================
   end do ! Optimization_Step do while
   !========================================================================================

   deallocate( Position_X_Source )
   deallocate( Position_Y_Source )

   deallocate( Relative_Thermal_Conductivity )
   deallocate( Relative_Electrical_Conductivity )

   !==========================================================================================================================
   if( Flag_ObjectiveFunction==11 .and. ( Flag_Physics==0 .or. Flag_Physics==1 ) )then
   !==========================================================================================================================
       deallocate( Temperature_Reference )
   end if
   !==========================================================================================================================
   if( Flag_ObjectiveFunction==11 .and. ( Flag_Physics==0 .or. Flag_Physics==2 ) )then
   !==========================================================================================================================
       deallocate( Electric_Potential_Reference )
   end if

   deallocate( Index_DV_Symmetric_2_GP_All )
   deallocate( Index_GP_All_2_DV_Symmetric )

   !==========================================================================================================================
   if( Flag_Physics==0 .or. Flag_Physics==1 )then
   !==========================================================================================================================
      deallocate( Value_Obj_Func_Normalize_Thermal )
   end if
   !==========================================================================================================================
   if( Flag_Physics==0 .or. Flag_Physics==2 )then
   !==========================================================================================================================
      deallocate( Value_Obj_Func_Normalize_DC )
   end if

   deallocate( Index_Grid_2_GridPoint )
   deallocate( Index_Grid_total )
   deallocate( plot_Number )
   deallocate( Position_Node_Preserve )
   deallocate( Index_Element_2_Node_Triangle_Preserve )
   deallocate( Class_Element_Triangle_Preserve )
   deallocate( Value_Plot_Preserve_Thermal )
   deallocate( Value_Plot_Preserve_DC )
   deallocate( Class_Node_Preserve )
   deallocate( Max_Position_Node_Preserve )
   deallocate( Element_and_LocalNode_Num_on_Electrical_Insulation_BC_Preserve )

   deallocate( Position_Node_Output )
   deallocate( Index_Element_2_Node_Triangle_Output )
   deallocate( Class_Element_Triangle_Output )
   deallocate( Value_Plot_Output_Thermal )
   deallocate( Value_Plot_Output_DC )
   deallocate( Class_Node_Output )
   deallocate( Max_Position_Node_Output )
   deallocate( Element_and_LocalNode_Number_on_Electrical_Insulation_BC_Output )
 
 
   deallocate( Obj_Func_All ) 
   deallocate( Obj_Func_Single_Thermal ) 
   deallocate( Obj_Func_Single_DC ) 
   
   deallocate( Class_GridPoint_Scattering )
   deallocate( Position_GridPoint_Scattering )
   deallocate( LSF_FixedDomain_GridPoint )
   deallocate( LSF_DesignDomain_GridPoint )
   deallocate( LSF_ExteriorDomain_GridPoint )
   deallocate( Design_Variable )
   deallocate( Design_Variable_Blackbox )
   deallocate( Design_Variable_Output )
   deallocate( LSF_Piecewise_Constant )
   deallocate( Diff_LSF_PC )
   deallocate( Average_Diff_LSF_PC )

   !deallocate( Class_GridPoint_PML )
   !deallocate( Position_GridPoint_PML )
   !deallocate( Index_Grid_PML_2_GridPoint_PML )

   deallocate( Correspondence_Number )
  
   !! PML
   !deallocate( Position_Node_OneGrid_PML )
   !deallocate( Number_Node_OneGrid_PML )
   !deallocate( Number_Element_OneGrid_PML )
   !deallocate( Class_Grid_PML )
   !deallocate( Class_Element_OneGrid_PML )
   !deallocate( Index_Element_PML_2_Node_PML )
  
   deallocate( Fitness_CMAES )
   deallocate( Fitness_Convergence )
   deallocate( Obj_Func_Real )
   deallocate( Fractal_Structure_All )
   deallocate( Perimeter_Structure_All )
   if( Flag_Perimeter_Constraint==2 ) deallocate( Perimeter_Implicit_All )
   deallocate( Volume_Structure_All )
   deallocate( Penalty_Volume_Constraint )
   deallocate( Individual_Optimal )
   deallocate( Obj_Func_Multi )
   deallocate( Termination )
   deallocate( Convergence_Ratio )
   deallocate( Convergence_Ratio_Average )
   deallocate( Convergence_Ratio_xmean )
   deallocate( Threshold_Convergence )
   deallocate( Threshold_Convergence_Average )
   deallocate( Threshold_Convergence_xmean )
 
   if( Flag_CPU_Time==1 )then
      deallocate( Time_CPU )
      deallocate( Generation_CPU )
      deallocate( Fitness_CPU_Threshold )
      deallocate( Obj_Func_CPU_Threshold )
   end if

   if( MyRank==0 )then
      write(*,*)'=============================================================='
      write(*,*)'=============================================================='
      write(*,*)'====      =   ==    ===  ==   ===    ==  ====  ==============='
      write(*,*)'====   ====   ==   =  =  ==   ==   ====  ====  ==============='
      write(*,*)'====     ==   ==   =     ==   ===   ===        ==============='
      write(*,*)'====   ====   ==   ==    ==   =====  ==  ====  ==============='
      write(*,*)'====   ====   ==   ===   ==   ==     ==  ====  ==============='
      write(*,*)'=============================================================='
      write(*,*)'=============================================================='
   end if
 
end program Main
 


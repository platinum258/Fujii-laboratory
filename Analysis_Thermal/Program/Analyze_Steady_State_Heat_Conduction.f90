
subroutine Analyze_Steady_State_Heat_Conduction &
       ( Flag_Structure, &
         Temperature_Lower_Side_Boundary, Temperature_Higher_Side_Boundary, &
         Flag_Thermal_Insulation_BC, &
         Thermal_Conductivity_Material, Thermal_Conductivity_Material_DesignDomain, & 
         Thermal_Conductivity_OuterDomain, Thermal_Conductivity_FixedDomain_Variable, &
         ID_Element_Material, ID_Element_OuterDomain, ID_Element_FixedDomain, ID_Element_Base_Material, & 
         ID_Element_Material_Exterior, ID_Element_Base_Material_Exterior, & 
         Number_Node, Number_Element_Triangle, &
         Position_Node, &
         Width_Matrix_LHS, Max_Position_Node, &
         Index_Element_2_Node_Triangle, Class_Element_Triangle, Thermal_conductivity_by_density, &
         Position_Source_X, Position_Source_Y, & 
         high, thickness, penalty_density, &
         !============================================================================================
         Temperature_Solution, Thermal_Conductivity_Triangle )



!subroutine Analyze_Steady_State_Heat_Conduction &
!       ( Temperature_Lower_Side_Boundary, Temperature_Higher_Side_Boundary, &
!         Number_Node, Number_Element_Triangle, &
!         Position_Node, &
!         Width_Matrix_LHS, Max_Position_Node, &
!         Index_Element_2_Node_Triangle, Class_Element_Triangle, &
!         !============================================================================================
!         Temperature_Solution )

   !$ use omp_lib
   use Parameters
   implicit none

   ! intent(in)
   double precision, intent(in) :: Temperature_Lower_Side_Boundary, Temperature_Higher_Side_Boundary 
   integer, intent(in) :: Flag_Thermal_Insulation_BC , Flag_Structure

   double precision, intent(in) :: Thermal_Conductivity_Material, Thermal_Conductivity_Material_DesignDomain 
   double precision, intent(in) :: Thermal_Conductivity_OuterDomain, Thermal_Conductivity_FixedDomain_Variable 

   integer, intent(in) :: ID_Element_Material, ID_Element_OuterDomain, ID_Element_FixedDomain, ID_Element_Base_Material  
   integer, intent(in) :: ID_Element_Material_Exterior, ID_Element_Base_Material_Exterior  

   integer, intent(in) :: Number_Node
   double precision, intent(in) :: Position_Node( 2, Number_Node ) 
 
   integer, intent(in) :: Number_Element_Triangle 
   integer, intent(in) :: Index_Element_2_Node_Triangle( 3, Number_Element_Triangle ) 
   integer, intent(in) :: Class_Element_Triangle( Number_Element_Triangle )
 
   double precision, intent(in) :: Max_Position_Node( 2, 2 ) 
   integer, intent(in) :: Width_Matrix_LHS

   double precision, intent(in) :: Position_Source_X, Position_Source_Y
 
   ! intent(out)
   !complex( kind( 0d0 ) ), intent(out) :: Temperature_Solution( Number_Node )
   double precision, intent(out) :: Temperature_Solution( Number_Node )
   double precision, intent(out) :: Thermal_Conductivity_Triangle( Number_Element_Triangle )
 
   integer :: i, j, k, e
   integer :: Counter
   integer :: NumThread_OMP, NumThread_MKL, NumThread_MPI

   integer, allocatable, dimension(:,:) :: Index_Element_2_Node_Triangle_Check
   double precision, allocatable, dimension(:) :: Area_Element_Triangle

   ! Model Data
   !complex( kind( 0d0 ) ), allocatable, dimension(:) :: Thermal_Conductivity_Triangle
   !double precision, allocatable, dimension(:) :: Thermal_Conductivity_Triangle
 
   ! Finite Element Method 
   double precision, allocatable, dimension(:,:,:) :: Position_Node_Element_Triangle
   double precision, allocatable, dimension(:,:) :: BasisFunction_A, BasisFunction_B, BasisFunction_C
   double precision, allocatable, dimension(:,:,:,:) :: Difference_Position_Node_Element_Triangle
   double precision, allocatable, dimension(:,:,:) :: BasisFunctionBC_Times_Difference_Position_Node
   double precision, allocatable, dimension(:,:) :: Basis_Function_Times_Position_Node
   !complex( kind( 0d0 ) ), allocatable, dimension(:,:,:) :: Coefficient_C, Coefficient_Xi, Coefficient_Eta
   !complex( kind( 0d0 ) ), allocatable, dimension(:,:,:) :: Coefficient_Xi_Eta, Coefficient_Xi_Xi, Coefficient_Eta_Eta
   double precision, allocatable, dimension(:,:,:) :: Coefficient_C, Coefficient_Xi, Coefficient_Eta
   double precision, allocatable, dimension(:,:,:) :: Coefficient_Xi_Eta, Coefficient_Xi_Xi, Coefficient_Eta_Eta
   
   !complex( kind( 0d0 ) ), allocatable, dimension(:,:,:) :: LocalMatrix_Triangle 
   double precision, allocatable, dimension(:,:,:) :: LocalMatrix_Triangle 

   !PEC
   integer, allocatable, dimension(:) :: Flag_Element_Triangle_on_Thermal_Insulator_BC
   integer :: Number_Outer_Element_on_PEC_BC

   ! Flat_PEC for Thermal cloak
   integer, allocatable, dimension(:) :: Element_Number_on_Flat_PEC
   integer :: Number_Element_on_Flat_PEC, Number_Element_on_Flat_PEC_tmp

   integer, allocatable, dimension(:,:) :: Node_Number_on_Flat_PEC_in_Element

   double precision, allocatable, dimension(:,:,:) :: Position_Node_PEC_BC, Position_Node_PEC_BC_tmp
   double precision, allocatable, dimension(:,:) :: Unit_Normal_Vector_on_PEC_BC, Unit_Normal_Vector_on_PEC_BC_tmp
   double precision, allocatable, dimension(:) :: Unit_Normal_Vector_on_Flat_PEC_BC 
   !complex( kind( 0d0 ) ), allocatable, dimension(:) :: Nodal_Value_Fixed_in_PEC 
   double precision, allocatable, dimension(:) :: Nodal_Value_Fixed_in_PEC 
   double precision, allocatable, dimension(:) :: Length_Edge_on_PEC_BC_in_Element 
   double precision :: Length_Edge_on_Flat_PEC_BC_in_Element 

   integer :: Number_NonZero  
   !complex( kind( 0d0 ) ), allocatable, dimension(:,:) :: GlobalMatrix
   double precision, allocatable, dimension(:,:) :: GlobalMatrix
   integer, allocatable, dimension(:,:) :: J_GlobalMatrix
   !complex( kind( 0d0 ) ), allocatable, dimension(:) ::  Global_Vector_RHS 
   double precision, allocatable, dimension(:) ::  Global_Vector_RHS 

   integer, allocatable, dimension(:) :: Flag_Node_Dirichlet_BC
   !complex( kind( 0d0 ) ), allocatable, dimension(:) :: Nodal_Value_on_Dirichlet_BC 
   double precision, allocatable, dimension(:) :: Nodal_Value_on_Dirichlet_BC 
 
   ! MUMPS
   !complex( kind( 0d0 ) ), allocatable, dimension(:) :: aa_zmumps
   double precision, allocatable, dimension(:) :: aa_zmumps
   integer, allocatable, dimension(:) :: ja_zmumps, ia_zmumps
   integer, allocatable, dimension(:) :: Flag_Element_Triangle_Removed

   integer :: Flag_Remove_Element, ID_Element_Thermal_Insulator 
   integer, parameter :: Flag_Symmetric= Flag_On
   integer, parameter :: Flag_Matrix_Upper= Flag_On

   character(len=1) :: match 

  
!****************************************************************************************************:
   double precision :: high,thickness
   integer :: penalty_density

   double precision :: Thermal_conductivity_by_density( Number_Element_Triangle )

   !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


   Flag_Remove_Element= Flag_Thermal_Insulation_BC
   ID_Element_Thermal_Insulator= ID_Element_FixedDomain

   !============================================================================================
   write(*,*)'   call Analyze_Steady_State_Heat_Conduction '
   !============================================================================================
    
   open(9999, file='Number_Thread')
   read(9999,*)NumThread_OMP, NumThread_MKL, NumThread_MPI
   write(*,*)'NumThread_OMP=', NumThread_OMP
   write(*,*)'NumThread_MKL=', NumThread_MKL
   write(*,*)'NumThread_MPI=', NumThread_MPI
   !$ call omp_set_num_threads( NumThread_OMP )
   !!!$ call mkl_set_num_threads( NumThread_MKL )
   close(9999)

   allocate( Index_Element_2_Node_Triangle_Check( 3, Number_Element_Triangle ) )

   do e= 1, Number_Element_Triangle
    do i= 1, 3
       Index_Element_2_Node_Triangle_Check( i, e )= Index_Element_2_Node_Triangle( i, e )
    end do
   end do
 
   !===============================================================================================
   write(*,*)'    =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= '
   write(*,*)'    Temperature_Lower_Side_Boundary=', Temperature_Lower_Side_Boundary
   write(*,*)'    Temperature_Higher_Side_Boundary=', Temperature_Higher_Side_Boundary
   write(*,*)'    =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= '
   !===============================================================================================

   !============================================================================================
   write(*,*)'    '
   write(*,*)'    <<<<<<<<<<<<<<<<<<<< Finite Element Method >>>>>>>>>>>>>>>>>>> '
   write(*,*)'    '
   !============================================================================================
 
   allocate( Position_Node_Element_Triangle( 2, 3, Number_Element_Triangle ) )
    
   !$omp parallel do default( none ) &
   !$omp private( e, i, j ) & 
   !$omp shared( Number_Element_Triangle, Position_Node, Position_Node_Element_Triangle ) & 
   !$omp shared( Index_Element_2_Node_Triangle ) 
   do e= 1, Number_Element_Triangle
    do i= 1, 3
       do j= 1, 2
        Position_Node_Element_Triangle( j, i, e )= Position_Node( j, Index_Element_2_Node_Triangle( i, e ) )
       end do
    end do
   end do

   allocate( Difference_Position_Node_Element_Triangle( 2, 3, 3, Number_Element_Triangle ) )

   do e= 1, Number_Element_Triangle
    do i= 1, 3
       do j= 1, 3
        do k= 1, 2
         Difference_Position_Node_Element_Triangle( k, j, i, e ) &
         =Position_Node_Element_Triangle( k, j, e ) -Position_Node_Element_Triangle( k, i, e )
        end do
       end do
    end do
   end do

   allocate( Area_Element_Triangle( Number_Element_Triangle ) )

   !$omp parallel do default( none ) &
   !$omp private( e ) & 
   !$omp shared( Number_Element_Triangle, Area_Element_Triangle, Difference_Position_Node_Element_Triangle ) 
   do e= 1, Number_Element_Triangle
    Area_Element_Triangle( e ) &
    =abs( Difference_Position_Node_Element_Triangle( 1, 1, 3, e ) & 
       *Difference_Position_Node_Element_Triangle( 2, 2, 3, e ) &
       -Difference_Position_Node_Element_Triangle( 1, 2, 3, e ) & 
       *Difference_Position_Node_Element_Triangle( 2, 1, 3, e ) )/2d0
   end do 
    
   !$omp parallel do default( none ) &
   !$omp private( e ) & 
   !$omp shared( Number_Element_Triangle, Area_Element_Triangle, Position_Node_Element_Triangle, Class_Element_Triangle, Index_Element_2_Node_Triangle ) 
   do e= 1, Number_Element_Triangle
    if( Area_Element_Triangle( e )==0.0d0 )then
       write(*,*)'Area_Element_Triangle( e )==0.0d0'
       write(*,*)'Element Number=', e, 'Class_Element=', Class_Element_Triangle( e )
       write(*,*)'Position Node 1=', Position_Node_Element_Triangle( 1, 1, e ), Position_Node_Element_Triangle( 2, 1, e ) 
       write(*,*)'Position Node 2=', Position_Node_Element_Triangle( 1, 2, e ), Position_Node_Element_Triangle( 2, 2, e ) 
       write(*,*)'Position Node 3=', Position_Node_Element_Triangle( 1, 3, e ), Position_Node_Element_Triangle( 2, 3, e ) 
       write(*,*)'Index_Element_2_Node_Triangle=', Index_Element_2_Node_Triangle( 1, e ), &
                  Index_Element_2_Node_Triangle( 2, e ), Index_Element_2_Node_Triangle( 3, e ) 
       call Output_Error( 'Analyze_Steady_State_Heat_Conduction', 166 )
    end if
   end do 
    
   allocate( BasisFunction_A( 3, Number_Element_Triangle ) )
   allocate( BasisFunction_B( 3, Number_Element_Triangle ) )
   allocate( BasisFunction_C( 3, Number_Element_Triangle ) )
   
   !$omp parallel do default( none )   &
   !$omp private( e ) &
   !$omp shared( Number_Element_Triangle, Position_Node_Element_Triangle ) &
   !$omp shared( Area_Element_Triangle, BasisFunction_A )
   do e= 1, Number_Element_Triangle
  
    BasisFunction_A( 1, e ) &
    = ( Position_Node_Element_Triangle( 1, 2, e ) *Position_Node_Element_Triangle( 2, 3, e ) &
       -Position_Node_Element_Triangle( 1, 3, e ) *Position_Node_Element_Triangle( 2, 2, e ) ) &
     /( 2d0 *Area_Element_Triangle( e ) )

    BasisFunction_A( 2, e ) &
    = ( Position_Node_Element_Triangle( 1, 3, e ) *Position_Node_Element_Triangle( 2, 1, e ) & 
       -Position_Node_Element_Triangle( 1, 1, e ) *Position_Node_Element_Triangle( 2, 3, e ) ) &
     /( 2d0 *Area_Element_Triangle( e ) )

    BasisFunction_A( 3, e ) & 
    = ( Position_Node_Element_Triangle( 1, 1, e ) *Position_Node_Element_Triangle( 2, 2, e ) & 
       -Position_Node_Element_Triangle( 1, 2, e ) *Position_Node_Element_Triangle( 2, 1, e ) ) & 
     /( 2d0 *Area_Element_Triangle( e ) )
    
   end do
 
   !$omp parallel do default( none )   &
   !$omp private( e ) &
   !$omp shared( Number_Element_Triangle, Difference_Position_Node_Element_Triangle ) &
   !$omp shared( Area_Element_Triangle, BasisFunction_B )
   do e= 1, Number_Element_Triangle

    BasisFunction_B( 1, e ) &
    = Difference_Position_Node_Element_Triangle( 2, 2, 3, e )/( 2d0 *Area_Element_Triangle( e ) )

    BasisFunction_B( 2, e ) & 
    = Difference_Position_Node_Element_Triangle( 2, 3, 1, e )/( 2d0 *Area_Element_Triangle( e ) )

    BasisFunction_B( 3, e ) & 
    = Difference_Position_Node_Element_Triangle( 2, 1, 2, e )/( 2d0 *Area_Element_Triangle( e ) )
    
   end do
  
   !$omp parallel do default( none )   &
   !$omp private( e ) &
   !$omp shared( Number_Element_Triangle, Difference_Position_Node_Element_Triangle ) &
   !$omp shared( Area_Element_Triangle, BasisFunction_C )
   do e= 1, Number_Element_Triangle

    BasisFunction_C( 1, e ) & 
    = Difference_Position_Node_Element_Triangle( 1, 3, 2, e )/( 2d0 *Area_Element_Triangle( e ) )

    BasisFunction_C( 2, e ) & 
    = Difference_Position_Node_Element_Triangle( 1, 1, 3, e )/( 2d0 *Area_Element_Triangle( e ) )

    BasisFunction_C( 3, e ) & 
    = Difference_Position_Node_Element_Triangle( 1, 2, 1, e )/( 2d0 *Area_Element_Triangle( e ) )
    
   end do
   
   allocate( Basis_Function_Times_Position_Node( 3, Number_Element_Triangle ) ) 

   !$omp parallel do default( none )   &
   !$omp private( e, i ) &
   !$omp shared( Number_Element_Triangle, Basis_Function_Times_Position_Node ) &
   !$omp shared( BasisFunction_A, BasisFunction_B, BasisFunction_C, Position_Node_Element_Triangle )
   do e= 1, Number_Element_Triangle
    do i= 1, 3
       Basis_Function_Times_Position_Node( i, e ) & 
       = BasisFunction_A( i, e ) & 
      +BasisFunction_B( i, e ) *Position_Node_Element_Triangle( 1, 3, e ) & 
      +BasisFunction_C( i, e ) *Position_Node_Element_Triangle( 2, 3, e ) 
    end do
   end do

   allocate( BasisFunctionBC_Times_Difference_Position_Node( 2, 3, Number_Element_Triangle ) )

   !$omp parallel do default( none )   &
   !$omp private( e, i, j ) &
   !$omp shared( Number_Element_Triangle, BasisFunctionBC_Times_Difference_Position_Node ) &
   !$omp shared( BasisFunction_A, BasisFunction_B, BasisFunction_C, Difference_Position_Node_Element_Triangle )
   do e= 1, Number_Element_Triangle
    do i= 1, 3
       do j= 1, 2
        BasisFunctionBC_Times_Difference_Position_Node( j, i, e ) &
        = BasisFunction_B( i, e ) *Difference_Position_Node_Element_Triangle( 1, j, 3, e ) &
         +BasisFunction_C( i, e ) *Difference_Position_Node_Element_Triangle( 2, j, 3, e )
       end do
    end do
   end do
 
   deallocate( Difference_Position_Node_Element_Triangle )

   allocate( Coefficient_C( 3, 3, Number_Element_Triangle ) )
   allocate( Coefficient_Xi( 3, 3, Number_Element_Triangle ) )
   allocate( Coefficient_Eta( 3, 3, Number_Element_Triangle ) )
   allocate( Coefficient_Xi_Eta( 3, 3, Number_Element_Triangle ) )
   allocate( Coefficient_Xi_Xi( 3, 3, Number_Element_Triangle ) )
   allocate( Coefficient_Eta_Eta( 3, 3, Number_Element_Triangle ) )
  
   !$omp parallel do default( none )   &
   !$omp private( e, i, j ) &
   !$omp shared( Number_Element_Triangle ) &
   !$omp shared( Basis_Function_Times_Position_Node, Coefficient_C ) 
   do e= 1, Number_Element_Triangle
    do i= 1, 3
       do j= 1, 3
        Coefficient_C( j, i, e ) & 
        = Basis_Function_Times_Position_Node( i, e )*Basis_Function_Times_Position_Node( j, e ) 
       end do
    end do
   end do

   !$omp parallel do default( none )   &
   !$omp private( e, i, j ) &
   !$omp shared( Number_Element_Triangle, Coefficient_Xi ) & 
   !$omp shared( Basis_Function_Times_Position_Node ) &
   !$omp shared( BasisFunctionBC_Times_Difference_Position_Node ) 
   do e= 1, Number_Element_Triangle
    do i= 1, 3
       do j= 1, 3
        Coefficient_Xi( j, i, e ) & 
        = Basis_Function_Times_Position_Node( j, e ) &
         *BasisFunctionBC_Times_Difference_Position_Node( 1, i, e ) &
         +Basis_Function_Times_Position_Node( i, e ) &
         *BasisFunctionBC_Times_Difference_Position_Node( 1, j, e ) 
       end do
    end do
   end do

   !$omp parallel do default( none )   &
   !$omp private( e, i, j ) &
   !$omp shared( Number_Element_Triangle, Coefficient_Eta ) & 
   !$omp shared( Basis_Function_Times_Position_Node ) &
   !$omp shared( BasisFunctionBC_Times_Difference_Position_Node ) 
   do e= 1, Number_Element_Triangle
    do i= 1, 3
       do j= 1, 3
        Coefficient_Eta( j, i, e ) &
        = Basis_Function_Times_Position_Node( j, e ) &
         *BasisFunctionBC_Times_Difference_Position_Node( 2, i, e ) &
         +Basis_Function_Times_Position_Node( i, e ) &
         *BasisFunctionBC_Times_Difference_Position_Node( 2, j, e ) 
       end do
    end do
   end do

   deallocate( Basis_Function_Times_Position_Node ) 

   !$omp parallel do default( none )   &
   !$omp private( e, i, j ) &
   !$omp shared( Number_Element_Triangle, Coefficient_Xi_Eta ) & 
   !$omp shared( BasisFunctionBC_Times_Difference_Position_Node ) 
   do e= 1, Number_Element_Triangle
    do i= 1, 3
       do j= 1, 3
        Coefficient_Xi_Eta( j, i, e ) & 
        = BasisFunctionBC_Times_Difference_Position_Node( 2, i, e ) &
         *BasisFunctionBC_Times_Difference_Position_Node( 1, j, e ) &
         +BasisFunctionBC_Times_Difference_Position_Node( 1, i, e ) &
         *BasisFunctionBC_Times_Difference_Position_Node( 2, j, e ) 
       end do
    end do
   end do

   !$omp parallel do default( none )   &
   !$omp private( e, i, j ) &
   !$omp shared( Number_Element_Triangle, Coefficient_Xi_Xi ) & 
   !$omp shared( BasisFunctionBC_Times_Difference_Position_Node ) 
   do e= 1, Number_Element_Triangle
    do i= 1, 3
       do j= 1, 3
        Coefficient_Xi_Xi( j, i, e ) & 
        = BasisFunctionBC_Times_Difference_Position_Node( 1, i, e ) &
         *BasisFunctionBC_Times_Difference_Position_Node( 1, j, e ) 
       end do
    end do
   end do

   !$omp parallel do default( none )   &
   !$omp private( e, i, j ) &
   !$omp shared( Number_Element_Triangle, Coefficient_Eta_Eta ) & 
   !$omp shared( BasisFunctionBC_Times_Difference_Position_Node ) 
   do e= 1, Number_Element_Triangle
    do i= 1, 3
       do j= 1, 3
        Coefficient_Eta_Eta( j, i, e ) & 
        = BasisFunctionBC_Times_Difference_Position_Node( 2, i, e ) &
         *BasisFunctionBC_Times_Difference_Position_Node( 2, j, e ) 
       end do
    end do
   end do
   
   deallocate( BasisFunctionBC_Times_Difference_Position_Node )

   !==============================================================================================================
   write(*,*)'        Class_Element --> Thermal_Conductivity '
   !==============================================================================================================

   !allocate( Thermal_Conductivity_Triangle( Number_Element_Triangle ) )

   do i= 1, Number_Element_Triangle
!    if( Class_Element_Triangle( i ) == ID_Element_Material )then
!       Thermal_Conductivity_Triangle( i )= Thermal_Conductivity_Material

!    else if( Class_Element_Triangle( i ) == ID_Element_OuterDomain )then
!       Thermal_Conductivity_Triangle( i )= Thermal_Conductivity_OuterDomain !Thermal_Conductivity_Material_DesignDomain

!    else if( Class_Element_Triangle( i ) == ID_Element_FixedDomain )then
!       Thermal_Conductivity_Triangle( i )= Thermal_Conductivity_FixedDomain_Variable

!    else if( Class_Element_Triangle( i ) == ID_Element_Base_Material )then
!       Thermal_Conductivity_Triangle( i )= Thermal_Conductivity_Material_DesignDomain

!    else if( Class_Element_Triangle( i ) == ID_Element_Material_Exterior )then
!       Thermal_Conductivity_Triangle( i )= Thermal_Conductivity_OuterDomain !Thermal_Conductivity_Material

!    else if( Class_Element_Triangle( i ) == ID_Element_Base_Material_Exterior )then
!       Thermal_Conductivity_Triangle( i )= Thermal_Conductivity_OuterDomain !Thermal_Conductivity_Material_DesignDomain

!    else
!       write(*,*)'ERROR'
!       write(*,*)'  Class_Element_Triangle( i )=', Class_Element_Triangle( i )
!       write(*,*)'  Analyze_Steady_State_Heat_Conduction.f90 152'
!    end if
!   end do 

   !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!    else
      Thermal_Conductivity_Triangle( i ) = 1.0d0 + ( high/thickness ) *Thermal_conductivity_by_density( i )**penalty_density
!    end if

      if ( Flag_Structure == 2 ) then
         Thermal_Conductivity_Triangle( i ) = 1.0d0
      else if( Flag_Structure == 1 ) then 
         Thermal_Conductivity_Triangle( i ) = 1.0d0 
         if ( Class_Element_Triangle( i ) == 2 ) then
              Thermal_Conductivity_Triangle( i ) = 0.0d0
         end if
      end if
          
   end do

   if ( Flag_Structure == 1 ) then
      open( 909, file='ahohogeaho.dat' ) 
         do i = 1, Number_Element_Triangle   
            write( 909,* ) Class_Element_Triangle( i ), Thermal_Conductivity_Triangle( i ) 
         end do
      close( 909 )
   end if

   !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 

   allocate( Flag_Element_Triangle_Removed( Number_Element_Triangle ) ) 

   allocate( LocalMatrix_Triangle( 3, 3, Number_Element_Triangle ) )

   LocalMatrix_Triangle= Zero 

   !==============================================================================================================
   if( Flag_Thermal_Insulation_BC==Flag_On )then
   !==============================================================================================================

    !==============================================================================================================
    ! Local Matrix in Triangle Element 
    !==============================================================================================================
   
    !$omp parallel do default( none )   &
    !$omp private( e, i, j ) &
    !$omp shared( Number_Element_Triangle, LocalMatrix_Triangle, ID_Element_Thermal_Insulator ) & 
    !$omp shared( Area_Element_Triangle, Thermal_Conductivity_Triangle ) &
    !$omp shared( BasisFunction_B, BasisFunction_C, Class_Element_Triangle )
    do e= 1, Number_Element_Triangle
       if( Class_Element_Triangle( e ) /= ID_Element_Thermal_Insulator )then
        do i= 1, 3
         do j= 1, 3
            LocalMatrix_Triangle( j, i, e ) & 
            = -Thermal_Conductivity_Triangle( e ) &
             *Area_Element_Triangle( e ) & 
             *( BasisFunction_B( i, e )*BasisFunction_B( j, e ) &
             +BasisFunction_C( i, e )*BasisFunction_C( j, e ) ) 
         end do
        end do
       end if
    end do

    !==============================================================================================================
    write(*,*)'    =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-='
    write(*,*)'    Implement Neumann Boundary Condition for Thermal Insulation of Cloaked Region'
    write(*,*)'    =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-='
    !==============================================================================================================

    allocate( Flag_Element_Triangle_on_Thermal_Insulator_BC( Number_Element_Triangle ) )
    allocate( Position_Node_PEC_BC( 2, 2, Number_Element_Triangle ) )

    call Detect_Node_and_Element_Triangle_on_Thermal_Insulator_BC &
       ( Number_Node, Number_Element_Triangle, &
       Index_Element_2_Node_Triangle, &
       ID_Element_Thermal_Insulator, Class_Element_Triangle, &
       Position_Node, &
       !================================================================
       Flag_Element_Triangle_on_Thermal_Insulator_BC, Number_Outer_Element_on_PEC_BC, &
       Position_Node_PEC_BC )

    allocate( Length_Edge_on_PEC_BC_in_Element( Number_Element_Triangle ) )

    do e= 1, Number_Element_Triangle
       if( Flag_Element_Triangle_on_Thermal_Insulator_BC( e )==Flag_On )then 
        Length_Edge_on_PEC_BC_in_Element( e ) &
        = sqrt( ( Position_Node_PEC_BC( 1, 2, e ) -Position_Node_PEC_BC( 1, 1, e ) )**2 &
           +( Position_Node_PEC_BC( 2, 2, e ) -Position_Node_PEC_BC( 2, 1, e ) )**2 )
       else if( Flag_Element_Triangle_on_Thermal_Insulator_BC( e )==Flag_Off )then 
        Length_Edge_on_PEC_BC_in_Element( e )= -1d10
       end if
    end do

    allocate( Position_Node_PEC_BC_tmp( 2, 2, Number_Outer_Element_on_PEC_BC ) )

    Counter= 0
    do e= 1, Number_Element_Triangle
       if( Flag_Element_Triangle_on_Thermal_Insulator_BC( e )==Flag_On )then 
        Counter= Counter +1
        do i= 1, 2
         do j= 1, 2
            Position_Node_PEC_BC_tmp( j, i, Counter )= Position_Node_PEC_BC( j, i, e )
         end do
        end do
       end if
    end do

    if( Counter/=Number_Outer_Element_on_PEC_BC )then
       write(*,*)'Number_Outer_Element_on_PEC_BC=', Number_Outer_Element_on_PEC_BC
       write(*,*)'Counter=', Counter 
       call Output_Error( 'Analyze_Steady_State_Heat_Conduction', 868 )
    end if

    call Change_Direction_Vector_for_Unit_Normal_Vector &
       ( Number_Outer_Element_on_PEC_BC, Position_Node_PEC_BC_tmp, 1 )

    allocate( Unit_Normal_Vector_on_PEC_BC_tmp( 2, Number_Outer_Element_on_PEC_BC ) )

    call Compute_Unit_Normal_Vector &
       ( Number_Outer_Element_on_PEC_BC, Position_Node_PEC_BC_tmp, Unit_Normal_Vector_on_PEC_BC_tmp, 0 )

    deallocate( Position_Node_PEC_BC_tmp )

    allocate( Unit_Normal_Vector_on_PEC_BC( 2, Number_Element_Triangle ) )

    Counter= 0
    do e= 1, Number_Element_Triangle
       if( Flag_Element_Triangle_on_Thermal_Insulator_BC( e )==Flag_On )then 
        Counter= Counter +1
        do i= 1, 2
         Unit_Normal_Vector_on_PEC_BC( i, e )= Unit_Normal_Vector_on_PEC_BC_tmp( i, Counter )
        end do
       else
        Unit_Normal_Vector_on_PEC_BC( 1:2, e )= 0d0 
       end if
    end do

    if( Counter/=Number_Outer_Element_on_PEC_BC )then
       write(*,*)'Number_Outer_Element_on_PEC_BC=', Number_Outer_Element_on_PEC_BC
       write(*,*)'Counter=', Counter 
       call Output_Error( 'Analyze_Steady_State_Heat_Conduction', 898 )
    end if

    deallocate( Unit_Normal_Vector_on_PEC_BC_tmp )

    !!Nabla T \cdot n = 0
    !do e= 1, Number_Element_Triangle
    !   if( Flag_Element_Triangle_on_Thermal_Insulator_BC( e )==Flag_On )then
    !    do i= 1, 3
    !       do j= 1, 3
    !        LocalMatrix_Triangle( j, i, e ) &
    !        = LocalMatrix_Triangle( j, i, e ) &
    !         +Thermal_Conductivity_Triangle( e ) &
    !         *Length_Edge_on_PEC_BC_in_Element( e ) &
    !         *( BasisFunction_B( j, e ) *Unit_Normal_Vector_on_PEC_BC( 1, e ) &
    !          +BasisFunction_C( j, e ) *Unit_Normal_Vector_on_PEC_BC( 2, e ) ) &
    !         *( BasisFunction_A( i, e ) &
    !         +BasisFunction_B( i, e )/2d0 &
    !          *( Position_Node_PEC_BC( 1, 1, e ) +Position_Node_PEC_BC( 1, 2, e ) ) &
    !         +BasisFunction_C( i, e )/2d0 &
    !          *( Position_Node_PEC_BC( 2, 1, e ) +Position_Node_PEC_BC( 2, 2, e ) ) )
    !       end do
    !    end do
    !   end if
    !end do

    deallocate( Length_Edge_on_PEC_BC_in_Element )
    deallocate( Unit_Normal_Vector_on_PEC_BC )
    deallocate( Position_Node_PEC_BC )
    deallocate( Flag_Element_Triangle_on_Thermal_Insulator_BC )

    !==============================================================================================================
    if( Flag_Thermal_Device==0 .or. Flag_Thermal_Device==10 .or. &
        Flag_Thermal_Device==20 .or. Flag_Thermal_Device==21 .or. Flag_Thermal_Device==30 )then
       write(*,*)'    =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-='
       !write(*,*)'    Implement Neumann Boundary Condition for Thermal cloak'
       write(*,*)'    Compute Boundary Integral of Direcrit Boundary'
       write(*,*)'    =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-='
    !==============================================================================================================

       Number_Element_on_Flat_PEC_tmp = Number_Grid_X_Scattering*4

       allocate( Element_Number_on_Flat_PEC( Number_Element_on_Flat_PEC_tmp ) )
       allocate( Node_Number_on_Flat_PEC_in_Element( 2, Number_Element_on_Flat_PEC_tmp ) )
   
       call Detect_Node_and_Element_Triangle_for_Boundary_Integral &
        ( Number_Node, Position_Node, &
          Number_Element_Triangle, Index_Element_2_Node_Triangle, &
          ID_Element_Thermal_Insulator, Class_Element_Triangle, &
          Number_Element_on_Flat_PEC_tmp, &
          Max_Position_Node, &
        !================================================================
          Element_Number_on_Flat_PEC, Number_Element_on_Flat_PEC, &
          Node_Number_on_Flat_PEC_in_Element )

       write(*,*)'Number_Element_on_Flat_PEC=', Number_Element_on_Flat_PEC

       allocate( Position_Node_PEC_BC( 2, 2, Number_Element_on_Flat_PEC ) )

       do e= 1, Number_Element_on_Flat_PEC
        do i= 1, 2
         do j= 1, 2
            Position_Node_PEC_BC( j, i, e )= Position_Node( j, Node_Number_on_Flat_PEC_in_Element( i, e ) )
         end do
        end do
       end do

       deallocate( Node_Number_on_Flat_PEC_in_Element )
       allocate( Unit_Normal_Vector_on_Flat_PEC_BC( 2 ) )

       Unit_Normal_Vector_on_Flat_PEC_BC( 1 )= 0d0
       Unit_Normal_Vector_on_Flat_PEC_BC( 2 )= -1d0

       Length_Edge_on_Flat_PEC_BC_in_Element= 1d0/Length_Characteristic

       do e= 1, Number_Element_on_Flat_PEC
        do i= 1, 3
         do j= 1, 3
            LocalMatrix_Triangle( j, i, Element_Number_on_Flat_PEC( e ) ) & 
            = LocalMatrix_Triangle( j, i, Element_Number_on_Flat_PEC( e ) ) & 
             +Thermal_Conductivity_Triangle( e ) &
             *Length_Edge_on_Flat_PEC_BC_in_Element & 
             *( BasisFunction_B( j, Element_Number_on_Flat_PEC( e ) ) *Unit_Normal_Vector_on_Flat_PEC_BC( 1 ) &
              +BasisFunction_C( j, Element_Number_on_Flat_PEC( e ) ) *Unit_Normal_Vector_on_Flat_PEC_BC( 2 ) ) &
             *( BasisFunction_A( i, Element_Number_on_Flat_PEC( e ) ) &
             +BasisFunction_B( i, Element_Number_on_Flat_PEC( e ) )/2d0 &
              *( Position_Node_PEC_BC( 1, 1, e ) +Position_Node_PEC_BC( 1, 2, e ) ) & 
             +BasisFunction_C( i, Element_Number_on_Flat_PEC( e ) )/2d0 & 
              *( Position_Node_PEC_BC( 2, 1, e ) +Position_Node_PEC_BC( 2, 2, e ) ) )  
         end do
        end do
       end do

       deallocate( Position_Node_PEC_BC )
       deallocate( Unit_Normal_Vector_on_Flat_PEC_BC )
       deallocate( Element_Number_on_Flat_PEC )
    end if

    !==============================================================================================================
    ! Create Flag_Element_Removed 
    !==============================================================================================================
    
    Flag_Element_Triangle_Removed= Flag_Off
   
    do e= 1, Number_Element_Triangle
       if( Class_Element_Triangle( e ) == ID_Element_FixedDomain )then
        Flag_Element_Triangle_Removed( e )= Flag_On
       end if
    end do

   !==============================================================================================================
   else if( Flag_Thermal_Insulation_BC==Flag_Off )then
   !==============================================================================================================
 
    !==============================================================================================================
    ! Local Matrix in Triangle Element 
    !==============================================================================================================
    !$omp parallel do default( none )   &
    !$omp private( e, i, j ) &
    !$omp shared( Number_Element_Triangle, LocalMatrix_Triangle ) & 
    !$omp shared( Area_Element_Triangle, Thermal_Conductivity_Triangle ) &
    !$omp shared( BasisFunction_B, BasisFunction_C, Class_Element_Triangle )
    do e= 1, Number_Element_Triangle
       do i= 1, 3
        do j= 1, 3
         LocalMatrix_Triangle( j, i, e ) & 
         = -Thermal_Conductivity_Triangle( e ) &
            *Area_Element_Triangle( e ) & 
            *( BasisFunction_B( i, e )*BasisFunction_B( j, e ) &
            +BasisFunction_C( i, e )*BasisFunction_C( j, e ) ) 
        end do
       end do
    end do

    Flag_Element_Triangle_Removed= Flag_Off
   !==============================================================================================================
   else ! Flag_Thermal_Insulation_BC /= 0, 1
   !==============================================================================================================
    call Output_Error( 'Analyze_Steady_State_Heat_Conduction', 838 )
   end if

   deallocate( Area_Element_Triangle )
   
   !==============================================================================================================
   ! Local Matrix --> Global Matrix 
   !==============================================================================================================

   allocate( GlobalMatrix( Number_Node, Width_Matrix_LHS ) ) 
   allocate( J_GlobalMatrix( Number_Node, Width_Matrix_LHS ) )
 
   !call Format_Local_Global_ZMUMPS&
   call Format_Local_Global_DMUMPS&
    ( LocalMatrix_Triangle, Index_Element_2_Node_Triangle, &
      3, Number_Node, Number_Element_Triangle, &
      1, Width_Matrix_LHS, Flag_Symmetric, &
      Class_Element_Triangle, Flag_Element_Triangle_Removed, Flag_Remove_Element, &
      !=====================================================
      GlobalMatrix, J_GlobalMatrix, Number_NonZero )

   deallocate( LocalMatrix_Triangle )
   deallocate( Flag_Element_Triangle_Removed ) 

   !==================================================================================================
   write(*,*)'       '
   write(*,*)'       =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-='
   write(*,*)'       Dirichlet Boundary Condition of Higher & Lower Temperatures'
   write(*,*)'       =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-='
   !==================================================================================================

   allocate( Global_Vector_RHS( Number_Node ) )
   Global_Vector_RHS= Zero

   allocate( Nodal_Value_on_Dirichlet_BC( Number_Node ) )
   allocate( Flag_Node_Dirichlet_BC( Number_Node ) )

   Nodal_Value_on_Dirichlet_BC= Zero 
   Flag_Node_Dirichlet_BC= Flag_Off

   if( Flag_Thermal_Device==0 .or. Flag_Thermal_Device==20 .or. Flag_Thermal_Device==30 )then
      !Lower Temperatures
     do i= 1, Number_Node 
      if( Max_Position_Node( 1, 1 ) -1d-8 <= Position_Node( 1, i ) .and. &
        Position_Node( 1, i ) <= Max_Position_Node( 1, 1 ) +1d-8  )then
         Flag_Node_Dirichlet_BC( i )= Flag_On
         Nodal_Value_on_Dirichlet_BC( i )= Temperature_Lower_Side_Boundary
      end if
  
      !Higher Temperatures
      if( Max_Position_Node( 1, 2 ) -1d-8 <= Position_Node( 1, i ) .and. &
        Position_Node( 1, i ) <= Max_Position_Node( 1, 2 ) +1d-8  )then
         Flag_Node_Dirichlet_BC( i )= Flag_On
         Nodal_Value_on_Dirichlet_BC( i )= Temperature_Higher_Side_Boundary
      end if
     end do
   else if(Flag_Thermal_Device==21 )then
      !Lower Temperatures
     do i= 1, Number_Node 
      if( Max_Position_Node( 1, 1 ) -1d-8 <= Position_Node( 1, i ) .and. &
        Position_Node( 1, i ) <= Max_Position_Node( 1, 1 ) +1d-8  )then
         Flag_Node_Dirichlet_BC( i )= Flag_On
         Nodal_Value_on_Dirichlet_BC( i )= Temperature_Lower_Side_Boundary
      end if

      if( Max_Position_Node( 1, 2 ) -1d-8 <= Position_Node( 1, i ) .and. &
        Position_Node( 1, i ) <= Max_Position_Node( 1, 2 ) +1d-8  )then
         Flag_Node_Dirichlet_BC( i )= Flag_On
         Nodal_Value_on_Dirichlet_BC( i )= Temperature_Lower_Side_Boundary
      end if
   
      if( Max_Position_Node( 2, 1 ) -1d-8 <= Position_Node( 2, i ) .and. &
        Position_Node( 2, i ) <= Max_Position_Node( 2, 1 ) +1d-8  )then
         Flag_Node_Dirichlet_BC( i )= Flag_On
         Nodal_Value_on_Dirichlet_BC( i )= Temperature_Lower_Side_Boundary
      end if

      if( Max_Position_Node( 2, 2 ) -1d-8 <= Position_Node( 2, i ) .and. &
        Position_Node( 2, i ) <= Max_Position_Node( 2, 2 ) +1d-8  )then
         Flag_Node_Dirichlet_BC( i )= Flag_On
         Nodal_Value_on_Dirichlet_BC( i )= Temperature_Lower_Side_Boundary
      end if

       ! Source 
       if( match( Position_Source_X, Position_Node( 1, i ) )=='y'.and. &
           match( Position_Source_Y, Position_Node( 2, i ) )=='y' )then
         Flag_Node_Dirichlet_BC( i )= Flag_On
         Nodal_Value_on_Dirichlet_BC( i )= Temperature_Higher_Side_Boundary
       end if

     end do
   else if( Flag_Thermal_Device==10 )then
     do i= 1, Number_Node
       !Higher Temperatures
       if( Max_Position_Node( 2, 1 ) -1d-8 <= Position_Node( 2, i ) .and. &
         Position_Node( 2, i ) <= Max_Position_Node( 2, 1 ) +1d-8  )then
         Flag_Node_Dirichlet_BC( i )= Flag_On
         Nodal_Value_on_Dirichlet_BC( i )= Temperature_Lower_Side_Boundary
       end if

       ! Source 
       if( Position_Source_X -1d-8 <= Position_Node( 1, i ) .and. &
           Position_Node( 1, i ) <= Position_Source_X +1d-8 .and. &
           Position_Source_Y -1d-8 <= Position_Node( 2, i ) .and. &
           Position_Node( 2, i ) <= Position_Source_Y +1d-8 )then
         Flag_Node_Dirichlet_BC( i )= Flag_On
         Nodal_Value_on_Dirichlet_BC( i )= Temperature_Higher_Side_Boundary
       end if
     end do
   end if



   !call Implement_Dirichlet_BC_complex_CSRF &
   call Implement_Dirichlet_BC_real_CSRF &
    ( Nodal_Value_on_Dirichlet_BC, Flag_Node_Dirichlet_BC, &
      Number_Node, Width_Matrix_LHS, &
      !======================================================
      GlobalMatrix, J_GlobalMatrix, Global_Vector_RHS )

   deallocate( Flag_Node_Dirichlet_BC )
   deallocate( Nodal_Value_on_Dirichlet_BC )

   !==============================================================================================================
   if( Flag_Thermal_Insulation_BC==Flag_On )then
   !==============================================================================================================
      !==================================================================================================
      write(*,*)'       '
      write(*,*)'       =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-='
      write(*,*)'       Fixed Nodal Vlaue in Thremal Insulator'
      write(*,*)'       =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-='
      !==================================================================================================
   
      allocate( Flag_Node_Dirichlet_BC( Number_Node ) )
   
      do i= 1, Number_Node
       Flag_Node_Dirichlet_BC( i )= Flag_Off
      end do
   
      do e= 1, Number_Element_Triangle
       if( Class_Element_Triangle( e ) == ID_Element_Thermal_Insulator )then
          do i= 1, 3
           Flag_Node_Dirichlet_BC( Index_Element_2_Node_Triangle( i, e ) )= Flag_On
          end do
       end if
      end do
   
      if( Flag_Thermal_Device==0 )then
        do e= 1, Number_Element_Triangle
         if( Class_Element_Triangle( e ) /= ID_Element_Thermal_Insulator )then
            do i= 1, 3
             Flag_Node_Dirichlet_BC( Index_Element_2_Node_Triangle( i, e ) )= Flag_Off
            end do
         end if
        end do
      end if
   
      Counter= 0
      do i= 1, Number_Node
       if( Flag_Node_Dirichlet_BC( i )==Flag_On )then
          Counter= Counter +1
       end if
      end do
   
      allocate( Nodal_Value_Fixed_in_PEC( Number_Node ) )
    
      do i= 1, Number_Node 
       Nodal_Value_Fixed_in_PEC( i )= 0d0 
      end do
   
      !if( Flag_Thermal_Device/=20 )then
         if( Counter > 0 )then
          !call Implement_Dirichlet_BC_complex_CSRF &
          call Implement_Dirichlet_BC_real_CSRF &
             ( Nodal_Value_Fixed_in_PEC, Flag_Node_Dirichlet_BC, &
             Number_Node, Width_Matrix_LHS, &
            !========================================================
             GlobalMatrix, J_GlobalMatrix, Global_Vector_RHS )
         end if
      !end if
   
      deallocate( Nodal_Value_Fixed_in_PEC )
      deallocate( Flag_Node_Dirichlet_BC )

   end if
   !==================================================================================================
   write(*,*)'       '
   write(*,*)'       =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-='
   write(*,*)'       Count Number Non-Zero '
   write(*,*)'       =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-='
   !==================================================================================================

   Number_NonZero= 0

   do j= 1, Width_Matrix_LHS
    do i= 1, Number_Node
       if( Flag_Symmetric==Flag_On )then
        if( i >= J_GlobalMatrix( i, j ) .and. J_GlobalMatrix( i, j ) > 0 )then
         Number_NonZero= Number_NonZero +1
        end if
       else if( Flag_Symmetric==Flag_Off )then
        if( 0 < J_GlobalMatrix( i, j ) .and. J_GlobalMatrix( i, j ) <= Number_Node  )then
         Number_NonZero= Number_NonZero +1
        end if
       else
        write(*,*)'Flag_Symmetric=', Flag_Symmetric
        call Output_Error( 'Analyze_Steady_State_Heat_Conduction', 1118 )
       end if
    end do
   end do

   write(*,*)'       Number_NonZero=', Number_NonZero
   write(*,*)'       Number_Node=',  Number_Node 
   write(*,*)'       Width_Matrix_LHS=', Width_Matrix_LHS

   !==================================================================================================
   write(*,*)'       Global Matrix --> MUMPS Format '
   !==================================================================================================

   allocate( aa_zmumps( Number_NonZero ) )
   allocate( ia_zmumps( Number_NonZero ) )
   allocate( ja_zmumps( Number_NonZero ) )

   !call Format_Global_aa_ZMUMPS&
   call Format_Global_aa_DMUMPS&
    ( GlobalMatrix, J_GlobalMatrix, & 
      Number_Node, Width_Matrix_LHS, Number_NonZero, Flag_Symmetric, Flag_Matrix_Upper, & 
      !==================================================== 
      aa_zmumps, ia_zmumps, ja_zmumps ) 
 
   write(*,*)'Flag_Symmetric=', Flag_Symmetric, 'at Main.f90' 

   if( Flag_Thermal_Insulation_BC==Flag_On )then 
 
    !call Check_ZMUMPS_Format & 
    call Check_DMUMPS_Format & 
       ( aa_zmumps, ia_zmumps, ja_zmumps, &  
       Number_Node, Number_NonZero, Flag_Symmetric, Flag_Matrix_Upper ) 
 
   else if( Flag_Thermal_Insulation_BC==Flag_Off )then 
 
    !call Check_ZMUMPS_Format & 
    call Check_DMUMPS_Format & 
       ( aa_zmumps, ia_zmumps, ja_zmumps, &  
       Number_Node, Number_NonZero, Flag_Symmetric, Flag_Matrix_Upper ) 
 
   end if 

   deallocate( GlobalMatrix )
   deallocate( J_GlobalMatrix )

   !==================================================================================================
   write(*,*)'       Solve Linear Equation'
   write(*,*)'        Number_Node=', Number_Node
   write(*,*)'        Number_NonZero=', Number_NonZero
   !==================================================================================================

   !call ME57( Number_Node, Number_NonZero, aa_zmumps, ia_zmumps, ja_zmumps, Global_Vector_RHS, Temperature_Solution )
   call MA57( Number_Node, Number_NonZero, aa_zmumps, ia_zmumps, ja_zmumps, Global_Vector_RHS, Temperature_Solution )

!   do i= 1, Number_Node
!      write(910,'(1x,es15.7)') Temperature_Solution( i )
!   end do
!
!   do e= 1, Number_Element_Triangle
!      if( Class_Element_Triangle( e ) == ID_Element_FixedDomain )then
!         do i= 1, 3
!            write(911,'(1x,es15.7)')Temperature_Solution( Index_Element_2_Node_Triangle( i, e ) )
!         end do
!      end if 
!   end do

   deallocate( aa_zmumps )
   deallocate( ia_zmumps )
   deallocate( ja_zmumps )

   deallocate( Global_Vector_RHS )

   !deallocate( Thermal_Conductivity_Triangle )
   
   deallocate( Position_Node_Element_Triangle )
   
   deallocate( BasisFunction_A )
   deallocate( BasisFunction_B )
   deallocate( BasisFunction_C )
    
   deallocate( Coefficient_C )
   deallocate( Coefficient_Xi )
   deallocate( Coefficient_Eta )
   deallocate( Coefficient_Xi_Eta )
   deallocate( Coefficient_Xi_Xi )
   deallocate( Coefficient_Eta_Eta ) 

   deallocate( Index_Element_2_Node_Triangle_Check )
   
   !==================================================================================================
   write(*,*)'       end subroutine Analyze_Steady_State_Heat_Conduction'
   !==================================================================================================
   return
end subroutine Analyze_Steady_State_Heat_Conduction



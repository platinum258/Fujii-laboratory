

subroutine Analyze_Light_Scattering_TE &
       ( Flag_Incident_Wave, Flag_PEC_Boundary_Condition, Flag_Material, &
         Frequency_Normalized, Incident_Angle, Frequency_Plasma_Normalized, Coefficient_Dumping_Drude, &
         Dielectric_Constant_Material, Dielectric_Constant_Air, & 
         Dielectric_Constant_OpenRegion, Dielectric_Constant_FixedDomain, &
         Dielectric_Constant_PML_Material, Dielectric_Constant_PML_Air, &
         Position_Dipole_X, Position_Dipole_Y, &
         ID_Element_Material, ID_Element_OpenRegion, ID_Element_FixedDomain, ID_Element_Air, ID_Element_Active, & 
         ID_Element_Material_Exterior, ID_Element_Air_Exterior, & 
         ID_Element_PML_Xp, ID_Element_PML_Xm, ID_Element_PML_Yp, ID_Element_PML_Ym, &
         ID_Element_PML_XpYp,  ID_Element_PML_XmYp, ID_Element_PML_XpYm,  ID_Element_PML_XmYm, &
         ID_Element_PML_Material, ID_Element_PML_Air, &
         !
         Number_Node, Number_Element_Triangle, Number_Element_Square,  &
         Position_Node, Class_Node, &
         Number_Node_Dirichlet_Boundary_PML, &
         Width_Matrix_LHS, Max_Position_Node, &
         Index_Element_2_Node_Triangle, Class_Element_Triangle, &
         Index_Element_2_Node_Square, Class_Element_Square, Class_Element_PML_Material, &
         ID_Node_Material, ID_Node_FixedDomain, ID_Node_Active, ID_Node_Air, ID_Node_OpenRegion, &
         ID_Node_PML_X, ID_Node_PML_Y, ID_Node_PML_XY, &
         !============================================================================================
         Field_Magnetic_All, Field_Magnetic_Incident, Field_Magnetic_Scattering )

   !$ use omp_lib
   use Parameters
   implicit none

   ! intent(in)
   double precision, intent(in) :: Frequency_Normalized, Incident_Angle 
   double precision, intent(in) :: Frequency_Plasma_Normalized, Coefficient_Dumping_Drude

   integer, intent(in) :: Flag_Incident_Wave, Flag_PEC_Boundary_Condition, Flag_Material 
   integer, intent(in) :: Number_Node, Number_Node_Dirichlet_Boundary_PML 

   double precision, intent(in) :: Position_Node( 2, Number_Node ) 
   integer, intent(in) :: Class_Node( Number_Node ) 
 
   complex( kind(0d0) ), intent(in) :: Dielectric_Constant_Material, Dielectric_Constant_Air  
   complex( kind(0d0) ), intent(in) :: Dielectric_Constant_OpenRegion, Dielectric_Constant_FixedDomain  

   double precision :: Dielectric_Constant_PML_Material, Dielectric_Constant_PML_Air

   double precision, intent(in) :: Position_Dipole_X, Position_Dipole_Y 

   integer, intent(in) :: Number_Element_Triangle, Number_Element_Square 
   integer, intent(in) :: Index_Element_2_Node_Triangle( 3, Number_Element_Triangle ) 
   integer, intent(in) :: Index_Element_2_Node_Square( 4, Number_Element_Square )
   integer, intent(in) :: Class_Element_Triangle( Number_Element_Triangle )
   integer, intent(in) :: Class_Element_Square( Number_Element_Square )
   integer, intent(in) :: Class_Element_PML_Material( Number_Element_Square )
 
   integer, intent(in) :: ID_Element_Material, ID_Element_OpenRegion, ID_Element_FixedDomain, ID_Element_Air, ID_Element_Active  
   integer, intent(in) :: ID_Element_Material_Exterior, ID_Element_Air_Exterior
   integer, intent(in) :: ID_Element_PML_Xp, ID_Element_PML_Xm, ID_Element_PML_Yp, ID_Element_PML_Ym
   integer, intent(in) :: ID_Element_PML_XpYp,  ID_Element_PML_XmYp, ID_Element_PML_XpYm,  ID_Element_PML_XmYm
   integer, intent(in) :: ID_Element_PML_Material, ID_Element_PML_Air
   integer, intent(in) :: ID_Node_Material, ID_Node_FixedDomain, ID_Node_Active, ID_Node_Air, ID_Node_OpenRegion 
   integer, intent(in) :: ID_Node_PML_X, ID_Node_PML_Y, ID_Node_PML_XY 
 
   double precision, intent(in) :: Max_Position_Node( 2, 2 ) 
   integer, intent(in) :: Width_Matrix_LHS
 
   ! intent(out)
   complex( kind( 0d0 ) ), intent(out) :: Field_Magnetic_Incident( Number_Node ) 
   complex( kind( 0d0 ) ), intent(out) :: Field_Magnetic_Scattering( Number_Node )
   complex( kind( 0d0 ) ), intent(out) :: Field_Magnetic_All( Number_Node )
 
   integer :: i, j, k, e
   integer :: Counter
   integer :: NumThread_OMP, NumThread_MKL, NumThread_MPI

   integer, allocatable, dimension(:,:) :: Index_Element_2_Node_Triangle_Check
   double precision, allocatable, dimension(:) :: Area_Element_Triangle

   ! Model Data

   complex( kind( 0d0 ) ), allocatable, dimension(:) :: Dielectric_Constant_Triangle
   double precision, allocatable, dimension(:) :: Dielectric_Constant_PML
 
   ! Finite Element Method 
   double precision, allocatable, dimension(:,:,:) :: Position_Node_Element_Triangle
   double precision, allocatable, dimension(:,:) :: BasisFunction_A, BasisFunction_B, BasisFunction_C
   double precision, allocatable, dimension(:,:,:,:) :: Difference_Position_Node_Element_Triangle
   double precision, allocatable, dimension(:,:,:) :: BasisFunctionBC_Times_Difference_Position_Node_TE
   double precision, allocatable, dimension(:,:) :: Basis_Function_Times_Position_Node_TE
   complex( kind( 0d0 ) ), allocatable, dimension(:,:,:) :: Coefficient_C, Coefficient_Xi, Coefficient_Eta
   complex( kind( 0d0 ) ), allocatable, dimension(:,:,:) :: Coefficient_Xi_Eta, Coefficient_Xi_Xi, Coefficient_Eta_Eta
   
   double precision, allocatable, dimension(:,:,:) :: Position_Node_Element_Square
   double precision, allocatable, dimension(:,:) :: x4Ele, y4Ele
   double precision, allocatable, dimension(:,:) :: Ni
   double precision, allocatable, dimension(:) :: Xii, Itai
   double precision, allocatable, dimension(:,:) :: Jaco
   double precision, allocatable, dimension(:) :: Xig, Itag, Wg
   double precision, allocatable, dimension(:,:) :: PxPXi, PxPeta, PyPXi, PyPeta
   double precision, allocatable, dimension(:,:) :: PNPXi, PNPeta
   double precision, allocatable, dimension(:,:,:) :: PNPx, PNPy
   complex( kind( 0d0 ) ), allocatable, dimension(:,:) :: GammaX, GammaY

   complex( kind( 0d0 ) ), allocatable, dimension(:,:,:) :: LocalMatrix_Triangle_Scattering 
   complex( kind( 0d0 ) ), allocatable, dimension(:,:,:) :: LocalMatrix_Triangle_Scattering_RHS 
   complex( kind( 0d0 ) ), allocatable, dimension(:,:,:) :: LocalMatrix_Triangle_Scattering_tmp

   complex( kind( 0d0 ) ), allocatable, dimension(:,:,:) :: LocalMatrix_Square_PML
   complex( kind( 0d0 ) ), allocatable, dimension(:,:,:,:) :: LocalMatrix_Square_PML_tmp

   !PEC
   integer, allocatable, dimension(:) :: Flag_Node_PEC_BC, Flag_Node_PEC_BC_Fixed, Flag_Element_Triangle_PEC_BC
   integer :: Number_Outer_Element_on_PEC_BC

   ! Flat_PEC for Carpet cloak
   integer, allocatable, dimension(:) :: Element_Number_on_Flat_PEC
   integer :: Number_Element_on_Flat_PEC

   integer, allocatable, dimension(:,:) :: Node_Number_on_Flat_PEC_in_Element

   double precision, allocatable, dimension(:,:,:) :: Position_Node_PEC_BC, Position_Node_PEC_BC_tmp
   double precision, allocatable, dimension(:,:) :: Unit_Normal_Vector_on_PEC_BC, Unit_Normal_Vector_on_PEC_BC_tmp
   double precision, allocatable, dimension(:) :: Unit_Normal_Vector_on_Flat_PEC_BC 
   complex( kind( 0d0 ) ), allocatable, dimension(:) :: Nodal_Value_Fixed_in_PEC 
   double precision, allocatable, dimension(:) :: Length_Edge_on_PEC_BC_in_Element 
   double precision :: Length_Edge_on_Flat_PEC_BC_in_Element 

   ! Linear equation
   complex( kind( 0d0 ) ), allocatable, dimension(:,:) :: GlobalMatrix_RHS 
   integer, allocatable, dimension(:,:) :: J_GlobalMatrix_RHS
   
   integer :: Number_NonZero, Number_NonZero_RHS 
   complex( kind( 0d0 ) ), allocatable, dimension(:,:) :: GlobalMatrix
   integer, allocatable, dimension(:,:) :: J_GlobalMatrix
   complex( kind( 0d0 ) ), allocatable, dimension(:) ::  Global_Vector_RHS 

   integer, allocatable, dimension(:) :: NodeNumber_Dirichlet_BC, Flag_Node_Dirichlet_BC
   complex( kind( 0d0 ) ), allocatable, dimension(:) :: Nodal_Value_on_Dirichlet_BC 
 
   ! MUMPS
   complex( kind( 0d0 ) ), allocatable, dimension(:) :: aa_zmumps
   integer, allocatable, dimension(:) :: ja_zmumps, ia_zmumps
   integer, allocatable, dimension(:) :: Flag_Element_Triangle_Removed, Flag_Element_Square_Removed

   ! Number of Integral Points of Gauss Legendre in Square Element
   integer, parameter :: Number_Point_Gauss_Legendre_Element_Square = 4
 
   !integer, parameter :: ID_Element_PEC= ID_Element_FixedDomain
   !integer, parameter :: ID_Element_Flat_PEC=ID_Element_PML_Ym 
   !integer, parameter :: Flag_Remove_Element= Flag_PEC_Boundary_Condition
   integer :: ID_Element_PEC, ID_Element_Flat_PEC, Flag_Remove_Element
   integer, parameter :: Flag_Symmetric= Flag_On
   integer, parameter :: Flag_Matrix_Upper= Flag_On

   ID_Element_PEC= ID_Element_FixedDomain
   ID_Element_Flat_PEC=ID_Element_PML_Ym 
   Flag_Remove_Element= Flag_PEC_Boundary_Condition

   !============================================================================================
   write(*,*)'   call Analyze_Light_Scattering_TE '
   !============================================================================================
    
   open(9999, file='Number_Thread')
   read(9999,*)NumThread_OMP!, NumThread_MKL, NumThread_MPI
   write(*,*)'NumThread_OMP=', NumThread_OMP
   !write(*,*)'NumThread_MKL=', NumThread_MKL
   !write(*,*)'NumThread_MPI=', NumThread_MPI
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
   write(*,*)'      =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= '
   write(*,*)'      Frequency_Normalized=', Frequency_Normalized
   write(*,*)'      Incident_Angle=', -Incident_Angle*360d0/( 2d0*Pi )
   write(*,*)'      =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= '
   !===============================================================================================

   !============================================================================================
   write(*,*)'      '
   write(*,*)'      <<<<<<<<<<<<<<<<<<<< Finite Element Method >>>>>>>>>>>>>>>>>>> '
   write(*,*)'      '
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
   !$omp shared( Number_Element_Triangle, Area_Element_Triangle, Position_Node_Element_Triangle ) 
   do e= 1, Number_Element_Triangle
      if( Area_Element_Triangle( e )==0.0d0 )then
         write(*,*)'Area_Element_Triangle( e )==0.0d0'
         write(*,*)'Element Number=', e
         call Output_Error( 'Analyze_Light_Scattering_TE', 166 )
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
   
   allocate( Basis_Function_Times_Position_Node_TE( 3, Number_Element_Triangle ) ) 

   !$omp parallel do default( none )   &
   !$omp private( e, i ) &
   !$omp shared( Number_Element_Triangle, Basis_Function_Times_Position_Node_TE ) &
   !$omp shared( BasisFunction_A, BasisFunction_B, BasisFunction_C, Position_Node_Element_Triangle )
   do e= 1, Number_Element_Triangle
      do i= 1, 3
         Basis_Function_Times_Position_Node_TE( i, e ) & 
         = BasisFunction_A( i, e ) & 
          +BasisFunction_B( i, e ) *Position_Node_Element_Triangle( 1, 3, e ) & 
          +BasisFunction_C( i, e ) *Position_Node_Element_Triangle( 2, 3, e ) 
      end do
   end do

   allocate( BasisFunctionBC_Times_Difference_Position_Node_TE( 2, 3, Number_Element_Triangle ) )

   !$omp parallel do default( none )   &
   !$omp private( e, i, j ) &
   !$omp shared( Number_Element_Triangle, BasisFunctionBC_Times_Difference_Position_Node_TE ) &
   !$omp shared( BasisFunction_A, BasisFunction_B, BasisFunction_C, Difference_Position_Node_Element_Triangle )
   do e= 1, Number_Element_Triangle
      do i= 1, 3
         do j= 1, 2
            BasisFunctionBC_Times_Difference_Position_Node_TE( j, i, e ) &
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
   !$omp shared( Basis_Function_Times_Position_Node_TE, Coefficient_C ) 
   do e= 1, Number_Element_Triangle
      do i= 1, 3
         do j= 1, 3
            Coefficient_C( j, i, e ) & 
            = Basis_Function_Times_Position_Node_TE( i, e )*Basis_Function_Times_Position_Node_TE( j, e ) 
         end do
      end do
   end do

   !$omp parallel do default( none )   &
   !$omp private( e, i, j ) &
   !$omp shared( Number_Element_Triangle, Coefficient_Xi ) & 
   !$omp shared( Basis_Function_Times_Position_Node_TE ) &
   !$omp shared( BasisFunctionBC_Times_Difference_Position_Node_TE ) 
   do e= 1, Number_Element_Triangle
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
   !$omp shared( Number_Element_Triangle, Coefficient_Eta ) & 
   !$omp shared( Basis_Function_Times_Position_Node_TE ) &
   !$omp shared( BasisFunctionBC_Times_Difference_Position_Node_TE ) 
   do e= 1, Number_Element_Triangle
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
   !$omp shared( Number_Element_Triangle, Coefficient_Xi_Eta ) & 
   !$omp shared( BasisFunctionBC_Times_Difference_Position_Node_TE ) 
   do e= 1, Number_Element_Triangle
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
   !$omp shared( Number_Element_Triangle, Coefficient_Xi_Xi ) & 
   !$omp shared( BasisFunctionBC_Times_Difference_Position_Node_TE ) 
   do e= 1, Number_Element_Triangle
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
   !$omp shared( Number_Element_Triangle, Coefficient_Eta_Eta ) & 
   !$omp shared( BasisFunctionBC_Times_Difference_Position_Node_TE ) 
   do e= 1, Number_Element_Triangle
      do i= 1, 3
         do j= 1, 3
            Coefficient_Eta_Eta( j, i, e ) & 
            = BasisFunctionBC_Times_Difference_Position_Node_TE( 2, i, e ) &
             *BasisFunctionBC_Times_Difference_Position_Node_TE( 2, j, e ) 
         end do
      end do
   end do
   
   deallocate( BasisFunctionBC_Times_Difference_Position_Node_TE )
 
   !==============================================================================================================
   ! Coordinate of Basis Function in Rectangle Element in PML
   !==============================================================================================================
   allocate( Xii( 4 ) )
   allocate( Itai( 4 ) )
   
   Xii( 1 )= -1d0
   Xii( 2 )=  1d0
   Xii( 3 )=  1d0
   Xii( 4 )= -1d0
   
   Itai( 1 )= -1d0
   Itai( 2 )= -1d0
   Itai( 3 )=  1d0
   Itai( 4 )=  1d0
   
   !==============================================================================================================
   ! Coordinate of Gauss Point for Gauss-Legendre Integration in Rectangle Element in PML
   !==============================================================================================================
   allocate( Xig( Number_Point_Gauss_Legendre_Element_Square ) )
   allocate( Itag( Number_Point_Gauss_Legendre_Element_Square ) )
   allocate( Wg( Number_Point_Gauss_Legendre_Element_Square ) )
   
   Xig( 1 )= -1d0/dsqrt( 3d0 )
   Xig( 2 )=  1d0/dsqrt( 3d0 )
   Xig( 3 )=  1d0/dsqrt( 3d0 )
   Xig( 4 )= -1d0/dsqrt( 3d0 )
   
   Itag( 1 )= -1d0/dsqrt( 3d0 )
   Itag( 2 )= -1d0/dsqrt( 3d0 )
   Itag( 3 )=  1d0/dsqrt( 3d0 )
   Itag( 4 )=  1d0/dsqrt( 3d0 )
   
   !==============================================================================================================
   ! Weight for Gauss-Legendre Integration in Rectangle Element in PML
   !==============================================================================================================
   do i= 1, Number_Point_Gauss_Legendre_Element_Square
      Wg( i )=1d0
   end do
   
   !==============================================================================================================
   ! (X,Y)_node --> (X,Y)_RectangleElement 
   !==============================================================================================================

   allocate( Position_Node_Element_Square( 2, 4, Number_Element_Square ) )
    
   !$omp parallel do default( none )   &
   !$omp private( E, I ) &
   !$omp shared(  Number_Element_Square, Position_Node, Position_Node_Element_Square, Index_Element_2_Node_Square ) 
   do e= 1, Number_Element_Square
      do i= 1, 4
         Position_Node_Element_Square( 1, i, e )= Position_Node( 1, Index_Element_2_Node_Square( i, e ) )
         Position_Node_Element_Square( 2, i, e )= Position_Node( 2, Index_Element_2_Node_Square( i, e ) )
      end do
   end do
    
   !==============================================================================================================
   ! Basis Function in Rectangle Element 
   !==============================================================================================================
    
   allocate( Ni( Number_Point_Gauss_Legendre_Element_Square, 4 ) )

   !$omp parallel do default( none )   &
   !$omp private(I, J) &
   !$omp shared( Xii, Xig, Itai, Itag, Ni) 
   do i= 1, 4
      do j= 1, Number_Point_Gauss_Legendre_Element_Square
         Ni( j, i )= ( 1d0 +Xii( i )*Xig( j ) )*( 1d0 +Itai(i)*Itag( j ) )/4d0
      end do
   end do
   
   allocate( x4Ele( Number_Point_Gauss_Legendre_Element_Square, Number_Element_Square ) )
   allocate( y4Ele( Number_Point_Gauss_Legendre_Element_Square, Number_Element_Square ) )
    
   !$omp parallel do default( none )   &
   !$omp private(E, J) &
   !$omp shared( Number_Element_Square, x4Ele, y4Ele)
   do e= 1, Number_Element_Square
      do j= 1, Number_Point_Gauss_Legendre_Element_Square
         x4Ele( j, e )=0d0
         y4Ele( j, e )=0d0
      end do
   end do
    
   !IMPOSSIBLE TO BE PARALLELED
   do e= 1, Number_Element_Square
      do i= 1, 4
         do j= 1, Number_Point_Gauss_Legendre_Element_Square
            x4Ele( j, e )=x4Ele( j, e )+Ni( j, i )*Position_Node_Element_Square( 1, i, e )
            y4Ele( j, e )=y4Ele( j, e )+Ni( j, i )*Position_Node_Element_Square( 2, i, e )
         end do
      end do
   end do
   
   allocate( PNPXi( Number_Point_Gauss_Legendre_Element_Square, 4 ) )
   allocate( PNPeta( Number_Point_Gauss_Legendre_Element_Square, 4 ) )
   
   !$omp parallel do default( none )   &
   !$omp private(I, J) &
   !$omp shared( Itag, Itai, Xig, Xii, PNPXi, PNPeta) 
   do i= 1, 4
      do j= 1, Number_Point_Gauss_Legendre_Element_Square 
         PNPXi( j, i )=Xii(I)*(1d0+Itai(I)*Itag(J))/4d0
         PNPeta( j, i )=(1d0+Xii(I)*Xig(J))*Itai(I)/4d0
      end do
   end do
    
   deallocate( Xig )
   deallocate( Itag )
   deallocate( Xii )
   deallocate( Itai )
    
   allocate( PxPXi( Number_Point_Gauss_Legendre_Element_Square, Number_Element_Square ) )
   allocate( PxPeta( Number_Point_Gauss_Legendre_Element_Square, Number_Element_Square ) )
   allocate( PyPXi( Number_Point_Gauss_Legendre_Element_Square, Number_Element_Square ) )
   allocate( PyPeta( Number_Point_Gauss_Legendre_Element_Square, Number_Element_Square ) )
    
   !$omp parallel do default( none )   &
   !$omp private(E,I) &
   !$omp shared( Number_Element_Square, PxPXi, PxPeta, PyPXi, PyPeta) 
   do e= 1, Number_Element_Square
      do i= 1, Number_Point_Gauss_Legendre_Element_Square 
         PxPXi( i, e )=  0d0
         PxPeta( i, e )= 0d0
         PyPXi( i, e )=  0d0
         PyPeta( i, e )= 0d0
      end do
   end do
   
   !impossible parallel
   do e= 1, Number_Element_Square
      do i= 1, 4
         do j= 1, Number_Point_Gauss_Legendre_Element_Square 
            PxPXi( j, e )=PxPXi( j, e )+ PNPXi( j, i )*Position_Node_Element_Square( 1, i, e )
            PxPeta( j, e )=PxPeta( j, e )+ PNPeta( j, i )*Position_Node_Element_Square( 1, i, e )
            PyPXi( j, e )=PyPXi( j, e )+ PNPXi( j, i )*Position_Node_Element_Square( 2, i, e )
            PyPeta( j, e )=PyPeta( j, e )+ PNPeta( j, i )*Position_Node_Element_Square( 2, i, e )
         end do
      end do
   end do
   
   deallocate( Position_Node_Element_Square )
   
   !==============================================================================================================
   ! Jacobian in Square Element 
   !==============================================================================================================
   allocate( Jaco( Number_Point_Gauss_Legendre_Element_Square, Number_Element_Square ) )
    
   !$omp parallel do default( none )   &
   !$omp private( e, i ) &
   !$omp shared( Number_Element_Square, PxPXi, PxPeta, PyPXi, PyPeta, Jaco) 
   do e= 1, Number_Element_Square
      do i= 1, Number_Point_Gauss_Legendre_Element_Square
         Jaco( i, e )= PxPXi( i, e )*PyPeta( i, e ) -PyPXi( i, e )*PxPeta( i, e )
      end do
   end do
   
   !==============================================================================================================
   ! PNPx: \frac{\partial N}{\partial x} 
   !==============================================================================================================
   allocate( PNPx( Number_Point_Gauss_Legendre_Element_Square, 4, Number_Element_Square ) )
   allocate( PNPy( Number_Point_Gauss_Legendre_Element_Square, 4, Number_Element_Square ) )
   
   !$omp parallel do default( none )   &
   !$omp private(E, I, J) &
   !$omp shared( Number_Element_Square, Jaco, PNPeta, PxPXi, PNPXi, PxPeta, PyPXi, PyPeta, PNPx, PNPy)  
   do e= 1, Number_Element_Square
      do i= 1, 4
         do j= 1, Number_Point_Gauss_Legendre_Element_Square 
            PNPx( j, i, e )=( PyPeta( j, e )*PNPXi( j, i )- PyPXi( j, e )*PNPeta( j, i ) )/Jaco( j, e )
            PNPy( j, i, e )=( -PxPeta( j, e )*PNPXi( j, i )+PxPXi( j, e )*PNPeta( j, i ) )/Jaco( j, e )
         end do
      end do
   end do
   
   deallocate( PNPXi )
   deallocate( PNPeta )
   deallocate( PxPXi )
   deallocate( PxPeta )
   deallocate( PyPXi )
   deallocate( PyPeta )
  
   !==============================================================================================================
   write(*,*)'            Class_Element --> Dielectric_Constant '
   !==============================================================================================================

   allocate( Dielectric_Constant_Triangle( Number_Element_Triangle ) )

   !===============================
   if( Flag_Material==0 )then
   !===============================
      do i= 1, Number_Element_Triangle
         if( Class_Element_Triangle( i ) == ID_Element_Material )then
            Dielectric_Constant_Triangle( i )= Dielectric_Constant_Material

         else if( Class_Element_Triangle( i ) == ID_Element_OpenRegion )then
            Dielectric_Constant_Triangle( i )= Dielectric_Constant_Air

         else if( Class_Element_Triangle( i ) == ID_Element_FixedDomain )then
            Dielectric_Constant_Triangle( i )= Dielectric_Constant_FixedDomain

         else if( Class_Element_Triangle( i ) == ID_Element_Air )then
            Dielectric_Constant_Triangle( i )= Dielectric_Constant_Air

         else if( Class_Element_Triangle( i ) == ID_Element_Active )then
            Dielectric_Constant_Triangle( i )= Dielectric_Constant_Air

         else if( Class_Element_Triangle( i ) == ID_Element_Material_Exterior )then
            Dielectric_Constant_Triangle( i )= Dielectric_Constant_Material

         else if( Class_Element_Triangle( i ) == ID_Element_Air_Exterior )then
            Dielectric_Constant_Triangle( i )= Dielectric_Constant_Air

         else
            write(*,*)'ERROR'
            write(*,*)'  Class_Element_Triangle( i )=', Class_Element_Triangle( i )
            write(*,*)'  Analyze_Light_Scattering_TE.f90 152'
         end if
      end do 
   !===============================
   else if( Flag_Material==1 )then
   !===============================
      do i= 1, Number_Element_Triangle
         if( Class_Element_Triangle( i ) == ID_Element_Material )then
            Dielectric_Constant_Triangle( i ) &
            = 1.0d0 -Frequency_Plasma_Normalized**2 &
              /( Frequency_Normalized*( Frequency_Normalized +Imaginary_Unit*Coefficient_Dumping_Drude ) )

         else if( Class_Element_Triangle( i ) == ID_Element_OpenRegion )then
            Dielectric_Constant_Triangle( i )= Dielectric_Constant_Air

         else if( Class_Element_Triangle( i ) == ID_Element_FixedDomain )then
            Dielectric_Constant_Triangle( i )= Dielectric_Constant_FixedDomain

         else if( Class_Element_Triangle( i ) == ID_Element_Air )then
            Dielectric_Constant_Triangle( i )= Dielectric_Constant_Air

         else if( Class_Element_Triangle( i ) == ID_Element_Active )then
            Dielectric_Constant_Triangle( i )= Dielectric_Constant_Air

         else if( Class_Element_Triangle( i ) == ID_Element_Material_Exterior )then
            Dielectric_Constant_Triangle( i )= Dielectric_Constant_Material

         else if( Class_Element_Triangle( i ) == ID_Element_Air_Exterior )then
            Dielectric_Constant_Triangle( i )= Dielectric_Constant_Air

         else
            write(*,*)'ERROR'
            write(*,*)'  Class_Element_Triangle( i )=', Class_Element_Triangle( i )
            write(*,*)'  Analyze_Light_Scattering_TE.f90 174'
         end if
      end do
   !===============================
   else if( Flag_Material==2 )then
   !===============================
      do i= 1, Number_Element_Triangle
         if( Class_Element_Triangle( i ) == ID_Element_Material )then
            Dielectric_Constant_Triangle( i )= Dielectric_Constant_Material

         else if( Class_Element_Triangle( i ) == ID_Element_OpenRegion )then
            Dielectric_Constant_Triangle( i )= Dielectric_Constant_Air

         else if( Class_Element_Triangle( i ) == ID_Element_FixedDomain )then
            Dielectric_Constant_Triangle( i ) &
            = 1.0d0 -Frequency_Plasma_Normalized**2 &
              /( Frequency_Normalized*( Frequency_Normalized +Imaginary_Unit*Coefficient_Dumping_Drude ) )

         else if( Class_Element_Triangle( i ) == ID_Element_Air )then
            Dielectric_Constant_Triangle( i )= Dielectric_Constant_Air

         else if( Class_Element_Triangle( i ) == ID_Element_Active )then
            Dielectric_Constant_Triangle( i )= Dielectric_Constant_Air

         else if( Class_Element_Triangle( i ) == ID_Element_Material_Exterior )then
            Dielectric_Constant_Triangle( i )= Dielectric_Constant_Material

         else if( Class_Element_Triangle( i ) == ID_Element_Air_Exterior )then
            Dielectric_Constant_Triangle( i )= Dielectric_Constant_Air

         else
            write(*,*)'ERROR'
            write(*,*)'  Class_Element_Triangle( i )=', Class_Element_Triangle( i )
            write(*,*)'  Analyze_Light_Scattering_TE.f90 174'
         end if
      end do
   !===============================
   else
   !===============================
      call Output_Error( 'Analyze_Light_Scattering_TE', 598 )
   end if

   !==============================================================================================================
   ! Absorbing Function of PML 
   !==============================================================================================================
   allocate( GammaX( Number_Point_Gauss_Legendre_Element_Square, Number_Element_Square ) )
   allocate( GammaY( Number_Point_Gauss_Legendre_Element_Square, Number_Element_Square ) )   
    
   !$omp parallel do default( none )   &
   !$omp private(E, J) &
   !$omp shared( Number_Element_Square, y4Ele, x4Ele, Max_Position_Node, Frequency_Normalized, GammaX, GammaY ) & 
   !$omp shared( Class_Element_Square, ID_Element_PML_Xp, ID_Element_PML_Xm, ID_Element_PML_Yp, ID_Element_PML_Ym )&
   !$omp shared( ID_Element_PML_XpYp, ID_Element_PML_XpYm, ID_Element_PML_XmYp, ID_Element_PML_XmYm ) 
   do e= 1, Number_Element_Square
      do j= 1, Number_Point_Gauss_Legendre_Element_Square 
         if( Class_Element_Square(E)==ID_Element_PML_Xp .or. Class_Element_Square(E)==ID_Element_PML_Xm )THEN
            GammaX( j, e ) &
            = 1d0 +Imaginary_Unit/( 2d0*Pi*Frequency_Normalized*(Max_Position_Node( 1, 2 )-abs(x4Ele( j, e ))) )
            GammaY( j, e )=One
         else if( Class_Element_Square(E)==ID_Element_PML_Yp .or. Class_Element_Square(E)==ID_Element_PML_Ym )THEN
            GammaX( j, e )=One
            GammaY( j, e ) &
            = 1d0 +Imaginary_Unit/( 2d0*Pi*Frequency_Normalized*(Max_Position_Node( 2, 2 )-abs(y4Ele( j, e ))) )
         else if( Class_Element_Square(E)==ID_Element_PML_XpYp .or. &
                Class_Element_Square(E)==ID_Element_PML_XpYm .or. &
                Class_Element_Square(E)==ID_Element_PML_XmYp .or. &
                Class_Element_Square(E)==ID_Element_PML_XmYm )then
            GammaX( j, e ) &
            = 1d0 +Imaginary_Unit/( 2d0*Pi*Frequency_Normalized*(Max_Position_Node( 1, 2 )-abs(x4Ele( j, e ))) )
            GammaY( j, e ) &
            = 1d0 +Imaginary_Unit/( 2d0*Pi*Frequency_Normalized*(Max_Position_Node( 2, 2 )-abs(y4Ele( j, e ))) )
         end if
      end do
   end do
   
   allocate( LocalMatrix_Triangle_Scattering_tmp( 3, 3, Number_Element_Triangle ) )
   allocate( LocalMatrix_Triangle_Scattering( 3, 3, Number_Element_Triangle ) )
   allocate( LocalMatrix_Triangle_Scattering_RHS( 3, 3, Number_Element_Triangle ) )

   do e= 1, Number_Element_Triangle
      do i= 1, 3
         do j= 1, 3
            LocalMatrix_Triangle_Scattering_tmp( j, i, e )= Zero 
            LocalMatrix_Triangle_Scattering( j, i, e )= Zero 
            LocalMatrix_Triangle_Scattering_RHS( j, i, e )= Zero 
         end do
      end do
   end do
 
   allocate( Flag_Element_Triangle_Removed( Number_Element_Triangle ) ) 
   allocate( Flag_Element_Square_Removed( Number_Element_Square ) ) 

   !==============================================================================================================
   if( Flag_PEC_Boundary_Condition==Flag_On )then
   !==============================================================================================================

      !==============================================================================================================
      ! Local Matrix in Triangle Element 
      !==============================================================================================================
     
      !$omp parallel do default( none )   &
      !$omp private( e, i, j ) &
      !$omp shared( Number_Element_Triangle, LocalMatrix_Triangle_Scattering_tmp ) & 
      !$omp shared( Frequency_Normalized, Area_Element_Triangle ) &
      !$omp shared( Coefficient_C, Coefficient_Xi, Coefficient_Eta ) &
      !$omp shared( Coefficient_Xi_Xi, Coefficient_Eta_Eta, Coefficient_Xi_Eta ) &
      !$omp shared( Class_Element_Triangle, ID_Element_PEC )
      do e= 1, Number_Element_Triangle
         if( Class_Element_Triangle( e ) /= ID_Element_PEC )then
            do i= 1, 3
               do j= 1, 3
                  LocalMatrix_Triangle_Scattering_tmp( j, i, e ) & 
                  = 1d0/2d0*Coefficient_C( j, i, e ) +1d0/6d0*Coefficient_Xi( j, i, e ) & 
                   +1d0/6d0*Coefficient_Eta( j, i, e ) +1d0/12d0*Coefficient_Xi_Xi( j, i, e ) & 
                   +1d0/12d0*Coefficient_Eta_Eta( j, i, e ) +1d0/24d0*Coefficient_Xi_Eta( j, i, e ) 
               end do
            end do
         end if
      end do
 
      !$omp parallel do default( none )   &
      !$omp private( e, i, j ) &
      !$omp shared( Number_Element_Triangle, LocalMatrix_Triangle_Scattering ) & 
      !$omp shared( Area_Element_Triangle, Dielectric_Constant_Triangle, LocalMatrix_Triangle_Scattering_tmp ) &
      !$omp shared( BasisFunction_B, BasisFunction_C, Frequency_Normalized ) &
      !$omp shared( Class_Element_Triangle, ID_Element_PEC )
      do e= 1, Number_Element_Triangle
         if( Class_Element_Triangle( e ) /= ID_Element_PEC )then
            do i= 1, 3
               do j= 1, 3
                  LocalMatrix_Triangle_Scattering( j, i, e ) & 
                  = -1d0/Dielectric_Constant_Triangle( e ) &
                     *Area_Element_Triangle( e ) & 
                     *( BasisFunction_B( i, e )*BasisFunction_B( j, e ) &
                     +BasisFunction_C( i, e )*BasisFunction_C( j, e ) ) &
                    +4d0*Pi*Pi*Frequency_Normalized*Frequency_Normalized &
                     *2d0*Area_Element_Triangle( e ) & 
                     *LocalMatrix_Triangle_Scattering_tmp( j, i, e )
               end do
            end do
         end if
      end do
 
      !!$omp parallel do default( none )   &
      !!$omp private( e, i, j ) &
      !!$omp shared( Number_Element_Triangle, LocalMatrix_Triangle_Scattering_RHS ) & 
      !!$omp shared( Area_Element_Triangle, Dielectric_Constant_Triangle, LocalMatrix_Triangle_Scattering_tmp ) &
      !!$omp shared( BasisFunction_B, BasisFunction_C, Frequency_Normalized ) & 
      !!$omp shared( Class_Element_Triangle )
      do e= 1, Number_Element_Triangle
         if( Class_Element_Triangle( e ) /= ID_Element_PEC )then
            do i= 1, 3
               do j= 1, 3
                  LocalMatrix_Triangle_Scattering_RHS( j, i, e ) & 
                  = ( 1d0/Dielectric_Constant_Triangle( e ) -1d0/Dielectric_Constant_OpenRegion ) & 
                    *Area_Element_Triangle( e ) & 
                    *( BasisFunction_B( i, e )*BasisFunction_B( j, e ) &
                      +BasisFunction_C( i, e )*BasisFunction_C( j, e ) ) 
               end do
            end do
         end if
      end do

      !==============================================================================================================
      write(*,*)'      =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-='
      write(*,*)'      Implement PEC Boundary Condition'
      write(*,*)'      =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-='
      !==============================================================================================================

      allocate( Flag_Element_Triangle_PEC_BC( Number_Element_Triangle ) )
      allocate( Position_Node_PEC_BC( 2, 2, Number_Element_Triangle ) )

      call Detect_Node_and_Element_Triangle_on_PEC_BC_TE&
         ( Number_Node, Number_Element_Triangle, &
           Index_Element_2_Node_Triangle, &
           ID_Element_PEC, Class_Element_Triangle, &
           Position_Node, &
         !================================================================
           Flag_Element_Triangle_PEC_BC, Number_Outer_Element_on_PEC_BC, &
           Position_Node_PEC_BC )

      allocate( Length_Edge_on_PEC_BC_in_Element( Number_Element_Triangle ) )

      do e= 1, Number_Element_Triangle
         if( Flag_Element_Triangle_PEC_BC( e )==Flag_On )then 
            Length_Edge_on_PEC_BC_in_Element( e ) &
            = sqrt( ( Position_Node_PEC_BC( 1, 2, e ) -Position_Node_PEC_BC( 1, 1, e ) )**2 &
                 +( Position_Node_PEC_BC( 2, 2, e ) -Position_Node_PEC_BC( 2, 1, e ) )**2 )
         else if( Flag_Element_Triangle_PEC_BC( e )==Flag_Off )then 
            Length_Edge_on_PEC_BC_in_Element( e )= -1d10
         end if
      end do

      allocate( Position_Node_PEC_BC_tmp( 2, 2, Number_Outer_Element_on_PEC_BC ) )

      Counter= 0
      do e= 1, Number_Element_Triangle
         if( Flag_Element_Triangle_PEC_BC( e )==Flag_On )then 
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
         call Output_Error( 'Analyze_Light_Scattering_TE', 868 )
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
         if( Flag_Element_Triangle_PEC_BC( e )==Flag_On )then 
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
         call Output_Error( 'Analyze_Light_Scattering_TE', 898 )
      end if

      deallocate( Unit_Normal_Vector_on_PEC_BC_tmp )

      do e= 1, Number_Element_Triangle
         if( Flag_Element_Triangle_PEC_BC( e )==Flag_On )then 
            do i= 1, 3
               do j= 1, 3
                  LocalMatrix_Triangle_Scattering_RHS( j, i, e ) & 
                  = LocalMatrix_Triangle_Scattering_RHS( j, i, e ) & 
                   +1d0/Dielectric_Constant_OpenRegion & 
                   *Length_Edge_on_PEC_BC_in_Element( e ) & 
                   *( BasisFunction_B( j, e ) *Unit_Normal_Vector_on_PEC_BC( 1, e ) &
                      +BasisFunction_C( j, e ) *Unit_Normal_Vector_on_PEC_BC( 2, e ) ) &
                   *( BasisFunction_A( i, e ) &
                     +BasisFunction_B( i, e )/2d0 &
                      *( Position_Node_PEC_BC( 1, 1, e ) +Position_Node_PEC_BC( 1, 2, e ) ) & 
                     +BasisFunction_C( i, e )/2d0 & 
                      *( Position_Node_PEC_BC( 2, 1, e ) +Position_Node_PEC_BC( 2, 2, e ) ) )  
               end do
            end do
         end if
      end do

      deallocate( Length_Edge_on_PEC_BC_in_Element )
      deallocate( Unit_Normal_Vector_on_PEC_BC )
      deallocate( Position_Node_PEC_BC )
      deallocate( Flag_Element_Triangle_PEC_BC )

      !==============================================================================================================
      if( Flag_Thermal_Device==10 )then
         write(*,*)'      =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-='
         write(*,*)'      Implement PEC Boundary Condition for Carpet cloak'
         write(*,*)'      =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-='
      !==============================================================================================================

         allocate( Element_Number_on_Flat_PEC( Number_Grid_X_Scattering*2 ) )
         allocate( Node_Number_on_Flat_PEC_in_Element( 2, Number_Grid_X_Scattering*2 ) )
   
         call Detect_Node_and_Element_Triangle_on_Flat_PEC_TE&
            ( Number_Node, Number_Element_Triangle, Number_Element_Square, &
              Index_Element_2_Node_Triangle, Index_Element_2_Node_Square, &
              ID_Element_PEC, Class_Element_Triangle, &
              ID_Element_Flat_PEC, Class_Element_Square, &
              Number_Grid_X_Scattering*2, &
            !================================================================
              Element_Number_on_Flat_PEC, Number_Element_on_Flat_PEC, &
              Node_Number_on_Flat_PEC_in_Element )

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
                  LocalMatrix_Triangle_Scattering_RHS( j, i, Element_Number_on_Flat_PEC( e ) ) & 
                  = LocalMatrix_Triangle_Scattering_RHS( j, i, Element_Number_on_Flat_PEC( e ) ) & 
                   +1d0/Dielectric_Constant_OpenRegion & 
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
      ! Create Flag_Element_Square_Removed 
      !==============================================================================================================
      
      Flag_Element_Triangle_Removed= Flag_Off
      Flag_Element_Square_Removed= Flag_Off
   
      do e= 1, Number_Element_Triangle
         if( Class_Element_Triangle( e ) == ID_Element_FixedDomain )then
            Flag_Element_Triangle_Removed( e )= Flag_On
         end if
      end do
   
      if( Flag_Thermal_Device==10 )then
         do e= 1, Number_Element_Square 
            !if( Class_Element_Square( e ) == ID_Element_PML_Ym )then
            if( Class_Element_Square( e ) == ID_Element_PML_Ym .or. &
                Class_Element_Square( e ) == ID_Element_PML_XmYm  .or. & 
                Class_Element_Square( e ) == ID_Element_PML_XpYm )then

               Flag_Element_Square_Removed( e )= Flag_On
            end if
         end do
      end if

   !==============================================================================================================
   else if( Flag_PEC_Boundary_Condition==Flag_Off )then
   !==============================================================================================================
 
      !==============================================================================================================
      ! Local Matrix in Triangle Element 
      !==============================================================================================================
     
      !$omp parallel do default( none )   &
      !$omp private( e, i, j ) &
      !$omp shared( Number_Element_Triangle, LocalMatrix_Triangle_Scattering_tmp ) & 
      !$omp shared( Frequency_Normalized, Area_Element_Triangle ) &
      !$omp shared( Coefficient_C, Coefficient_Xi, Coefficient_Eta ) &
      !$omp shared( Coefficient_Xi_Xi, Coefficient_Eta_Eta, Coefficient_Xi_Eta )
      do e= 1, Number_Element_Triangle
         do i= 1, 3
            do j= 1, 3
               LocalMatrix_Triangle_Scattering_tmp( j, i, e ) & 
               = 1d0/2d0*Coefficient_C( j, i, e ) +1d0/6d0*Coefficient_Xi( j, i, e ) & 
                +1d0/6d0*Coefficient_Eta( j, i, e ) +1d0/12d0*Coefficient_Xi_Xi( j, i, e ) & 
                +1d0/12d0*Coefficient_Eta_Eta( j, i, e ) +1d0/24d0*Coefficient_Xi_Eta( j, i, e ) 
            end do
         end do
      end do
     
      !$omp parallel do default( none )   &
      !$omp private( e, i, j ) &
      !$omp shared( Number_Element_Triangle, LocalMatrix_Triangle_Scattering ) & 
      !$omp shared( Area_Element_Triangle, Dielectric_Constant_Triangle, LocalMatrix_Triangle_Scattering_tmp ) &
      !$omp shared( BasisFunction_B, BasisFunction_C, Frequency_Normalized ) 
      do e= 1, Number_Element_Triangle
         do i= 1, 3
            do j= 1, 3
               LocalMatrix_Triangle_Scattering( j, i, e ) & 
               = -Area_Element_Triangle( e ) & 
                  *( BasisFunction_B( i, e )*BasisFunction_B( j, e )+ BasisFunction_C( i, e )*BasisFunction_C( j, e ) ) &
                 +Dielectric_Constant_Triangle( e ) & 
                  *4d0*Pi*Pi*Frequency_Normalized*Frequency_Normalized*2d0*Area_Element_Triangle( e ) & 
                  *LocalMatrix_Triangle_Scattering_tmp( j, i, e )
            end do
         end do
      end do
    
      !$omp parallel do default( none )   &
      !$omp private( e, i, j ) &
      !$omp shared( Number_Element_Triangle, LocalMatrix_Triangle_Scattering_RHS ) & 
      !$omp shared( Area_Element_Triangle, Dielectric_Constant_Triangle, LocalMatrix_Triangle_Scattering_tmp ) &
      !$omp shared( BasisFunction_B, BasisFunction_C, Frequency_Normalized, Dielectric_Constant_OpenRegion ) 
      do e= 1, Number_Element_Triangle
         do i= 1, 3
            do j= 1, 3
               LocalMatrix_Triangle_Scattering_RHS( j, i, e ) & 
               = -( Dielectric_Constant_Triangle( e ) -Dielectric_Constant_OpenRegion ) & 
                  *4d0*Pi*Pi*Frequency_Normalized*Frequency_Normalized*2d0*Area_Element_Triangle( e ) & 
                  *LocalMatrix_Triangle_Scattering_tmp( j, i, e )
            end do
         end do
      end do

      Flag_Element_Triangle_Removed= Flag_Off
      Flag_Element_Square_Removed= Flag_Off
   !==============================================================================================================
   else ! Flag_PEC_Boundary_Condition /= 0, 1
   !==============================================================================================================
      call Output_Error( 'Analyze_Light_Scattering_TE', 838 )
   end if

   deallocate( Area_Element_Triangle )

   !==============================================================================================================
   ! Relative Permittivity in PML ( Square Element )
   !==============================================================================================================

   allocate( Dielectric_Constant_PML( Number_Element_Square ) )

   do e= 1, Number_Element_Square
      if( Class_Element_PML_Material( e )==ID_Element_PML_Air )then
         Dielectric_Constant_PML( e )= Dielectric_Constant_PML_Air

      else if( Class_Element_PML_Material( e )==ID_Element_PML_Material )then
         Dielectric_Constant_PML( e )= Dielectric_Constant_PML_Material

      else
         call Output_Error( 'Analyze_Light_Scattering_TE', 751 )
      end if
   end do
 
   !==============================================================================================================
   ! Local Matrix in Square Element 
   !==============================================================================================================
   allocate( LocalMatrix_Square_PML( 4, 4, Number_Element_Square) )
   allocate( LocalMatrix_Square_PML_tmp(4, 4, 4, Number_Element_Square) )
    
   !$omp parallel do default( none )   &
   !$omp private( e, i, j ) &
   !$omp shared( Number_Element_Square, LocalMatrix_Square_PML) 
   do e= 1, Number_Element_Square
      do i= 1, 4
         do j= 1, 4
            LocalMatrix_Square_PML( j, i, e )= ZERO
         end do
      end do
   end do
  
   !$omp parallel do default( none )   &
   !$omp private( e, i, j ) &
   !$omp shared( Number_Element_Square, LocalMatrix_Square_PML_tmp, Wg ) &
   !$omp shared( Dielectric_Constant_PML, GammaY, GammaX, Jaco, PNPy, PNPx, Frequency_Normalized, Ni )
   do e= 1, Number_Element_Square
      do i= 1, 4
         do j= 1, 4
            do k= 1, Number_Point_Gauss_Legendre_Element_Square 
               LocalMatrix_Square_PML_tmp( k, j, i, e ) &
               = -Wg( k )*GammaY( k, e )/GammaX( k, e ) &
                  *Jaco( k, e )*PNPx( k, i, e )*PNPx( k, j, e ) &
                 -Wg( k )*GammaX( k, e )/GammaY( k, e ) & 
                  *Jaco( k, e )*PNPy( k, i, e )*PNPy( k, j, e ) &
                 +Wg( k )*( 4d0*Pi*Pi*Frequency_Normalized*Frequency_Normalized ) & 
                  *Dielectric_Constant_PML( e )*GammaX( k, e )*GammaY( k, e )*Jaco( k, e )*Ni( k, i )*Ni( k, j )
            end do
         end do
      end do
   end do
 
   deallocate( Dielectric_Constant_PML )
 
   !$omp parallel do default( none )   &
   !$omp private( e, i, j, k ) &
   !$omp shared( Number_Element_Square, LocalMatrix_Square_PML, LocalMatrix_Square_PML_tmp )
   do e= 1, Number_Element_Square
      do i= 1, 4
         do j= 1, 4
            do k= 1, Number_Point_Gauss_Legendre_Element_Square 
               LocalMatrix_Square_PML( j, i, e ) &
               = LocalMatrix_Square_PML( j, i, e ) +LocalMatrix_Square_PML_tmp( k, j, i, e ) 
            end do
         end do
      end do
   end do
  
   deallocate( LocalMatrix_Square_PML_tmp )
   deallocate( GammaX )
   deallocate( GammaY )
   
   !==============================================================================================================
   ! Local Matrix --> Global Matrix 
   !==============================================================================================================

   allocate( GlobalMatrix( Number_Node, Width_Matrix_LHS ) ) 
   allocate( J_GlobalMatrix( Number_Node, Width_Matrix_LHS ) )
 
   call Format_Local_Global_ZMUMPS&
      ( LocalMatrix_Triangle_Scattering, Index_Element_2_Node_Triangle, &
        3, Number_Node, Number_Element_Triangle, &
        1, Width_Matrix_LHS, Flag_Symmetric, &
        Class_Element_Triangle, Flag_Element_Triangle_Removed, Flag_Remove_Element, &
        !=====================================================
        GlobalMatrix, J_GlobalMatrix, Number_NonZero )

   call Format_Local_Global_ZMUMPS&
      ( LocalMatrix_Square_PML, Index_Element_2_Node_Square, &
        4, Number_Node, Number_Element_Square, &
        0, Width_Matrix_LHS, Flag_Symmetric, &
        Class_Element_Square, Flag_Element_Square_Removed, Flag_Remove_Element, &
        !=====================================================
        GlobalMatrix, J_GlobalMatrix, Number_NonZero )

   deallocate( LocalMatrix_Triangle_Scattering )
   deallocate( LocalMatrix_Square_PML )

   allocate( GlobalMatrix_RHS( Number_Node, Width_Matrix_LHS ) )
   allocate( J_GlobalMatrix_RHS( Number_Node, Width_Matrix_LHS ) )

   call Format_Local_Global_ZMUMPS&
      ( LocalMatrix_Triangle_Scattering_RHS, Index_Element_2_Node_Triangle, &
        3, Number_Node, Number_Element_Triangle, &
        1, Width_Matrix_LHS, Flag_Symmetric, &
        Class_Element_Triangle, Flag_Element_Triangle_Removed, Flag_Remove_Element, &
        !=====================================================
        GlobalMatrix_RHS, J_GlobalMatrix_RHS, Number_NonZero_RHS )
 
   deallocate( Flag_Element_Triangle_Removed ) 
   deallocate( Flag_Element_Square_Removed ) 
   deallocate( LocalMatrix_Triangle_Scattering_RHS ) 
   deallocate( LocalMatrix_Triangle_Scattering_tmp )

   !==============================================================================================================
   ! Incident Wave 
   !==============================================================================================================

   do i= 1, Number_Node
      Field_Magnetic_Incident( i )= Zero
   end do

   call Compute_Incident_Wave &
      ( Flag_Incident_Wave, Position_Node, Class_Node, Number_Node, & 
        Frequency_Normalized, Incident_Angle, &
        ID_Node_Material, ID_Node_Air, ID_Node_OpenRegion, ID_Node_Active, ID_Node_FixedDomain, &
        ID_Node_PML_X, ID_Node_PML_Y, ID_Node_PML_XY, &
        real( Dielectric_Constant_OpenRegion ), &
        Position_Dipole_x, Position_Dipole_Y, &
        !==================================================================
        Field_Magnetic_Incident )

   if( Flag_PEC_Boundary_Condition==Flag_On )then
      allocate( Flag_Node_PEC_BC_Fixed( Number_Node ) )
 
      Flag_Node_PEC_BC_Fixed= Flag_On
   
      do e= 1, Number_Element_Triangle
         if( Class_Element_Triangle( e ) /= ID_Element_PEC )then
            do i= 1, 3
               Flag_Node_PEC_BC_Fixed( Index_Element_2_Node_Triangle( i, e ) )= Flag_Off
            end do
         end if
      end do
      do i= 1, Number_Node
         if( Flag_Node_PEC_BC_Fixed( i )==Flag_On )then
            Field_Magnetic_Incident( i )= Zero
         end if
      end do
      deallocate( Flag_Node_PEC_BC_Fixed )
   end if
 
   allocate( Global_Vector_RHS( Number_Node ) )

   !$omp parallel do default( none )   &
   !$omp private( i ) &
   !$omp shared( Number_Node, Global_Vector_RHS )
   do i= 1, Number_Node
      Global_Vector_RHS( i )= Zero
   end do

   !!$omp parallel do default( none )   &
   !!$omp private( i, j ) &
   !!$omp shared( Number_Node, Width_Matrix_LHS, Global_Vector_RHS ) &
   !!$omp shared( GlobalMatrix_RHS, J_GlobalMatrix_RHS, Field_Magnetic_Incident )
   do i= 1, Number_Node
      do j= 1, Width_Matrix_LHS
         if( J_GlobalMatrix_RHS( i, j )/=0 )then
            Global_Vector_RHS( i )&
            = Global_Vector_RHS( i ) &
             +GlobalMatrix_RHS( i, j ) *Field_Magnetic_Incident( J_GlobalMatrix_RHS( i, j ) )
         end if
      end do
   end do

   deallocate( GlobalMatrix_RHS )
   deallocate( J_GlobalMatrix_RHS )

   !==================================================================================================
   write(*,*)'         '
   write(*,*)'         =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-='
   write(*,*)'         Dirichlet Boundary Condition of PML'
   write(*,*)'         =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-='
   !==================================================================================================

   allocate( NodeNumber_Dirichlet_BC( 2 ) )

   NodeNumber_Dirichlet_BC( 1 )= Number_Node -Number_Node_Dirichlet_Boundary_PML +1 
   NodeNumber_Dirichlet_BC( 2 )= Number_Node

   allocate( Nodal_Value_on_Dirichlet_BC( Number_Node ) )
   allocate( Flag_Node_Dirichlet_BC( Number_Node ) )

   do i= 1, Number_Node 
      Nodal_Value_on_Dirichlet_BC( i )= Zero
      Flag_Node_Dirichlet_BC( i )= Flag_Off
   end do

   do i= NodeNumber_Dirichlet_BC( 1 ), NodeNumber_Dirichlet_BC( 2 ) 
      Flag_Node_Dirichlet_BC( i )= Flag_On
   end do

   deallocate( NodeNumber_Dirichlet_BC )

   call Implement_Dirichlet_BC_complex_CSRF &
      ( Nodal_Value_on_Dirichlet_BC, Flag_Node_Dirichlet_BC, &
        Number_Node, Width_Matrix_LHS, &
        !======================================================
        GlobalMatrix, J_GlobalMatrix, Global_Vector_RHS )

   deallocate( Flag_Node_Dirichlet_BC )
   deallocate( Nodal_Value_on_Dirichlet_BC )

   !==================================================================================================
   write(*,*)'         '
   write(*,*)'         =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-='
   write(*,*)'         Fixed Nodal Vlaue in PEC'
   write(*,*)'         =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-='
   !==================================================================================================

   allocate( Flag_Node_PEC_BC_Fixed( Number_Node ) )

   do i= 1, Number_Node
      Flag_Node_PEC_BC_Fixed( i )= Flag_Off
   end do

   do e= 1, Number_Element_Triangle
      if( Class_Element_Triangle( e ) == ID_Element_PEC )then
         do i= 1, 3
            Flag_Node_PEC_BC_Fixed( Index_Element_2_Node_Triangle( i, e ) )= Flag_On
         end do
      end if
   end do

   do e= 1, Number_Element_Triangle
      if( Class_Element_Triangle( e ) /= ID_Element_PEC )then
         do i= 1, 3
            Flag_Node_PEC_BC_Fixed( Index_Element_2_Node_Triangle( i, e ) )= Flag_Off
         end do
      end if
   end do

   if( Flag_Thermal_Device==10 )then
      do e= 1, Number_Element_Square 
         if( Class_Element_Square( e ) == ID_Element_PML_Ym .or. &
             Class_Element_Square( e ) == ID_Element_PML_XmYm  .or. & 
             Class_Element_Square( e ) == ID_Element_PML_XpYm )then

            do i= 1, 4
               Flag_Node_PEC_BC_Fixed( Index_Element_2_Node_Square( i, e ) )= Flag_On
            end do
         end if
      end do
      do e= 1, Number_Element_Triangle
         if( Class_Element_Triangle( e ) /= ID_Element_PEC )then
            do i= 1, 3
               Flag_Node_PEC_BC_Fixed( Index_Element_2_Node_Triangle( i, e ) )= Flag_Off
            end do
         end if
      end do
   end if

   allocate( Nodal_Value_Fixed_in_PEC( Number_Node ) )
 
   do i= 1, Number_Node 
      Nodal_Value_Fixed_in_PEC( i )= -Field_Magnetic_Incident( i )
   end do

   call Implement_Dirichlet_BC_complex_CSRF &
      ( Nodal_Value_Fixed_in_PEC, Flag_Node_PEC_BC_Fixed, &
        Number_Node, Width_Matrix_LHS, &
       !========================================================
        GlobalMatrix, J_GlobalMatrix, Global_Vector_RHS )

   deallocate( Nodal_Value_Fixed_in_PEC )


   deallocate( Flag_Node_PEC_BC_Fixed )

   !==================================================================================================
   write(*,*)'         '
   write(*,*)'         =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-='
   write(*,*)'         Count Number Non-Zero '
   write(*,*)'         =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-='
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
            call Output_Error( 'Analyze_Light_Scattering_TE', 1118 )
         end if
      end do
   end do

   write(*,*)'         Number_NonZero=', Number_NonZero
   write(*,*)'         Number_Node=',  Number_Node 
   write(*,*)'         Width_Matrix_LHS=', Width_Matrix_LHS

   !==================================================================================================
   write(*,*)'         Global Matrix --> MUMPS Format '
   !==================================================================================================

   allocate( aa_zmumps( Number_NonZero ) )
   allocate( ia_zmumps( Number_NonZero ) )
   allocate( ja_zmumps( Number_NonZero ) )

   call Format_Global_aa_ZMUMPS&
      ( GlobalMatrix, J_GlobalMatrix, & 
        Number_Node, Width_Matrix_LHS, Number_NonZero, Flag_Symmetric, Flag_Matrix_Upper, & 
        !==================================================== 
        aa_zmumps, ia_zmumps, ja_zmumps ) 
 
   write(*,*)'Flag_Symmetric=', Flag_Symmetric, 'at Main.f90' 
 
   if( Flag_PEC_Boundary_Condition==Flag_On )then 
 
      call Check_ZMUMPS_Format & 
         ( aa_zmumps, ia_zmumps, ja_zmumps, &  
           Number_Node, Number_NonZero, Flag_Symmetric, Flag_Matrix_Upper ) 
 
   else if( Flag_PEC_Boundary_Condition==Flag_Off )then 
 
      call Check_ZMUMPS_Format & 
         ( aa_zmumps, ia_zmumps, ja_zmumps, &  
           Number_Node, Number_NonZero, Flag_Symmetric, Flag_Matrix_Upper ) 
 
   end if 

   deallocate( GlobalMatrix )
   deallocate( J_GlobalMatrix )


   !==================================================================================================
   write(*,*)'         Solve Linear Equation'
   !==================================================================================================

   call ME57( Number_Node, Number_NonZero, aa_zmumps, ia_zmumps, ja_zmumps, Global_Vector_RHS, Field_Magnetic_Scattering )

   deallocate( aa_zmumps )
   deallocate( ia_zmumps )
   deallocate( ja_zmumps )

   !$omp parallel do default( none )   &
   !$omp private( i ) &
   !$omp shared( Number_Node, Field_Magnetic_All, Field_Magnetic_Scattering, Field_Magnetic_Incident )
   do i= 1, Number_Node 
      Field_Magnetic_All( i )= Field_Magnetic_Scattering( i ) +Field_Magnetic_Incident( i )
   end do

   deallocate( Global_Vector_RHS )

   deallocate( Dielectric_Constant_Triangle )
   
   deallocate( Position_Node_Element_Triangle )
    
   deallocate( x4Ele ) 
   deallocate( y4Ele )
   
   deallocate( PNPx )
   deallocate( PNPy )
   
   deallocate( Jaco )
   deallocate( Wg )
   deallocate( Ni )
   
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
   
   return
end subroutine Analyze_Light_Scattering_TE



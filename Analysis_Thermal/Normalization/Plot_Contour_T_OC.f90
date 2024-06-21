
subroutine Plot_Contour_T_OC  &
       ( File_Number, &
         Flag_Thermal_Insulation_BC, Value_Plotted, Number_Level, Optimization_Step,  &
         Position_Node, Number_Node, &
         Index_Element_2_Node, Number_Element, &
         Class_Element, ID_Element_Structure, ID_Element_FixedDomain, &
         Number_PV, Position_Arrow_PV, Length_Maximum_Vector, Level_PV, &
         Objective_Function ) 

   !use omp_lib
   use Parameters
   implicit none

   integer, intent(in) :: File_Number, Flag_Thermal_Insulation_BC
   integer, intent(in) :: Number_Node, Number_Element
   double precision, intent(inout) :: Value_Plotted( Number_Node )

   integer, intent(in) :: Number_Level
   integer, intent(in) :: Optimization_Step

   double precision, intent(in) :: Position_Node( 2, Number_Node )
   integer, intent(in) :: Index_Element_2_Node( 4, Number_Element )
   integer, intent(in) :: Class_Element( Number_Element )
   integer, intent(in) :: ID_Element_Structure, ID_Element_FixedDomain
   integer, intent(in) :: Number_PV 
   double precision, intent(in) :: Position_Arrow_PV( 2, 2, Number_PV ) ! ( xy, from to, Number_PV )
   !double precision, intent(inout) :: Length_Maximum_Vector 
   double precision, intent(in) :: Length_Maximum_Vector 
   
   double precision, intent(in) :: Level_PV( Number_PV )
   double precision, intent(in) :: Objective_Function( 0:Number_Optimization_Step ) 
 
   character(len=Length_Character_Optimization_Step) :: Optimization_Step_Character
   character(len=8) :: Format_Filenumber
   
   integer :: e, i, j, k, l 
   
   character(len=256) :: Filename_Result

   double precision :: Length_Maximum_Vector_Position 

   integer, allocatable, dimension(:,:) :: Index_Element_2_Node_tmp
   double precision, allocatable, dimension(:) :: Edge_Length 
   integer :: Counter_Edge 
   
   double precision, allocatable, dimension(:) :: Position_Minimum 
   double precision, allocatable, dimension(:) :: Position_Maximum 
   double precision, allocatable, dimension(:) :: Translation_Position, Width_PostScript 
   double precision, allocatable, dimension(:,:,:) :: Position_Maximum_Color_Var 
   double precision, allocatable, dimension(:) :: Translate_Color_Var 
   
   double precision, allocatable, dimension(:,:) :: Position_Node_Plot_Mesh 

   double precision, allocatable, dimension(:,:) :: Difference_Position_1_Minus_NewPoint
   double precision, allocatable, dimension(:,:) :: Difference_Position_2_Minus_NewPoint 
   
   double precision :: Minimum_Value_Plotted, Maximum_Value_Plotted
   integer, allocatable, dimension(:) :: Level_Node, Counter_Node_Level
   double precision, allocatable, dimension(:,:) :: Value_Boundary_Level 
   double precision, allocatable, dimension(:,:) :: RGB
   double precision, allocatable, dimension(:) :: RGB_Lowest, RGB_Highest
   integer, allocatable, dimension(:) :: Level_Boundary 

   integer, allocatable, dimension(:) :: Number_Node_Difference_Level
   integer, allocatable, dimension(:) :: Level_Low, Level_High
   integer, allocatable, dimension(:,:) :: Local_Node_Number_Edge_Line
   integer, allocatable, dimension(:,:,:) :: Local_Node_Number_Edge_Line_Plot
   double precision, allocatable, dimension(:,:) :: Value_Boundary_Level_Line
   double precision, allocatable, dimension(:,:,:) :: Position_Node_Boundary_Level_Plot_1
   double precision, allocatable, dimension(:,:,:) :: Position_Node_Boundary_Level_Plot_2

   integer, allocatable, dimension(:) :: Local_Node_Number_Different_Level

   integer, allocatable, dimension(:,:) :: Local_Node_Number_Level_Max_2_Min
   integer, allocatable, dimension(:,:) :: Number_Node_Difference_Level_3_Node
   double precision, allocatable, dimension(:,:,:,:) :: Position_Node_Boundary_Level_Plot
   double precision, allocatable, dimension(:,:,:) :: Position_Node_Boundary_Level_Plot_tmp

   integer, allocatable, dimension(:) :: Number_Level_Element
   integer, allocatable, dimension(:,:) :: Number_Node_Level_Element
   double precision, allocatable, dimension(:,:,:,:) :: Position_Node_Plot_Mesh_Element
   integer, allocatable, dimension(:,:) :: Level_Element

   integer, allocatable, dimension(:,:,:) :: Counter_Level_Element_Memory_Check
   integer, allocatable, dimension(:) :: Max_Number_Level_Element_tmp 
   integer :: Max_Number_Level_Element

   integer, allocatable, dimension(:) :: Flag_Plot_Element

   integer :: Counter_Percent

   integer :: Number_Edge_Plot_Boundary 
   integer, allocatable, dimension(:) :: Element_Number_Plot_Boundary, Element_Number_Plot_Boundary_tmp
   integer, allocatable, dimension(:,:) :: Local_Node_Number_Plot_Boundary, Local_Node_Number_Plot_Boundary_tmp
   double precision, allocatable, dimension(:,:,:) :: Position_Boundary_Plot, Position_Boundary_Plot_tmp

   !double precision, allocatable, dimension(:,:,:) :: Position_Arrow_PV_Plot
   !double precision, allocatable, dimension(:) :: Length_Arrow_PV_Plot, Stemthick_Plot, &
   !                                               Headthick_Plot, Headlength_Plot
   !double precision :: Stemthick, Headthick, Headlength

   double precision, allocatable, dimension(:) :: Level_Flux_tmp 

   integer, allocatable, dimension(:) :: Level_Contour 
   integer :: Number_Level_Contour, StepSize_Contour, Intercept_Contour
   integer, allocatable, dimension(:) :: Counter_Contour
   integer  :: Counter_Contour_All, Counter_tmp
   double precision, allocatable, dimension(:,:,:) :: Position_Contour_tmp, Position_Contour
   double precision, allocatable, dimension(:,:) :: Value_Contour_Element, Value_Node_in_Element
   double precision, allocatable, dimension(:) :: Value_Contour
   double precision :: Position_Interpolated( 2 ), Position_Node_Line( 2, 2, 2 )
   double precision :: Position_Node_1_tmp( 2 ), Position_Node_2_tmp( 2 ) 
   double precision :: Value_Node_Line( 2, 2 )
   double precision :: Value_Node_Line_1_tmp, Value_Node_Line_2_tmp
   integer :: LN( 3, 3 ) 
   !integer, allocatable, dimension(:,:,:) :: Local_Node_Number_Edge_Line_Plot

   !================================================
   write(*,*)'===================================================='
   write(*,*)'      call Plot_Contour_T_OC'
   write(*,*)'===================================================='
   !================================================

   allocate( Position_Minimum( 2 ) )
   allocate( Position_Maximum( 2 ) )
 
   !Minimum
   Position_Minimum( 1 )= Integer_Initialization_Plus 
   Position_Minimum( 2 )= Integer_Initialization_Plus
   !Maximum
   Position_Maximum( 1 )= Integer_Initialization_Minus
   Position_Maximum( 2 )= Integer_Initialization_Minus
 
   ! No parallel 
   do e= 1, Number_Element
      do i= 1, 4
         if( Position_Minimum( 1 ) > Position_Node( 1, Index_Element_2_Node( i, e ) ) )then
            Position_Minimum( 1 ) = Position_Node( 1, Index_Element_2_Node( i, e ) )
         end if

         if( Position_Maximum( 1 ) < Position_Node( 1, Index_Element_2_Node( i, e ) ) )then
            Position_Maximum( 1 ) = Position_Node( 1, Index_Element_2_Node( i, e ) )
         end if

         if( Position_Minimum( 2 ) > Position_Node( 2, Index_Element_2_Node( i, e ) ) )then
            Position_Minimum( 2 ) = Position_Node( 2, Index_Element_2_Node( i, e ) )
         end if

         if( Position_Maximum( 2 ) < Position_Node( 2, Index_Element_2_Node( i, e ) ) )then
            Position_Maximum( 2 ) = Position_Node( 2, Index_Element_2_Node( i, e ) )
         end if
      end do
   end do
 
   allocate( Translation_Position( 2 ) )

   Translation_Position( 1 )= Translation_X
   Translation_Position( 2 )= Translation_Y

   allocate( Width_PostScript( 2 ) )

   Width_PostScript( 1 )= WidthPS_X 
   Width_PostScript( 2 )= WidthPS_Y

   !======================================================================================================================
   ! Color Var
   !======================================================================================================================
   if( Flag_Color_Var==1 )then
   
      allocate( Translate_Color_Var( 2 ) )
   
      Translate_Color_Var( 1 )= ( Position_Maximum( 1 ) -Position_Minimum( 1 ) )*0.1d0
      Translate_Color_Var( 2 )= ( Position_Maximum( 2 ) -Position_Minimum( 2 ) )*0.15d0
   
      allocate( Position_Maximum_Color_Var( 2, 2, 2 ) )
   
      Position_Maximum_Color_Var( 1, 1, 1 )= Position_Maximum( 1 ) +Translate_Color_Var( 1 ) 
      Position_Maximum_Color_Var( 2, 1, 1 )= Position_Minimum( 2 ) +Translate_Color_Var( 2 ) 
   
      Position_Maximum_Color_Var( 1, 1, 2 )= Position_Maximum_Color_Var( 1, 1, 1 ) 
      Position_Maximum_Color_Var( 2, 1, 2 )= Position_Maximum_Color_Var( 2, 1, 1 ) 
   
      Position_Maximum_Color_Var( 1, 2, 1 )= Position_Maximum_Color_Var( 1, 1, 1 )  
      Position_Maximum_Color_Var( 2, 2, 1 )= Position_Maximum( 2 ) -Translate_Color_Var( 2 ) 
   
      Position_Maximum_Color_Var( 1, 2, 2 )= Position_Maximum_Color_Var( 1, 2, 1 )  
      Position_Maximum_Color_Var( 2, 2, 2 )= Position_Maximum_Color_Var( 2, 2, 1 ) 
   
      deallocate( Translate_Color_Var )
       
      ! No Parallel 
      do i= 1, 2
         do j= 1, 2
            do k= 1, 2
               Position_Maximum_Color_Var( k, j, i ) &
               = ( Position_Maximum_Color_Var( k, j, i ) -Position_Minimum( k ) )  &
                /( Position_Maximum( k ) -Position_Minimum( k ) )*dble( Width_PostScript( k ) ) 
            end do
         end do
      end do
   
      do i= 1, 2
         do j= 1, 2
            do k= 1, 2
               Position_Maximum_Color_Var( k, j, i )= Position_Maximum_Color_Var( k, j, i ) +Translation_Position( k ) 
            end do
         end do
      end do
   
      Position_Maximum_Color_Var( 1, 1, 2 )= Position_Maximum_Color_Var( 1, 1, 1 ) +Width_Color_Var
      Position_Maximum_Color_Var( 1, 2, 2 )= Position_Maximum_Color_Var( 1, 2, 1 ) +Width_Color_Var
   
   else if( Flag_Color_Var==2 )then
   
      allocate( Translate_Color_Var( 2 ) )
   
      Translate_Color_Var( 1 )= ( Position_Maximum( 1 ) -Position_Minimum( 1 ) )*0.25d0
      Translate_Color_Var( 2 )= ( Position_Maximum( 2 ) -Position_Minimum( 2 ) )*0.02d0
   
      allocate( Position_Maximum_Color_Var( 2, 2, 2 ) )
   
      Position_Maximum_Color_Var( 1, 1, 1 )= Position_Minimum( 1 ) +Translate_Color_Var( 1 ) 
      Position_Maximum_Color_Var( 2, 1, 1 )= Position_Maximum( 2 ) +Translate_Color_Var( 2 ) 
   
      Position_Maximum_Color_Var( 1, 1, 2 )= Position_Maximum( 1 ) -Translate_Color_Var( 1 ) 
      Position_Maximum_Color_Var( 2, 1, 2 )= Position_Maximum_Color_Var( 2, 1, 1 ) 
   
      Position_Maximum_Color_Var( 1, 2, 1 )= Position_Maximum_Color_Var( 1, 1, 1 )  
      Position_Maximum_Color_Var( 2, 2, 1 )= Position_Maximum_Color_Var( 2, 1, 1 ) 
   
      Position_Maximum_Color_Var( 1, 2, 2 )= Position_Maximum_Color_Var( 1, 1, 2 ) 
      Position_Maximum_Color_Var( 2, 2, 2 )= Position_Maximum_Color_Var( 2, 2, 1 ) 
   
      deallocate( Translate_Color_Var )
      
      do i= 1, 2
         do j= 1, 2
            do k= 1, 2
               Position_Maximum_Color_Var( k, j, i ) &
               = ( Position_Maximum_Color_Var( k, j, i ) -Position_Minimum( k ) )  &
                /( Position_Maximum( k ) -Position_Minimum( k ) )*dble( Width_PostScript( k ) ) 
            end do
         end do
      end do
   
      do i= 1, 2
         do j= 1, 2
            do k= 1, 2
               Position_Maximum_Color_Var( k, j, i )= Position_Maximum_Color_Var( k, j, i ) +Translation_Position( k ) 
            end do
         end do
      end do
   
      Position_Maximum_Color_Var( 2, 2, 1 )= Position_Maximum_Color_Var( 2, 1, 1 ) +Width_Color_Var
      Position_Maximum_Color_Var( 2, 2, 2 )= Position_Maximum_Color_Var( 2, 1, 2 ) +Width_Color_Var

   end if

   !======================================================================================================================
   ! Position_Node --> Position_Node_Plot 
   !======================================================================================================================
   allocate( Position_Node_Plot_Mesh( 2, Number_Node ) )
  
   !$omp parallel do default( none ) &
   !$omp private( i, j ) &
   !$omp shared( Number_Node, Position_Node_Plot_Mesh, Position_Node, Position_Minimum, Position_Maximum, Width_PostScript ) 
   do i= 1, Number_Node
      do j= 1, 2 
         Position_Node_Plot_Mesh( j, i )= ( Position_Node( j, i ) -Position_Minimum( j ) )  &
                             /( Position_Maximum( j ) -Position_Minimum( j ) )*dble( Width_PostScript( j ) ) 
      end do
   end do
 
   !$omp parallel do default( none ) &
   !$omp private( i, j ) &
   !$omp shared( Number_Node, Position_Node_Plot_Mesh, Translation_Position ) 
   do i= 1, Number_Node
      do j= 1, 2 
         Position_Node_Plot_Mesh( j, i )= Position_Node_Plot_Mesh( j, i ) +Translation_Position( j ) 
      end do
   end do

   !================================================
   ! Level
   !================================================

   if( Flag_Range_Value_Plot==0 )then
      Minimum_Value_Plotted= Double_Precision_Initialization_Plus 
      Maximum_Value_Plotted= Double_Precision_Initialization_Minus
   
      do i= 1, Number_Node
         if( Value_Plotted( i ) > Maximum_Value_Plotted )then
            Maximum_Value_Plotted= Value_Plotted( i ) 
         else if( Value_Plotted( i ) < Minimum_Value_Plotted )then
            Minimum_Value_Plotted= Value_Plotted( i ) 
         end if
      end do

   else if( Flag_Range_Value_Plot==1 )then
      Maximum_Value_Plotted= Maximum_Value_Plotted_Fixed 
      Minimum_Value_Plotted= Minimum_Value_Plotted_Fixed
   end if

   if( Maximum_Value_Plotted < Minimum_Value_Plotted )then
      write(*,*)'         Minimum_Value_Plotted=', Minimum_Value_Plotted
      write(*,*)'         Maximum_Value_Plotted=', Maximum_Value_Plotted
      call Output_Error( 'Plot_Contour_T_OC', 270 )
   end if

   write(*,*)'         Maximum_Value_Plotted=', Maximum_Value_Plotted
   write(*,*)'         Minimum_Value_Plotted=', Minimum_Value_Plotted

   allocate( Value_Boundary_Level( 2, Number_Level ) )
 
   !$omp parallel do default( none ) &
   !$omp private( i ) &
   !$omp shared( Number_Level, Value_Boundary_Level, Maximum_Value_Plotted, Minimum_Value_Plotted ) 
   do i= 1, Number_Level
      Value_Boundary_Level( 1, i ) &
      = ( Maximum_Value_Plotted -Minimum_Value_Plotted )/ Number_Level*( i-1 ) !+Minimum_Value_Plotted
      Value_Boundary_Level( 2, i ) &
      = ( Maximum_Value_Plotted -Minimum_Value_Plotted )/ Number_Level*( i )  !+Minimum_Value_Plotted 
   end do 

   !$omp parallel do default( none ) &
   !$omp private( i ) &
   !$omp shared( Number_Node, Value_Plotted, Minimum_Value_Plotted ) 
   do i= 1, Number_Node 
      Value_Plotted( i )= Value_Plotted( i ) -Minimum_Value_Plotted
   end do

   allocate( Level_Node( Number_Node ) )
   allocate( Counter_Node_Level( Number_Level ) )

   !$omp parallel do default( none ) &
   !$omp private( i ) &
   !$omp shared( Number_Level, Counter_Node_Level ) 
   do i= 1, Number_Level
      Counter_Node_Level( i )= 0
   end do

   do j= 1, Number_Node
      Level_Node( j )= Integer_Initialization

      if( Value_Plotted( j ) < 0.0D0 )then
         Level_Node( j )= 0
      else if( Value_Plotted( j ) > Maximum_Value_Plotted -Minimum_Value_Plotted )then
         Level_Node( j )= Number_Level +1
      else
         do i= 1, Number_Level
            if( Value_Plotted( j ) >= Value_Boundary_Level( 1, i )*Tolerance_Digit_Error_Minus .and. & 
                Value_Plotted( j ) <= Value_Boundary_Level( 2, i )*Tolerance_Digit_Error_Plus )then
    
               Level_Node( j )= i  
               Counter_Node_Level( i )= Counter_Node_Level( i ) +1
            end if 
         end do
      end if
   end do

   deallocate( Counter_Node_Level )


   !===============================================================================================================
   ! Count Node Level for Preventing Memory Corruption
   ! OK : Counter_Level_Element_Memory_Check > Max_Number_Level_Element 
   !===============================================================================================================
   allocate( Counter_Level_Element_Memory_Check( 4, 4, Number_Element ) )
   allocate( Max_Number_Level_Element_tmp( Number_Element ) )

   !$omp parallel do default( none ) &
   !$omp private( e, i, j ) &
   !$omp shared( Number_Element, Max_Number_Level_Element_tmp, Counter_Level_Element_Memory_Check ) & 
   !$omp shared( Level_Node, Index_Element_2_Node ) 
   do e= 1, Number_Element 
      Max_Number_Level_Element_tmp( e )= 0   

      do i= 1, 4 
         do j= 1, 4 
            Counter_Level_Element_Memory_Check( j, i, e ) &
            = abs( Level_Node( Index_Element_2_Node( i, e ) ) -Level_Node( Index_Element_2_Node( j, e ) ) ) +1 
         end do
      end do

      do i= 1, 4 
         do j= 1, 4 
            if( Counter_Level_Element_Memory_Check( j, i, e ) > Max_Number_Level_Element_tmp( e ) )then 
               Max_Number_Level_Element_tmp( e )= Counter_Level_Element_Memory_Check( j, i, e )
            end if
         end do
      end do
   end do

   deallocate( Counter_Level_Element_Memory_Check )

   Max_Number_Level_Element= 0

   do e= 1, Number_Element
      if( Max_Number_Level_Element_tmp( e ) > Max_Number_Level_Element )then
         Max_Number_Level_Element= Max_Number_Level_Element_tmp( e )
      end if
   end do

   write(*,*)'==============================================================='
   write(*,*)'         Max_Number_Level_Element=', Max_Number_Level_Element
   write(*,*)'         Max_Number_Level_Element_Limit=', Max_Number_Level_Element_Limit
   write(*,*)'==============================================================='

   if( Max_Number_Level_Element > Max_Number_Level_Element_Limit )then
      write(*,*)'         Max_Number_Level_Element > Max_Number_Level_Element_Limit '
     !go to 2272 
      ! go to end  of this subroutine 
      !call Output_Error( 'Plot_Contour_T_OC', 353 ) 
   end if

   deallocate( Max_Number_Level_Element_tmp )


   if( Flag_Range_Value_Plot==0 )then
      !$omp parallel do default( none ) &
      !$omp private( i ) &
      !$omp shared( Number_Node, Level_Node, Number_Level, Value_Plotted, Value_Boundary_Level ) 
      do i= 1, Number_Node 
         if( Level_Node( i ) <= 0 )then
            write(*,*)'Level_Node( i )=', Level_Node( i )
            write(*,*)'Number_Level=', Number_Level
            write(*,*)'Value_Plotted( i )=', Value_Plotted( i )
            write(*,*)'Value_Boundary_Level( 1, 1 )', Value_Boundary_Level( 1, 1 )
            write(*,*)'Value_Boundary_Level( 2, Number_Level )', Value_Boundary_Level( 2, Number_Level )
            write(*,*)'i=', i
            call Output_Error( 'Plot_Contour_T_OC', 334 ) 
         else if( Level_Node( i ) > Number_Level )then
            write(*,*)'Level_Node( i )=', Level_Node( i )
            write(*,*)'Number_Level=', Number_Level
            write(*,*)'Value_Plotted( i )=', Value_Plotted( i )
            write(*,*)'Value_Boundary_Level( 1, 1 )', Value_Boundary_Level( 1, 1 )
            write(*,*)'Value_Boundary_Level( 2, Number_Level )', Value_Boundary_Level( 2, Number_Level )
            write(*,*)'i=', i
            call Output_Error( 'Plot_Contour_T_OC', 343 ) 
         end if
      end do
   end if

   !======================================================================================================================
   ! RGB 
   !======================================================================================================================
   allocate( RGB( 3, Number_Level ) )
   allocate( RGB_Lowest( 3 ) )
   allocate( RGB_Highest( 3 ) )

   call Set_RGB( Flag_RGB_Color, Number_Level, RGB, RGB_Lowest, RGB_Highest )

   !do i= 1, Number_Level
   !   write(901,*) i, RGB( 1, i ) 
   !   write(902,*) i, RGB( 2, i )
   !   write(903,*) i, RGB( 3, i )
   !end do

   !================================================
   !write(*,*)'         Plot Configuration'
   !================================================

   allocate( Number_Level_Element( Number_Element ) )
   allocate( Number_Node_Level_Element( Max_Number_Level_Element, Number_Element ) )
   
   allocate( Position_Node_Plot_Mesh_Element( 2, 4, Max_Number_Level_Element, Number_Element ) )
   allocate( Level_Element( Max_Number_Level_Element, Number_Element ) )
  
   allocate( Local_Node_Number_Different_Level( Number_Element ) ) 

   allocate( Level_High( Number_Element ) ) 
   allocate( Level_Low( Number_Element ) ) 
   allocate( Number_Node_Difference_Level( Number_Element ) ) 
 
   allocate( Local_Node_Number_Edge_Line_Plot( 2, 2, Number_Element ) )

   allocate( Local_Node_Number_Edge_Line( 2, Number_Element ) )

   allocate( Value_Boundary_Level_Line( Max_Number_Level_Element, Number_Element ) )

   allocate( Number_Node_Difference_Level_3_Node( 3, Number_Element ) ) 

   allocate( Difference_Position_1_Minus_NewPoint( 2, Number_Element ) ) 
   allocate( Difference_Position_2_Minus_NewPoint( 2, Number_Element ) ) 
 
   allocate( Local_Node_Number_Level_Max_2_Min( 3, Number_Element ) )

   allocate( Position_Node_Boundary_Level_Plot( 2, Max_Number_Level_Element, 3, Number_Element ) ) 
   allocate( Position_Node_Boundary_Level_Plot_tmp( 2, Max_Number_Level_Element, Number_Element ) )

   allocate( Position_Node_Boundary_Level_Plot_1( 2, Max_Number_Level_Element, Number_Element ) )
   allocate( Position_Node_Boundary_Level_Plot_2( 2, Max_Number_Level_Element, Number_Element ) )



   write( Format_Filenumber, '(A2,I1,A1,I1,A1)' ) & 
   '(i', Length_Character_Optimization_Step,'.',Length_Character_Optimization_Step,')'
  
   write( Optimization_Step_Character, Format_Filenumber ) Optimization_Step
 
   Filename_Result= trim( "Contour_T_OC_"//trim(Optimization_Step_Character)//".ps" )

   !==================================================================================================
   write(*,*) '         Create Element Data'
   !==================================================================================================

   allocate( Flag_Plot_Element( Number_Element ) )

      do e= 1, Number_Element

         Flag_Plot_Element( e )=Integer_Initialization_Plus 

         !=======================================================================================
         if( Level_Node( Index_Element_2_Node( 1, e ) ) <= 0 .or. & 
             Level_Node( Index_Element_2_Node( 2, e ) ) <= 0 .or. & 
             Level_Node( Index_Element_2_Node( 3, e ) ) <= 0 )then
         !=======================================================================================

            call Create_PlotData_Element_Homogeneous &
               ( Position_Node_Plot_Mesh, Index_Element_2_Node, Level_Node, & 
                 e, Number_Node, Number_Element, Number_Level, Max_Number_Level_Element, & 
                 !=================================================================
                 Number_Level_Element, Number_Node_Level_Element, Position_Node_Plot_Mesh_Element, Level_Element )

            Flag_Plot_Element( e )= 0

            !call Plot_Gradation_Element_Monocromatic &
            !   ( File_Number, e, Number_Element, Max_Number_Level_Element, & 
            !     Number_Level_Element, Number_Node_Level_Element, Position_Node_Plot_Mesh_Element, 0 )
 
         !=======================================================================================
         ! Node Value Exceed the Value of Color Range
         else if( Level_Node( Index_Element_2_Node( 1, e ) ) >= Number_Level +1 .or. & 
                Level_Node( Index_Element_2_Node( 2, e ) ) >= Number_Level +1 .or. & 
                Level_Node( Index_Element_2_Node( 3, e ) ) >= Number_Level +1 )then
         !=======================================================================================

            call Create_PlotData_Element_Homogeneous &
               ( Position_Node_Plot_Mesh, Index_Element_2_Node, Level_Node, & 
                 e, Number_Node, Number_Element, Number_Level, Max_Number_Level_Element, & 
                 !=================================================================
                 Number_Level_Element, Number_Node_Level_Element, Position_Node_Plot_Mesh_Element, Level_Element )

            Flag_Plot_Element( e )= 1

            !call Plot_Gradation_Element_Monocromatic &
            !   ( File_Number, e, Number_Element, Max_Number_Level_Element, & 
            !     Number_Level_Element, Number_Node_Level_Element, Position_Node_Plot_Mesh_Element, 1 )
 
         !=======================================================================================
         else if( Level_Node( Index_Element_2_Node( 1, e ) ) == Level_Node( Index_Element_2_Node( 2, e ) ) .and. & 
                Level_Node( Index_Element_2_Node( 2, e ) ) == Level_Node( Index_Element_2_Node( 3, e ) ) .and. & 
                Level_Node( Index_Element_2_Node( 3, e ) ) == Level_Node( Index_Element_2_Node( 1, e ) ) )then 
         !=======================================================================================
    
            call Create_PlotData_Element_Homogeneous &
               ( Position_Node_Plot_Mesh, Index_Element_2_Node, Level_Node, & 
                 e, Number_Node, Number_Element, Number_Level, Max_Number_Level_Element, & 
                 !=================================================================
                 Number_Level_Element, Number_Node_Level_Element, Position_Node_Plot_Mesh_Element, Level_Element )

            Flag_Plot_Element( e )= 2

            !call Plot_Gradation_Element &
            !   ( File_Number, e, Number_Element, Number_Level, Max_Number_Level_Element, & 
            !     Number_Level_Element, Number_Node_Level_Element, &
            !     Position_Node_Plot_Mesh_Element, RGB, Level_Element )
 
         !=======================================================================================
         else if( Level_Node( Index_Element_2_Node( 1, e ) ) == Level_Node( Index_Element_2_Node( 2, e ) ) .and. & 
                Level_Node( Index_Element_2_Node( 1, e ) ) > Level_Node( Index_Element_2_Node( 3, e ) ) .and. &
                Level_Node( Index_Element_2_Node( 2, e ) ) > Level_Node( Index_Element_2_Node( 3, e ) ) )then 
         !=======================================================================================

            Local_Node_Number_Different_Level( e )= 3 

            Local_Node_Number_Edge_Line_Plot( 1, 1, e )= Local_Node_Number_Different_Level( e ) 
            Local_Node_Number_Edge_Line_Plot( 2, 1, e )= 1
            Local_Node_Number_Edge_Line_Plot( 1, 2, e )= Local_Node_Number_Different_Level( e )
            Local_Node_Number_Edge_Line_Plot( 2, 2, e )= 2

            Level_High( e )= Level_Node( Index_Element_2_Node( Local_Node_Number_Edge_Line_Plot( 2, 1, e ), e ) ) 
            Level_Low( e )=  Level_Node( Index_Element_2_Node( Local_Node_Number_Different_Level( e ), e ) )

            Number_Node_Difference_Level( e )= Level_High( e ) -Level_Low( e ) 

            if( Number_Node_Difference_Level( e )<=0 )then
               write(*,*)'Number_Node_Difference_Level( e )==0'
               !stop
              !go to 2272
            end if

            do i= 1, Number_Node_Difference_Level( e )
               Value_Boundary_Level_Line( i, e )& 
               = Value_Boundary_Level&
                 ( 1, i +Level_Node( Index_Element_2_Node( Local_Node_Number_Different_Level( e ), e ) ) )  
            end do

            if( Level_Node( Index_Element_2_Node( Local_Node_Number_Different_Level( e ), e ) ) &
                +Number_Node_Difference_Level( e ) & 
                /= Level_Node( Index_Element_2_Node( Local_Node_Number_Edge_Line_Plot( 2, 1, e ), e ) ) .or. &
               Level_Node( Index_Element_2_Node( Local_Node_Number_Different_Level( e ), e ) ) &
                +Number_Node_Difference_Level( e ) & 
                /= Level_Node( Index_Element_2_Node( Local_Node_Number_Edge_Line_Plot( 2, 2, e ), e ) ) )then
 
               write(*,*)'Level_Node( Index_Element_2_Node( Local_Node_Number_Different_Level( e ), e ) )=', & 
               Level_Node( Index_Element_2_Node( Local_Node_Number_Different_Level( e ), e ) )  
               write(*,*)'Number_Node_Difference_Level( e )=', Number_Node_Difference_Level( e )
               write(*,*)'Level_Node( Index_Element_2_Node( Local_Node_Number_Edge_Line_Plot( 2, 1, e ), e ) )=', & 
                      Level_Node( Index_Element_2_Node( Local_Node_Number_Edge_Line_Plot( 2, 1, e ), e ) ) 
               write(*,*)'Level_Node( Index_Element_2_Node( Local_Node_Number_Edge_Line_Plot( 2, 2, e ), e ) )=', &
                      Level_Node( Index_Element_2_Node( Local_Node_Number_Edge_Line_Plot( 2, 2, e ), e ) ) 
               write(*,*)'Plot_Contour_T_OC.f90 264'
               !stop
              !go to 2272
            end if

            !=================================================================================
            ! 3 --> 1
            !=================================================================================

            Local_Node_Number_Edge_Line( 1, e )= Local_Node_Number_Different_Level( e )
            Local_Node_Number_Edge_Line( 2, e )= Local_Node_Number_Edge_Line_Plot( 2, 1, e )

            do i= 1, Number_Node_Difference_Level( e )

               call Compute_Position_Point_Plot_Result &
                  ( Value_Plotted, Value_Boundary_Level_Line, &
                    Local_Node_Number_Edge_Line, &
                    Position_Node_Plot_Mesh, Number_Node, Number_Node_Difference_Level( e ), &
                    Index_Element_2_Node, Number_Element, e, Max_Number_Level_Element, &
                    !======================================================================
                    Position_Node_Boundary_Level_Plot_1 )

            end do

            do i= 1, Number_Node_Difference_Level( e )
               do k= 1, 2
                  Difference_Position_1_Minus_NewPoint( k, e ) &
                  = Position_Node_Plot_Mesh( k, Index_Element_2_Node( Local_Node_Number_Edge_Line( 1, e ) ,e ) ) & 
                   -Position_Node_Boundary_Level_Plot_1( k, i, e )

                  Difference_Position_2_Minus_NewPoint( k, e ) &
                  = Position_Node_Plot_Mesh( k, Index_Element_2_Node( Local_Node_Number_Edge_Line( 2, e ) ,e ) ) & 
                   -Position_Node_Boundary_Level_Plot_1( k, i, e )
               end do

               if( Difference_Position_1_Minus_NewPoint( 1, e )*Difference_Position_2_Minus_NewPoint( 1, e ) &
                   > Tolerance_Plot .or. &
                   Difference_Position_1_Minus_NewPoint( 2, e )*Difference_Position_2_Minus_NewPoint( 2, e ) &
                   > Tolerance_Plot )then    

                  write(287,*)'e=', e 
                  write(287,*)'Number_Node_Difference_Level( e )=', Number_Node_Difference_Level( e ) 
                  do k= 1, Number_Node_Difference_Level( e )
                     write(287,*) k, Position_Node_Boundary_Level_Plot_1( 1, k, e ), & 
                               Position_Node_Boundary_Level_Plot_1( 2, k, e ) 
                  end do
   
     
                  write(287,*)' ' 
                  write(*,*)'Plot_Contour_T_OC.f90 304' 
                  !stop
                 !go to 2272
               end if
            end do
    
            !=================================================================================
            ! 3 --> 2
            !=================================================================================

            Local_Node_Number_Edge_Line( 1, e )= Local_Node_Number_Different_Level( e )
            Local_Node_Number_Edge_Line( 2, e )= Local_Node_Number_Edge_Line_Plot( 2, 2, e )

            do i= 1, Number_Node_Difference_Level( e )

               call Compute_Position_Point_Plot_Result &
                  ( Value_Plotted, Value_Boundary_Level_Line, &
                    Local_Node_Number_Edge_Line, &
                    Position_Node_Plot_Mesh, Number_Node, Number_Node_Difference_Level( e ), &
                    Index_Element_2_Node, Number_Element, e, Max_Number_Level_Element, &
                    !======================================================================
                    Position_Node_Boundary_Level_Plot_2 )

            end do

            do i= 1, Number_Node_Difference_Level( e )
               do k= 1, 2
                  Difference_Position_1_Minus_NewPoint( k, e ) &
                  = Position_Node_Plot_Mesh( k, Index_Element_2_Node( Local_Node_Number_Edge_Line( 1, e ) ,e ) ) & 
                   -Position_Node_Boundary_Level_Plot_2( k, i, e )

                  Difference_Position_2_Minus_NewPoint( k, e ) &
                  = Position_Node_Plot_Mesh( k, Index_Element_2_Node( Local_Node_Number_Edge_Line( 2, e ) ,e ) ) & 
                   -Position_Node_Boundary_Level_Plot_2( k, i, e )
               end do

               if( Difference_Position_1_Minus_NewPoint( 1, e )*Difference_Position_2_Minus_NewPoint( 1, e ) &
                   > Tolerance_Plot .or. &
                   Difference_Position_1_Minus_NewPoint( 2, e )*Difference_Position_2_Minus_NewPoint( 2, e ) & 
                   > Tolerance_Plot )then    

                  write(288,*)'e=', e 
                  write(288,*)'Number_Node_Difference_Level( e )=', Number_Node_Difference_Level( e ) 

                  do k= 1, Number_Node_Difference_Level( e )
                     write(288,*) k, Position_Node_Boundary_Level_Plot_2( 1, k, e ), & 
                               Position_Node_Boundary_Level_Plot_2( 2, k, e ) 
                  end do
   
                  write(288,*)' ' 
                  write(*,*)'Plot_Contour_T_OC.f90 381' 
                  !stop
                 !go to 2272
               end if
            end do
 
            Flag_Plot_Element( e )= 3

            !call Plot_Element_InHomogeneous &
            !   ( Position_Node_Plot_Mesh, Index_Element_2_Node, Level_Node, RGB, & 
            !     File_Number, e, Number_Node, Number_Element, Number_Level, Max_Number_Level_Element, & 
            !     Number_Node_Difference_Level( e ), & 
            !     Local_Node_Number_Different_Level( e ), Local_Node_Number_Edge_Line_Plot,  &
            !     Position_Node_Boundary_Level_Plot_1, Position_Node_Boundary_Level_Plot_2 )


 
            !call Create_PlotData_Element_InHomogeneous &
            !   ( Position_Node_Plot_Mesh, Index_Element_2_Node, Level_Node, & 
            !     e, Number_Node, Number_Element, & 
            !     Number_Node_Difference_Level( e ), &
            !     Local_Node_Number_Different_Level( e ), Local_Node_Number_Edge_Line_Plot,  &
            !     Position_Node_Boundary_Level_Plot_1, Position_Node_Boundary_Level_Plot_2, &
            !     !=================================================================
            !     Number_Level_Element, Number_Node_Level_Element, Position_Node_Plot_Mesh_Element, Level_Element )

            !call Plot_Gradation_Element &
            !   ( File_Number, e, Number_Element, Number_Level, & 
            !     Number_Level_Element, Number_Node_Level_Element, Position_Node_Plot_Mesh_Element, RGB, Level_Element )
 
         !=======================================================================================
         else if( Level_Node( Index_Element_2_Node( 1, e ) ) == Level_Node( Index_Element_2_Node( 2, e ) ) .and. & 
                Level_Node( Index_Element_2_Node( 1, e ) ) < Level_Node( Index_Element_2_Node( 3, e ) ) .and. & 
                Level_Node( Index_Element_2_Node( 2, e ) ) < Level_Node( Index_Element_2_Node( 3, e ) ) )then 
         !=======================================================================================

            Local_Node_Number_Different_Level( e )= 3 

            Local_Node_Number_Edge_Line_Plot( 1, 1, e )= Local_Node_Number_Different_Level( e ) 
            Local_Node_Number_Edge_Line_Plot( 2, 1, e )= 1
            Local_Node_Number_Edge_Line_Plot( 1, 2, e )= Local_Node_Number_Different_Level( e )
            Local_Node_Number_Edge_Line_Plot( 2, 2, e )= 2

            Level_High( e )=  Level_Node( Index_Element_2_Node( Local_Node_Number_Different_Level( e ), e ) )
            Level_Low( e )= Level_Node( Index_Element_2_Node( Local_Node_Number_Edge_Line_Plot( 2, 1, e ), e ) ) 

            Number_Node_Difference_Level( e )= Level_High( e ) -Level_Low( e ) 

            if( Number_Node_Difference_Level( e )<=0 )then
               write(*,*)'Number_Node_Difference_Level( e )==0'
               !stop
              !go to 2272
            end if

            do i= 1, Number_Node_Difference_Level( e )
               Value_Boundary_Level_Line( i, e )& 
               = Value_Boundary_Level&
                 ( 2, -i +Level_Node( Index_Element_2_Node( Local_Node_Number_Different_Level( e ), e ) ) )  
            end do

            if( Level_Node( Index_Element_2_Node( Local_Node_Number_Different_Level( e ), e ) ) & 
                -Number_Node_Difference_Level( e ) & 
                /= Level_Node( Index_Element_2_Node( Local_Node_Number_Edge_Line_Plot( 2, 1, e ), e ) ) .or. &
               Level_Node( Index_Element_2_Node( Local_Node_Number_Different_Level( e ), e ) ) & 
                -Number_Node_Difference_Level( e ) & 
                /= Level_Node( Index_Element_2_Node( Local_Node_Number_Edge_Line_Plot( 2, 2, e ), e ) ) )then
 
               write(*,*)'Level_Node( Index_Element_2_Node( Local_Node_Number_Different_Level( e ), e ) )=', & 
               Level_Node( Index_Element_2_Node( Local_Node_Number_Different_Level( e ), e ) )  
               write(*,*)'Number_Node_Difference_Level( e )=', Number_Node_Difference_Level( e )
               write(*,*)'Level_Node( Index_Element_2_Node( Local_Node_Number_Edge_Line_Plot( 2, 1, e ), e ) )=', & 
                      Level_Node( Index_Element_2_Node( Local_Node_Number_Edge_Line_Plot( 2, 1, e ), e ) ) 
               write(*,*)'Level_Node( Index_Element_2_Node( Local_Node_Number_Edge_Line_Plot( 2, 2, e ), e ) )=', &
                      Level_Node( Index_Element_2_Node( Local_Node_Number_Edge_Line_Plot( 2, 2, e ), e ) ) 
               write(*,*)'Plot_Contour_T_OC.f90 448'
               !stop
              !go to 2272
            end if

            !=================================================================================
            ! 3 --> 1
            !=================================================================================

            Local_Node_Number_Edge_Line( 1, e )= Local_Node_Number_Different_Level( e )
            Local_Node_Number_Edge_Line( 2, e )= Local_Node_Number_Edge_Line_Plot( 2, 1, e )

            do i= 1, Number_Node_Difference_Level( e )

               call Compute_Position_Point_Plot_Result &
                  ( Value_Plotted, Value_Boundary_Level_Line, &
                    Local_Node_Number_Edge_Line, &
                    Position_Node_Plot_Mesh, Number_Node, Number_Node_Difference_Level( e ), &
                    Index_Element_2_Node, Number_Element, e, Max_Number_Level_Element, &
                    !======================================================================
                    Position_Node_Boundary_Level_Plot_1 )

            end do

            do i= 1, Number_Node_Difference_Level( e )
               do k= 1, 2
                  Difference_Position_1_Minus_NewPoint( k, e ) &
                  = Position_Node_Plot_Mesh( k, Index_Element_2_Node( Local_Node_Number_Edge_Line( 1, e ) ,e ) ) & 
                   -Position_Node_Boundary_Level_Plot_1( k, i, e )

                  Difference_Position_2_Minus_NewPoint( k, e ) &
                  = Position_Node_Plot_Mesh( k, Index_Element_2_Node( Local_Node_Number_Edge_Line( 2, e ) ,e ) ) & 
                   -Position_Node_Boundary_Level_Plot_1( k, i, e )
               end do

               if( Difference_Position_1_Minus_NewPoint( 1, e )*Difference_Position_2_Minus_NewPoint( 1, e ) &
                   > Tolerance_Plot .or. &
                   Difference_Position_1_Minus_NewPoint( 2, e )*Difference_Position_2_Minus_NewPoint( 2, e ) &
                   > Tolerance_Plot )then    

                  write(287,*)'e=', e 
                  write(287,*)'Number_Node_Difference_Level( e )=', Number_Node_Difference_Level( e ) 
                  write(287,*)'Level_Node( Index_Element_2_Node( Local_Node_Number_Edge_Line( 1, e ), e ) )=', & 
                           Level_Node( Index_Element_2_Node( Local_Node_Number_Edge_Line( 1, e ), e ) ) 
                  write(287,*)'Level_Node( Index_Element_2_Node( Local_Node_Number_Edge_Line( 2, e ), e ) )=', &
                           Level_Node( Index_Element_2_Node( Local_Node_Number_Edge_Line( 2, e ), e ) ) 
                  write(287,*) Position_Node_Plot_Mesh &
                           ( 1, Index_Element_2_Node(Local_Node_Number_Edge_Line( 1, e ) ,e ) ), &
                           Position_Node_Plot_Mesh &
                           ( 2, Index_Element_2_Node(Local_Node_Number_Edge_Line( 1, e ) ,e ) )
                  write(287,*)' ' 
                  do k= 1, Number_Node_Difference_Level( e )
                     write(287,*) Position_Node_Boundary_Level_Plot_1( 1, k, e ), & 
                              Position_Node_Boundary_Level_Plot_1( 2, k, e ) 
                  end do
   
                  write(287,*)' ' 
                  write(287,*) Position_Node_Plot_Mesh &
                           ( 1, Index_Element_2_Node(Local_Node_Number_Edge_Line( 2, e ) ,e ) ), &
                           Position_Node_Plot_Mesh & 
                           ( 2, Index_Element_2_Node(Local_Node_Number_Edge_Line( 2, e ) ,e ) ) 
     
                  write(287,*)' ' 
                  write(287,*)'=====' 
                  write(287,*)' ' 
                  write(*,*)'Plot_Contour_T_OC.f90 509' 
                  !stop
                 !go to 2272
               end if
            end do
    
            !=================================================================================
            ! 3 --> 2
            !=================================================================================

            Local_Node_Number_Edge_Line( 1, e )= Local_Node_Number_Different_Level( e )
            Local_Node_Number_Edge_Line( 2, e )= Local_Node_Number_Edge_Line_Plot( 2, 2, e )

            do i= 1, Number_Node_Difference_Level( e )

               call Compute_Position_Point_Plot_Result &
                  ( Value_Plotted, Value_Boundary_Level_Line, &
                    Local_Node_Number_Edge_Line, &
                    Position_Node_Plot_Mesh, Number_Node, Number_Node_Difference_Level( e ), &
                    Index_Element_2_Node, Number_Element, e, Max_Number_Level_Element, &
                    !======================================================================
                    Position_Node_Boundary_Level_Plot_2 )

            end do

            do i= 1, Number_Node_Difference_Level( e )
               do k= 1, 2
                  Difference_Position_1_Minus_NewPoint( k, e ) &
                  = Position_Node_Plot_Mesh( k, Index_Element_2_Node( Local_Node_Number_Edge_Line( 1, e ) ,e ) ) & 
                   -Position_Node_Boundary_Level_Plot_2( k, i, e )

                  Difference_Position_2_Minus_NewPoint( k, e ) &
                  = Position_Node_Plot_Mesh( k, Index_Element_2_Node( Local_Node_Number_Edge_Line( 2, e ) ,e ) ) & 
                   -Position_Node_Boundary_Level_Plot_2( k, i, e )
               end do

               if( Difference_Position_1_Minus_NewPoint( 1, e )*Difference_Position_2_Minus_NewPoint( 1, e ) &
                   > Tolerance_Plot .or. &
                   Difference_Position_1_Minus_NewPoint( 2, e )*Difference_Position_2_Minus_NewPoint( 2, e ) &
                   > Tolerance_Plot )then    

                  write(288,*)'e=', e 
                  write(288,*)'Number_Node_Difference_Level( e )=', Number_Node_Difference_Level( e ) 
                  do k= 1, Number_Node_Difference_Level( e )
                     write(288,*) Position_Node_Boundary_Level_Plot_2( 1, k, e ), & 
                              Position_Node_Boundary_Level_Plot_2( 2, k, e ) 
                  end do
   
                  write(288,*)' ' 
                  write(*,*)'Plot_Contour_T_OC.f90 573' 
                  !stop
                 !go to 2272
               end if
            end do
    
            Flag_Plot_Element( e )= 3

            !call Plot_Element_InHomogeneous &
            !   ( Position_Node_Plot_Mesh, Index_Element_2_Node, Level_Node, RGB, & 
            !     File_Number, e, Number_Node, Number_Element, Number_Level, Max_Number_Level_Element, & 
            !     Number_Node_Difference_Level( e ), &
            !     Local_Node_Number_Different_Level( e ), Local_Node_Number_Edge_Line_Plot,  &
            !     Position_Node_Boundary_Level_Plot_1, Position_Node_Boundary_Level_Plot_2 )
! 
!            call Create_PlotData_Element_InHomogeneous &
!               ( Position_Node_Plot_Mesh, Index_Element_2_Node, Level_Node, & 
!                 e, Number_Node, Number_Element, & 
!                 Number_Node_Difference_Level( e ), &
!                 Local_Node_Number_Different_Level( e ), Local_Node_Number_Edge_Line_Plot,  &
!                 Position_Node_Boundary_Level_Plot_1, Position_Node_Boundary_Level_Plot_2, &
!                 !=================================================================
!                 Number_Level_Element, Number_Node_Level_Element, Position_Node_Plot_Mesh_Element, Level_Element )
!
!            call Plot_Gradation_Element &
!               ( File_Number, e, Number_Element, Number_Level, & 
!                 Number_Level_Element, Number_Node_Level_Element, Position_Node_Plot_Mesh_Element, RGB, Level_Element )

         !=======================================================================================
         else if( Level_Node( Index_Element_2_Node( 2, e ) ) == Level_Node( Index_Element_2_Node( 3, e ) ) .and. & 
                Level_Node( Index_Element_2_Node( 2, e ) ) > Level_Node( Index_Element_2_Node( 1, e ) ) .and. & 
                Level_Node( Index_Element_2_Node( 3, e ) ) > Level_Node( Index_Element_2_Node( 1, e ) ) )then 
         !=======================================================================================

            Local_Node_Number_Different_Level( e )= 1 


            Local_Node_Number_Edge_Line_Plot( 1, 1, e )= Local_Node_Number_Different_Level( e ) 
            Local_Node_Number_Edge_Line_Plot( 2, 1, e )= 2
            Local_Node_Number_Edge_Line_Plot( 1, 2, e )= Local_Node_Number_Different_Level( e )
            Local_Node_Number_Edge_Line_Plot( 2, 2, e )= 3

            Level_High( e )= Level_Node( Index_Element_2_Node( Local_Node_Number_Edge_Line_Plot( 2, 1, e ), e ) ) 
            Level_Low( e )=  Level_Node( Index_Element_2_Node( Local_Node_Number_Different_Level( e ), e ) )

            Number_Node_Difference_Level( e )= Level_High( e ) -Level_Low( e ) 

            if( Number_Node_Difference_Level( e )<=0 )then
               write(*,*)'Number_Node_Difference_Level( e )==0'
               !stop
              !go to 2272
            end if


            do i= 1, Number_Node_Difference_Level( e )
               Value_Boundary_Level_Line( i, e )& 
               = Value_Boundary_Level&
                 ( 1, i +Level_Node( Index_Element_2_Node( Local_Node_Number_Different_Level( e ), e ) ) )  
            end do

            if( Level_Node( Index_Element_2_Node( Local_Node_Number_Different_Level( e ), e ) ) &
               +Number_Node_Difference_Level( e ) & 
                /= Level_Node( Index_Element_2_Node( Local_Node_Number_Edge_Line_Plot( 2, 1, e ), e ) ) .or. &
               Level_Node( Index_Element_2_Node( Local_Node_Number_Different_Level( e ), e ) ) &
               +Number_Node_Difference_Level( e ) & 
                /= Level_Node( Index_Element_2_Node( Local_Node_Number_Edge_Line_Plot( 2, 2, e ), e ) ) )then
 
               write(*,*)'Level_Node( Index_Element_2_Node( Local_Node_Number_Different_Level( e ), e ) )=', & 
               Level_Node( Index_Element_2_Node( Local_Node_Number_Different_Level( e ), e ) )  
               write(*,*)'Number_Node_Difference_Level( e )=', Number_Node_Difference_Level( e )
               write(*,*)'Level_Node( Index_Element_2_Node( Local_Node_Number_Edge_Line_Plot( 2, 1, e ), e ) )=', & 
                      Level_Node( Index_Element_2_Node( Local_Node_Number_Edge_Line_Plot( 2, 1, e ), e ) ) 
               write(*,*)'Level_Node( Index_Element_2_Node( Local_Node_Number_Edge_Line_Plot( 2, 2, e ), e ) )=', &
                      Level_Node( Index_Element_2_Node( Local_Node_Number_Edge_Line_Plot( 2, 2, e ), e ) ) 
               write(*,*)'Plot_Contour_T_OC.f90 634'
               !stop
              !go to 2272
            end if

            !=================================================================================
            ! 3 --> 1
            !=================================================================================

            Local_Node_Number_Edge_Line( 1, e )= Local_Node_Number_Different_Level( e )
            Local_Node_Number_Edge_Line( 2, e )= Local_Node_Number_Edge_Line_Plot( 2, 1, e )

            do i= 1, Number_Node_Difference_Level( e )

               call Compute_Position_Point_Plot_Result &
                  ( Value_Plotted, Value_Boundary_Level_Line, &
                    Local_Node_Number_Edge_Line, &
                    Position_Node_Plot_Mesh, Number_Node, Number_Node_Difference_Level( e ), &
                    Index_Element_2_Node, Number_Element, e, Max_Number_Level_Element, &
                    !======================================================================
                    Position_Node_Boundary_Level_Plot_1 )

            end do

            do i= 1, Number_Node_Difference_Level( e )
               do k= 1, 2
                  Difference_Position_1_Minus_NewPoint( k, e ) &
                  = Position_Node_Plot_Mesh( k, Index_Element_2_Node( Local_Node_Number_Edge_Line( 1, e ) ,e ) ) & 
                   -Position_Node_Boundary_Level_Plot_1( k, i, e )

                  Difference_Position_2_Minus_NewPoint( k, e ) &
                  = Position_Node_Plot_Mesh( k, Index_Element_2_Node( Local_Node_Number_Edge_Line( 2, e ) ,e ) ) & 
                   -Position_Node_Boundary_Level_Plot_1( k, i, e )
               end do

               if( Difference_Position_1_Minus_NewPoint( 1, e )*Difference_Position_2_Minus_NewPoint( 1, e ) &
                   > Tolerance_Plot .or. &
                   Difference_Position_1_Minus_NewPoint( 2, e )*Difference_Position_2_Minus_NewPoint( 2, e ) &
                   > Tolerance_Plot )then    

                  write(287,*)'e=', e 
                  write(287,*)'Number_Node_Difference_Level( e )=', Number_Node_Difference_Level( e ) 
                  write(287,*)'Level_Node( Index_Element_2_Node( Local_Node_Number_Edge_Line( 1, e ), e ) )=', & 
                           Level_Node( Index_Element_2_Node( Local_Node_Number_Edge_Line( 1, e ), e ) ) 
                  write(287,*)'Level_Node( Index_Element_2_Node( Local_Node_Number_Edge_Line( 2, e ), e ) )=', &
                           Level_Node( Index_Element_2_Node( Local_Node_Number_Edge_Line( 2, e ), e ) ) 
                  do k= 1, Number_Node_Difference_Level( e )
                     write(287,*) Position_Node_Boundary_Level_Plot_1( 1, k, e ), & 
                              Position_Node_Boundary_Level_Plot_1( 2, k, e ) 
                  end do
   
                  write(287,*)' ' 
                  write(*,*)'Plot_Contour_T_OC.f90 702' 
                  !stop
                 !go to 2272
               end if
            end do
    
            !=================================================================================
            ! 3 --> 2
            !=================================================================================

            Local_Node_Number_Edge_Line( 1, e )= Local_Node_Number_Different_Level( e )
            Local_Node_Number_Edge_Line( 2, e )= Local_Node_Number_Edge_Line_Plot( 2, 2, e )

            do i= 1, Number_Node_Difference_Level( e )

               call Compute_Position_Point_Plot_Result &
                  ( Value_Plotted, Value_Boundary_Level_Line, &
                    Local_Node_Number_Edge_Line, &
                    Position_Node_Plot_Mesh, Number_Node, Number_Node_Difference_Level( e ), &
                    Index_Element_2_Node, Number_Element, e, Max_Number_Level_Element, &
                    !======================================================================
                    Position_Node_Boundary_Level_Plot_2 )

            end do

            do i= 1, Number_Node_Difference_Level( e )

               do k= 1, 2
                  Difference_Position_1_Minus_NewPoint( k, e ) &
                  = Position_Node_Plot_Mesh( k, Index_Element_2_Node( Local_Node_Number_Edge_Line( 1, e ) ,e ) ) & 
                   -Position_Node_Boundary_Level_Plot_2( k, i, e )

                  Difference_Position_2_Minus_NewPoint( k, e ) &
                  = Position_Node_Plot_Mesh( k, Index_Element_2_Node( Local_Node_Number_Edge_Line( 2, e ) ,e ) ) & 
                   -Position_Node_Boundary_Level_Plot_2( k, i, e )
               end do

               if( Difference_Position_1_Minus_NewPoint( 1, e )*Difference_Position_2_Minus_NewPoint( 1, e ) &
                   > Tolerance_Plot .or. &
                   Difference_Position_1_Minus_NewPoint( 2, e )*Difference_Position_2_Minus_NewPoint( 2, e ) &
                   > Tolerance_Plot )then  

                  write(288,*)'e=', e 
                  write(288,*)'Number_Node_Difference_Level( e )=', Number_Node_Difference_Level( e ) 
                  do k= 1, Number_Node_Difference_Level( e )
                     write(288,*) Position_Node_Boundary_Level_Plot_2( 1, k, e ), & 
                              Position_Node_Boundary_Level_Plot_2( 2, k, e ) 
                  end do
   
                  write(288,*)' ' 
                  write(*,*)'Plot_Contour_T_OC.f90 763' 
                  !stop
                 !go to 2272
               end if
            end do
    
            Flag_Plot_Element( e )= 3

            !call Plot_Element_InHomogeneous &
            !   ( Position_Node_Plot_Mesh, Index_Element_2_Node, Level_Node, RGB, & 
            !     File_Number, e, Number_Node, Number_Element, Number_Level, Max_Number_Level_Element, & 
            !     Number_Node_Difference_Level( e ), &
            !     Local_Node_Number_Different_Level( e ), Local_Node_Number_Edge_Line_Plot,  &
            !     Position_Node_Boundary_Level_Plot_1, Position_Node_Boundary_Level_Plot_2 )
! 
!            call Create_PlotData_Element_InHomogeneous &
!               ( Position_Node_Plot_Mesh, Index_Element_2_Node, Level_Node, & 
!                 e, Number_Node, Number_Element, & 
!                 Number_Node_Difference_Level( e ), &
!                 Local_Node_Number_Different_Level( e ), Local_Node_Number_Edge_Line_Plot,  &
!                 Position_Node_Boundary_Level_Plot_1, Position_Node_Boundary_Level_Plot_2, &
!                 !=================================================================
!                 Number_Level_Element, Number_Node_Level_Element, Position_Node_Plot_Mesh_Element, Level_Element )
!
!            call Plot_Gradation_Element &
!               ( File_Number, e, Number_Element, Number_Level, & 
!                 Number_Level_Element, Number_Node_Level_Element, Position_Node_Plot_Mesh_Element, RGB, Level_Element )


         !=======================================================================================
         else if( Level_Node( Index_Element_2_Node( 2, e ) ) == Level_Node( Index_Element_2_Node( 3, e ) ) .and. & 
                Level_Node( Index_Element_2_Node( 2, e ) ) < Level_Node( Index_Element_2_Node( 1, e ) ) .and. & 
                Level_Node( Index_Element_2_Node( 3, e ) ) < Level_Node( Index_Element_2_Node( 1, e ) ) )then 
         !=======================================================================================

            Local_Node_Number_Different_Level( e )= 1 


            Local_Node_Number_Edge_Line_Plot( 1, 1, e )= Local_Node_Number_Different_Level( e ) 
            Local_Node_Number_Edge_Line_Plot( 2, 1, e )= 2
            Local_Node_Number_Edge_Line_Plot( 1, 2, e )= Local_Node_Number_Different_Level( e )
            Local_Node_Number_Edge_Line_Plot( 2, 2, e )= 3

            Level_High( e )=  Level_Node( Index_Element_2_Node( Local_Node_Number_Different_Level( e ), e ) )
            Level_Low( e )= Level_Node( Index_Element_2_Node( Local_Node_Number_Edge_Line_Plot( 2, 1, e ), e ) ) 

            Number_Node_Difference_Level( e )= Level_High( e ) -Level_Low( e ) 

            if( Number_Node_Difference_Level( e )<=0 )then
               write(*,*)'Number_Node_Difference_Level( e )==0'
               !stop
              !go to 2272
            end if


            do i= 1, Number_Node_Difference_Level( e )
               Value_Boundary_Level_Line( i, e )& 
               = Value_Boundary_Level&
                ( 2, -i +Level_Node( Index_Element_2_Node( Local_Node_Number_Different_Level( e ), e ) ) )  
            end do

            if( Level_Node( Index_Element_2_Node( Local_Node_Number_Different_Level( e ), e ) ) & 
                -Number_Node_Difference_Level( e ) & 
                /= Level_Node( Index_Element_2_Node( Local_Node_Number_Edge_Line_Plot( 2, 1, e ), e ) ) .or. &
               Level_Node( Index_Element_2_Node( Local_Node_Number_Different_Level( e ), e ) ) & 
                -Number_Node_Difference_Level( e ) & 
                /= Level_Node( Index_Element_2_Node( Local_Node_Number_Edge_Line_Plot( 2, 2, e ), e ) ) )then
 
               write(*,*)'Level_Node( Index_Element_2_Node( Local_Node_Number_Different_Level( e ), e ) )=', & 
               Level_Node( Index_Element_2_Node( Local_Node_Number_Different_Level( e ), e ) )  
               write(*,*)'Number_Node_Difference_Level( e )=', Number_Node_Difference_Level( e )
               write(*,*)'Level_Node( Index_Element_2_Node( Local_Node_Number_Edge_Line_Plot( 2, 1, e ), e ) )=', & 
                      Level_Node( Index_Element_2_Node( Local_Node_Number_Edge_Line_Plot( 2, 1, e ), e ) ) 
               write(*,*)'Level_Node( Index_Element_2_Node( Local_Node_Number_Edge_Line_Plot( 2, 2, e ), e ) )=', &
                      Level_Node( Index_Element_2_Node( Local_Node_Number_Edge_Line_Plot( 2, 2, e ), e ) ) 
               write(*,*)'Plot_Contour_T_OC.f90 448'
               !stop
              !go to 2272
            end if

            !=================================================================================
            ! 3 --> 1
            !=================================================================================

            Local_Node_Number_Edge_Line( 1, e )= Local_Node_Number_Different_Level( e )
            Local_Node_Number_Edge_Line( 2, e )= Local_Node_Number_Edge_Line_Plot( 2, 1, e )

            do i= 1, Number_Node_Difference_Level( e )

               call Compute_Position_Point_Plot_Result &
                  ( Value_Plotted, Value_Boundary_Level_Line, &
                    Local_Node_Number_Edge_Line, &
                    Position_Node_Plot_Mesh, Number_Node, Number_Node_Difference_Level( e ), &
                    Index_Element_2_Node, Number_Element, e, Max_Number_Level_Element, &
                    !======================================================================
                    Position_Node_Boundary_Level_Plot_1 )

            end do

            do i= 1, Number_Node_Difference_Level( e )
               do k= 1, 2
                  Difference_Position_1_Minus_NewPoint( k, e ) &
                  = Position_Node_Plot_Mesh( k, Index_Element_2_Node( Local_Node_Number_Edge_Line( 1, e ) ,e ) ) & 
                   -Position_Node_Boundary_Level_Plot_1( k, i, e )

                  Difference_Position_2_Minus_NewPoint( k, e ) &
                  = Position_Node_Plot_Mesh( k, Index_Element_2_Node( Local_Node_Number_Edge_Line( 2, e ) ,e ) ) & 
                   -Position_Node_Boundary_Level_Plot_1( k, i, e )
               end do

               if( Difference_Position_1_Minus_NewPoint( 1, e )*Difference_Position_2_Minus_NewPoint( 1, e ) & 
                   > Tolerance_Plot .or. &
                   Difference_Position_1_Minus_NewPoint( 2, e )*Difference_Position_2_Minus_NewPoint( 2, e ) &
                   > Tolerance_Plot )then    

                  write(287,*)'e=', e 
                  write(287,*)'Number_Node_Difference_Level( e )=', Number_Node_Difference_Level( e ) 
                  do k= 1, Number_Node_Difference_Level( e )
                     write(287,*) Position_Node_Boundary_Level_Plot_1( 1, k, e ), & 
                              Position_Node_Boundary_Level_Plot_1( 2, k, e ) 
                  end do
   
                  write(287,*)' ' 
                  write(*,*)'Plot_Contour_T_OC.f90 888' 
                  !stop
                 !go to 2272
               end if
            end do
    
            !=================================================================================
            ! 3 --> 2
            !=================================================================================

            Local_Node_Number_Edge_Line( 1, e )= Local_Node_Number_Different_Level( e )
            Local_Node_Number_Edge_Line( 2, e )= Local_Node_Number_Edge_Line_Plot( 2, 2, e )

            do i= 1, Number_Node_Difference_Level( e )

               call Compute_Position_Point_Plot_Result &
                  ( Value_Plotted, Value_Boundary_Level_Line, &
                    Local_Node_Number_Edge_Line, &
                    Position_Node_Plot_Mesh, Number_Node, Number_Node_Difference_Level( e ), &
                    Index_Element_2_Node, Number_Element, e, Max_Number_Level_Element, &
                    !======================================================================
                    Position_Node_Boundary_Level_Plot_2 )

            end do

            do i= 1, Number_Node_Difference_Level( e )
               do k= 1, 2
                  Difference_Position_1_Minus_NewPoint( k, e ) &
                  = Position_Node_Plot_Mesh( k, Index_Element_2_Node( Local_Node_Number_Edge_Line( 1, e ) ,e ) ) & 
                   -Position_Node_Boundary_Level_Plot_2( k, i, e )

                  Difference_Position_2_Minus_NewPoint( k, e ) &
                  = Position_Node_Plot_Mesh( k, Index_Element_2_Node( Local_Node_Number_Edge_Line( 2, e ) ,e ) ) & 
                   -Position_Node_Boundary_Level_Plot_2( k, i, e )
               end do

               if( Difference_Position_1_Minus_NewPoint( 1, e )*Difference_Position_2_Minus_NewPoint( 1, e ) &
                   > Tolerance_Plot  .or. &
                   Difference_Position_1_Minus_NewPoint( 2, e )*Difference_Position_2_Minus_NewPoint( 2, e ) &
                   > Tolerance_Plot )then    

                  write(288,*)'e=', e 
                  write(288,*)'Number_Node_Difference_Level( e )=', Number_Node_Difference_Level( e ) 
                  do k= 1, Number_Node_Difference_Level( e )
                     write(288,*) Position_Node_Boundary_Level_Plot_2( 1, k, e ), & 
                              Position_Node_Boundary_Level_Plot_2( 2, k, e ) 
                  end do
   
                  write(288,*)' ' 
                  write(*,*)'Plot_Contour_T_OC.f90 949' 
                  !stop
                 !go to 2272
               end if
            end do
    
            Flag_Plot_Element( e )= 3

            !call Plot_Element_InHomogeneous &
            !   ( Position_Node_Plot_Mesh, Index_Element_2_Node, Level_Node, RGB, & 
            !     File_Number, e, Number_Node, Number_Element, Number_Level, Max_Number_Level_Element, & 
            !     Number_Node_Difference_Level( e ), &
            !     Local_Node_Number_Different_Level( e ), Local_Node_Number_Edge_Line_Plot,  &
            !     Position_Node_Boundary_Level_Plot_1, Position_Node_Boundary_Level_Plot_2 )

! 
!            call Create_PlotData_Element_InHomogeneous &
!               ( Position_Node_Plot_Mesh, Index_Element_2_Node, Level_Node, & 
!                 e, Number_Node, Number_Element, & 
!                 Number_Node_Difference_Level( e ), & 
!                 Local_Node_Number_Different_Level( e ), Local_Node_Number_Edge_Line_Plot,  &
!                 Position_Node_Boundary_Level_Plot_1, Position_Node_Boundary_Level_Plot_2, &
!                 !=================================================================
!                 Number_Level_Element, Number_Node_Level_Element, Position_Node_Plot_Mesh_Element, Level_Element )
!
!            call Plot_Gradation_Element &
!               ( File_Number, e, Number_Element, Number_Level, & 
!                 Number_Level_Element, Number_Node_Level_Element, Position_Node_Plot_Mesh_Element, RGB, Level_Element )



         !=======================================================================================
         else if( Level_Node( Index_Element_2_Node( 3, e ) ) == Level_Node( Index_Element_2_Node( 1, e ) ) .and. & 
                Level_Node( Index_Element_2_Node( 1, e ) ) > Level_Node( Index_Element_2_Node( 2, e ) ) .and. & 
                Level_Node( Index_Element_2_Node( 3, e ) ) > Level_Node( Index_Element_2_Node( 2, e ) ) )then 
         !=======================================================================================

            Local_Node_Number_Different_Level( e )= 2 


            Local_Node_Number_Edge_Line_Plot( 1, 1, e )= Local_Node_Number_Different_Level( e ) 
            Local_Node_Number_Edge_Line_Plot( 2, 1, e )= 3
            Local_Node_Number_Edge_Line_Plot( 1, 2, e )= Local_Node_Number_Different_Level( e )
            Local_Node_Number_Edge_Line_Plot( 2, 2, e )= 1

            Level_High( e )= Level_Node( Index_Element_2_Node( Local_Node_Number_Edge_Line_Plot( 2, 1, e ), e ) ) 
            Level_Low( e )=  Level_Node( Index_Element_2_Node( Local_Node_Number_Different_Level( e ), e ) )

            Number_Node_Difference_Level( e )= Level_High( e ) -Level_Low( e ) 

            if( Number_Node_Difference_Level( e )<=0 )then
               write(*,*)'Number_Node_Difference_Level( e )==0'
               !stop
              !go to 2272
            end if


            do i= 1, Number_Node_Difference_Level( e )
               Value_Boundary_Level_Line( i, e )& 
               = Value_Boundary_Level&
                 ( 1, i +Level_Node( Index_Element_2_Node( Local_Node_Number_Different_Level( e ), e ) ) )  
            end do

            if( Level_Node( Index_Element_2_Node( Local_Node_Number_Different_Level( e ), e ) ) &
                +Number_Node_Difference_Level( e ) & 
                /= Level_Node( Index_Element_2_Node( Local_Node_Number_Edge_Line_Plot( 2, 1, e ), e ) ) .or. &
               Level_Node( Index_Element_2_Node( Local_Node_Number_Different_Level( e ), e ) ) & 
                +Number_Node_Difference_Level( e ) & 
                /= Level_Node( Index_Element_2_Node( Local_Node_Number_Edge_Line_Plot( 2, 2, e ), e ) ) )then
 
               write(*,*)'Level_Node( Index_Element_2_Node( Local_Node_Number_Different_Level( e ), e ) )=', & 
               Level_Node( Index_Element_2_Node( Local_Node_Number_Different_Level( e ), e ) )  
               write(*,*)'Number_Node_Difference_Level( e )=', Number_Node_Difference_Level( e )
               write(*,*)'Level_Node( Index_Element_2_Node( Local_Node_Number_Edge_Line_Plot( 2, 1, e ), e ) )=', & 
                      Level_Node( Index_Element_2_Node( Local_Node_Number_Edge_Line_Plot( 2, 1, e ), e ) ) 
               write(*,*)'Level_Node( Index_Element_2_Node( Local_Node_Number_Edge_Line_Plot( 2, 2, e ), e ) )=', &
                      Level_Node( Index_Element_2_Node( Local_Node_Number_Edge_Line_Plot( 2, 2, e ), e ) ) 
               write(*,*)'Plot_Contour_T_OC.f90 826'
               !stop
              !go to 2272
            end if

            !=================================================================================
            ! 3 --> 1
            !=================================================================================

            Local_Node_Number_Edge_Line( 1, e )= Local_Node_Number_Different_Level( e )
            Local_Node_Number_Edge_Line( 2, e )= Local_Node_Number_Edge_Line_Plot( 2, 1, e )

            do i= 1, Number_Node_Difference_Level( e )

               call Compute_Position_Point_Plot_Result &
                  ( Value_Plotted, Value_Boundary_Level_Line, &
                    Local_Node_Number_Edge_Line, &
                    Position_Node_Plot_Mesh, Number_Node, Number_Node_Difference_Level( e ), &
                    Index_Element_2_Node, Number_Element, e, Max_Number_Level_Element, &
                    !======================================================================
                    Position_Node_Boundary_Level_Plot_1 )

            end do

            do i= 1, Number_Node_Difference_Level( e )
               do k= 1, 2
                  Difference_Position_1_Minus_NewPoint( k, e ) &
                  = Position_Node_Plot_Mesh( k, Index_Element_2_Node( Local_Node_Number_Edge_Line( 1, e ) ,e ) ) & 
                   -Position_Node_Boundary_Level_Plot_1( k, i, e )

                  Difference_Position_2_Minus_NewPoint( k, e ) &
                  = Position_Node_Plot_Mesh( k, Index_Element_2_Node( Local_Node_Number_Edge_Line( 2, e ) ,e ) ) & 
                   -Position_Node_Boundary_Level_Plot_1( k, i, e )
               end do

               if( Difference_Position_1_Minus_NewPoint( 1, e )*Difference_Position_2_Minus_NewPoint( 1, e ) &
                   > Tolerance_Plot .or. &
                   Difference_Position_1_Minus_NewPoint( 2, e )*Difference_Position_2_Minus_NewPoint( 2, e ) & 
                   > Tolerance_Plot )then    

                  write(287,*)'e=', e 
                  write(287,*)'Number_Node_Difference_Level( e )=', Number_Node_Difference_Level( e ) 
                  do k= 1, Number_Node_Difference_Level( e )
                     write(287,*) Position_Node_Boundary_Level_Plot_1( 1, k, e ), & 
                              Position_Node_Boundary_Level_Plot_1( 2, k, e ) 
                  end do
   
                  write(287,*)' ' 
                  write(*,*)'Plot_Contour_T_OC.f90 1073' 
                  !stop
                 !go to 2272
               end if
            end do
    
            !=================================================================================
            ! 3 --> 2
            !=================================================================================

            Local_Node_Number_Edge_Line( 1, e )= Local_Node_Number_Different_Level( e )
            Local_Node_Number_Edge_Line( 2, e )= Local_Node_Number_Edge_Line_Plot( 2, 2, e )

            do i= 1, Number_Node_Difference_Level( e )

               call Compute_Position_Point_Plot_Result &
                  ( Value_Plotted, Value_Boundary_Level_Line, &
                    Local_Node_Number_Edge_Line, &
                    Position_Node_Plot_Mesh, Number_Node, Number_Node_Difference_Level( e ), &
                    Index_Element_2_Node, Number_Element, e, Max_Number_Level_Element, &
                    !======================================================================
                    Position_Node_Boundary_Level_Plot_2 )

            end do

            do i= 1, Number_Node_Difference_Level( e )
               do k= 1, 2
                  Difference_Position_1_Minus_NewPoint( k, e ) &
                  = Position_Node_Plot_Mesh( k, Index_Element_2_Node( Local_Node_Number_Edge_Line( 1, e ) ,e ) ) & 
                   -Position_Node_Boundary_Level_Plot_2( k, i, e )

                  Difference_Position_2_Minus_NewPoint( k, e ) &
                  = Position_Node_Plot_Mesh( k, Index_Element_2_Node( Local_Node_Number_Edge_Line( 2, e ) ,e ) ) & 
                   -Position_Node_Boundary_Level_Plot_2( k, i, e )
               end do

               if( Difference_Position_1_Minus_NewPoint( 1, e )*Difference_Position_2_Minus_NewPoint( 1, e ) &
                   > Tolerance_Plot .or. &
                   Difference_Position_1_Minus_NewPoint( 2, e )*Difference_Position_2_Minus_NewPoint( 2, e ) &
                   > Tolerance_Plot )then    

                  write(288,*)'e=', e 
                  write(288,*)'Number_Node_Difference_Level( e )=', Number_Node_Difference_Level( e ) 
                  do k= 1, Number_Node_Difference_Level( e )
                     write(288,*) Position_Node_Boundary_Level_Plot_2( 1, k, e ), & 
                              Position_Node_Boundary_Level_Plot_2( 2, k, e ) 
                  end do
   
                  write(288,*)' ' 
                  write(*,*)'Plot_Contour_T_OC.f90 1134' 
                  !stop
                 !go to 2272
               end if
            end do
    
            Flag_Plot_Element( e )= 3

            !call Plot_Element_InHomogeneous &
            !   ( Position_Node_Plot_Mesh, Index_Element_2_Node, Level_Node, RGB, & 
            !     File_Number, e, Number_Node, Number_Element, Number_Level, Max_Number_Level_Element, & 
            !     Number_Node_Difference_Level( e ), &
            !     Local_Node_Number_Different_Level( e ), Local_Node_Number_Edge_Line_Plot,  &
            !     Position_Node_Boundary_Level_Plot_1, Position_Node_Boundary_Level_Plot_2 )
 
            !call Create_PlotData_Element_InHomogeneous &
            !   ( Position_Node_Plot_Mesh, Index_Element_2_Node, Level_Node, & 
            !     e, Number_Node, Number_Element, & 
            !     Number_Node_Difference_Level( e ), & 
            !     Local_Node_Number_Different_Level( e ), Local_Node_Number_Edge_Line_Plot,  &
            !     Position_Node_Boundary_Level_Plot_1, Position_Node_Boundary_Level_Plot_2, &
            !     !=================================================================
            !     Number_Level_Element, Number_Node_Level_Element, Position_Node_Plot_Mesh_Element, Level_Element )

            !call Plot_Gradation_Element &
            !   ( File_Number, e, Number_Element, Number_Level, & 
            !     Number_Level_Element, Number_Node_Level_Element, Position_Node_Plot_Mesh_Element, RGB, Level_Element )

         !=======================================================================================
         else if( Level_Node( Index_Element_2_Node( 3, e ) ) == Level_Node( Index_Element_2_Node( 1, e ) ) .and. & 
                Level_Node( Index_Element_2_Node( 1, e ) ) < Level_Node( Index_Element_2_Node( 2, e ) ) .and. & 
                Level_Node( Index_Element_2_Node( 3, e ) ) < Level_Node( Index_Element_2_Node( 2, e ) ) )then 
         !=======================================================================================

            Local_Node_Number_Different_Level( e )= 2 


            Local_Node_Number_Edge_Line_Plot( 1, 1, e )= Local_Node_Number_Different_Level( e ) 
            Local_Node_Number_Edge_Line_Plot( 2, 1, e )= 1
            Local_Node_Number_Edge_Line_Plot( 1, 2, e )= Local_Node_Number_Different_Level( e )
            Local_Node_Number_Edge_Line_Plot( 2, 2, e )= 3

            Level_High( e )=  Level_Node( Index_Element_2_Node( Local_Node_Number_Different_Level( e ), e ) )
            Level_Low( e )= Level_Node( Index_Element_2_Node( Local_Node_Number_Edge_Line_Plot( 2, 1, e ), e ) ) 

            Number_Node_Difference_Level( e )= Level_High( e ) -Level_Low( e ) 

            if( Number_Node_Difference_Level( e )<=0 )then
               write(*,*)'Number_Node_Difference_Level( e )==0'
               !stop
              !go to 2272
            end if


            do i= 1, Number_Node_Difference_Level( e )
               Value_Boundary_Level_Line( i, e )& 
               = Value_Boundary_Level&
                ( 2, -i +Level_Node( Index_Element_2_Node( Local_Node_Number_Different_Level( e ), e ) ) )  
            end do

            if( Level_Node( Index_Element_2_Node( Local_Node_Number_Different_Level( e ), e ) ) & 
                -Number_Node_Difference_Level( e ) & 
                /= Level_Node( Index_Element_2_Node( Local_Node_Number_Edge_Line_Plot( 2, 1, e ), e ) ) .or. &
               Level_Node( Index_Element_2_Node( Local_Node_Number_Different_Level( e ), e ) ) & 
                -Number_Node_Difference_Level( e ) & 
                /= Level_Node( Index_Element_2_Node( Local_Node_Number_Edge_Line_Plot( 2, 2, e ), e ) ) )then
 
               write(*,*)'Level_Node( Index_Element_2_Node( Local_Node_Number_Different_Level( e ), e ) )=', & 
               Level_Node( Index_Element_2_Node( Local_Node_Number_Different_Level( e ), e ) )  
               write(*,*)'Number_Node_Difference_Level( e )=', Number_Node_Difference_Level( e )
               write(*,*)'Level_Node( Index_Element_2_Node( Local_Node_Number_Edge_Line_Plot( 2, 1, e ), e ) )=', & 
                      Level_Node( Index_Element_2_Node( Local_Node_Number_Edge_Line_Plot( 2, 1, e ), e ) ) 
               write(*,*)'Level_Node( Index_Element_2_Node( Local_Node_Number_Edge_Line_Plot( 2, 2, e ), e ) )=', &
                      Level_Node( Index_Element_2_Node( Local_Node_Number_Edge_Line_Plot( 2, 2, e ), e ) ) 
               write(*,*)'Plot_Contour_T_OC.f90 448'
               !stop
              !go to 2272
            end if

            !=================================================================================
            ! 3 --> 1
            !=================================================================================

            Local_Node_Number_Edge_Line( 1, e )= Local_Node_Number_Different_Level( e )
            Local_Node_Number_Edge_Line( 2, e )= Local_Node_Number_Edge_Line_Plot( 2, 1, e )

            do i= 1, Number_Node_Difference_Level( e )

               call Compute_Position_Point_Plot_Result &
                  ( Value_Plotted, Value_Boundary_Level_Line, &
                    Local_Node_Number_Edge_Line, &
                    Position_Node_Plot_Mesh, Number_Node, Number_Node_Difference_Level( e ), &
                    Index_Element_2_Node, Number_Element, e, Max_Number_Level_Element, &
                    !======================================================================
                    Position_Node_Boundary_Level_Plot_1 )

            end do

            do i= 1, Number_Node_Difference_Level( e )
               do k= 1, 2
                  Difference_Position_1_Minus_NewPoint( k, e ) &
                  = Position_Node_Plot_Mesh( k, Index_Element_2_Node( Local_Node_Number_Edge_Line( 1, e ) ,e ) ) & 
                   -Position_Node_Boundary_Level_Plot_1( k, i, e )

                  Difference_Position_2_Minus_NewPoint( k, e ) &
                  = Position_Node_Plot_Mesh( k, Index_Element_2_Node( Local_Node_Number_Edge_Line( 2, e ) ,e ) ) & 
                   -Position_Node_Boundary_Level_Plot_1( k, i, e )
               end do

               if( Difference_Position_1_Minus_NewPoint( 1, e )*Difference_Position_2_Minus_NewPoint( 1, e ) &
                   > Tolerance_Plot .or. &
                   Difference_Position_1_Minus_NewPoint( 2, e )*Difference_Position_2_Minus_NewPoint( 2, e ) &
                   > Tolerance_Plot )then    

                  write(287,*)'e=', e 
                  write(287,*)'Number_Node_Difference_Level( e )=', Number_Node_Difference_Level( e ) 
                  do k= 1, Number_Node_Difference_Level( e )
                     write(287,*) Position_Node_Boundary_Level_Plot_1( 1, k, e ), & 
                              Position_Node_Boundary_Level_Plot_1( 2, k, e ) 
                  end do
   
                  write(287,*)' ' 
                  write(*,*)'Plot_Contour_T_OC.f90 1259' 
                  !stop
                 !go to 2272
               end if
            end do
    
            !=================================================================================
            ! 3 --> 2
            !=================================================================================

            Local_Node_Number_Edge_Line( 1, e )= Local_Node_Number_Different_Level( e )
            Local_Node_Number_Edge_Line( 2, e )= Local_Node_Number_Edge_Line_Plot( 2, 2, e )

            do i= 1, Number_Node_Difference_Level( e )

               call Compute_Position_Point_Plot_Result &
                  ( Value_Plotted, Value_Boundary_Level_Line, &
                    Local_Node_Number_Edge_Line, &
                    Position_Node_Plot_Mesh, Number_Node, Number_Node_Difference_Level( e ), &
                    Index_Element_2_Node, Number_Element, e, Max_Number_Level_Element, &
                    !======================================================================
                    Position_Node_Boundary_Level_Plot_2 )

            end do

            do i= 1, Number_Node_Difference_Level( e )
               do k= 1, 2
                  Difference_Position_1_Minus_NewPoint( k, e ) &
                  = Position_Node_Plot_Mesh( k, Index_Element_2_Node( Local_Node_Number_Edge_Line( 1, e ) ,e ) ) & 
                   -Position_Node_Boundary_Level_Plot_2( k, i, e )

                  Difference_Position_2_Minus_NewPoint( k, e ) &
                  = Position_Node_Plot_Mesh( k, Index_Element_2_Node( Local_Node_Number_Edge_Line( 2, e ) ,e ) ) & 
                   -Position_Node_Boundary_Level_Plot_2( k, i, e )
               end do

               if( Difference_Position_1_Minus_NewPoint( 1, e )*Difference_Position_2_Minus_NewPoint( 1, e ) &
                   > Tolerance_Plot .or. &
                   Difference_Position_1_Minus_NewPoint( 2, e )*Difference_Position_2_Minus_NewPoint( 2, e ) &
                   > Tolerance_Plot )then    

                  write(288,*)'e=', e 
                  write(288,*)'Number_Node_Difference_Level( e )=', Number_Node_Difference_Level( e ) 
                  do k= 1, Number_Node_Difference_Level( e )
                     write(288,*) Position_Node_Boundary_Level_Plot_2( 1, k, e ), & 
                              Position_Node_Boundary_Level_Plot_2( 2, k, e ) 
                  end do
   
                  write(288,*)' ' 
                  write(*,*)'Plot_Contour_T_OC.f90 1320' 
                  !stop
                 !go to 2272
               end if
            end do
    
            Flag_Plot_Element( e )= 3

            !call Plot_Element_InHomogeneous &
            !   ( Position_Node_Plot_Mesh, Index_Element_2_Node, Level_Node, RGB, & 
            !     File_Number, e, Number_Node, Number_Element, Number_Level, Max_Number_Level_Element, & 
            !     Number_Node_Difference_Level( e ), &
            !     Local_Node_Number_Different_Level( e ), Local_Node_Number_Edge_Line_Plot,  &
            !     Position_Node_Boundary_Level_Plot_1, Position_Node_Boundary_Level_Plot_2 )

! 
!            call Create_PlotData_Element_InHomogeneous &
!               ( Position_Node_Plot_Mesh, Index_Element_2_Node, Level_Node, & 
!                 e, Number_Node, Number_Element, & 
!                 Number_Node_Difference_Level( e ), & 
!                 Local_Node_Number_Different_Level( e ), Local_Node_Number_Edge_Line_Plot,  &
!                 Position_Node_Boundary_Level_Plot_1, Position_Node_Boundary_Level_Plot_2, &
!                 !=================================================================
!                 Number_Level_Element, Number_Node_Level_Element, Position_Node_Plot_Mesh_Element, Level_Element )
!
!            call Plot_Gradation_Element &
!               ( File_Number, e, Number_Element, Number_Level, & 
!                 Number_Level_Element, Number_Node_Level_Element, Position_Node_Plot_Mesh_Element, RGB, Level_Element )



         !=======================================================================================
         else if( Level_Node( Index_Element_2_Node( 1, e ) ) > Level_Node( Index_Element_2_Node( 2, e ) ) .and. & 
                Level_Node( Index_Element_2_Node( 1, e ) ) > Level_Node( Index_Element_2_Node( 3, e ) ) .and. & 
                Level_Node( Index_Element_2_Node( 2, e ) ) > Level_Node( Index_Element_2_Node( 3, e ) ) )then 
         !=======================================================================================

            Local_Node_Number_Level_Max_2_Min( 1, e )= 1 
            Local_Node_Number_Level_Max_2_Min( 2, e )= 2
            Local_Node_Number_Level_Max_2_Min( 3, e )= 3

            Number_Node_Difference_Level_3_Node( 1, e ) &
            = Level_Node( Index_Element_2_Node( Local_Node_Number_Level_Max_2_Min( 1, e ), e ) ) & 
             -Level_Node( Index_Element_2_Node( Local_Node_Number_Level_Max_2_Min( 2, e ), e ) ) 

            Number_Node_Difference_Level_3_Node( 2, e ) &
            = Level_Node( Index_Element_2_Node( Local_Node_Number_Level_Max_2_Min( 1, e ), e ) ) & 
             -Level_Node( Index_Element_2_Node( Local_Node_Number_Level_Max_2_Min( 3, e ), e ) ) 

            Number_Node_Difference_Level_3_Node( 3, e ) &
            = Level_Node( Index_Element_2_Node( Local_Node_Number_Level_Max_2_Min( 2, e ), e ) ) & 
             -Level_Node( Index_Element_2_Node( Local_Node_Number_Level_Max_2_Min( 3, e ), e ) ) 

            do i= 1, 3
               do k= 1, Number_Node_Difference_Level_3_Node( 2, e )
                  do l= 1, 2 
                      Position_Node_Boundary_Level_Plot( l, k, i, e )= 0.0d0
                  end do
               end do
            end do

            !=====================================================================================
            ! Compute Position Boundary
            !=====================================================================================
            do i= 1, 3

               if( i==3 )then
                  do k= 1, Number_Node_Difference_Level_3_Node( i, e )
                     Value_Boundary_Level_Line( k, e ) &
                     = Value_Boundary_Level&
                       ( 2, -k +Level_Node( Index_Element_2_Node( Local_Node_Number_Level_Max_2_Min( 2, e ), e ) ) )  
                  end do
                  Local_Node_Number_Edge_Line( 1, e )= Local_Node_Number_Level_Max_2_Min( 2, e )
                  Local_Node_Number_Edge_Line( 2, e )= Local_Node_Number_Level_Max_2_Min( 3, e ) 
               else
                  do k= 1, Number_Node_Difference_Level_3_Node( i, e )
                     Value_Boundary_Level_Line( k, e ) &
                     = Value_Boundary_Level&
                       ( 2, -k +Level_Node( Index_Element_2_Node( Local_Node_Number_Level_Max_2_Min( 1, e ), e ) ) )  
                  end do
                  Local_Node_Number_Edge_Line( 1, e )= Local_Node_Number_Level_Max_2_Min( 1, e )

                  if( i==1 )then
                     Local_Node_Number_Edge_Line( 2, e )= Local_Node_Number_Level_Max_2_Min( 2, e ) 
                  else if( i==2 )then
                     Local_Node_Number_Edge_Line( 2, e )= Local_Node_Number_Level_Max_2_Min( 3, e )
                  end if 
               end if


               call Compute_Position_Point_Plot_Result &
                  ( Value_Plotted, Value_Boundary_Level_Line, &
                    Local_Node_Number_Edge_Line, &
                    Position_Node_Plot_Mesh, Number_Node, Number_Node_Difference_Level_3_Node( i, e ),  &
                    Index_Element_2_Node, Number_Element, e, Max_Number_Level_Element, &
                    !======================================================================
                    Position_Node_Boundary_Level_Plot_tmp )

               do k= 1, Number_Node_Difference_Level_3_Node( i, e )
                  do l= 1, 2
                     Position_Node_Boundary_Level_Plot( l, k, i, e ) & 
                     = Position_Node_Boundary_Level_Plot_tmp( l, k, e )
                  end do
               end do

            end do

            Flag_Plot_Element( e )= 4

            !call Plot_Element_InHomogeneous_3_Node &
            !   ( Position_Node_Plot_Mesh, Index_Element_2_Node, Level_Node, RGB, & 
            !     File_Number, e, Number_Node, Number_Element, Number_Level, Max_Number_Level_Element, & 
            !     Number_Node_Difference_Level_3_Node, Local_Node_Number_Level_Max_2_Min, &
            !     Position_Node_Boundary_Level_Plot )

         !=======================================================================================
         else if( Level_Node( Index_Element_2_Node( 1, e ) ) > Level_Node( Index_Element_2_Node( 2, e ) ) .and. & 
                Level_Node( Index_Element_2_Node( 1, e ) ) > Level_Node( Index_Element_2_Node( 3, e ) ) .and. & 
                Level_Node( Index_Element_2_Node( 3, e ) ) > Level_Node( Index_Element_2_Node( 2, e ) ) )then 
         !=======================================================================================

            Local_Node_Number_Level_Max_2_Min( 1, e )= 1 
            Local_Node_Number_Level_Max_2_Min( 2, e )= 3
            Local_Node_Number_Level_Max_2_Min( 3, e )= 2

            Number_Node_Difference_Level_3_Node( 1, e ) &
            = Level_Node( Index_Element_2_Node( Local_Node_Number_Level_Max_2_Min( 1, e ), e ) ) & 
             -Level_Node( Index_Element_2_Node( Local_Node_Number_Level_Max_2_Min( 2, e ), e ) ) 

            Number_Node_Difference_Level_3_Node( 2, e ) &
            = Level_Node( Index_Element_2_Node( Local_Node_Number_Level_Max_2_Min( 1, e ), e ) ) & 
             -Level_Node( Index_Element_2_Node( Local_Node_Number_Level_Max_2_Min( 3, e ), e ) ) 

            Number_Node_Difference_Level_3_Node( 3, e ) &
            = Level_Node( Index_Element_2_Node( Local_Node_Number_Level_Max_2_Min( 2, e ), e ) ) & 
             -Level_Node( Index_Element_2_Node( Local_Node_Number_Level_Max_2_Min( 3, e ), e ) ) 

            do i= 1, 3
               do k= 1, Number_Node_Difference_Level_3_Node( 2, e )
                  do l= 1, 2 
                      Position_Node_Boundary_Level_Plot( l, k, i, e )= 0.0d0
                  end do
               end do
            end do

            !=====================================================================================
            ! Compute Position Boundary
            !=====================================================================================
            do i= 1, 3

               if( i==3 )then
                  do k= 1, Number_Node_Difference_Level_3_Node( i, e )
                     Value_Boundary_Level_Line( k, e ) &
                     = Value_Boundary_Level&
                       ( 2, -k +Level_Node( Index_Element_2_Node( Local_Node_Number_Level_Max_2_Min( 2, e ), e ) ) )  
                  end do
                  Local_Node_Number_Edge_Line( 1, e )= Local_Node_Number_Level_Max_2_Min( 2, e )
                  Local_Node_Number_Edge_Line( 2, e )= Local_Node_Number_Level_Max_2_Min( 3, e ) 
               else
                  do k= 1, Number_Node_Difference_Level_3_Node( i, e )
                     Value_Boundary_Level_Line( k, e ) &
                     = Value_Boundary_Level&
                       ( 2, -k +Level_Node( Index_Element_2_Node( Local_Node_Number_Level_Max_2_Min( 1, e ), e ) ) )  
                  end do
                  Local_Node_Number_Edge_Line( 1, e )= Local_Node_Number_Level_Max_2_Min( 1, e )

                  if( i==1 )then
                     Local_Node_Number_Edge_Line( 2, e )= Local_Node_Number_Level_Max_2_Min( 2, e ) 
                  else if( i==2 )then
                     Local_Node_Number_Edge_Line( 2, e )= Local_Node_Number_Level_Max_2_Min( 3, e )
                  end if 
               end if


               call Compute_Position_Point_Plot_Result &
                  ( Value_Plotted, Value_Boundary_Level_Line, &
                    Local_Node_Number_Edge_Line, &
                    Position_Node_Plot_Mesh, Number_Node, Number_Node_Difference_Level_3_Node( i, e ),  &
                    Index_Element_2_Node, Number_Element, e, Max_Number_Level_Element, &
                    !======================================================================
                    Position_Node_Boundary_Level_Plot_tmp )

               do k= 1, Number_Node_Difference_Level_3_Node( i, e )
                  do l= 1, 2
                     Position_Node_Boundary_Level_Plot( l, k, i, e ) & 
                     = Position_Node_Boundary_Level_Plot_tmp( l, k, e )
                  end do
               end do

            end do

            Flag_Plot_Element( e )= 4

            !call Plot_Element_InHomogeneous_3_Node &
            !   ( Position_Node_Plot_Mesh, Index_Element_2_Node, Level_Node, RGB, & 
            !     File_Number, e, Number_Node, Number_Element, Number_Level, Max_Number_Level_Element, & 
            !     Number_Node_Difference_Level_3_Node, Local_Node_Number_Level_Max_2_Min, &
            !     Position_Node_Boundary_Level_Plot )


         !=======================================================================================
         else if( Level_Node( Index_Element_2_Node( 2, e ) ) > Level_Node( Index_Element_2_Node( 1, e ) ) .and. & 
                Level_Node( Index_Element_2_Node( 2, e ) ) > Level_Node( Index_Element_2_Node( 3, e ) ) .and. & 
                Level_Node( Index_Element_2_Node( 1, e ) ) > Level_Node( Index_Element_2_Node( 3, e ) ) )then 
         !=======================================================================================

            Local_Node_Number_Level_Max_2_Min( 1, e )= 2 
            Local_Node_Number_Level_Max_2_Min( 2, e )= 1
            Local_Node_Number_Level_Max_2_Min( 3, e )= 3

            Number_Node_Difference_Level_3_Node( 1, e ) &
            = Level_Node( Index_Element_2_Node( Local_Node_Number_Level_Max_2_Min( 1, e ), e ) ) & 
             -Level_Node( Index_Element_2_Node( Local_Node_Number_Level_Max_2_Min( 2, e ), e ) ) 

            Number_Node_Difference_Level_3_Node( 2, e ) &
            = Level_Node( Index_Element_2_Node( Local_Node_Number_Level_Max_2_Min( 1, e ), e ) ) & 
             -Level_Node( Index_Element_2_Node( Local_Node_Number_Level_Max_2_Min( 3, e ), e ) ) 

            Number_Node_Difference_Level_3_Node( 3, e ) &
            = Level_Node( Index_Element_2_Node( Local_Node_Number_Level_Max_2_Min( 2, e ), e ) ) & 
             -Level_Node( Index_Element_2_Node( Local_Node_Number_Level_Max_2_Min( 3, e ), e ) ) 

            do i= 1, 3
               do k= 1, Number_Node_Difference_Level_3_Node( 2, e )
                  do l= 1, 2 
                      Position_Node_Boundary_Level_Plot( l, k, i, e )= 0.0d0
                  end do
               end do
            end do

            !=====================================================================================
            ! Compute Position Boundary
            !=====================================================================================
            do i= 1, 3

               if( i==3 )then
                  do k= 1, Number_Node_Difference_Level_3_Node( i, e )
                     Value_Boundary_Level_Line( k, e ) &
                     = Value_Boundary_Level&
                       ( 2, -k +Level_Node( Index_Element_2_Node( Local_Node_Number_Level_Max_2_Min( 2, e ), e ) ) )  
                  end do
                  Local_Node_Number_Edge_Line( 1, e )= Local_Node_Number_Level_Max_2_Min( 2, e )
                  Local_Node_Number_Edge_Line( 2, e )= Local_Node_Number_Level_Max_2_Min( 3, e ) 
               else
                  do k= 1, Number_Node_Difference_Level_3_Node( i, e )
                     Value_Boundary_Level_Line( k, e ) &
                     = Value_Boundary_Level&
                       ( 2, -k +Level_Node( Index_Element_2_Node( Local_Node_Number_Level_Max_2_Min( 1, e ), e ) ) )  
                  end do
                  Local_Node_Number_Edge_Line( 1, e )= Local_Node_Number_Level_Max_2_Min( 1, e )

                  if( i==1 )then
                     Local_Node_Number_Edge_Line( 2, e )= Local_Node_Number_Level_Max_2_Min( 2, e ) 
                  else if( i==2 )then
                     Local_Node_Number_Edge_Line( 2, e )= Local_Node_Number_Level_Max_2_Min( 3, e )
                  end if 
               end if


               call Compute_Position_Point_Plot_Result &
                  ( Value_Plotted, Value_Boundary_Level_Line, &
                    Local_Node_Number_Edge_Line, &
                    Position_Node_Plot_Mesh, Number_Node, Number_Node_Difference_Level_3_Node( i, e ),  &
                    Index_Element_2_Node, Number_Element, e, Max_Number_Level_Element, &
                    !======================================================================
                    Position_Node_Boundary_Level_Plot_tmp )

               do k= 1, Number_Node_Difference_Level_3_Node( i, e )
                  do l= 1, 2
                     Position_Node_Boundary_Level_Plot( l, k, i, e ) & 
                     = Position_Node_Boundary_Level_Plot_tmp( l, k, e )
                  end do
               end do

            end do

            Flag_Plot_Element( e )= 4

            !call Plot_Element_InHomogeneous_3_Node &
            !   ( Position_Node_Plot_Mesh, Index_Element_2_Node, Level_Node, RGB, & 
            !     File_Number, e, Number_Node, Number_Element, Number_Level, Max_Number_Level_Element, & 
            !     Number_Node_Difference_Level_3_Node, Local_Node_Number_Level_Max_2_Min, &
            !     Position_Node_Boundary_Level_Plot )

         !=======================================================================================
         else if( Level_Node( Index_Element_2_Node( 2, e ) ) > Level_Node( Index_Element_2_Node( 1, e ) ) .and. & 
                Level_Node( Index_Element_2_Node( 2, e ) ) > Level_Node( Index_Element_2_Node( 3, e ) ) .and. & 
                Level_Node( Index_Element_2_Node( 3, e ) ) > Level_Node( Index_Element_2_Node( 1, e ) ) )then 
         !=======================================================================================

            Local_Node_Number_Level_Max_2_Min( 1, e )= 2 
            Local_Node_Number_Level_Max_2_Min( 2, e )= 3
            Local_Node_Number_Level_Max_2_Min( 3, e )= 1

            Number_Node_Difference_Level_3_Node( 1, e ) &
            = Level_Node( Index_Element_2_Node( Local_Node_Number_Level_Max_2_Min( 1, e ), e ) ) & 
             -Level_Node( Index_Element_2_Node( Local_Node_Number_Level_Max_2_Min( 2, e ), e ) ) 

            Number_Node_Difference_Level_3_Node( 2, e ) &
            = Level_Node( Index_Element_2_Node( Local_Node_Number_Level_Max_2_Min( 1, e ), e ) ) & 
             -Level_Node( Index_Element_2_Node( Local_Node_Number_Level_Max_2_Min( 3, e ), e ) ) 

            Number_Node_Difference_Level_3_Node( 3, e ) &
            = Level_Node( Index_Element_2_Node( Local_Node_Number_Level_Max_2_Min( 2, e ), e ) ) & 
             -Level_Node( Index_Element_2_Node( Local_Node_Number_Level_Max_2_Min( 3, e ), e ) ) 

            do i= 1, 3
               do k= 1, Number_Node_Difference_Level_3_Node( 2, e )
                  do l= 1, 2 
                      Position_Node_Boundary_Level_Plot( l, k, i, e )= 0.0d0
                  end do
               end do
            end do

            !=====================================================================================
            ! Compute Position Boundary
            !=====================================================================================
            do i= 1, 3

               if( i==3 )then
                  do k= 1, Number_Node_Difference_Level_3_Node( i, e )
                     Value_Boundary_Level_Line( k, e ) &
                     = Value_Boundary_Level&
                       ( 2, -k +Level_Node( Index_Element_2_Node( Local_Node_Number_Level_Max_2_Min( 2, e ), e ) ) )  
                  end do
                  Local_Node_Number_Edge_Line( 1, e )= Local_Node_Number_Level_Max_2_Min( 2, e )
                  Local_Node_Number_Edge_Line( 2, e )= Local_Node_Number_Level_Max_2_Min( 3, e ) 
               else
                  do k= 1, Number_Node_Difference_Level_3_Node( i, e )
                     Value_Boundary_Level_Line( k, e ) &
                     = Value_Boundary_Level&
                       ( 2, -k +Level_Node( Index_Element_2_Node( Local_Node_Number_Level_Max_2_Min( 1, e ), e ) ) )  
                  end do
                  Local_Node_Number_Edge_Line( 1, e )= Local_Node_Number_Level_Max_2_Min( 1, e )

                  if( i==1 )then
                     Local_Node_Number_Edge_Line( 2, e )= Local_Node_Number_Level_Max_2_Min( 2, e ) 
                  else if( i==2 )then
                     Local_Node_Number_Edge_Line( 2, e )= Local_Node_Number_Level_Max_2_Min( 3, e )
                  end if 
               end if


               call Compute_Position_Point_Plot_Result &
                  ( Value_Plotted, Value_Boundary_Level_Line, &
                    Local_Node_Number_Edge_Line, &
                    Position_Node_Plot_Mesh, Number_Node, Number_Node_Difference_Level_3_Node( i, e ),  &
                    Index_Element_2_Node, Number_Element, e, Max_Number_Level_Element, &
                    !======================================================================
                    Position_Node_Boundary_Level_Plot_tmp )

               do k= 1, Number_Node_Difference_Level_3_Node( i, e )
                  do l= 1, 2
                     Position_Node_Boundary_Level_Plot( l, k, i, e ) & 
                     = Position_Node_Boundary_Level_Plot_tmp( l, k, e )
                  end do
               end do

            end do

            Flag_Plot_Element( e )= 4

            !call Plot_Element_InHomogeneous_3_Node &
            !   ( Position_Node_Plot_Mesh, Index_Element_2_Node, Level_Node, RGB, & 
            !     File_Number, e, Number_Node, Number_Element, Number_Level, Max_Number_Level_Element, & 
            !     Number_Node_Difference_Level_3_Node, Local_Node_Number_Level_Max_2_Min, &
            !     Position_Node_Boundary_Level_Plot )

         !=======================================================================================
         else if( Level_Node( Index_Element_2_Node( 3, e ) ) > Level_Node( Index_Element_2_Node( 1, e ) ) .and. & 
                Level_Node( Index_Element_2_Node( 3, e ) ) > Level_Node( Index_Element_2_Node( 2, e ) ) .and. & 
                Level_Node( Index_Element_2_Node( 1, e ) ) > Level_Node( Index_Element_2_Node( 2, e ) ) )then 
         !=======================================================================================

            Local_Node_Number_Level_Max_2_Min( 1, e )= 3 
            Local_Node_Number_Level_Max_2_Min( 2, e )= 1
            Local_Node_Number_Level_Max_2_Min( 3, e )= 2


            Number_Node_Difference_Level_3_Node( 1, e ) &
            = Level_Node( Index_Element_2_Node( Local_Node_Number_Level_Max_2_Min( 1, e ), e ) ) & 
             -Level_Node( Index_Element_2_Node( Local_Node_Number_Level_Max_2_Min( 2, e ), e ) ) 

            Number_Node_Difference_Level_3_Node( 2, e ) &
            = Level_Node( Index_Element_2_Node( Local_Node_Number_Level_Max_2_Min( 1, e ), e ) ) & 
             -Level_Node( Index_Element_2_Node( Local_Node_Number_Level_Max_2_Min( 3, e ), e ) ) 

            Number_Node_Difference_Level_3_Node( 3, e ) &
            = Level_Node( Index_Element_2_Node( Local_Node_Number_Level_Max_2_Min( 2, e ), e ) ) & 
             -Level_Node( Index_Element_2_Node( Local_Node_Number_Level_Max_2_Min( 3, e ), e ) ) 

            do i= 1, 3
               do k= 1, Number_Node_Difference_Level_3_Node( 2, e )
                  do l= 1, 2 
                      Position_Node_Boundary_Level_Plot( l, k, i, e )= 0.0d0
                  end do
               end do
            end do

            !=====================================================================================
            ! Compute Position Boundary
            !=====================================================================================
            do i= 1, 3

               if( i==3 )then
                  do k= 1, Number_Node_Difference_Level_3_Node( i, e )
                     Value_Boundary_Level_Line( k, e ) &
                     = Value_Boundary_Level&
                       ( 2, -k +Level_Node( Index_Element_2_Node( Local_Node_Number_Level_Max_2_Min( 2, e ), e ) ) )  
                  end do
                  Local_Node_Number_Edge_Line( 1, e )= Local_Node_Number_Level_Max_2_Min( 2, e )
                  Local_Node_Number_Edge_Line( 2, e )= Local_Node_Number_Level_Max_2_Min( 3, e ) 
               else
                  do k= 1, Number_Node_Difference_Level_3_Node( i, e )
                     Value_Boundary_Level_Line( k, e ) &
                     = Value_Boundary_Level&
                       ( 2, -k +Level_Node( Index_Element_2_Node( Local_Node_Number_Level_Max_2_Min( 1, e ), e ) ) )  
                  end do
                  Local_Node_Number_Edge_Line( 1, e )= Local_Node_Number_Level_Max_2_Min( 1, e )

                  if( i==1 )then
                     Local_Node_Number_Edge_Line( 2, e )= Local_Node_Number_Level_Max_2_Min( 2, e ) 
                  else if( i==2 )then
                     Local_Node_Number_Edge_Line( 2, e )= Local_Node_Number_Level_Max_2_Min( 3, e )
                  end if 
               end if


               call Compute_Position_Point_Plot_Result &
                  ( Value_Plotted, Value_Boundary_Level_Line, &
                    Local_Node_Number_Edge_Line, &
                    Position_Node_Plot_Mesh, Number_Node, Number_Node_Difference_Level_3_Node( i, e ),  &
                    Index_Element_2_Node, Number_Element, e, Max_Number_Level_Element, &
                    !======================================================================
                    Position_Node_Boundary_Level_Plot_tmp )

               do k= 1, Number_Node_Difference_Level_3_Node( i, e )
                  do l= 1, 2
                     Position_Node_Boundary_Level_Plot( l, k, i, e ) & 
                     = Position_Node_Boundary_Level_Plot_tmp( l, k, e )
                  end do
               end do

            end do

            Flag_Plot_Element( e )= 4

            !call Plot_Element_InHomogeneous_3_Node &
            !   ( Position_Node_Plot_Mesh, Index_Element_2_Node, Level_Node, RGB, & 
            !     File_Number, e, Number_Node, Number_Element, Number_Level, Max_Number_Level_Element, & 
            !     Number_Node_Difference_Level_3_Node, Local_Node_Number_Level_Max_2_Min, &
            !     Position_Node_Boundary_Level_Plot )

         !=======================================================================================
         else if( Level_Node( Index_Element_2_Node( 3, e ) ) > Level_Node( Index_Element_2_Node( 1, e ) ) .and. & 
                Level_Node( Index_Element_2_Node( 3, e ) ) > Level_Node( Index_Element_2_Node( 2, e ) ) .and. & 
                Level_Node( Index_Element_2_Node( 2, e ) ) > Level_Node( Index_Element_2_Node( 1, e ) ) )then 
         !=======================================================================================

            Local_Node_Number_Level_Max_2_Min( 1, e )= 3 
            Local_Node_Number_Level_Max_2_Min( 2, e )= 2
            Local_Node_Number_Level_Max_2_Min( 3, e )= 1


            Number_Node_Difference_Level_3_Node( 1, e ) &
            = Level_Node( Index_Element_2_Node( Local_Node_Number_Level_Max_2_Min( 1, e ), e ) ) & 
             -Level_Node( Index_Element_2_Node( Local_Node_Number_Level_Max_2_Min( 2, e ), e ) ) 

            Number_Node_Difference_Level_3_Node( 2, e ) &
            = Level_Node( Index_Element_2_Node( Local_Node_Number_Level_Max_2_Min( 1, e ), e ) ) & 
             -Level_Node( Index_Element_2_Node( Local_Node_Number_Level_Max_2_Min( 3, e ), e ) ) 

            Number_Node_Difference_Level_3_Node( 3, e ) &
            = Level_Node( Index_Element_2_Node( Local_Node_Number_Level_Max_2_Min( 2, e ), e ) ) & 
             -Level_Node( Index_Element_2_Node( Local_Node_Number_Level_Max_2_Min( 3, e ), e ) ) 

            do i= 1, 3
               do k= 1, Number_Node_Difference_Level_3_Node( 2, e )
                  do l= 1, 2 
                      Position_Node_Boundary_Level_Plot( l, k, i, e )= 0.0d0
                  end do
               end do
            end do

            !=====================================================================================
            ! Compute Position Boundary
            !=====================================================================================
            do i= 1, 3

               if( i==3 )then
                  do k= 1, Number_Node_Difference_Level_3_Node( i, e )
                     Value_Boundary_Level_Line( k, e ) &
                     = Value_Boundary_Level&
                       ( 2, -k +Level_Node( Index_Element_2_Node( Local_Node_Number_Level_Max_2_Min( 2, e ), e ) ) )  
                  end do
                  Local_Node_Number_Edge_Line( 1, e )= Local_Node_Number_Level_Max_2_Min( 2, e )
                  Local_Node_Number_Edge_Line( 2, e )= Local_Node_Number_Level_Max_2_Min( 3, e ) 
               else
                  do k= 1, Number_Node_Difference_Level_3_Node( i, e )
                     Value_Boundary_Level_Line( k, e ) &
                     = Value_Boundary_Level&
                       ( 2, -k +Level_Node( Index_Element_2_Node( Local_Node_Number_Level_Max_2_Min( 1, e ), e ) ) )  
                  end do
                  Local_Node_Number_Edge_Line( 1, e )= Local_Node_Number_Level_Max_2_Min( 1, e )

                  if( i==1 )then
                     Local_Node_Number_Edge_Line( 2, e )= Local_Node_Number_Level_Max_2_Min( 2, e ) 
                  else if( i==2 )then
                     Local_Node_Number_Edge_Line( 2, e )= Local_Node_Number_Level_Max_2_Min( 3, e )
                  end if 
               end if


               call Compute_Position_Point_Plot_Result &
                  ( Value_Plotted, Value_Boundary_Level_Line, &
                    Local_Node_Number_Edge_Line, &
                    Position_Node_Plot_Mesh, Number_Node, Number_Node_Difference_Level_3_Node( i, e ),  &
                    Index_Element_2_Node, Number_Element, e, Max_Number_Level_Element, & 
                    !======================================================================
                    Position_Node_Boundary_Level_Plot_tmp )

               do k= 1, Number_Node_Difference_Level_3_Node( i, e )
                  do l= 1, 2
                     Position_Node_Boundary_Level_Plot( l, k, i, e ) & 
                     = Position_Node_Boundary_Level_Plot_tmp( l, k, e )
                  end do
               end do

            end do

            Flag_Plot_Element( e )= 4

            !call Plot_Element_InHomogeneous_3_Node &
            !   ( Position_Node_Plot_Mesh, Index_Element_2_Node, Level_Node, RGB, & 
            !     File_Number, e, Number_Node, Number_Element, Number_Level, Max_Number_Level_Element, & 
            !     Number_Node_Difference_Level_3_Node, Local_Node_Number_Level_Max_2_Min, &
            !     Position_Node_Boundary_Level_Plot )

         !=======================================================================================
         else
         !=======================================================================================
            write(*,*)'ERROR'
            write(*,*)'  Level_Node( Index_Element_2_Node( 1, e ) )=', Level_Node( Index_Element_2_Node( 1, e ) )
            write(*,*)'  Level_Node( Index_Element_2_Node( 2, e ) )=', Level_Node( Index_Element_2_Node( 2, e ) )
            write(*,*)'  Level_Node( Index_Element_2_Node( 3, e ) )=', Level_Node( Index_Element_2_Node( 3, e ) )
            !stop
           !go to 2272
         end if
       
      end do ! do e= 1, Number_Element

      !==================================================================================================
      write(*,*) '         Make Contour Data' 
      !==================================================================================================

      LN( 1, 1 )= 1
      LN( 2, 1 )= 2
      LN( 3, 1 )= 3
      LN( 1, 2 )= 2
      LN( 2, 2 )= 3
      LN( 3, 2 )= 1
      LN( 1, 3 )= 3
      LN( 2, 3 )= 1
      LN( 3, 3 )= 2

      Number_Level_Contour = 25 !10
      allocate( Level_Contour( Number_Level_Contour ) )

      if( mod( Number_Level_Plot_EM,Number_Level_Contour )/=0 )then
         write(*,*)'==================================================================================='
         write(*,*)'Number_Level_Plot_EM=', Number_Level_Plot_EM
         write(*,*)'Number_Level_Contour=', Number_Level_Contour
         write(*,*)'mod( Number_Level_Plot_EM,Number_Level_Contour )=', mod( Number_Level_Plot_EM,Number_Level_Contour ), '.=0'
         write(*,*)'Plot_Contour_T_OC.f90 2164'
         write(*,*)'==================================================================================='
         stop
      end if

      StepSize_Contour= Number_Level_Plot_EM / Number_Level_Contour
      Intercept_Contour= StepSize_Contour/2 

      do i= 1, Number_Level_Contour
         Level_Contour( i )= Intercept_Contour+StepSize_Contour*( i-1 )
      end do

      allocate( Value_Contour( Number_Level_Contour ) )

      do i= 1, Number_Level_Contour
         Value_Contour( i )= minval( Value_Plotted ) +( maxval( Value_Plotted ) -minval( Value_Plotted ) )/ Number_Level *Level_Contour( i ) 
      end do

      allocate( Value_Node_in_Element( 4, Number_Element ) )

      do e= 1, Number_Element
         do i= 1, 4 
            Value_Node_in_Element( i, e ) = Value_Plotted( Index_Element_2_Node( i, e ) )
         end do
      end do

      !write(2147,*) Level_Contour

      allocate( Counter_Contour( Number_Element ) )
      Counter_Contour=0

      allocate( Position_Contour_tmp( 2, Max_Number_Level_Element, Number_Element ) )
      allocate( Value_Contour_Element( Max_Number_Level_Element, Number_Element ) )

      Position_Contour_tmp= 0.0d0

      ! Level_High( e ) -Level_Low( e )
      do e= 1, Number_Element
         do i= 1, Number_Level_Contour
            !if( Level_Low( e ) <= Value_Contour( i ) .and. Value_Contour( i ) <= Level_High( e ) )then
            do j= 1, 3
               if( Value_Node_in_Element( LN( 1, j ), e ) <= Value_Contour( i ) .and. &
                   Value_Node_in_Element( LN( 2, j ), e ) >= Value_Contour( i ) .and. &
                   Value_Node_in_Element( LN( 3, j ), e ) >= Value_Contour( i )  )then 

                  Counter_Contour( e ) = Counter_Contour( e ) +1
               else if( Value_Node_in_Element( LN( 1, j ), e ) >= Value_Contour( i ) .and. &
                        Value_Node_in_Element( LN( 2, j ), e ) <= Value_Contour( i ) .and. &
                        Value_Node_in_Element( LN( 3, j ), e ) <= Value_Contour( i )  )then 

                  Counter_Contour( e ) = Counter_Contour( e ) +1
               end if
            end do
         end do
      end do

      Position_Node_Line=0d0
      Value_Node_Line=0d0
      Counter_Contour_All = 0

      do e= 1, Number_Element
         !if( Counter_Contour( e ) >= 1 )then
         if( Class_Element( e )==ID_Element_FixedDomain .and. Flag_Thermal_Insulation_BC==1 )then
            Counter_Contour( e ) = 0
         else if( Counter_Contour( e ) >= 1 )then
            do i= 1, Number_Level_Contour
               do j= 1, 3
                  if( Value_Node_in_Element( LN( 1, j ), e ) <= Value_Contour( i ) .and. &
                      Value_Node_in_Element( LN( 2, j ), e ) >= Value_Contour( i ) .and. &
                      Value_Node_in_Element( LN( 3, j ), e ) >= Value_Contour( i ) )then 

                     Position_Node_Line( 1, 1, 1 )=Position_Node( 1, Index_Element_2_Node( LN( 1, j ), e ) )
                     Position_Node_Line( 2, 1, 1 )=Position_Node( 2, Index_Element_2_Node( LN( 1, j ), e ) )

                     Position_Node_Line( 1, 2, 1 )=Position_Node( 1, Index_Element_2_Node( LN( 2, j ), e ) )
                     Position_Node_Line( 2, 2, 1 )=Position_Node( 2, Index_Element_2_Node( LN( 2, j ), e ) )

                     Value_Node_Line( 1, 1 )= Value_Node_in_Element( LN( 1, j ), e ) 
                     Value_Node_Line( 2, 1 )= Value_Node_in_Element( LN( 2, j ), e )   

                     Position_Node_Line( 1, 1, 2 )=Position_Node( 1, Index_Element_2_Node( LN( 1, j ), e ) )
                     Position_Node_Line( 2, 1, 2 )=Position_Node( 2, Index_Element_2_Node( LN( 1, j ), e ) )

                     Position_Node_Line( 1, 2, 2 )=Position_Node( 1, Index_Element_2_Node( LN( 3, j ), e ) )
                     Position_Node_Line( 2, 2, 2 )=Position_Node( 2, Index_Element_2_Node( LN( 3, j ), e ) )

                     Value_Node_Line( 1, 2 )= Value_Node_in_Element( LN( 1, j ), e )  
                     Value_Node_Line( 2, 2 )= Value_Node_in_Element( LN( 3, j ), e )   

                     !Value_Contour_Element( Counter_Contour( e ), e ) = Value_Boundary_Level( 1, Value_Contour( i ) )  
                     !Value_Contour_Element( Counter_Contour( e ), e ) = Value_Contour( i ) 
                     Value_Contour_Element( 1, e ) = Value_Contour( i ) 

                     goto 2213
                  else if( Value_Node_in_Element( LN( 1, j ), e ) >= Value_Contour( i ) .and. &
                           Value_Node_in_Element( LN( 2, j ), e ) <= Value_Contour( i ) .and. &
                           Value_Node_in_Element( LN( 3, j ), e ) <= Value_Contour( i )  )then 

                     Position_Node_Line( 1, 1, 1 )=Position_Node( 1, Index_Element_2_Node( LN( 1, j ), e ) )
                     Position_Node_Line( 2, 1, 1 )=Position_Node( 2, Index_Element_2_Node( LN( 1, j ), e ) )

                     Position_Node_Line( 1, 2, 1 )=Position_Node( 1, Index_Element_2_Node( LN( 2, j ), e ) )
                     Position_Node_Line( 2, 2, 1 )=Position_Node( 2, Index_Element_2_Node( LN( 2, j ), e ) )

                     Value_Node_Line( 1, 1 )= Value_Node_in_Element( LN( 1, j ), e )  
                     Value_Node_Line( 2, 1 )= Value_Node_in_Element( LN( 2, j ), e )   

                     Position_Node_Line( 1, 1, 2 )=Position_Node( 1, Index_Element_2_Node( LN( 1, j ), e ) )
                     Position_Node_Line( 2, 1, 2 )=Position_Node( 2, Index_Element_2_Node( LN( 1, j ), e ) )

                     Position_Node_Line( 1, 2, 2 )=Position_Node( 1, Index_Element_2_Node( LN( 3, j ), e ) )
                     Position_Node_Line( 2, 2, 2 )=Position_Node( 2, Index_Element_2_Node( LN( 3, j ), e ) )

                     Value_Node_Line( 1, 2 )= Value_Node_in_Element( LN( 1, j ), e )  
                     Value_Node_Line( 2, 2 )= Value_Node_in_Element( LN( 3, j ), e )   

                     !Value_Contour_Element( Counter_Contour( e ), e ) = Value_Boundary_Level( 1, Value_Contour( i ) )  
                     !Value_Contour_Element( Counter_Contour( e ), e ) = Value_Contour( i ) 
                     Value_Contour_Element( 1, e ) = Value_Contour( i ) 

                     goto 2213
                  end if
               end do
            end do

            2213 continue

            do j= 1, 2

               Position_Node_1_tmp( : )= Position_Node_Line( :, 1, j )
               Position_Node_2_tmp( : )= Position_Node_Line( :, 2, j )
               Value_Node_Line_1_tmp= Value_Node_Line( 1, j )
               Value_Node_Line_2_tmp= Value_Node_Line( 2, j )

               call Compute_Position_Point_Plot_Scalor &
                    ( Value_Contour_Element( 1, e ), &
                      Position_Node_1_tmp, Position_Node_2_tmp, &
                      Value_Node_Line_1_tmp, Value_Node_Line_2_tmp, &
                      !======================================================================
                      Position_Interpolated )

               Position_Contour_tmp( :, j, e ) = Position_Interpolated( : )

            end do
            Counter_Contour_All = Counter_Contour_All +1
         end if
      end do

      allocate( Position_Contour( 2, 2, Counter_Contour_All ) )

      Counter_tmp = 0
      do e= 1, Number_Element 
         if( Counter_Contour( e ) >= 1 )then
            Counter_tmp = Counter_tmp +1
            do i= 1, 2
               do j= 1, 2
                  Position_Contour( j, i, Counter_tmp ) = Position_Contour_tmp( j, i, e )
               end do
            end do
         end if 
      end do

      !do e= 1, Counter_Contour_All 
      !   write(2325,*) Position_Contour( 1, 1, e ), Position_Contour( 2, 1, e )
      !   write(2325,*) Position_Contour( 1, 2, e ), Position_Contour( 2, 2, e )
      !   write(2325,*)' ' 
      !end do

      do e= 1, Counter_Contour_All 
         do i= 1, 2
            do j= 1, 2
               Position_Contour( j, i, e )&
               = ( Position_Contour( j, i, e ) -Position_Minimum( j ) )  &
                /( Position_Maximum( j ) -Position_Minimum( j ) )*dble( Width_PostScript( j ) ) &
                +Translation_Position( j ) 
            end do
         end do
      end do

      deallocate( Value_Node_in_Element )
      deallocate( Position_Contour_tmp )
      deallocate( Value_Contour_Element )
      deallocate( Value_Contour )
      deallocate( Level_Contour )
      deallocate( Counter_Contour )

      !==================================================================================================
      write(*,*) '         Detect Boundary of Configuration'
      !==================================================================================================

      allocate( Index_Element_2_Node_tmp( 3, Number_Element ) )

      do e= 1, Number_Element
         do i= 1, 3
            Index_Element_2_Node_tmp( i, e )= Index_Element_2_Node( i, e )
         end do
      end do

      allocate( Element_Number_Plot_Boundary_tmp( Number_Element ) )
      allocate( Local_Node_Number_Plot_Boundary_tmp( 2, Number_Element ) )

      call Detect_Boundary_Data &
         ( ID_Element_Structure, &
           Number_Node, Position_Node, &
           Number_Element, 3, Index_Element_2_Node_tmp, Class_Element, &
           0, &
         !====================================================================================================
           Number_Edge_Plot_Boundary, Element_Number_Plot_Boundary_tmp, Local_Node_Number_Plot_Boundary_tmp ) 

      write(*,*)'         Number_Edge_Plot_Boundary=', Number_Edge_Plot_Boundary

      deallocate( Index_Element_2_Node_tmp )

      if( Number_Edge_Plot_Boundary==0 )then
         !call Output_Error( 'Plot_Contour_T_OC', 2281 ) 
      end if

      do e= 1, Number_Edge_Plot_Boundary
         if( Class_Element( Element_Number_Plot_Boundary_tmp( e ) ) /= ID_Element_Structure )then 
            write(*,*)'Class_Element( Element_Number_Plot_Boundary_tmp( e ) )=', Class_Element( Element_Number_Plot_Boundary_tmp( e ) )
            call Output_Error( 'Plot_Contour_T_OC', 2298 ) 
         end if
      end do

      !==================================================================================================
      write(*,*) '         Boundary Data --> Node Data'
      !==================================================================================================

      allocate( Position_Boundary_Plot_tmp( 2, 2, Number_Edge_Plot_Boundary ) )
 
      do e= 1, Number_Edge_Plot_Boundary
         do i= 1, 2
            do j= 1, 2
               Position_Boundary_Plot_tmp( j, i, e )&
               = Position_Node( j, Index_Element_2_Node( Local_Node_Number_Plot_Boundary_tmp( i, e ), Element_Number_Plot_Boundary_tmp( e ) ) ) 
            end do
         end do
      end do

      !!==================================================================================================
      !write(*,*) '         Output Boundary Data .......... fort.998 '
      !!==================================================================================================
      !do e= 1, Number_Edge_Plot_Boundary
      !   write( 998, * ) Position_Boundary_Plot_tmp( 1, 1, e ), Position_Boundary_Plot_tmp( 2, 1, e )
      !   write( 998, * ) Position_Boundary_Plot_tmp( 1, 2, e ), Position_Boundary_Plot_tmp( 2, 2, e )
      !   write( 998, * ) ' '
      !end do

      !==================================================================================================
      write(*,*) '         Check Edge Length'
      !==================================================================================================
      
      allocate( Edge_Length( Number_Edge_Plot_Boundary ) )

      do e= 1, Number_Edge_Plot_Boundary
         Edge_Length( e ) &
         = sqrt( ( Position_Boundary_Plot_tmp( 1, 1, e ) -Position_Boundary_Plot_tmp( 1, 2, e ) )**2 &
              +( Position_Boundary_Plot_tmp( 2, 1, e ) -Position_Boundary_Plot_tmp( 2, 2, e ) )**2 ) 

      end do

      Counter_Edge= 0
      do e= 1, Number_Edge_Plot_Boundary
         if( Edge_Length( e ) > 1d0/200d0*sqrt(2d0) +1d-4 )then
            Counter_Edge= Counter_Edge +1
         end if
      end do
      write(*,*)'Counter_Edge=', Counter_Edge 
      write(*,*)'Number_Edge_Plot_Boundary=', Number_Edge_Plot_Boundary

      deallocate( Edge_Length )

      !==================================================================================================
      write(*,*) '         Sizing : Position_Boundary_Plot_tmp --> Position_Boundary_Plot'
      !==================================================================================================

      allocate( Position_Boundary_Plot( 2, 2, Number_Edge_Plot_Boundary ) )

      !$omp parallel do default( none ) &
      !$omp private( e, i, j ) &
      !$omp shared( Number_Edge_Plot_Boundary, Position_Boundary_Plot, Position_Boundary_Plot_tmp ) & 
      !$omp shared( Position_Minimum, Position_Maximum, Translation_Position, Width_PostScript ) 
      do e= 1, Number_Edge_Plot_Boundary
         do i= 1, 2
            do j= 1, 2
               Position_Boundary_Plot( j, i, e )&
               = ( Position_Boundary_Plot_tmp( j, i, e ) -Position_Minimum( j ) )  &
                /( Position_Maximum( j ) -Position_Minimum( j ) )*dble( Width_PostScript( j ) ) &
                +Translation_Position( j ) 
            end do
         end do
      end do

      deallocate( Position_Boundary_Plot_tmp )
      deallocate( Element_Number_Plot_Boundary_tmp )
      deallocate( Local_Node_Number_Plot_Boundary_tmp )

      !==================================================================================================
      write(*,*) '         Output Boundary Data .......... fort.999 '
      !==================================================================================================
      do e= 1, Number_Edge_Plot_Boundary
         write( 999, * ) Position_Boundary_Plot( 1, 1, e ), Position_Boundary_Plot( 2, 1, e )
         write( 999, * ) Position_Boundary_Plot( 1, 2, e ), Position_Boundary_Plot( 2, 2, e )
         write( 999, * ) ' '
      end do

      !!==================================================================================================
      !write(*,*) '         Sizing : Position_Arrow_PV --> Position_Arrow_PV_Plot'
      !!==================================================================================================

      !allocate( Position_Arrow_PV_Plot( 2, 2, Number_PV ) )
      !allocate( Length_Arrow_PV_Plot( Number_PV ) )

      !do e= 1, Number_PV
      !   do i= 1, 2
      !      do j= 1, 2
      !         Position_Arrow_PV_Plot( j, i, e )&
      !         = ( Position_Arrow_PV( j, i, e ) -Position_Minimum( j ) )  &
      !          /( Position_Maximum( j ) -Position_Minimum( j ) )*dble( Width_PostScript( j ) ) &
      !          +Translation_Position( j ) 
      !      end do
      !   end do
      !end do

      !do e= 1, Number_PV
      !   do i= 1, 2
      !      do j= 1, 2
      !         Length_Arrow_PV_Plot( e )&
      !         = sqrt( ( Position_Arrow_PV_Plot( 1, 2, e ) -Position_Arrow_PV_Plot( 1, 1, e ) )**2  &
      !                +( Position_Arrow_PV_Plot( 2, 2, e ) -Position_Arrow_PV_Plot( 2, 1, e ) )**2 )  
      !      end do
      !   end do
      !end do

      !Length_Maximum_Vector_Position &
      != Length_Maximum_Vector /( Position_Maximum( 1 ) -Position_Minimum( 1 ) )*dble( Width_PostScript( 1 ) )

      !==================================================================================================
      write(*,*) '         Plot Element Data ..........'
      !==================================================================================================

      !hoge 
      do i= 1, 2 
         Position_Minimum( i )= Translation_Position( i ) 
         Position_Maximum( i )= Translation_Position( i ) +Width_PostScript( i )
      end do

      open( File_Number, file=Filename_Result, status='replace' )
      write(File_Number,*)'%!PS-Adobe-2.0'
      write(File_Number,*)'%BoundingBox: 0 0', Width_PostScript( 1 ), Width_PostScript( 2 )
      write(File_Number,*)'%EndHeader'
      write(File_Number,*)'stroke'
      write(File_Number,*)'0.0001 setlinewidth'
      write(File_Number,*)'0.3 setgray'

      Counter_Percent= 1 

      do e= 1, Number_Element
         
         if( e== Number_Element/100*Counter_Percent )then
            !write(*,'(a,i3,a)', advance='no')'[', Counter_Percent, '%]'
            write(*,'(a,i3,a,$)')'[', Counter_Percent, '%]'
            Counter_Percent= Counter_Percent +1
         end if 

         2244 continue

         if( Class_Element( e )==ID_Element_FixedDomain .and. Flag_Thermal_Insulation_BC==1 )then

         else if( Flag_Plot_Element( e )==3 )then

            call Plot_Element_InHomogeneous &
               ( Position_Node_Plot_Mesh, Index_Element_2_Node, Level_Node, RGB, & 
                 File_Number, e, Number_Node, Number_Element, Number_Level, Max_Number_Level_Element, & 
                 Number_Node_Difference_Level( e ), & 
                 Local_Node_Number_Different_Level( e ), Local_Node_Number_Edge_Line_Plot,  &
                 Position_Node_Boundary_Level_Plot_1, Position_Node_Boundary_Level_Plot_2 )

         else if( Flag_Plot_Element( e )==4 )then

            call Plot_Element_InHomogeneous_3_Node &
               ( Position_Node_Plot_Mesh, Index_Element_2_Node, Level_Node, RGB, & 
                 File_Number, e, Number_Node, Number_Element, Number_Level, Max_Number_Level_Element, & 
                 Number_Node_Difference_Level_3_Node, Local_Node_Number_Level_Max_2_Min, &
                 Position_Node_Boundary_Level_Plot )

         else if( Flag_Plot_Element( e )==2 )then

            call Plot_Gradation_Element &
               ( File_Number, e, Number_Element, Number_Level, Max_Number_Level_Element, & 
                 Number_Level_Element, Number_Node_Level_Element, &
                 Position_Node_Plot_Mesh_Element, RGB, Level_Element )

         else if( Flag_Plot_Element( e )==1 )then

            call Plot_Gradation_Element_Monocromatic &
               ( File_Number, e, Number_Element, Max_Number_Level_Element, & 
                 Number_Level_Element, Number_Node_Level_Element, Position_Node_Plot_Mesh_Element, 1, & 
                 RGB_Lowest, RGB_Highest )

         else if( Flag_Plot_Element( e )==0 )then

            call Plot_Gradation_Element_Monocromatic &
               ( File_Number, e, Number_Element, Max_Number_Level_Element, & 
                 Number_Level_Element, Number_Node_Level_Element, Position_Node_Plot_Mesh_Element, 0, & 
                 RGB_Lowest, RGB_Highest )

         else
            call Output_Error( 'Plot_Contour_T_OC', 2220 ) 
         end if

      end do

      deallocate( Flag_Plot_Element )

      !==================================================================================================
      ! Plot contour 
      !==================================================================================================

      write( File_Number, * )'1 1 1 setrgbcolor'
      do e= 1, Counter_Contour_All 
         write( File_Number, * )'1 setlinewidth'
         write( File_Number, * ) 'newpath'
         write( File_Number, * ) Position_Contour( 1, 1, e ), Position_Contour( 2, 1, e ), 'moveto' 
         write( File_Number, * ) Position_Contour( 1, 2, e ), Position_Contour( 2, 2, e ), 'lineto' 
         write( File_Number, * )'stroke'
      end do

      !==================================================================================================
      ! Plot Fixed Domain
      !==================================================================================================

      if( Flag_FixedDomain==1 )then

         if( Flag_Thermal_Device==20 .or. Flag_Thermal_Device==21 .or. Flag_Thermal_Device==30 )then
            write( File_Number, * )'1 setlinewidth'
            write( File_Number, * )'0 0 0 setrgbcolor'
            write( File_Number, * ) 'newpath'
            write( File_Number, * ) Position_Center_FixedDomain_X_PS, Position_Center_FixedDomain_Y_PS, &
                         Radius_FixedDomain_Plot_PS,' 0 360 arc'
            write( File_Number, * )'closepath'
            write( File_Number, * )'stroke'
         else 
            write( File_Number, * )'1 setlinewidth'
            write( File_Number, * ) 'newpath'
            write( File_Number, * ) Position_Center_FixedDomain_X_PS, Position_Center_FixedDomain_Y_PS, &
                         Radius_FixedDomain_Plot_PS,' 0 360 arc'
            write( File_Number, * )'gsave'
            write( File_Number, * )'0 0 0 setrgbcolor'
            write( File_Number, * )'fill'
            write( File_Number, * )'grestore'
            write( File_Number, * )'0 0 0 setrgbcolor'
            write( File_Number, * )'stroke'
         end if

      else if( Flag_FixedDomain==2 .and. Flag_Thermal_Device/=210 )then
         write( File_Number, * )'1 setlinewidth'
         write( File_Number, * ) Position_Corner_FixedDomain_X1_PS, Position_Corner_FixedDomain_Y1_PS, 'moveto' 
         write( File_Number, * ) Position_Corner_FixedDomain_X2_PS, Position_Corner_FixedDomain_Y1_PS, 'lineto' 
         write( File_Number, * ) Position_Corner_FixedDomain_X2_PS, Position_Corner_FixedDomain_Y2_PS, 'lineto' 
         write( File_Number, * ) Position_Corner_FixedDomain_X1_PS, Position_Corner_FixedDomain_Y2_PS, 'lineto' 
         write( File_Number, * ) Position_Corner_FixedDomain_X1_PS, Position_Corner_FixedDomain_Y1_PS, 'lineto' 
         write( File_Number, * )'stroke'

         if( Flag_Thermal_Device==50 .or. Flag_Thermal_Device==51 )then
             write( File_Number, * )'1 setlinewidth'
             write( File_Number, * ) 'newpath'
             write( File_Number, * ) ( Position_Corner_FixedDomain_X1_PS +Position_Corner_FixedDomain_X2_PS )/2d0, &
                           Position_Corner_FixedDomain_Y2_PS, & 
                           Radius_FixedDomain_Plot_PS,' 0 180 arc'
             write( File_Number, * )'closepath'
             write( File_Number, * )'stroke'
         end if

      end if

      !==================================================================================================
      ! Plot Design Domain
      !==================================================================================================

      if( Flag_DesignDomain==1 )then

         write( File_Number, * )'0 0 0 setrgbcolor'
         if( Flag_Thermal_Device==50 .or. Flag_Thermal_Device==51 )then
            write( File_Number, * )'1 setlinewidth'
            write( File_Number, * ) 'newpath'
            write( File_Number, * ) Position_Center_DesignDomain_X_PS, &
                            Position_Center_DesignDomain_Y_PS, &
                            Radius_DesignDomain_Plot_PS,' 0 180 arc'
            write( File_Number, * )'closepath'
            write( File_Number, * )'stroke'
         else
            write( File_Number, * )'1 setlinewidth'
            write( File_Number, * ) 'newpath'
            write( File_Number, * ) Position_Center_DesignDomain_X_PS, Position_Center_DesignDomain_Y_PS, &
                            Radius_DesignDomain_Plot_PS,' 0 360 arc'
            write( File_Number, * )'closepath'
            write( File_Number, * )'stroke'
         end if

      else if( Flag_DesignDomain==2 )then
         write( File_Number, * )'1 setlinewidth'
         write( File_Number, * ) Position_Corner_DesignDomain_X1_PS, Position_Corner_DesignDomain_Y1_PS, 'moveto' 
         write( File_Number, * ) Position_Corner_DesignDomain_X2_PS, Position_Corner_DesignDomain_Y1_PS, 'lineto' 
         write( File_Number, * ) Position_Corner_DesignDomain_X2_PS, Position_Corner_DesignDomain_Y2_PS, 'lineto' 
         write( File_Number, * ) Position_Corner_DesignDomain_X1_PS, Position_Corner_DesignDomain_Y2_PS, 'lineto' 
         write( File_Number, * ) Position_Corner_DesignDomain_X1_PS, Position_Corner_DesignDomain_Y1_PS, 'lineto' 
         write( File_Number, * )'stroke'
      end if


      !==================================================================================================
      ! Plot Outer Domain 
      !==================================================================================================

         write( File_Number, * )'0.5 setlinewidth'
         write( File_Number, * )'0.0 setgray' 
         write( File_Number, * ) Position_Corner_WholeDomain_X1_PS, Position_Corner_WholeDomain_Y1_PS, 'moveto' 
         write( File_Number, * ) Position_Corner_WholeDomain_X2_PS, Position_Corner_WholeDomain_Y1_PS, 'lineto' 
         write( File_Number, * ) Position_Corner_WholeDomain_X2_PS, Position_Corner_WholeDomain_Y2_PS, 'lineto' 
         write( File_Number, * ) Position_Corner_WholeDomain_X1_PS, Position_Corner_WholeDomain_Y2_PS, 'lineto' 
         write( File_Number, * ) Position_Corner_WholeDomain_X1_PS, Position_Corner_WholeDomain_Y1_PS, 'lineto' 
         write( File_Number, * )'stroke'
      !==================================================================================================
      ! Plot Configuration
      !==================================================================================================

      open( 2536, file='./Edge.ps', status='replace' )
         write(2536,*)'%!PS-Adobe-2.0'
         write(2536,*)'%BoundingBox: 0 0', Width_PostScript( 1 ), Width_PostScript( 2 )
         write(2536,*)'%EndHeader'
         write(2536,*)'stroke'

         do e= 1, Number_Edge_Plot_Boundary
            write( 2536, * )'0.5 setlinewidth'
            write( 2536, * ) 'newpath'
            write( 2536, * ) Position_Boundary_Plot( 1, 1, e ), Position_Boundary_Plot( 2, 1, e ), 'moveto' 
            write( 2536, * ) Position_Boundary_Plot( 1, 2, e ), Position_Boundary_Plot( 2, 2, e ), 'lineto' 
            write( 2536, * )'stroke'
         end do
      close( 2536 )

      write( File_Number, * )'0 0 0 setrgbcolor'
      do e= 1, Number_Edge_Plot_Boundary
         write( File_Number, * )'0.5 setlinewidth'
         write( File_Number, * ) 'newpath'
         write( File_Number, * ) Position_Boundary_Plot( 1, 1, e ), Position_Boundary_Plot( 2, 1, e ), 'moveto' 
         write( File_Number, * ) Position_Boundary_Plot( 1, 2, e ), Position_Boundary_Plot( 2, 2, e ), 'lineto' 
         write( File_Number, * )'stroke'
      end do

      !!==================================================================================================
      !! Plot Poynting Vector 
      !!==================================================================================================

      !!Stemthick= 0.2d0!*Length_Maximum_Vector_Position
      !!Headthick= Length_Maximum_Vector_Position/7d0 
      !!Headlength= Length_Maximum_Vector_Position/7d0
      !Stemthick= 0.1d0*Length_Maximum_Vector_Position
      !Headthick= Length_Maximum_Vector_Position/7d0 
      !Headlength= Length_Maximum_Vector_Position/5d0

      !allocate( Stemthick_Plot( Number_PV ) )
      !allocate( Headthick_Plot( Number_PV ) )
      !allocate( Headlength_Plot( Number_PV ) )

      !Stemthick_Plot(:)= Length_Arrow_PV_Plot(:)*0.1
      !Headthick_Plot(:)= Length_Arrow_PV_Plot(:)*0.2
      !Headlength_Plot(:)= Length_Arrow_PV_Plot(:)*0.5

      !allocate( Level_Flux_tmp( Number_PV ) )

      !do e= 1, Number_PV
      !   Level_Flux_tmp( e )= Level_PV( e )*0.9d0 +0.1d0
      !end do

      !!write( File_Number, * ) '# plot heat flux arrow'
      !do e= 1, Number_PV
      !   call Plot_Arrow_Postscript&
      !      ( File_Number, 2, 'w', &
      !        Position_Arrow_PV_Plot( 1, 1, e ), Position_Arrow_PV_Plot( 2, 1, e ), &
      !        Position_Arrow_PV_Plot( 1, 2, e ), Position_Arrow_PV_Plot( 2, 2, e ), & 
      !        Stemthick_Plot( e ), Headthick_Plot( e ), Headlength_Plot( e ) )
      !        !Stemthick*Level_Flux_tmp( e ), Headthick*Level_Flux_tmp( e ), Headlength*Level_Flux_tmp( e ) )
      !        !Stemthick, Headthick, Headlength )
      !end do
      !!write( File_Number, * ) '# end plot heat flux arrow'

      !deallocate( Level_Flux_tmp )

      !==================================================================================================
      ! Plot Color Var 
      !==================================================================================================
      if( Flag_Color_Var==1 )then
         call Plot_Color_Var_Right_Side & 
            ( File_Number, RGB, Number_Level, Position_Maximum_Color_Var, Minimum_Value_Plotted, Maximum_Value_Plotted )
      else if( Flag_Color_Var==2 )then
         call Plot_Color_Var_Upper_Side & 
            ( File_Number, RGB, Number_Level, Position_Maximum_Color_Var, Minimum_Value_Plotted, Maximum_Value_Plotted )
      end if


      !write( File_Number, * ) '/Helvetica findfont 15 scalefont setfont'
      !write( File_Number, * ) Position_Minimum( 1 ), Position_Minimum( 2 )-20d0, 'moveto'
      !write( File_Number, * ) '(Objective Function=) show'
 
      !write( File_Number, * ) '/Helvetica findfont 15 scalefont setfont'
      !write( File_Number, * ) ( Position_Minimum( 1 )+Position_Maximum( 1 ) )*1d0/2d0, &
      !               Position_Minimum( 2 )-20d0, 'moveto'
      !write( File_Number, * ) '(Optimization Step=) show'
      !write( File_Number, * ) '/Helvetica findfont 18 scalefont setfont'
      !write( File_Number, * ) Position_Maximum( 1 )-50d0, Position_Minimum( 2 )-20d0, 'moveto'
      !write( File_Number, 2127 )'(', Optimization_Step, ') show'
      !2127 format(A1,I5,A6) 

      !==================================================================================================

      write( File_Number, * ) '/Times-italic findfont 18 scalefont setfont'
      write( File_Number, * ) Position_Minimum( 1 ) +15d0, Position_Minimum( 2 )-20d0, 'moveto'
      write( File_Number, * ) '(F=) show'
      write( File_Number, * ) '/Times-Roman findfont 18 scalefont setfont'
      write( File_Number, * ) ( Position_Minimum( 1 )+Position_Maximum( 1 ) )*2d0/16d0, &
                      Position_Minimum( 2 )-20d0, 'moveto'
      write( File_Number, 2126 )'(', Objective_Function( Optimization_Step ), ') show'
      2126 format(A1,Es11.4,A6) 

      write( File_Number, * ) '/Symbol findfont 18 scalefont setfont'
      write( File_Number, * ) ( Position_Minimum( 1 )+Position_Maximum( 1 ) )*1d0/2d0, & 
                      Position_Minimum( 2 )-20d0, 'moveto'
      write( File_Number, * ) '(\164=) show'
      write( File_Number, * ) '/Times-Roman findfont 18 scalefont setfont'
      write( File_Number, * ) ( Position_Minimum( 1 )+Position_Maximum( 1 ) )*1d0/2d0+20d0, &
                      Position_Minimum( 2 )-20d0, 'moveto'
      !write( File_Number, 2127 )'(', Coefficient_Complexity_Tau_Normal, ') show'
      2127 format(A1,ES8.1,A6) 

      close( File_Number )

      deallocate( Position_Minimum )
      deallocate( Position_Maximum )

      write(*,*)'         ==> HF_T_OC__', Optimization_Step_Character, '.ps' 
     
   if( Flag_Color_Var==1 .or. Flag_Color_Var==2 )then
      deallocate( Position_Maximum_Color_Var )
   end if

   deallocate( Translation_Position )
   deallocate( Width_PostScript )

   deallocate( Position_Node_Plot_Mesh )
 
   deallocate( RGB )
   deallocate( RGB_Lowest )
   deallocate( RGB_Highest )
   deallocate( Level_Node )
   deallocate( Value_Boundary_Level )

   deallocate( Number_Level_Element )
   deallocate( Position_Node_Plot_Mesh_Element )
   deallocate( Level_Element )
   deallocate( Local_Node_Number_Different_Level ) 
   deallocate( Level_High ) 
   deallocate( Level_Low ) 
   deallocate( Number_Node_Difference_Level ) 
 
   deallocate( Local_Node_Number_Edge_Line_Plot )

   deallocate( Local_Node_Number_Edge_Line )
   deallocate( Value_Boundary_Level_Line )
   deallocate( Number_Node_Difference_Level_3_Node ) 
   deallocate( Difference_Position_1_Minus_NewPoint )            
   deallocate( Difference_Position_2_Minus_NewPoint )            
 
   deallocate( Local_Node_Number_Level_Max_2_Min )

   deallocate( Position_Node_Boundary_Level_Plot ) 
   deallocate( Position_Node_Boundary_Level_Plot_tmp )
   deallocate( Position_Node_Boundary_Level_Plot_1 )
   deallocate( Position_Node_Boundary_Level_Plot_2 )

   deallocate( Position_Boundary_Plot )
   !deallocate( Position_Arrow_PV_Plot )
   !deallocate( Length_Arrow_PV_Plot )
   !deallocate( Stemthick_Plot )
   !deallocate( Headthick_Plot )
   !deallocate( Headlength_Plot )
   deallocate( Position_Contour )


   2272 continue

   return
end subroutine Plot_Contour_T_OC


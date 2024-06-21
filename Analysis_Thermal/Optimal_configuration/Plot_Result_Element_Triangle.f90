
subroutine Plot_Result_Element_Triangle &
       ( File_Number, &
         Value_Plotted, Number_Level, Optimization_Step,  &
         Position_Node, Number_Node, &
         Index_Element_2_Node, Number_Element, &
         Objective_Function ) 

   !use omp_lib
   use Parameters
   implicit none

   integer, intent(in) :: File_Number
   integer, intent(in) :: Number_Node, Number_Element
   double precision, intent(inout) :: Value_Plotted( Number_Node )

   integer, intent(in) :: Number_Level
   integer, intent(in) :: Optimization_Step

   double precision, intent(in) :: Position_Node( 2, Number_Node )
   integer, intent(in) :: Index_Element_2_Node( 4, Number_Element )
 
   double precision, intent(in) :: Objective_Function( 0:Number_Optimization_Step ) 
 
   character(len=Length_Character_Optimization_Step) :: Optimization_Step_Character
   character(len=8) :: Format_Filenumber
   
   integer :: e, i, j, k, l 
   
   character(len=256) :: Filename_Result
   
   double precision, allocatable, dimension(:) :: Position_Minimum 
   double precision, allocatable, dimension(:) :: Position_Maximum 
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

   !================================================
   write(*,*)'      call Plot_Result_Element_Triangle'
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
            Position_Maximum_Color_Var( 1, j, i ) &
            = ( Position_Maximum_Color_Var( 1, j, i ) -Position_Minimum( 1 ) )  &
             /( Position_Maximum( 1 ) -Position_Minimum( 1 ) )*dble( WidthPS_X ) 
            Position_Maximum_Color_Var( 2, j, i ) &
            = ( Position_Maximum_Color_Var( 2, j, i ) -Position_Minimum( 2 ) )  &
             /( Position_Maximum( 2 ) -Position_Minimum( 2 ) )*dble( WidthPS_Y ) 
         end do
      end do
   
      do i= 1, 2
         do j= 1, 2
            Position_Maximum_Color_Var( 1, j, i )= Position_Maximum_Color_Var( 1, j, i ) +Translation_X 
            Position_Maximum_Color_Var( 2, j, i )= Position_Maximum_Color_Var( 2, j, i ) +Translation_Y
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
            Position_Maximum_Color_Var( 1, j, i ) &
            = ( Position_Maximum_Color_Var( 1, j, i ) -Position_Minimum( 1 ) )  &
             /( Position_Maximum( 1 ) -Position_Minimum( 1 ) )*dble( WidthPS_X ) 
            Position_Maximum_Color_Var( 2, j, i ) &
            = ( Position_Maximum_Color_Var( 2, j, i ) -Position_Minimum( 2 ) )  &
             /( Position_Maximum( 2 ) -Position_Minimum( 2 ) )*dble( WidthPS_Y ) 
         end do
      end do
   
      do i= 1, 2
         do j= 1, 2
            Position_Maximum_Color_Var( 1, j, i )= Position_Maximum_Color_Var( 1, j, i ) +Translation_X 
            Position_Maximum_Color_Var( 2, j, i )= Position_Maximum_Color_Var( 2, j, i ) +Translation_Y
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
   !$omp shared( Number_Node, Position_Node_Plot_Mesh, Position_Node, Position_Minimum, Position_Maximum ) 
   do i= 1, Number_Node
      Position_Node_Plot_Mesh( 1, i )= ( Position_Node( 1, i ) -Position_Minimum( 1 ) )  &
                          /( Position_Maximum( 1 ) -Position_Minimum( 1 ) )*dble( WidthPS_X ) 
      Position_Node_Plot_Mesh( 2, i )= ( Position_Node( 2, i ) -Position_Minimum( 2 ) )  &
                          /( Position_Maximum( 2 ) -Position_Minimum( 2 ) )*dble( WidthPS_Y ) 
   end do
 
   !$omp parallel do default( none ) &
   !$omp private( i ) &
   !$omp shared( Number_Node, Position_Node_Plot_Mesh ) 
   do i= 1, Number_Node
      Position_Node_Plot_Mesh( 1, i )= Position_Node_Plot_Mesh( 1, i ) +Translation_X 
      Position_Node_Plot_Mesh( 2, i )= Position_Node_Plot_Mesh( 2, i ) +Translation_Y
   end do

   Position_Minimum( 1 )= Translation_X 
   Position_Minimum( 2 )= Translation_Y 
   Position_Maximum( 1 )= Translation_X +WidthPS_X
   Position_Maximum( 2 )= Translation_Y +WidthPS_Y

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
      write(*,*)'Minimum_Value_Plotted=', Minimum_Value_Plotted
      write(*,*)'Maximum_Value_Plotted=', Maximum_Value_Plotted
      call Output_Error( 'Plot_Result_Element_Triangle', 270 )
   end if

   write(*,*)'      Maximum_Value_Plotted=', Maximum_Value_Plotted
   write(*,*)'      Minimum_Value_Plotted=', Minimum_Value_Plotted

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
   write(*,*)'      Max_Number_Level_Element=', Max_Number_Level_Element
   write(*,*)'      Max_Number_Level_Element_Limit=', Max_Number_Level_Element_Limit
   write(*,*)'==============================================================='

   if( Max_Number_Level_Element > Max_Number_Level_Element_Limit )then
      write(*,*)'      Max_Number_Level_Element > Max_Number_Level_Element_Limit '
     !go to 2272 
      ! go to end  of this subroutine 
      !call Output_Error( 'Plot_Result_Element_Triangle', 353 ) 
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
            call Output_Error( 'Plot_Result_Element_Triangle', 334 ) 
         else if( Level_Node( i ) > Number_Level )then
            write(*,*)'Level_Node( i )=', Level_Node( i )
            write(*,*)'Number_Level=', Number_Level
            write(*,*)'Value_Plotted( i )=', Value_Plotted( i )
            write(*,*)'Value_Boundary_Level( 1, 1 )', Value_Boundary_Level( 1, 1 )
            write(*,*)'Value_Boundary_Level( 2, Number_Level )', Value_Boundary_Level( 2, Number_Level )
            write(*,*)'i=', i
            call Output_Error( 'Plot_Result_Element_Triangle', 343 ) 
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
 
   Filename_Result= trim( "Temperature_Distribution_"//trim(Optimization_Step_Character)//".ps" )


   !==================================================================================================
   write(*,*) '      Create Element Data'
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
               write(*,*)'Plot_Result_Element_Triangle.f90 264'
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
                  write(*,*)'Plot_Result_Element_Triangle.f90 304' 
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
                  write(*,*)'Plot_Result_Element_Triangle.f90 381' 
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
               write(*,*)'Plot_Result_Element_Triangle.f90 448'
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
                  write(*,*)'Plot_Result_Element_Triangle.f90 509' 
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
                  write(*,*)'Plot_Result_Element_Triangle.f90 573' 
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
               write(*,*)'Plot_Result_Element_Triangle.f90 634'
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
                  write(*,*)'Plot_Result_Element_Triangle.f90 702' 
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
                  write(*,*)'Plot_Result_Element_Triangle.f90 763' 
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
               write(*,*)'Plot_Result_Element_Triangle.f90 448'
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
                  write(*,*)'Plot_Result_Element_Triangle.f90 888' 
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
                  write(*,*)'Plot_Result_Element_Triangle.f90 949' 
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
               write(*,*)'Plot_Result_Element_Triangle.f90 826'
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
                  write(*,*)'Plot_Result_Element_Triangle.f90 1073' 
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
                  write(*,*)'Plot_Result_Element_Triangle.f90 1134' 
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
               write(*,*)'Plot_Result_Element_Triangle.f90 448'
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
                  write(*,*)'Plot_Result_Element_Triangle.f90 1259' 
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
                  write(*,*)'Plot_Result_Element_Triangle.f90 1320' 
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
       
      end do

      !==================================================================================================
      write(*,*) '      Plot Element Data'
      !==================================================================================================

      open( File_Number, file=Filename_Result, status='replace' )
      write(File_Number,*)'%!PS-Adobe-2.0'
      write(File_Number,*)'%BoundingBox: 0 0', WidthPS_X, WidthPS_Y
      write(File_Number,*)'%EndHeader'
      write(File_Number,*)'stroke'
      write(File_Number,*)'0.001 setlinewidth'
      
      !write(File_Number,*)dble( Translation_X ), dble( Translation_Y ), 'moveto'
      !write(File_Number,*)dble( Translation_X ) +WidthPS_X, dble( Translation_Y ) +0D0, 'lineto'
      !write(File_Number,*)dble( Translation_X ) +WidthPS_X, dble( Translation_Y ) +WidthPS_Y, 'lineto'
      !write(File_Number,*)dble( Translation_X ) +0D0, dble( Translation_Y ) +WidthPS_Y, 'lineto'
      !write(File_Number,*)dble( Translation_X ), dble( Translation_Y ), 'lineto'
      !write(File_Number,*)'stroke'
    
      write(File_Number,*)'0.3 setgray'

      Counter_Percent= 1 

      do e= 1, Number_Element
         
         if( e== Number_Element/100*Counter_Percent )then
            !write(*,'(a,i3,a)', advance='no')'[', Counter_Percent, '%]'
            write(*,'(a,i3,a,$)')'[', Counter_Percent, '%]'
            Counter_Percent= Counter_Percent +1
         end if 

         2244 continue

         if( Flag_Plot_Element( e )==3 )then

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
            call Output_Error( 'Plot_Result_Element_Triangle', 2220 ) 
         end if

      end do

      deallocate( Flag_Plot_Element )

      !==================================================================================================
      ! Plot Fixed Domain
      !==================================================================================================

      !if( Flag_FixedDomain==1 .and. Flag_PEC_Boundary_Condition==0 )then
      if( Flag_FixedDomain==1 )then
         write( File_Number, * )'1 setlinewidth'
         write( File_Number, * ) 'newpath'
         write( File_Number, * ) Position_Center_FixedDomain_X_PS, Position_Center_FixedDomain_Y_PS, &
                         Radius_FixedDomain_Plot_PS,' 0 360 arc'
         write( File_Number, * )'closepath'
         write( File_Number, * )'stroke'

      else if( Flag_FixedDomain==2 .and. Flag_Thermal_Device/=210 )then
         write( File_Number, * )'1 setlinewidth'
         write( File_Number, * ) Position_Corner_FixedDomain_X1_PS, Position_Corner_FixedDomain_Y1_PS, 'moveto' 
         write( File_Number, * ) Position_Corner_FixedDomain_X2_PS, Position_Corner_FixedDomain_Y1_PS, 'lineto' 
         write( File_Number, * ) Position_Corner_FixedDomain_X2_PS, Position_Corner_FixedDomain_Y2_PS, 'lineto' 
         write( File_Number, * ) Position_Corner_FixedDomain_X1_PS, Position_Corner_FixedDomain_Y2_PS, 'lineto' 
         write( File_Number, * ) Position_Corner_FixedDomain_X1_PS, Position_Corner_FixedDomain_Y1_PS, 'lineto' 
         write( File_Number, * )'stroke'

         !if( Flag_Thermal_Device==50 .or. Flag_Thermal_Device==51 )then
         !    write( File_Number, * )'1 setlinewidth'
         !    write( File_Number, * ) 'newpath'
         !    write( File_Number, * ) ( Position_Corner_FixedDomain_X1_PS +Position_Corner_FixedDomain_X2_PS )/2d0, &
         !                  Position_Corner_FixedDomain_Y2_PS, & 
         !                  Radius_FixedDomain_Plot_PS,' 0 180 arc'
         !    write( File_Number, * )'closepath'
         !    write( File_Number, * )'stroke'
         !end if

      end if

      !==================================================================================================
      ! Plot Design Domain
      !==================================================================================================

      if( Flag_DesignDomain==1 )then

         if( Flag_Thermal_Device==0 )then 
            write( File_Number, * )'1 setlinewidth'
            write( File_Number, * ) 'newpath'
            write( File_Number, * ) Position_Center_DesignDomain_X_PS, Position_Center_DesignDomain_Y_PS, &
                            Radius_DesignDomain_Plot_PS,' 0 360 arc'
            write( File_Number, * )'closepath'
            write( File_Number, * )'stroke'
         end if

      !else if( Flag_DesignDomain==2 )then
      !   write( File_Number, * )'1 setlinewidth'
      !   write( File_Number, * ) Position_Corner_DesignDomain_X1_PS, Position_Corner_DesignDomain_Y1_PS, 'moveto' 
      !   write( File_Number, * ) Position_Corner_DesignDomain_X2_PS, Position_Corner_DesignDomain_Y1_PS, 'lineto' 
      !   write( File_Number, * ) Position_Corner_DesignDomain_X2_PS, Position_Corner_DesignDomain_Y2_PS, 'lineto' 
      !   write( File_Number, * ) Position_Corner_DesignDomain_X1_PS, Position_Corner_DesignDomain_Y2_PS, 'lineto' 
      !   write( File_Number, * ) Position_Corner_DesignDomain_X1_PS, Position_Corner_DesignDomain_Y1_PS, 'lineto' 
      !   write( File_Number, * )'stroke'
      end if

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

      write(*,*)' '
      write(*,*) File_Number
      write(*,*)'         ==> Temperature_Distribution_', Optimization_Step_Character, '.ps' 
     
   if( Flag_Color_Var==1 .or. Flag_Color_Var==2 )then
      deallocate( Position_Maximum_Color_Var )
   end if

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

   2272 continue

   return
end subroutine Plot_Result_Element_Triangle



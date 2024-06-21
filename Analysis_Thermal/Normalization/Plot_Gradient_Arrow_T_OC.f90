
subroutine Plot_Gradient_Arrow_T_OC  &
     ( FN, &
      Flag_Thermal_Insulation_BC, Physical_Quantity_Plotted, Thermal_Conductivity_Element, &
      Number_Level, Optimization_Step,  &
      Position_Node, Number_Node, &
      Index_Element_2_Node, Number_Element, &
      Class_Element, ID_Element_Structure, ID_Element_FixedDomain, &
      Objective_Function ) 

  !use omp_lib
  use Parameters
  implicit none

  integer, intent(in) :: FN, Flag_Thermal_Insulation_BC
  integer, intent(in) :: Number_Node, Number_Element
  double precision, intent(inout) :: Physical_Quantity_Plotted( Number_Node )
  double precision, intent(in) :: Thermal_Conductivity_Element( Number_Element )

  integer, intent(in) :: Number_Level
  integer, intent(in) :: Optimization_Step

  double precision, intent(in) :: Position_Node( 2, Number_Node )
  integer, intent(in) :: Index_Element_2_Node( 4, Number_Element )
  integer, intent(in) :: Class_Element( Number_Element )
  integer, intent(in) :: ID_Element_Structure, ID_Element_FixedDomain
  
  double precision, intent(in) :: Objective_Function( 0:Number_Optimization_Step ) 
 
  character(len=Length_Character_Optimization_Step) :: Optimization_Step_Character
  character(len=8) :: Format_Filenumber
  
  character(len=256) :: Filename_Result

  double precision, allocatable, dimension(:) :: Position_Minimum, Position_Minimum_Plot 
  double precision, allocatable, dimension(:) :: Position_Maximum, Position_Maximum_Plot
  double precision, allocatable, dimension(:) :: Translation_Position, Width_PostScript 
  double precision, allocatable, dimension(:) :: Translate_Color_Var 
  double precision, allocatable, dimension(:,:,:) :: Position_Maximum_Color_Var 

  double precision :: Minimum_Value_Plotted, Maximum_Value_Plotted
  double precision, allocatable, dimension(:,:) :: Value_Boundary_Level 
  double precision, allocatable, dimension(:,:) :: RGB
  double precision, allocatable, dimension(:) :: RGB_Lowest, RGB_Highest

  integer :: e, i, j, k, l, gpno, Element_Number_include_Arrow
  double precision, allocatable, dimension(:,:,:) :: Position_Node_in_Element, Position_Node_Plot_Element
  double precision, allocatable, dimension(:,:) :: Position_Node_Plot 
  double precision, allocatable, dimension(:,:) :: Value_OF_Element 

  double precision, allocatable, dimension(:,:,:,:) :: Difference_Position_Node_in_Element
  double precision, allocatable, dimension(:) :: Area_Element 
  double precision, allocatable, dimension(:,:) :: BasisFunction_B, BasisFunction_C

  double precision, allocatable, dimension(:,:) :: Gradient_in_Element 
  double precision, allocatable, dimension(:) :: Abs_Heat_Flux
  double precision, allocatable, dimension(:) :: Value_Plotted

  integer, allocatable, dimension(:) :: Level_Element

  integer :: Counter_Percent

  integer, allocatable, dimension(:,:) :: Element_Number_Arrow
  integer, allocatable, dimension(:) :: Number_Grid
  double precision, allocatable, dimension(:,:,:) :: Position_Arrow, Direction_Arrow
  double precision, allocatable, dimension(:) :: Position_Arrow_tmp

  double precision, allocatable, dimension(:) :: Number_Line, Counter_Line
  double precision :: Line_Length_X, Line_Length_X_Default, Ratio_Default, Vector_Length


  double precision :: Min_Size_Arrow, Max_Size_Arrow, Base_Size 
  integer, allocatable, dimension(:,:) :: Level_Arrow

  !if( Flag_Thermal_Device==30 ) Structure='asymmetry'

  write( Format_Filenumber, '(A2,I1,A1,I1,A1)' ) & 
  '(i', Length_Character_Optimization_Step,'.',Length_Character_Optimization_Step,')'
  
  write( Optimization_Step_Character, Format_Filenumber ) Optimization_Step
 
  Filename_Result= trim( "Heat_Flux_Arrow_OC_"//trim(Optimization_Step_Character)//".ps" )

  !================================================
  write(*,*)'===================================================='
  write(*,*)'    call Plot_Gradient_Arrow_T_OC'
  write(*,*)'===================================================='
  !================================================

  !===================================================================================
  write(*,*)'    Position_Node --> Position_Node_in_Element '
  !===================================================================================
  allocate( Position_Node_in_Element( 2, 3, Number_Element ) )
  
  do e= 1, Number_Element
    do i= 1, 3
      do j= 1, 2
        Position_Node_in_Element( j, i, e )= Position_Node( j, Index_Element_2_Node( i, e ) )
      end do
    end do
  end do

  !===================================================================================
  write(*,*)'    Physical_Quantity_Plotted --> Value_OF_Element '
  !===================================================================================
  allocate( Value_OF_Element( 3, Number_Element ) )

  do e= 1, Number_Element
    do i= 1, 3
      Value_OF_Element( i, e ) = Physical_Quantity_Plotted( Index_Element_2_Node( i, e ) )
    end do
  end do

  do e= 1, Number_Element
    do i= 1, 3
      if( Normal_Temperature_Lower_Side_Boundary > Value_OF_Element( i, e ) .or. &
          Normal_Temperature_Higher_Side_Boundary < Value_OF_Element( i, e ) )then 

        write(*,*)'====================================================='
        write(*,*)'Value_OF_Element( i, e )=', Value_OF_Element( i, e )
        write(*,*)'e=', e, 'i=', i
        call Output_Error( 'Plot_Gradient_Arrow_T_OC', 118 )
        write(*,*)'====================================================='
      end if
    end do
  end do

  allocate( Difference_Position_Node_in_Element( 2, 3, 3, Number_Element ) )

  do e= 1, Number_Element
    do i= 1, 3
      do j= 1, 3
        do k= 1, 2
          Difference_Position_Node_in_Element( k, j, i, e ) &
          =Position_Node_in_Element( k, j, e ) -Position_Node_in_Element( k, i, e )
        end do
      end do
    end do
  end do

  !===================================================================================
  write(*,*)'    Compute Basis Function '
  !===================================================================================
  allocate( Area_Element( Number_Element ) )

  Area_Element( : ) &
  =abs( Difference_Position_Node_in_Element( 1, 1, 3, : ) & 
       *Difference_Position_Node_in_Element( 2, 2, 3, : ) &
       -Difference_Position_Node_in_Element( 1, 2, 3, : ) & 
       *Difference_Position_Node_in_Element( 2, 1, 3, : ) )/2.0d0

  allocate( BasisFunction_B( 3, Number_Element ) )
  allocate( BasisFunction_C( 3, Number_Element ) )

  BasisFunction_B( 1, : ) &
  = Difference_Position_Node_in_Element( 2, 2, 3, : )/( 2.0d0 *Area_Element( : ) )

  BasisFunction_B( 2, : ) & 
  = Difference_Position_Node_in_Element( 2, 3, 1, : )/( 2.0d0 *Area_Element( : ) )

  BasisFunction_B( 3, : ) & 
  = Difference_Position_Node_in_Element( 2, 1, 2, : )/( 2.0d0 *Area_Element( : ) )


  BasisFunction_C( 1, : ) & 
  = Difference_Position_Node_in_Element( 1, 3, 2, : )/( 2.0d0 *Area_Element( : ) )

  BasisFunction_C( 2, : ) & 
  = Difference_Position_Node_in_Element( 1, 1, 3, : )/( 2.0d0 *Area_Element( : ) )

  BasisFunction_C( 3, : ) & 
  = Difference_Position_Node_in_Element( 1, 2, 1, : )/( 2.0d0 *Area_Element( : ) )

  deallocate( Difference_Position_Node_in_Element )
  deallocate( Area_Element )

  do e = 1, Number_Element
     write(173,*) BasisFunction_B( 1, e ), BasisFunction_B( 2, e ), BasisFunction_B( 3, e )
     write(174,*) BasisFunction_C( 1, e ), BasisFunction_C( 2, e ), BasisFunction_C( 3, e )
  enddo

  !===================================================================================
  write(*,*)'    Compute Gradient of Physical Quantity in Element '
  !===================================================================================
  allocate( Gradient_in_Element( 2, Number_Element ) )
  Gradient_in_Element= 0.0d0

  do e = 1, Number_Element
    do i = 1, 3 
      Gradient_in_Element( 1, e ) = Gradient_in_Element( 1, e ) +BasisFunction_B( i, e ) *Value_OF_Element( i, e )
      Gradient_in_Element( 2, e ) = Gradient_in_Element( 2, e ) +BasisFunction_C( i, e ) *Value_OF_Element( i, e )
    end do
  end do

  !===================================================================================
  write(*,*)'    Compute Absolute Value of Heat Flux '
  !===================================================================================
  allocate( Abs_Heat_Flux( Number_Element ) )

  do e = 1, Number_Element
    Abs_Heat_Flux( e ) = sqrt( ( -Gradient_in_Element( 1, e ) )**2 +( -Gradient_in_Element( 2, e ) )**2 )
  end do

  deallocate( BasisFunction_B )
  deallocate( BasisFunction_C )


  allocate( Value_Plotted( Number_Element ) )

  do e = 1, Number_Element
    Value_Plotted( e ) = Abs_Heat_Flux( e ) *Thermal_Conductivity_Element( e )&
                         /( Thermal_Conductivity_OpenRegion*0.25d0 ) ! Abs_Heat_Flux_Fe = 0.25
    !Value_Plotted( e ) = Thermal_Conductivity_Element( e )
    write(193,*) Thermal_Conductivity_Element( e ), Abs_Heat_Flux( e ), Value_Plotted( e )
    write(194,*) Gradient_in_Element( 1, e ), Gradient_in_Element( 2, e )
    write(195,*) Value_OF_Element( 1, e ), Value_OF_Element( 2, e ), Value_OF_Element( 3, e ) 
  end do
  write(196,*) minval( Abs_Heat_Flux ), maxval( Abs_Heat_Flux )

  deallocate( Value_OF_Element )
  deallocate( Abs_Heat_Flux )

  !======================================================================================================================
  write(*,*)'    Decide the Frame of Position ' 
  !======================================================================================================================
  allocate( Position_Minimum( 2 ) )
  allocate( Position_Maximum( 2 ) )
 
  !Minimum
  Position_Minimum( 1 )= minval( Position_Node( 1, : ) ) 
  Position_Minimum( 2 )= minval( Position_Node( 2, : ) )
  !Maximum
  Position_Maximum( 1 )= maxval( Position_Node( 1, : ) ) 
  Position_Maximum( 2 )= maxval( Position_Node( 2, : ) ) 
 
  allocate( Translation_Position( 2 ) )

  Translation_Position( 1 )= Translation_X
  Translation_Position( 2 )= Translation_Y

  allocate( Width_PostScript( 2 ) )

  Width_PostScript( 1 )= WidthPS_X 
  Width_PostScript( 2 )= WidthPS_Y

  !======================================================================================================================
  write(*,*)'    Position_Node --> Position_Node_Plot' 
  !======================================================================================================================
  allocate( Position_Node_Plot( 2, Number_Node ) )
  
  do i= 1, Number_Node
    do j= 1, 2 
      Position_Node_Plot( j, i ) &
      =  ( Position_Node( j, i ) -Position_Minimum( j ) )  &
        /( Position_Maximum( j ) -Position_Minimum( j ) )*dble( Width_PostScript( j ) ) 
    end do
  end do
 
  do i= 1, Number_Node
    do j= 1, 2 
      Position_Node_Plot( j, i )= Position_Node_Plot( j, i ) +Translation_Position( j ) 
    end do
  end do

  !===================================================================================
  write(*,*)'   Position_Node_Plot  --> Position_Node_Plot_Element '
  !===================================================================================
  allocate( Position_Node_Plot_Element( 2, 3, Number_Element ) )
  
  do e= 1, Number_Element
    do i= 1, 3
      do j= 1, 2
        Position_Node_Plot_Element( j, i, e )= Position_Node_Plot( j, Index_Element_2_Node( i, e ) )
      end do
    end do
  end do

  !===================================================================================
  write(*,*)'   Position_Node_Plot  --> Position_Minimum_Plot & Position_Maximum_Plot '
  !===================================================================================
  allocate( Position_Minimum_Plot( 2 ) )
  allocate( Position_Maximum_Plot( 2 ) )
 
  !Minimum
  Position_Minimum_Plot( 1 )= minval( Position_Node_Plot( 1, : ) ) 
  Position_Minimum_Plot( 2 )= minval( Position_Node_Plot( 2, : ) )
  !Maximum
  Position_Maximum_Plot( 1 )= maxval( Position_Node_Plot( 1, : ) ) 
  Position_Maximum_Plot( 2 )= maxval( Position_Node_Plot( 2, : ) ) 

  deallocate( Position_Node_Plot )

  !================================================
  ! Level
  !================================================

  if( Flag_Range_Value_Plot_Gradient==0 )then
    Minimum_Value_Plotted= minval( Value_Plotted ) 
    Maximum_Value_Plotted= maxval( Value_Plotted )  
  else if( Flag_Range_Value_Plot_Gradient==1 )then
    Maximum_Value_Plotted= Maximum_Value_Plotted_Gradient 
    Minimum_Value_Plotted= Minimum_Value_Plotted_Gradient
  end if

  write(264,*) Minimum_Value_Plotted, Maximum_Value_Plotted

  allocate( Value_Boundary_Level( 2, Number_Level ) )
 
  do i= 1, Number_Level
    Value_Boundary_Level( 1, i ) &
    = ( Maximum_Value_Plotted -Minimum_Value_Plotted )/ Number_Level*( i-1 ) !+Minimum_Value_Plotted
    Value_Boundary_Level( 2, i ) &
    = ( Maximum_Value_Plotted -Minimum_Value_Plotted )/ Number_Level*( i )  !+Minimum_Value_Plotted 
  end do 

  Value_Plotted( : )= Value_Plotted( : ) -Minimum_Value_Plotted

  allocate( Level_Element( Number_Element ) )

  do j= 1, Number_Element 
    Level_Element( j )= Integer_Initialization

    if( Value_Plotted( j ) < 0.0d0 )then
      Level_Element( j )= 0
    else if( Value_Plotted( j ) >= Maximum_Value_Plotted -Minimum_Value_Plotted )then
      Level_Element( j )= Number_Level +1
    else
      do i= 1, Number_Level
        if( Value_Plotted( j ) >= Value_Boundary_Level( 1, i ) .and. & 
            Value_Plotted( j ) < Value_Boundary_Level( 2, i ) )then
   
          Level_Element( j )= i  
        end if 
      end do
    end if
  end do

  deallocate( Value_Boundary_Level )

  !======================================================================================================================
  ! RGB 
  !======================================================================================================================
  allocate( RGB( 3, Number_Level ) )
  allocate( RGB_Lowest( 3 ) )
  allocate( RGB_Highest( 3 ) )

  call Set_RGB( Flag_RGB_Color_Gradient, Number_Level, RGB, RGB_Lowest, RGB_Highest )

  !do i= 1, Number_Level
  !  write(901,*) i, RGB( 1, i ) 
  !  write(902,*) i, RGB( 2, i )
  !  write(903,*) i, RGB( 3, i )
  !end do



  !==================================================================================================
  write(*,*) '         Plot Element Data ..........'
  !==================================================================================================

  open( FN, file=Filename_Result, status='replace' )
  write(FN,*)'%!PS-Adobe-2.0'
  write(FN,*)'%BoundingBox: 0 0', Width_PostScript( 1 ), Width_PostScript( 2 )
  write(FN,*)'%EndHeader'
  write(FN,*)'stroke'
  write(FN,*)'0.0001 setlinewidth'
  write(FN,*)'0.3 setgray'

  Counter_Percent= 1 

  do e= 1, Number_Element
         
    if( e== Number_Element/100*Counter_Percent )then
      !write(*,'(a,i3,a)', advance='no')'[', Counter_Percent, '%]'
      write(*,'(a,i3,a,$)')'[', Counter_Percent, '%]'
      Counter_Percent= Counter_Percent +1
    end if 

    write(FN,*)'gsave'
    write(FN,*) Position_Node_Plot_Element( 1, 1, e ), Position_Node_Plot_Element( 2, 1, e ), 'moveto'
    write(FN,*) Position_Node_Plot_Element( 1, 2, e ), Position_Node_Plot_Element( 2, 2, e ), 'lineto'
    write(FN,*) Position_Node_Plot_Element( 1, 3, e ), Position_Node_Plot_Element( 2, 3, e ), 'lineto'
    write(FN,*) Position_Node_Plot_Element( 1, 1, e ), Position_Node_Plot_Element( 2, 1, e ), 'lineto'
    write(FN,*)'closepath'
    if( 1 <= Level_Element( e ) .and. Level_Element( e ) <= 100 )then
      write(FN,*) RGB( 1, Level_Element( e ) ), RGB( 2, Level_Element( e ) ), RGB( 3, Level_Element( e ) ), 'setrgbcolor' 
    else if( Level_Element( e )==0  )then
      write(FN,*) RGB_Lowest( 1 ), RGB_Lowest( 2 ), RGB_Lowest( 3 ), 'setrgbcolor' 
    else if( Level_Element( e )==101  )then
      write(FN,*) RGB_Highest( 1 ), RGB_Highest( 2 ), RGB_Highest( 3 ), 'setrgbcolor' 
    end if
    write(FN,*)'fill'
    write(FN,*)'grestore'
    write(FN,*)'newpath'

  end do

  deallocate( RGB_Lowest )
  deallocate( RGB_Highest )

   !======================================================================================================================
   ! Plot Fixed Domain & Design Domain 
   !======================================================================================================================
   write( FN, * )'0 0 0 setrgbcolor'
         !write( FN, * )'1 1 1 setrgbcolor'
   write( FN, * )'1 setlinewidth'
   write( FN, * ) 'newpath'
   write( FN, * ) Position_Center_FixedDomain_X_PS, Position_Center_FixedDomain_Y_PS, &
                         Radius_FixedDomain_Plot_PS,' 0 360 arc'
   write( FN, * )'closepath'
   write( FN, * )'stroke'

   write( FN, * )'0 0 0 setrgbcolor'
   write( FN, * )'1 setlinewidth'
   write( FN, * ) 'newpath'
   write( FN, * ) Position_Center_DesignDomain_X_PS, Position_Center_DesignDomain_Y_PS, &
                  Radius_DesignDomain_Plot_PS,' 0 360 arc'
   write( FN, * )'closepath'
   write( FN, * )'stroke'
   

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

  deallocate( Width_PostScript )

   !==================================================================================================
   write(*,*)'Plot Color Var'
   !==================================================================================================
   if( Flag_Color_Var==1 )then
      call Plot_Color_Var_Right_Side & 
         ( FN, RGB, Number_Level, Position_Maximum_Color_Var, Minimum_Value_Plotted, Maximum_Value_Plotted )
   else if( Flag_Color_Var==2 )then
      call Plot_Color_Var_Upper_Side & 
         ( FN, RGB, Number_Level, Position_Maximum_Color_Var, Minimum_Value_Plotted, Maximum_Value_Plotted )
   end if


  deallocate( RGB )

  deallocate( Position_Maximum_Color_Var )


  deallocate( Position_Minimum )
  deallocate( Position_Maximum )
  deallocate( Translation_Position )

   !==================================================================================================
   write(*,*)'call Plot Arrow '
   !==================================================================================================

   allocate( Number_Grid( 2 ) )
   Number_Grid( 1 )= 40
   Number_Grid( 2 )= Number_Grid( 1 )*1.5d0

   allocate( Element_Number_Arrow( Number_Grid( 1 ), Number_Grid( 2 ) ) )
   allocate( Position_Arrow( 2, Number_Grid( 1 ), Number_Grid( 2 ) ) )

   do i= 1, Number_Grid( 2 )
      do j= 1, Number_Grid( 1 )
         Position_Arrow( 1, j, i ) &
         = ( Position_Maximum_Plot( 1 ) -Position_Minimum_Plot( 1 ) )/Number_Grid( 1 ) *( j-0.5 ) &
            +Position_Minimum_Plot( 1 ) +1d-8
         Position_Arrow( 2, j, i ) &
         = ( Position_Maximum_Plot( 2 ) -Position_Minimum_Plot( 2 ) )/Number_Grid( 2 ) *( i-0.5 ) &
            +Position_Minimum_Plot( 2 ) -1d-8
         !Position_Arrow( 1, j, i ) &
         != ( Position_Maximum_Plot( 1 ) -Position_Minimum_Plot( 1 ) )/Number_Grid( 1 ) *( j-0.5 ) &
         !   +Position_Minimum_Plot( 1 ) 
         !Position_Arrow( 2, j, i ) &
         != ( Position_Maximum_Plot( 2 ) -Position_Minimum_Plot( 2 ) )/Number_Grid( 2 ) *( i-0.5 ) &
         !   +Position_Minimum_Plot( 2 )
      end do
   end do

   allocate( Position_Arrow_tmp( 2 ) )
   allocate( Direction_Arrow( 2, Number_Grid( 1 ), Number_Grid( 2 ) ) )

   do i= 1, Number_Grid( 2 )
      do j= 1, Number_Grid( 1 )
         do k= 1, 2
            Position_Arrow_tmp( k )= Position_Arrow( k, j, i )
         end do

         call Detect_Element_inwhich_Point_Exists &
         ( Number_Element, Position_Node_Plot_Element, Position_Arrow_tmp, &
           Element_Number_include_Arrow )

         Element_Number_Arrow( j, i )= Element_Number_include_Arrow

         do k= 1, 2
            Direction_Arrow( k, j, i )= -Gradient_in_Element( k, Element_Number_include_Arrow )
         end do

      end do
   end do

   allocate( Level_Arrow( Number_Grid( 1 ), Number_Grid( 2 ) ) )

   do i= 1, Number_Grid( 2 )
      do j= 1, Number_Grid( 1 )
         Level_Arrow( j, i )= Level_Element( Element_Number_Arrow( j, i ) )
      end do
   end do

   Base_Size= ( Position_Maximum_Plot( 1 ) -Position_Minimum_Plot( 1 ) )/Number_Grid( 1 )
   Max_Size_Arrow= Base_Size*2d0
   Min_Size_Arrow= Base_Size/10.0d0 

   do i= 1, Number_Grid( 2 )
      do j= 1, Number_Grid( 1 )
         call Plot_Single_Arrow&
              ( FN, 'w', 1, &
                Direction_Arrow( 1, j, i ), Direction_Arrow( 2, j, i ), &
                Position_Arrow( 1, j, i ), Position_Arrow( 2, j, i ), Level_Arrow( j, i ), 100, &
                Min_Size_Arrow, Max_Size_Arrow, 30 ) 
      end do
   end do

   deallocate( Direction_Arrow )
   deallocate( Level_Element )
   deallocate( Level_Arrow )
   deallocate( Element_Number_Arrow )
   deallocate( Position_Arrow )
   deallocate( Number_Grid )
   deallocate( Position_Arrow_tmp )

   !write( File_Number, * ) '# plot heat flux arrow'
   !do e= 1, Number_PV
   !   call Plot_Arrow_Postscript&
   !      ( File_Number, 2, 'w', &
   !        Position_Arrow_PV_Plot( 1, 1, e ), Position_Arrow_PV_Plot( 2, 1, e ), &
   !        Position_Arrow_PV_Plot( 1, 2, e ), Position_Arrow_PV_Plot( 2, 2, e ), & 
   !        Stemthick_Plot( e ), Headthick_Plot( e ), Headlength_Plot( e ) )
   !        !Stemthick*Level_Flux_tmp( e ), Headthick*Level_Flux_tmp( e ), Headlength*Level_Flux_tmp( e ) )
   !        !Stemthick, Headthick, Headlength )
   !end do



  !==================================================================================================
  write(*,*)'Plot Stream Line' 
  !==================================================================================================

!  Number_Line_X= 362 !362 242 
!  Number_Line_Y= 24 
!  !Number_Line_Y= 24 
!  Ratio_Default= 1d0
!  Structure='asymmetry'
!  !Structure='symmetry'
!
!  if( Structure=='symmetry' )then
!     call Plot_Stream_Line_Symmetry &
!        ( FN, Number_Element, Number_Line_X, Number_Line_Y, Ratio_Default, & 
!          Position_Minimum_Plot, Position_Maximum_Plot, &
!          Position_Node_Plot_Element, Gradient_in_Element, & 
!          Optimization_Step )
!
!  else if( Structure=='asymmetry' )then
!     call Plot_Stream_Line &
!        ( FN, Number_Element, Number_Line_X, Number_Line_Y, Ratio_Default, & 
!          Position_Minimum_Plot, Position_Maximum_Plot, &
!          Position_Node_Plot_Element, Gradient_in_Element, & 
!          Optimization_Step )
!  else
!     write(*,*)'Structure=', Structure
!     call Output_Error( 'Plot_Gradient_Arrow_T_OC', 558 )
!  end if

  deallocate( Gradient_in_Element )
  deallocate( Position_Node_Plot_Element )

  deallocate( Position_Minimum_Plot )
  deallocate( Position_Maximum_Plot )
 
  2272 continue

  write(*,*) '     end subroutine Plot_Gradient_Arrow_T_OC' 
  return
end subroutine Plot_Gradient_Arrow_T_OC



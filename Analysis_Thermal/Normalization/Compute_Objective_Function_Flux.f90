
subroutine Compute_Objective_Function_Flux &
       ( Flag_Structure, ID_Element, &
         Value_Objective_Function, Number_Node_Value_Objective_Function, &
         Position_Node, Number_Node, &
         Index_Element_2_Node, Class_Element, Number_Element, &
         Thermal_Conductivity_Evaluated, & 
         Value_ObjectiveFunction_Normalize, &
         !=================================================================
         Objective_Function, Objective_Function_Value )

  !$use onp_lib
  use Parameters
  implicit none

  integer, intent(in) :: Flag_Structure, ID_Element 
  integer, intent(in) :: Number_Node_Value_Objective_Function 
  integer, intent(in) :: Number_Node, Number_Element

  !complex( kind( 0d0 ) ), intent(in) :: Value_Objective_Function( Number_Node_Value_Objective_Function )
  double precision, intent(in) :: Value_Objective_Function( Number_Node_Value_Objective_Function )

  double precision, intent(in) :: Position_Node( 2, Number_Node )
  integer, intent(in) :: Index_Element_2_Node( 3, Number_Element ) 
  integer, intent(in) :: Class_Element( Number_Element ) 

  double precision, intent(in) :: Value_ObjectiveFunction_Normalize
  double precision, intent(in) :: Thermal_Conductivity_Evaluated 

  double precision, intent(out) :: Objective_Function, Objective_Function_Value

  integer :: e, i, j, k
  double precision, allocatable, dimension(:,:,:) :: Position_Node_Element_Triangle

  integer, allocatable, dimension(:) :: Element_Number_on_Integrated_Line
  integer, allocatable, dimension(:) :: Counter_Node_in_Element_on_Integrated_Line
  integer, allocatable, dimension(:) :: Flag_Element_Number_on_Integrated_Line
  integer, allocatable, dimension(:,:) :: Local_Node_tmp, Local_Node 
  integer :: Number_Element_on_Integrated_Line, Counter_tmp

  double precision, allocatable, dimension(:) :: Edge_Legnth_Integrated_Element 
  double precision, allocatable, dimension(:,:,:) :: Position_Node_Element_on_Line
  double precision, allocatable, dimension(:,:,:,:) :: Difference_Position_Node_Element_on_Line
  double precision, allocatable, dimension(:) :: Area_Element_on_Line 
  double precision, allocatable, dimension(:,:) :: BasisFunction_B_on_Line, BasisFunction_C_on_Line
  double precision, allocatable, dimension(:,:) :: Value_OF_Element_Line

  double precision :: Unit_Normal_Vector( 2 ) !=( -1.0d0, 0.0d0 ) ! for negative x direction

  if( Flag_Thermal_Device==20 .or. Flag_Thermal_Device==21 )then
    Unit_Normal_Vector( 1 ) = -1.0d0
    Unit_Normal_Vector( 2 ) = 0.0d0
  else if( Flag_Thermal_Device==30 )then
    Unit_Normal_Vector( 1 ) = 1.0d0
    Unit_Normal_Vector( 2 ) = 0.0d0
  end if

  !=======================================================================================================
  write(*,*)'   call Compute_Objective_Function_Flux ' 
  !=======================================================================================================

  !===================================================================================
  write(*,*)'      Position_Node --> Position_Node_Element_Triangle '
  !===================================================================================
  allocate( Position_Node_Element_Triangle( 2, 3, Number_Element ) )
  
  do e= 1, Number_Element
    do i= 1, 3
      do j= 1, 2
        Position_Node_Element_Triangle( j, i, e )= Position_Node( j, Index_Element_2_Node( i, e ) )
      end do
    end do
  end do

  !=======================================================================================================
  write(*,*)'      Detect Element on Integrated Line ' 
  !=======================================================================================================

  allocate( Counter_Node_in_Element_on_Integrated_Line( Number_Element ) ) 
  allocate( Flag_Element_Number_on_Integrated_Line( Number_Element ) ) 
  allocate( Local_Node_tmp( 2, Number_Element ) ) 
  Number_Element_on_Integrated_Line = 0
  Counter_Node_in_Element_on_Integrated_Line = 0
  Counter_tmp = 0

  do e= 1, Number_Element
    if( Class_Element( e )==ID_Element )then
      Counter_tmp = Counter_tmp +1

      do i= 1, 3
        ! Detect the Elements including 2 Nodes on the Integrated Line
        if( 0.0d0 -1d-8 <= Position_Node_Element_Triangle( 1, i, e ) .and. & 
            Position_Node_Element_Triangle( 1, i, e ) <= 0.0d0 +1d-8 )then 

          Counter_Node_in_Element_on_Integrated_Line( e ) = Counter_Node_in_Element_on_Integrated_Line( e ) +1

          Local_Node_tmp( Counter_Node_in_Element_on_Integrated_Line( e ), e ) = i
        end if
      end do

      if( Counter_Node_in_Element_on_Integrated_Line( e ) == 2 )then
        Number_Element_on_Integrated_Line = Number_Element_on_Integrated_Line +1
      else if( Counter_Node_in_Element_on_Integrated_Line( e ) > 2 )then
        write(*,*)'e=', e 
        write(*,*)'Counter_Node_in_Element_on_Integrated_Line( e )=', Counter_Node_in_Element_on_Integrated_Line( e ) 
        call Output_Error( 'Compute_Objective_Function_Flux', 103 )
      end if
    end if
  end do

  if( Counter_tmp  == 0 )then
    write(*,*)'======================================================================='
    write(*,*)'Counter_tmp=', Counter_tmp 
    call Output_Error( 'Compute_Objective_Function_Flux', 109 )
  end if

  do e= 1, Number_Element
    if( Counter_Node_in_Element_on_Integrated_Line( e ) == 2 )then
      do i= 1, 2
        if( Local_Node_tmp( i, e ) < 1 .or. 3 < Local_Node_tmp( i, e ) )then
          write(*,*)'======================================================================='
          write(*,*)'Local_Node_tmp( i, e )=', Local_Node_tmp( i, e )
          write(*,*)'e=', e, 'i=', i
          call Output_Error( 'Compute_Objective_Function_Flux', 117 ) 
        end if
      end do
    end if
  end do

  allocate( Element_Number_on_Integrated_Line( Number_Element_on_Integrated_Line ) ) 
  allocate( Local_Node( 2, Number_Element_on_Integrated_Line ) ) 
  Counter_tmp = 0

  do e= 1, Number_Element
    if( Counter_Node_in_Element_on_Integrated_Line( e ) == 2 )then
      Counter_tmp = Counter_tmp +1
      Element_Number_on_Integrated_Line( Counter_tmp ) = e
      do i= 1, 2
        Local_Node( i, Counter_tmp ) = Local_Node_tmp( i, e )
      end do
    end if
  end do

  if( Number_Element_on_Integrated_Line /= Counter_tmp )then
    write(*,*)'======================================================================='
    write(*,*)'Number_Element_on_Integrated_Line=', Number_Element_on_Integrated_Line
    write(*,*)'Counter_tmp=', Counter_tmp
    call Output_Error( 'Compute_Objective_Function_Flux', 123 )
  end if

  if( Number_Element_on_Integrated_Line == 0 )then
    write(*,*)'======================================================================='
    write(*,*)'Number_Element_on_Integrated_Line=', Number_Element_on_Integrated_Line
    write(*,*)'Counter_tmp=', Counter_tmp
    call Output_Error( 'Compute_Objective_Function_Flux', 139 )
  end if


  do e= 1, Number_Element_on_Integrated_Line
    do i= 1, 2
      if( Local_Node( i, e ) < 1 .or. 3 < Local_Node( i, e ) )then
        write(*,*)'======================================================================='
        write(*,*)'Local_Node( i, e )=', Local_Node( i, e )
        write(*,*)'e=', e, 'i=', i
        call Output_Error( 'Compute_Objective_Function_Flux', 137 ) 
      end if
    end do
  end do

!  do e= 1, Number_Element_on_Integrated_Line
!      do i= 1, 2
!        write(174,*) Position_Node_Element_Triangle( 1, Local_Node( i, e ), Element_Number_on_Integrated_Line( e ) ), &
!                     Position_Node_Element_Triangle( 2, Local_Node( i, e ), Element_Number_on_Integrated_Line( e ) )
!      end do 
!      write(174,*)' ' 
!  end do
!
!  do e= 1, Number_Element
!    if( Counter_Node_in_Element_on_Integrated_Line( e ) == 2 )then
!      do i= 1, 2
!        write(173,*) Position_Node_Element_Triangle( 1, Local_Node_tmp( i, e ), e ), &
!                     Position_Node_Element_Triangle( 2, Local_Node_tmp( i, e ), e )
!      end do 
!      write(173,*)' ' 
!    end if
!  end do
!
!  do e= 1, Number_Element_on_Integrated_Line
!    do i= 1, 2
!      write(172,*) Position_Node_Element_Triangle( 1, Local_Node( i, e ), Element_Number_on_Integrated_Line( e ) ), &
!                   Position_Node_Element_Triangle( 2, Local_Node( i, e ), Element_Number_on_Integrated_Line( e ) )
!    end do 
!    write(172,*)' ' 
!  end do
!
!  do e= 1, Number_Element_on_Integrated_Line
!      do i= 1, 3
!        write(171,*) Position_Node_Element_Triangle( 1, i, Element_Number_on_Integrated_Line( e ) ), &
!                     Position_Node_Element_Triangle( 2, i, Element_Number_on_Integrated_Line( e ) )
!      end do 
!      write(171,*) Position_Node_Element_Triangle( 1, 1, Element_Number_on_Integrated_Line( e ) ), &
!                   Position_Node_Element_Triangle( 2, 1, Element_Number_on_Integrated_Line( e ) )
!      write(171,*)' ' 
!  end do
!
!  do e= 1, Number_Element
!    if( Counter_Node_in_Element_on_Integrated_Line( e ) >= 2 )then
!      do i= 1, 3
!        write(170,*) Position_Node_Element_Triangle( 1, i, e ), Position_Node_Element_Triangle( 2, i, e )
!      end do 
!      write(170,*) Position_Node_Element_Triangle( 1, 1, e ), Position_Node_Element_Triangle( 2, 1, e )
!      write(170,*)' ' 
!    end if
!  end do

  allocate( Position_Node_Element_on_Line( 2, 3, Number_Element_on_Integrated_Line ) )

  do e= 1, Number_Element_on_Integrated_Line
    !do i= 1, 3 
    !  do j= 1, 2 
    !    Position_Node_Element_on_Line( j, i, e )= Position_Node_Element_Triangle( j, i, Element_Number_on_Integrated_Line( e ) ) 
    !  end do
    !end do
    Position_Node_Element_on_Line( :, :, e )= Position_Node_Element_Triangle( :, :, Element_Number_on_Integrated_Line( e ) ) 
  end do

!  do e= 1, Number_Element_on_Integrated_Line
!    write(223,*) Position_Node_Element_on_Line( 1, Local_Node( 1, e ), e ), Position_Node_Element_on_Line( 2, Local_Node( 1, e ), e )
!    write(223,*) Position_Node_Element_on_Line( 1, Local_Node( 2, e ), e ), Position_Node_Element_on_Line( 2, Local_Node( 2, e ), e )
!    write(223,*)' '
!  end do

  allocate( Edge_Legnth_Integrated_Element( Number_Element_on_Integrated_Line ) )

  do e= 1, Number_Element_on_Integrated_Line
    Edge_Legnth_Integrated_Element( e ) &
    = sqrt( &
      (  Position_Node_Element_on_Line( 1, Local_Node( 1, e ), e ) &
        -Position_Node_Element_on_Line( 1, Local_Node( 2, e ), e ) )**2 &
     +(  Position_Node_Element_on_Line( 2, Local_Node( 1, e ), e ) &
        -Position_Node_Element_on_Line( 2, Local_Node( 2, e ), e ) )**2 )

    !write(224,*) Edge_Legnth_Integrated_Element( e )
  end do

  allocate( Value_OF_Element_Line( 3, Number_Element_on_Integrated_Line ) )

  do e= 1, Number_Element_on_Integrated_Line
    do i= 1, 3
      Value_OF_Element_Line( i, e ) = Value_Objective_Function( Index_Element_2_Node( i, Element_Number_on_Integrated_Line( e ) ) )
    end do
  end do

  allocate( Difference_Position_Node_Element_on_Line( 2, 3, 3, Number_Element_on_Integrated_Line ) )

  do e= 1, Number_Element_on_Integrated_Line
    do i= 1, 3
      do j= 1, 3
        do k= 1, 2
          Difference_Position_Node_Element_on_Line( k, j, i, e ) &
          =Position_Node_Element_on_Line( k, j, e ) -Position_Node_Element_on_Line( k, i, e )
        end do
      end do
    end do
  end do

  allocate( Area_Element_on_Line( Number_Element_on_Integrated_Line ) )

  Area_Element_on_Line( : ) &
  =abs( Difference_Position_Node_Element_on_Line( 1, 1, 3, : ) & 
       *Difference_Position_Node_Element_on_Line( 2, 2, 3, : ) &
       -Difference_Position_Node_Element_on_Line( 1, 2, 3, : ) & 
       *Difference_Position_Node_Element_on_Line( 2, 1, 3, : ) )/2.0d0

  allocate( BasisFunction_B_on_Line( 3, Number_Element_on_Integrated_Line ) )
  allocate( BasisFunction_C_on_Line( 3, Number_Element_on_Integrated_Line ) )

  BasisFunction_B_on_Line( 1, : ) &
  = Difference_Position_Node_Element_on_Line( 2, 2, 3, : )/( 2.0d0 *Area_Element_on_Line( : ) )

  BasisFunction_B_on_Line( 2, : ) & 
  = Difference_Position_Node_Element_on_Line( 2, 3, 1, : )/( 2.0d0 *Area_Element_on_Line( : ) )

  BasisFunction_B_on_Line( 3, : ) & 
  = Difference_Position_Node_Element_on_Line( 2, 1, 2, : )/( 2.0d0 *Area_Element_on_Line( : ) )

      
  BasisFunction_C_on_Line( 1, : ) & 
  = Difference_Position_Node_Element_on_Line( 1, 3, 2, : )/( 2.0d0 *Area_Element_on_Line( : ) )

  BasisFunction_C_on_Line( 2, : ) & 
  = Difference_Position_Node_Element_on_Line( 1, 1, 3, : )/( 2.0d0 *Area_Element_on_Line( : ) )

  BasisFunction_C_on_Line( 3, : ) & 
  = Difference_Position_Node_Element_on_Line( 1, 2, 1, : )/( 2.0d0 *Area_Element_on_Line( : ) )






!  do e= 1, Number_Element_on_Integrated_Line
!    do i= 1, 3
!      write(912,'(1x,es15.8,1x,es15.8,1x,es15.8)')Value_OF_Element_Line( i, e ), BasisFunction_B_on_Line( i, e ), BasisFunction_C_on_Line( i, e ) 
!    end do
!    write(914,'(1x,es15.8)') Edge_Legnth_Integrated_Element( e )
!  end do






  Objective_Function_Value= 0.0d0 

  ! Unit Normal n= ( -1.0d0, 0.0d0 )

  do e = 1, Number_Element_on_Integrated_Line
    do i= 1, 3
      Objective_Function_Value &
      = Objective_Function_Value &
       -( BasisFunction_B_on_Line( i, e )*Unit_Normal_Vector( 1 ) &
         +BasisFunction_C_on_Line( i, e )*Unit_Normal_Vector( 2 ) ) &
       *Value_OF_Element_Line( i, e )*Edge_Legnth_Integrated_Element( e )
    end do 
  end do 

  ! Average of Values of 2 Element sandwitching an Integrated Line
  Objective_Function_Value = 0.5d0 *Thermal_Conductivity_Evaluated * Objective_Function_Value

  write(913,*)'Objective_Function_Value=', Objective_Function_Value
  write(913,*)'Value_ObjectiveFunction_Normalize=', Value_ObjectiveFunction_Normalize

   if( Flag_Structure==2 )then
      !Objective_Function= Objective_Function_Value/Objective_Function_Value
      Objective_Function= Objective_Function_Value
   else
      Objective_Function= Objective_Function_Value/Value_ObjectiveFunction_Normalize
   end if

  deallocate( Difference_Position_Node_Element_on_Line )
  deallocate( Area_Element_on_Line )
  deallocate( BasisFunction_B_on_Line )
  deallocate( BasisFunction_C_on_Line )

  deallocate( Position_Node_Element_Triangle )

  deallocate( Element_Number_on_Integrated_Line ) 
  deallocate( Counter_Node_in_Element_on_Integrated_Line ) 
  deallocate( Flag_Element_Number_on_Integrated_Line ) 
  deallocate( Local_Node_tmp ) 
  deallocate( Local_Node ) 

  deallocate( Edge_Legnth_Integrated_Element )
  deallocate( Position_Node_Element_on_Line )
  deallocate( Value_OF_Element_Line )

  !=======================================================================================================
  write(*,*)'   end Compute_Objective_Function_Flux ' 
  !=======================================================================================================
  return
end subroutine Compute_Objective_Function_Flux



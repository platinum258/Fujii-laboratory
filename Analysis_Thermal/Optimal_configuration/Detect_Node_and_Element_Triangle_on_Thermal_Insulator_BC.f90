
subroutine Detect_Node_and_Element_Triangle_on_Thermal_Insulator_BC&
       ( Number_Node, Number_Element_Triangle, &
         Index_Element_2_Node_Triangle, &
         ID_Element_PEC, Class_Element_Triangle, &
         Position_Node, & 
      !================================================================
         Flag_Element_Triangle_PEC_BC, Number_Outer_Element_on_PEC_BC, &
         Position_Node_PEC_BC )

   !$ use omp_lib
   use Parameters
   implicit none

   ! intent(in)
   integer, intent(in) :: Number_Node, Number_Element_Triangle

   integer, intent(in) :: ID_Element_PEC 
   integer, intent(in) :: Index_Element_2_Node_Triangle( 3, Number_Element_Triangle ) 
   integer, intent(in) :: Class_Element_Triangle( Number_Element_Triangle )
   double precision, intent(in) :: Position_Node( 2, Number_Node )

   integer, intent(out) :: Flag_Element_Triangle_PEC_BC( Number_Element_Triangle ) 
   integer, intent(out) :: Number_Outer_Element_on_PEC_BC
   double precision, intent(out) :: Position_Node_PEC_BC( 2, 2, Number_Element_Triangle )
   
   integer :: e, i, j 
   integer, allocatable, dimension(:) :: Flag_Node_Outside_PEC, Flag_Node_Inside_PEC
   integer, allocatable, dimension(:) :: Counter_Node_in_Element, Flag_Node_PEC_BC 
   integer, allocatable, dimension(:,:) :: Local_Node_Number_on_PEC_BC 
   integer :: Counter

   !======================================================================================== 
   write(*,*)'         call Detect_Node_and_Element_Triangle_on_Thermal_Insulator_BC'
   !======================================================================================== 

   allocate( Flag_Node_Outside_PEC( Number_Node ) )
   allocate( Flag_Node_Inside_PEC( Number_Node ) )

   do i= 1, Number_Node 
      Flag_Node_Outside_PEC( i )= 0 
      Flag_Node_Inside_PEC( i )= 0
   end do

   do e= 1, Number_Element_Triangle
      if( Class_Element_Triangle( e )/=ID_Element_PEC )then
         do i= 1, 3 
            Flag_Node_Outside_PEC( Index_Element_2_Node_Triangle( i, e ) )= 1
         end do
      end if
      if( Class_Element_Triangle( e )==ID_Element_PEC )then
         do i= 1, 3 
            Flag_Node_Inside_PEC( Index_Element_2_Node_Triangle( i, e ) )= 2
         end do
      end if
   end do

   allocate( Flag_Node_PEC_BC( Number_Node ) )

   do i= 1, Number_Node
      Flag_Node_PEC_BC( i )= Flag_Node_Outside_PEC( i ) +Flag_Node_Inside_PEC( i ) 
   end do

   do i= 1, Number_Node
      if( Flag_Node_PEC_BC( i ) < 0 .or. 3 < Flag_Node_PEC_BC( i ) )then
         write(*,*)'Flag_Node_PEC_BC( i )=', Flag_Node_PEC_BC( i )
         call Output_Error( 'Detect_Node_and_Element_Triangle_on_Thermal_Insulator_BC', 62 )
      end if
   end do

   deallocate( Flag_Node_Outside_PEC )
   deallocate( Flag_Node_Inside_PEC )

   allocate( Counter_Node_in_Element( Number_Element_Triangle ) )
   allocate( Local_Node_Number_on_PEC_BC( 2, Number_Element_Triangle ) )

   do e= 1, Number_Element_Triangle
      Flag_Element_Triangle_PEC_BC( e )= Flag_Off
      Counter_Node_in_Element( e )= 0
   end do

   Local_Node_Number_on_PEC_BC= -1d8

   do e= 1, Number_Element_Triangle
      do i= 1, 3 
         if( Flag_Node_PEC_BC( Index_Element_2_Node_Triangle( i, e ) )==3 )then
            Counter_Node_in_Element( e )= Counter_Node_in_Element( e ) +1
            Local_Node_Number_on_PEC_BC( Counter_Node_in_Element( e ), e )= i
         end if
      end do
   end do

   deallocate( Flag_Node_PEC_BC )

   Counter= 0

   do e= 1, Number_Element_Triangle
      if( Counter_Node_in_Element( e )==2 .and. Class_Element_Triangle( e )/=ID_Element_PEC )then
         Flag_Element_Triangle_PEC_BC( e )= Flag_On
         Counter= Counter +1

         do i= 1, 2
            do j= 1, 2
               Position_Node_PEC_BC( j, i, e )= Position_Node( j, Index_Element_2_Node_Triangle( Local_Node_Number_on_PEC_BC( i, e ), e ) )
            end do
         end do
      else
         Position_Node_PEC_BC( :, :, e )= -1d8
      end if
   end do

   deallocate( Counter_Node_in_Element )
   deallocate( Local_Node_Number_on_PEC_BC )

   Number_Outer_Element_on_PEC_BC= Counter

   !======================================================================================== 
   write(*,*)'         end Detect_Node_and_Element_Triangle_on_Thermal_Insulator_BC'
   !======================================================================================== 

   return
end subroutine Detect_Node_and_Element_Triangle_on_Thermal_Insulator_BC



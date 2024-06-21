
subroutine Detect_Node_and_Element_Triangle_for_Boundary_Integral&
       ( Number_Node, Position_Node, &
         Number_Element_Triangle, Index_Element_2_Node_Triangle, &
         ID_Element_Thermal_Insulator, Class_Element_Triangle, &
         Number_Element_on_Flat_PEC_tmp, &
         Max_Position_Node, &
       !================================================================
         Element_Number_on_Flat_PEC_tmp, Number_Element_on_Flat_PEC, &
         Node_Number_on_Flat_PEC_in_Element )

   !$ use omp_lib
   use Parameters
   implicit none

   ! intent(in)
   integer, intent(in) :: Number_Node, Number_Element_Triangle
   double precision, intent(in) :: Position_Node( 2, Number_Node ) 

   integer, intent(in) :: ID_Element_Thermal_Insulator 
   integer, intent(in) :: Index_Element_2_Node_Triangle( 3, Number_Element_Triangle ) 
   integer, intent(in) :: Class_Element_Triangle( Number_Element_Triangle )
   integer, intent(in) :: Number_Element_on_Flat_PEC_tmp
   double precision, intent(in) :: Max_Position_Node( 2, 2 )

   integer, intent(out) :: Element_Number_on_Flat_PEC_tmp( Number_Element_on_Flat_PEC_tmp )
   integer, intent(out) :: Number_Element_on_Flat_PEC
   integer, intent(out) :: Node_Number_on_Flat_PEC_in_Element( 2, Number_Element_on_Flat_PEC_tmp )
 
   integer :: e, i, Counter 
   integer, allocatable, dimension(:) :: Flag_Node_Outside_PEC, Flag_Node_Inside_Flat_PEC, Flag_Node_ALL
   integer, allocatable, dimension(:) :: Counter_Node_in_Element_on_Flat_PEC 
   integer, allocatable, dimension(:,:) :: Local_Node_on_Flat_PEC

   !======================================================================================== 
   write(*,*)'         call Detect_Node_and_Element_Triangle_for_Boundary_Integral'
   !======================================================================================== 

   allocate( Flag_Node_ALL( Number_Node ) )

   Flag_Node_ALL= 0

   if( Flag_Thermal_Device==0 .or. Flag_Thermal_Device==20 .or. Flag_Thermal_Device==30 )then
      do i= 1, Number_Node
         if( Max_Position_Node( 1, 1 ) -1d-8 <= Position_Node( 1, i ) .and. &
             Position_Node( 1, i ) <= Max_Position_Node( 1, 1 ) +1d-8 )then
            Flag_Node_ALL( i )= 3 
         end if
         if( Max_Position_Node( 1, 2 ) -1d-8 <= Position_Node( 1, i ) .and. &
             Position_Node( 1, i ) <= Max_Position_Node( 1, 2 ) +1d-8 )then
            Flag_Node_ALL( i )= 3 
         end if
      end do
   else if( Flag_Thermal_Device==21 )then 
      do i= 1, Number_Node
         if( Max_Position_Node( 1, 1 ) -1d-8 <= Position_Node( 1, i ) .and. &
             Position_Node( 1, i ) <= Max_Position_Node( 1, 1 ) +1d-8 )then
            Flag_Node_ALL( i )= 3 
         end if
         if( Max_Position_Node( 1, 2 ) -1d-8 <= Position_Node( 1, i ) .and. &
             Position_Node( 1, i ) <= Max_Position_Node( 1, 2 ) +1d-8 )then
            Flag_Node_ALL( i )= 3 
         end if
         if( Max_Position_Node( 2, 1 ) -1d-8 <= Position_Node( 2, i ) .and. &
             Position_Node( 2, i ) <= Max_Position_Node( 2, 1 ) +1d-8 )then
            Flag_Node_ALL( i )= 3 
         end if
         if( Max_Position_Node( 2, 2 ) -1d-8 <= Position_Node( 2, i ) .and. &
             Position_Node( 2, i ) <= Max_Position_Node( 2, 2 ) +1d-8 )then
            Flag_Node_ALL( i )= 3 
         end if
      end do
   else if( Flag_Thermal_Device==10 )then ! aho

      do e= 1, Number_Element_Triangle
         if( Class_Element_Triangle( e ) == ID_Element_Thermal_Insulator )then
            do i= 1, 3
               Flag_Node_ALL( Index_Element_2_Node_Triangle( i, e ) )= 3
            end do
         end if
      end do

      do i= 1, Number_Node 
         if( Max_Position_Node( 2, 1 ) -1d-8 <= Position_Node( 2, i ) .and. &
             Position_Node( 2, i ) <= Max_Position_Node( 2, 1 ) +1d-8 )then
            Flag_Node_ALL( i )= 3
         end if
      end do 
   else
      write(*,*)'Flag_Thermal_Device=', Flag_Thermal_Device
      call Output_Error( 'Detect_Node_and_Element_Triangle_for_Boundary_Integral', 91 )
   end if

   do i= 1, Number_Node
      if( Flag_Node_ALL( i ) < 0 .or. 3 < Flag_Node_ALL( i ) )then
         write(*,*)'Flag_Node_ALL( i )=', Flag_Node_ALL( i )
         call Output_Error( 'Detect_Node_and_Element_Triangle_for_Boundary_Integral', 90 )
      end if
   end do

   allocate( Counter_Node_in_Element_on_Flat_PEC( Number_Element_Triangle ) )

   do e= 1, Number_Element_Triangle
      Counter_Node_in_Element_on_Flat_PEC( e )= 0
   end do

   allocate( Local_Node_on_Flat_PEC( 2, Number_Element_Triangle ) )

   do e= 1, Number_Element_Triangle
      if( Class_Element_Triangle( e )/=ID_Element_Thermal_Insulator )then
         Counter= 0
         do i= 1, 3 
            if( Flag_Node_ALL( Index_Element_2_Node_Triangle( i, e ) )==3 )then
               Counter_Node_in_Element_on_Flat_PEC( e )= Counter_Node_in_Element_on_Flat_PEC( e ) +1

               Counter= Counter +1
               if( Counter>=3 )then
                  write(*,*)'Counter=', Counter
                  call Output_Error( 'Detect_Node_and_Element_Triangle_for_Boundary_Integral', 117 ) 
               end if
               Local_Node_on_Flat_PEC( Counter, e )= i
            end if
         end do
      end if
   end do

   deallocate( Flag_Node_ALL )

   Node_Number_on_Flat_PEC_in_Element = 0
   Counter= 0
   do e= 1, Number_Element_Triangle
      if( Counter_Node_in_Element_on_Flat_PEC( e )==2 )then
         Counter= Counter +1
         Element_Number_on_Flat_PEC_tmp( Counter )= e

         if( Counter > Number_Element_on_Flat_PEC_tmp )then
            write(*,*)'Number_Element_on_Flat_PEC_tmp=', Number_Element_on_Flat_PEC_tmp
            write(*,*)'Counter=', Counter
            call Output_Error( 'Detect_Node_and_Element_Triangle_for_Boundary_Integral', 144 ) 
         end if

         Node_Number_on_Flat_PEC_in_Element( 1, Counter )&
         = Index_Element_2_Node_Triangle( Local_Node_on_Flat_PEC( 1, e ), e )
         Node_Number_on_Flat_PEC_in_Element( 2, Counter )&
         = Index_Element_2_Node_Triangle( Local_Node_on_Flat_PEC( 2, e ), e )

      end if
   end do

   deallocate( Counter_Node_in_Element_on_Flat_PEC )
   deallocate( Local_Node_on_Flat_PEC )

   Number_Element_on_Flat_PEC= Counter

   if( Number_Element_on_Flat_PEC_tmp < Number_Element_on_Flat_PEC )then
      write(*,*)'Number_Element_on_Flat_PEC_tmp < Number_Element_on_Flat_PEC'
      write(*,*)'Number_Element_on_Flat_PEC_tmp =', Number_Element_on_Flat_PEC_tmp
      write(*,*)'Number_Element_on_Flat_PEC =', Number_Element_on_Flat_PEC
      call Output_Error( 'Detect_Node_and_Element_Triangle_for_Boundary_Integral', 125 )
   end if

   if( maxval( Node_Number_on_Flat_PEC_in_Element ) > Number_Node )then
      write(*,*) 'maxval( Node_Number_on_Flat_PEC_in_Element ) > Number_Node' 
      write(*,*) 'maxval( Node_Number_on_Flat_PEC_in_Element )', maxval( Node_Number_on_Flat_PEC_in_Element ) 
      write(*,*) 'Number_Node=', Number_Node 
      call Output_Error( 'Detect_Node_and_Element_Triangle_for_Boundary_Integral', 128 )
   end if

   !======================================================================================== 
   write(*,*)'         end Detect_Node_and_Element_Triangle_for_Boundary_Integral'
   !======================================================================================== 

   return
end subroutine Detect_Node_and_Element_Triangle_for_Boundary_Integral



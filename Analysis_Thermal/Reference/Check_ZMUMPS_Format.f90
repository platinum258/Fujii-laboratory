

subroutine Check_ZMUMPS_Format( aa_zmumps, ia_zmumps, ja_zmumps, Number_Node, Number_NonZero, Flag_Symmetric, Flag_Matrix_Upper )

   !use omp_lib
   use Parameters
   implicit none

   integer, intent(in) :: Number_Node, Number_NonZero 
   complex( kind( 0d0 ) ), intent(in) :: aa_zmumps( Number_NonZero )
   integer, intent(in) :: ia_zmumps( Number_NonZero )
   integer, intent(in) :: ja_zmumps( Number_NonZero )
   integer, intent(in) :: Flag_Symmetric, Flag_Matrix_Upper

   integer :: i 
   integer, allocatable, dimension(:) :: Counter_Row, Counter_Column
   integer :: I_Max_AVM, J_Max_AVM
   integer :: I_Min_AVM, J_Min_AVM

   !========================================================================================
   write(*,*)'         call Check_ZMUMPS_Format'
   !========================================================================================

   if( Number_Node <= 0 .or. Number_Node > Number_NonZero )then
       call Output_Error( 'Check_ZMUMPS_Format', 24 )
   end if

   if( Number_NonZero <= 0 .or. Number_NonZero > Integer_Initialization_Plus )then
       call Output_Error( 'Check_ZMUMPS_Format', 28 )
   end if

   I_Max_AVM= 0
   J_Max_AVM= 0
   I_Min_AVM= Number_Node
   J_Min_AVM= Number_Node

   do i= 1, Number_NonZero
      if( aa_zmumps( i )==Zero )then
         write(*,*)'i=', i
         write(*,*)'aa_zmumps( i )==Zero'
         write(*,*)'ia_zmumps( i )=', ia_zmumps( i )
         write(*,*)'ja_zmumps( i )=', ja_zmumps( i )
         call Output_Error( 'Check_ZMUMPS_Format', 21 )
      end if
 
      if( ia_zmumps( i ) > Number_Node .or. ia_zmumps( i ) <= 0 )then

         write(*,*)'ia_zmumps( i ) > Number_Node'
         write(*,*)'i=', i
         write(*,*)'ia_zmumps( i )=', ia_zmumps( i )
         write(*,*)'Number_Node=', Number_Node
         call Output_Error( 'Check_ZMUMPS_Format', 33 )
      end if

      if( ja_zmumps( i ) > Number_Node .or. ja_zmumps( i ) <= 0 )then

         write(*,*)'ja_zmumps( i ) > Number_Node'
         write(*,*)'i=', i
         write(*,*)'ja_zmumps( i )=', ja_zmumps( i )
         write(*,*)'Number_Node=', Number_Node
         call Output_Error( 'Check_ZMUMPS_Format', 42 )
      end if

      if( ia_zmumps( i ) < ja_zmumps( i ) .and. Flag_Symmetric==Flag_On .and. Flag_Matrix_Upper==Flag_Off )then
         write(*,*)'ia_zmumps( i )=', ia_zmumps( i )
         write(*,*)'ja_zmumps( i )=', ja_zmumps( i )
         write(*,*)'Flag_Symmetric=', Flag_Symmetric
         write(*,*)'Flag_Matrix_Upper=', Flag_Matrix_Upper
         call Output_Error( 'Check_ZMUMPS_Format', 61 )
      end if

      if( ia_zmumps( i ) > ja_zmumps( i ) .and. Flag_Symmetric==Flag_On .and. Flag_Matrix_Upper==Flag_On )then
         write(*,*)'ia_zmumps( i )=', ia_zmumps( i )
         write(*,*)'ja_zmumps( i )=', ja_zmumps( i )
         write(*,*)'Flag_Symmetric=', Flag_Symmetric
         write(*,*)'Flag_Matrix_Upper=', Flag_Matrix_Upper
         call Output_Error( 'Check_ZMUMPS_Format', 75 )
      end if

      if( I_Max_AVM < ia_zmumps( i ) )then
         I_Max_AVM= ia_zmumps( i )
      end if

      if( I_Min_AVM > ia_zmumps( i ) )then
         I_Min_AVM= ia_zmumps( i )
      end if

      if( J_Max_AVM < ja_zmumps( i ) )then
         J_Max_AVM= ja_zmumps( i )
      end if

      if( J_Min_AVM > ja_zmumps( i ) )then
         J_Min_AVM= ja_zmumps( i )
      end if

   end do

   if( I_Max_AVM/=Number_Node )then
      write(*,*)'I_Max_AVM/=Number_Node'
      write(*,*)'I_Max_AVM=', I_Max_AVM 
      write(*,*)'Number_Node=', Number_Node
      call Output_Error( 'Check_ZMUMPS_Format', 63 )
   end if 

   if( J_Max_AVM/=Number_Node )then
      write(*,*)'J_Max_AVM/=Number_Node'
      write(*,*)'J_Max_AVM=', J_Max_AVM 
      write(*,*)'Number_Node=', Number_Node
      call Output_Error( 'Check_ZMUMPS_Format', 71 )
   end if

   if( I_Min_AVM/=1 )then
      write(*,*)'I_Min_AVM/= 1'
      write(*,*)'I_Min_AVM=', I_Min_AVM 
      write(*,*)'Number_Node=', Number_Node
      call Output_Error( 'Check_ZMUMPS_Format', 88 )
   end if 

   if( J_Min_AVM/=1 )then
      write(*,*)'J_Min_AVM/= 1'
      write(*,*)'J_Min_AVM=', J_Min_AVM 
      write(*,*)'Number_Node=', Number_Node
      call Output_Error( 'Check_ZMUMPS_Format', 95 )
   end if

   write(*,*)'            (', I_Min_AVM, J_Min_AVM, ') --> (', I_Max_AVM, J_Max_AVM, ')'
   write(*,*)'            Number_Node=', Number_Node
   write(*,*)'            Number_NonZero=', Number_NonZero

   allocate( Counter_Row( Number_Node ) )
   allocate( Counter_Column( Number_Node ) )

   do i= 1, Number_Node
      Counter_Row( i )= 0
      Counter_Column( i )= 0
   end do

   do i= 1, Number_NonZero
      Counter_Row( ia_zmumps( i ) )= Counter_Row( ia_zmumps( i ) ) +1
      Counter_Column( ja_zmumps( i ) )= Counter_Column( ja_zmumps( i ) ) +1
   end do

   deallocate( Counter_Row )
   deallocate( Counter_Column )

   return
end subroutine Check_ZMUMPS_Format


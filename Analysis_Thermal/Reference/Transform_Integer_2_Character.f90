
subroutine Transform_Integer_2_Character &
       ( Integer_input, Length_Character_Integer_input, &
         Character_output )

   !$use omp_lib
   use Parameters
   implicit none

   integer, intent(in) :: Integer_input 
   integer, intent(in) :: Length_Character_Integer_input

   character(len=*), intent(out) :: Character_output

   character(len=8) :: Format_Filenumber
   character(len=Length_Character_Integer_input-1) :: Character_output_tmp

   !====================================================================================
   write(*,*)'    call Transform_Integer_2_Character'
   !====================================================================================
   write(*,*)'       Integer_input=', Integer_input
   write(*,*)'       Length_Character_Integer_input=', Length_Character_Integer_input

   if( Integer_input >= 0 )then

      write( Format_Filenumber, '(A2,I1,A1,I1,A1)' ) & 
      '(i', Length_Character_Integer_input,'.',Length_Character_Integer_input,')'
     
      write( Character_output, Format_Filenumber ) Integer_input

   else 

      write( Format_Filenumber, '(A2,I1,A1,I1,A1)' ) & 
      '(i', Length_Character_Integer_input-1,'.',Length_Character_Integer_input-1,')'
     
      write( Character_output_tmp, Format_Filenumber ) abs(Integer_input)

      Character_output= trim("-"//trim(Character_output_tmp))

   end if


   write(*,*) Character_output, Length_Character_Integer_input



   !====================================================================================
   write(*,*)'    end Transform_Integer_2_Character'
   !====================================================================================
   return
end subroutine Transform_Integer_2_Character


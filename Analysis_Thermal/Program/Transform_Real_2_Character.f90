
subroutine Transform_Real_2_Character &
       ( Real_input, Length_Character_Real_input, &
         Character_output )

   !$use omp_lib
   use Parameters
   implicit none

   double precision, intent(in) :: Real_input 
   integer, intent(in) :: Length_Character_Real_input

   character(len=*), intent(out) :: Character_output

   character(len=8) :: Format_Filenumber

   integer :: Integer_Part, Decimal_Part, Number_Power_Decimal_Part, Number_Power_Integer_Part

   integer :: Length_Character_Integer_Part, Length_Character_Decimal_Part
   integer :: Length_Character_Integer_Part_tmp

   double precision :: Decimal, Torelance_Dicimal

   character(:), allocatable :: Character_Integer_Part, Character_Decimal_Part

   !====================================================================================
   write(*,*)'   call Transform_Real_2_Character'
   !====================================================================================

   Torelance_Dicimal=1d-8

   Integer_Part= int( Real_input )
   if( Real_input >= 0d0 )then
      if( Integer_Part==0 )then
         Number_Power_Integer_Part= 1 
         Decimal= Real_input -dble( Integer_Part )
      else
         Number_Power_Integer_Part= Int( log10( dble( Integer_Part ) ) ) +1 
         Decimal= Real_input -dble( Integer_Part )
      end if
   else
      if( Integer_Part==0 )then
         Number_Power_Integer_Part= 2 
         Decimal= -( Real_input  -dble( Integer_Part ) ) +1d-6 
      else
         Number_Power_Integer_Part= Int( log10( dble( -Integer_Part ) ) ) +2 
         Decimal= -( Real_input  -dble( Integer_Part ) ) +1d-6 
      end if
   end if
   Length_Character_Integer_Part= Number_Power_Integer_Part

   call Obtain_Number_Power_Decimal_Part &
      ( Real_input, Number_Power_Decimal_Part )

   !call Obtain_Number_Power_Decimal_Part &
   !   ( Decimal, Number_Power_Decimal_Part )

   write(*,*) '      Real_input=', Real_input
   write(*,*) '      Decimal=', Decimal, 'Number_Power_Decimal_Part=', Number_Power_Decimal_Part
   write(*,*) '      Number_Power_Decimal_Part=', Number_Power_Decimal_Part

   Length_Character_Decimal_Part= Number_Power_Decimal_Part

   !Decimal_Part= Int( Decimal*10**(Number_Power_Decimal_Part) )
   Decimal_Part= Int( Decimal*( 1d0 +Torelance_Dicimal ) *10**(Number_Power_Decimal_Part) )

   write(*,*) '      Integer_Part=', Integer_Part, 'Decimal_Part=', Decimal_Part
   !Decimal, Number_Power_Decimal_Part
   !write(*,*) Real_input, Integer_Part, Decimal_Part, Decimal, Number_Power_Decimal_Part

   !if( Integer_Part >= 0 )then
   if( Real_input >= 0 )then

      write( Format_Filenumber, '(A2,I1,A1,I1,A1)' ) & 
      '(i', Length_Character_Integer_Part,'.',Length_Character_Integer_Part,')'
   
      !write(*,*)'Length_Character_Integer_Part=', Length_Character_Integer_Part
      allocate( character(Length_Character_Integer_Part)::Character_Integer_Part )
      write( Character_Integer_Part, Format_Filenumber ) Integer_Part
   
      !write(*,*)'Length_Character_Decimal_Part=', Length_Character_Decimal_Part
      write( Format_Filenumber, '(A2,I1,A1,I1,A1)' ) & 
      '(i', Length_Character_Decimal_Part,'.',Length_Character_Decimal_Part,')'
     
      allocate( character(Length_Character_Decimal_Part)::Character_Decimal_Part )
      write( Character_Decimal_Part, Format_Filenumber ) Decimal_Part
   
      Character_output= trim( trim(Character_Integer_Part)//"."//trim(Character_Decimal_Part) )

   else

      Length_Character_Integer_Part_tmp= Length_Character_Integer_Part -1

      write( Format_Filenumber, '(A2,I1,A1,I1,A1)' ) & 
      '(i', Length_Character_Integer_Part_tmp,'.', Length_Character_Integer_Part_tmp,')'
   
      !write(*,*)'Length_Character_Integer_Part=', Length_Character_Integer_Part
      allocate( character(Length_Character_Integer_Part_tmp)::Character_Integer_Part )
      write( Character_Integer_Part, Format_Filenumber ) -Integer_Part
   
      !write(*,*)'Length_Character_Decimal_Part=', Length_Character_Decimal_Part
      write( Format_Filenumber, '(A2,I1,A1,I1,A1)' ) & 
      '(i', Length_Character_Decimal_Part,'.',Length_Character_Decimal_Part,')'
     
      allocate( character(Length_Character_Decimal_Part)::Character_Decimal_Part )
      write( Character_Decimal_Part, Format_Filenumber ) Decimal_Part
   
      Character_output= trim( "-"//trim(Character_Integer_Part)//"."//trim(Character_Decimal_Part) )

   end if

   write(*,*) '      ', Character_Integer_Part, ".", Character_Decimal_Part 
   write(*,*) '      ', Character_output, Length_Character_Real_input

   deallocate( Character_Integer_Part )
   deallocate( Character_Decimal_Part )

   !====================================================================================
   write(*,*)'   end Transform_Real_2_Character'
   !====================================================================================
   return
end subroutine Transform_Real_2_Character


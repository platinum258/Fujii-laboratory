
subroutine Obtain_Number_Power_Decimal_Part &
           ( Real_input, Number_Power_Decimal_Part )

     !$use omp_lib
     use Parameters
     implicit none

     double precision, intent(in) :: Real_input 
     integer, intent(out) :: Number_Power_Decimal_Part 

     integer :: i

     integer :: Number_Power_Decimal_Part_tmp
     double precision :: Difference_Real_Integer, Difference_Real_Integer_Previous 
     double precision :: abs_Real_input 
     double precision :: Real_input_tmp
     double precision :: Power_operation 

     integer :: Number_limit_power 
     integer :: Number_limit_power_zero, Counter_power_zero
     double precision :: Tolerance_digit
     double precision :: Real_10fold 

     integer :: Integer_Part
     double precision :: Decimal_Part
     character(len=1) :: match

     !====================================================================================
     write(*,*)'          call Obtain_Number_Power_Decimal_Part'
     !====================================================================================

     Integer_Part= int( Real_input )
     if( Real_input > 0.0d0 )then
          Decimal_Part= Real_input -dble( Integer_Part )
     else if( match( Real_input, 0.0d0 )=='y' )then
          Number_Power_Decimal_Part=1 
          goto 680 
     else
          Decimal_Part= -( Real_input -dble( Integer_Part ) ) +1d-6 
     end if

     write(*,*) Real_input
     Number_limit_power= 8- int( log10( abs( Real_input ) ) )
     Power_operation= 10d0**(Number_limit_power)
     Tolerance_digit= 1d-1**( Number_limit_power +3 )
     Number_limit_power_zero= 8

     if( Real_input >= 0.0d0 )then
          Real_input_tmp= dble( int( ( Real_input +Tolerance_digit )*Power_operation ) )/Power_operation 
     else
          Real_input_tmp= dble( int( ( -Real_input +Tolerance_digit )*Power_operation ) )/Power_operation 
     end if

     abs_Real_input= abs( Real_input_tmp +Tolerance_digit )

     Number_Power_Decimal_Part_tmp= Power_operation

     !Real_10fold= int( Real_input_tmp )*10d0
     !Difference_Real_Integer_Previous= 0d0

     write(*,*)'               Real_input=', Real_input
     !do i= 0, Number_limit_power
     !     !Difference_Real_Integer= abs( abs_Real_input*10d0**i -Real_10fold ) 
     !     Difference_Real_Integer= abs( abs_Real_input*10d0**i -int( abs_Real_input*10d0**( i ) ) ) 
 
     !     write(*,*)'                i=', i, 'Difference_Real_Integer=', Difference_Real_Integer
     !     !if( Difference_Real_Integer <= Difference_Real_Integer_Previous*( 1.0d0 -Tolerance_digit ) )then
     !     if( i >= 1 .and. &
     !         Difference_Real_Integer <= Difference_Real_Integer_Previous )then
     !          Number_Power_Decimal_Part_tmp= i
     !     end if

     !     Difference_Real_Integer_Previous= Difference_Real_Integer
     !     !Real_10fold= dble( int( abs_Real_input*10d0**(i+1) ) )
     !end do
     
     Difference_Real_Integer_Previous= abs( abs_Real_input -int( abs_Real_input ) )
     do i= 1, Number_limit_power
          Difference_Real_Integer= abs( abs_Real_input*( 10.0d0**i ) -int( abs_Real_input*( 10.0d0**i ) ) ) 
 
          write(*,*)'                i=', i, 'Difference_Real_Integer=', Difference_Real_Integer
          if( Difference_Real_Integer <= Difference_Real_Integer_Previous )then
               Number_Power_Decimal_Part_tmp= i
          end if

          Difference_Real_Integer_Previous= Difference_Real_Integer
     end do


     Number_Power_Decimal_Part= Number_Power_Decimal_Part_tmp

     680 continue
     write(*,*)'               Number_Power_Decimal_Part=', Number_Power_Decimal_Part
     !====================================================================================
     write(*,*)'          end Obtain_Number_Power_Decimal_Part'
     !====================================================================================
     return
end subroutine Obtain_Number_Power_Decimal_Part


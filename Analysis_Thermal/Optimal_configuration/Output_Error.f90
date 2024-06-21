
subroutine Output_Error( Name_Subroutine, Line_Number )

   implicit none

   character(len=*), intent(in) :: Name_Subroutine
   integer, intent(in) :: Line_Number

   write( *, * )'=============================================='
   write( *, * )'ERROR'
   write( *, * ) Name_Subroutine, '.f90', Line_Number
   write( *, * )'=============================================='

   open( 100, file='./Error', status='replace' )
      write( 100, * )'ERROR'
      write( 100, * ) Name_Subroutine, '.f90', Line_Number
   close( 100 )

   stop
end subroutine Output_Error 


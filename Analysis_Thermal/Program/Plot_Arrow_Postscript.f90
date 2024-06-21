
subroutine Plot_Arrow_Postscript &
       ( File_Number, Flag_Type_Arrow, Color_Arrow, &
         Position_Arrow_Root_X, Position_Arrow_Root_Y, &
         Position_Arrow_Tip_X, Position_Arrow_Tip_Y, & 
         Stemthick, Headthick, Headlength )

   !use omp_lib
   use Parameters
   implicit none

   integer, intent(in) :: File_Number, Flag_Type_Arrow
   double precision, intent(in) :: Position_Arrow_Root_X, Position_Arrow_Root_Y
   double precision, intent(in) :: Position_Arrow_Tip_X, Position_Arrow_Tip_Y
   double precision, intent(in) :: Stemthick, Headthick, Headlength  
   character(len=1), intent(in) :: Color_Arrow

   integer :: e, i, j, k, l 
   double precision :: Ratio, Length_Arrow 
   double precision :: Position_Arrow_tmp( 2, 2 ), Position_Arrow( 2, 4 )
   double precision :: RGB( 3 )

   if( Color_Arrow=='w' )then
      RGB( 1 )=1d0
      RGB( 2 )=1d0
      RGB( 3 )=1d0
   else if( Color_Arrow=='r' )then
      RGB( 1 )=1.0d0
      RGB( 2 )=0.0d0
      RGB( 3 )=0.0d0
   else
      call Output_Error( 'Plot_Arrow_Postscript', 32 )
   end if

   Length_Arrow &
   = sqrt( ( Position_Arrow_Tip_X -Position_Arrow_Root_X )**2 +( Position_Arrow_Tip_Y -Position_Arrow_Root_Y )**2 )

   Position_Arrow_tmp( 1, 1 )= -( Position_Arrow_Tip_X-Position_Arrow_Root_X )/Length_Arrow *Headlength
   Position_Arrow_tmp( 2, 1 )= -( Position_Arrow_Tip_Y-Position_Arrow_Root_Y )/Length_Arrow *Headlength

   Ratio= Position_Arrow_tmp( 2, 1 )/Position_Arrow_tmp( 1, 1 )

   Position_Arrow_tmp( 2, 2 )= Headthick/sqrt( Ratio*Ratio +1d0 )
   Position_Arrow_tmp( 1, 2 )= -Ratio*Position_Arrow_tmp( 2, 2 ) 

   Position_Arrow( 1, 1 )= Position_Arrow_Tip_X
   Position_Arrow( 2, 1 )= Position_Arrow_Tip_Y
   Position_Arrow( 1, 2 )= Position_Arrow_Tip_X +Position_Arrow_tmp( 1, 1 ) +Position_Arrow_tmp( 1, 2 )
   Position_Arrow( 2, 2 )= Position_Arrow_Tip_Y +Position_Arrow_tmp( 2, 1 ) +Position_Arrow_tmp( 2, 2 )
   Position_Arrow( 1, 3 )= Position_Arrow_Tip_X +Position_Arrow_tmp( 1, 1 ) -Position_Arrow_tmp( 1, 2 )
   Position_Arrow( 2, 3 )= Position_Arrow_Tip_Y +Position_Arrow_tmp( 2, 1 ) -Position_Arrow_tmp( 2, 2 )
   Position_Arrow( 1, 4 )= Position_Arrow_Tip_X +Position_Arrow_tmp( 1, 1 )*0.8d0 
   Position_Arrow( 2, 4 )= Position_Arrow_Tip_Y +Position_Arrow_tmp( 2, 1 )*0.8d0 

   write( File_Number, * ) RGB( 1 ), RGB( 2 ), RGB( 3 ), 'setrgbcolor' 
   write( File_Number, * ) Stemthick, 'setlinewidth'
   write( File_Number, * ) 'newpath'
   write( File_Number, * ) Position_Arrow_Root_X, Position_Arrow_Root_Y, 'moveto' 
   write( File_Number, * ) Position_Arrow_Tip_X, Position_Arrow_Tip_Y, 'lineto' 
   write( File_Number, * )'stroke'

   if( Flag_Type_Arrow==0 )then
      write( File_Number, * ) Stemthick, 'setlinewidth'
      write( File_Number, * ) 'newpath'
      write( File_Number, * ) Position_Arrow( 1, 1 ), Position_Arrow( 2, 1 ), 'moveto' 
      write( File_Number, * ) Position_Arrow( 1, 2 ), Position_Arrow( 2, 2 ), 'lineto' 
      write( File_Number, * )'stroke'

      write( File_Number, * ) Stemthick, 'setlinewidth'
      write( File_Number, * ) 'newpath'
      write( File_Number, * ) Position_Arrow( 1, 1 ), Position_Arrow( 2, 1 ), 'moveto' 
      write( File_Number, * ) Position_Arrow( 1, 3 ), Position_Arrow( 2, 3 ), 'lineto' 
      write( File_Number, * )'stroke'
   else if( Flag_Type_Arrow==1 )then
      write( File_Number, * ) '0.00001 setlinewidth'
      write( File_Number, * ) 'gsave'
      write( File_Number, * ) 'newpath'
      write( File_Number, * ) Position_Arrow( 1, 1 ), Position_Arrow( 2, 1 ), 'moveto' 
      write( File_Number, * ) Position_Arrow( 1, 2 ), Position_Arrow( 2, 2 ), 'lineto' 
      write( File_Number, * ) Position_Arrow( 1, 3 ), Position_Arrow( 2, 3 ), 'lineto' 
      write( File_Number, * )'closepath'
      write( File_Number, * ) RGB( 1 ), RGB( 2 ), RGB( 3 ), 'setrgbcolor' 
      write( File_Number, * )'fill'
      write( File_Number, * )'grestore'
   else if( Flag_Type_Arrow==2 )then
      write( File_Number, * ) '0.00001 setlinewidth'
      write( File_Number, * ) 'gsave'
      write( File_Number, * ) 'newpath'
      write( File_Number, * ) Position_Arrow( 1, 1 ), Position_Arrow( 2, 1 ), 'moveto' 
      write( File_Number, * ) Position_Arrow( 1, 2 ), Position_Arrow( 2, 2 ), 'lineto' 
      write( File_Number, * ) Position_Arrow( 1, 4 ), Position_Arrow( 2, 4 ), 'lineto' 
      write( File_Number, * )'closepath'
      write( File_Number, * ) RGB( 1 ), RGB( 2 ), RGB( 3 ), 'setrgbcolor' 
      write( File_Number, * )'fill'
      write( File_Number, * )'grestore'
      write( File_Number, * ) 'gsave'
      write( File_Number, * ) 'newpath'
      write( File_Number, * ) Position_Arrow( 1, 1 ), Position_Arrow( 2, 1 ), 'moveto' 
      write( File_Number, * ) Position_Arrow( 1, 3 ), Position_Arrow( 2, 3 ), 'lineto' 
      write( File_Number, * ) Position_Arrow( 1, 4 ), Position_Arrow( 2, 4 ), 'lineto' 
      write( File_Number, * )'closepath'
      write( File_Number, * ) RGB( 1 ), RGB( 2 ), RGB( 3 ), 'setrgbcolor' 
      write( File_Number, * )'fill'
      write( File_Number, * )'grestore'
   end if

   return
end subroutine Plot_Arrow_Postscript 


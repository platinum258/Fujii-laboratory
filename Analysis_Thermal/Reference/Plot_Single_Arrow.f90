
subroutine Plot_Single_Arrow&
          ( File_Number, Color_Arrow, Flag_Type_Arrow, Direction_Arrow_X, Direction_Arrow_Y, &
            Position_Arrow_X, Position_Arrow_Y, Level_Arrow, Max_Level_Arrow, Min_Size_Arrow, Max_Size_Arrow, Angle_Arrow ) 

   implicit none
   integer, intent(in) :: File_Number, Flag_Type_Arrow, Level_Arrow, Max_Level_Arrow, Angle_Arrow
   character(len=1), intent(in) :: Color_Arrow
   double precision, intent(in) :: Position_Arrow_X, Position_Arrow_Y
   double precision, intent(in) :: Direction_Arrow_X, Direction_Arrow_Y
   double precision, intent(in) :: Min_Size_Arrow, Max_Size_Arrow

   double precision :: Position_Arrow_Root_X, Position_Arrow_Root_Y
   double precision :: Position_Arrow_Tip_X, Position_Arrow_Tip_Y
   double precision :: Position_Arrow_Tip_Stem_X, Position_Arrow_Tip_Stem_Y
   double precision :: Stemthick, Headthick, Headlength  
   double precision, allocatable, dimension(:) :: Direction_Normalize

   integer :: e, i, j, k, l 
   double precision :: Ratio, Length_Arrow 
   double precision :: Position_Arrow_tmp( 2, 2 ), Position_Arrow( 2, 4 )
   double precision :: RGB( 3 )

   allocate( Direction_Normalize( 2 ) )

   Direction_Normalize( 1 )&
   = Direction_Arrow_X/( sqrt( Direction_Arrow_X**2+ Direction_Arrow_Y**2 ) )
   Direction_Normalize( 2 )&
   = Direction_Arrow_Y/( sqrt( Direction_Arrow_X**2+ Direction_Arrow_Y**2 ) )

   Length_Arrow &
   = ( Max_Size_Arrow -Min_Size_Arrow )/dble( Max_Level_Arrow )*dble( Level_Arrow ) +Min_Size_Arrow

   !Stemthick= 0.1d0*Length_Arrow
   Stemthick= 0.01d0*Max_Size_Arrow
   if( Angle_Arrow==45 )then
      Headthick= Length_Arrow/7d0 
      Headlength= Length_Arrow/5d0
   else if( Angle_Arrow==30 )then
      Headthick= Length_Arrow/10d0 
      Headlength= Length_Arrow/4d0
      !Headthick= Max_Size_Arrow/7d0/2d0
      !Headlength= Max_Size_Arrow/4d0/2d0
   end if

   Position_Arrow_Root_X &
   =Position_Arrow_X -Length_Arrow/2.0d0*Direction_Normalize( 1 )  
   Position_Arrow_Root_Y &
   =Position_Arrow_Y -Length_Arrow/2.0d0*Direction_Normalize( 2 )  

   Position_Arrow_Tip_X &
   =Position_Arrow_X +Length_Arrow/2.0d0*Direction_Normalize( 1 )  
   Position_Arrow_Tip_Y &
   =Position_Arrow_Y +Length_Arrow/2.0d0*Direction_Normalize( 2 )  

   if( Flag_Type_Arrow==0 )then
      Position_Arrow_Tip_Stem_X &
      =Position_Arrow_X +( Length_Arrow/2.0d0 -Stemthick*0.9d0 )*Direction_Normalize( 1 )  
      Position_Arrow_Tip_Stem_Y &
      =Position_Arrow_Y +( Length_Arrow/2.0d0 -Stemthick*0.9d0 )*Direction_Normalize( 2 )  
   else
      Position_Arrow_Tip_Stem_X &
      =Position_Arrow_X +( Length_Arrow/2.0d0 -Headlength*0.9d0 )*Direction_Normalize( 1 )  
      Position_Arrow_Tip_Stem_Y &
      =Position_Arrow_Y +( Length_Arrow/2.0d0 -Headlength*0.9d0 )*Direction_Normalize( 2 )  
   end if

   if( Color_Arrow=='w' )then
      RGB( 1 )=1d0
      RGB( 2 )=1d0
      RGB( 3 )=1d0
   else if( Color_Arrow=='r' )then
      RGB( 1 )=1.0d0
      RGB( 2 )=0.0d0
      RGB( 3 )=0.0d0
   else
      call Output_Error( 'Plot_Single_Arrow', 32 )
   end if

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
   !write( File_Number, * ) Position_Arrow_Tip_X, Position_Arrow_Tip_Y, 'lineto' 
   write( File_Number, * ) Position_Arrow_Tip_Stem_X, Position_Arrow_Tip_Stem_Y, 'lineto' 
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
   else 
      write(*,*)'=================================='
      write(*,*)'Error Plot_Single_Arrow.f90 131'
      write(*,*)'=================================='
      stop
   end if

   deallocate( Direction_Normalize )

   return
end subroutine Plot_Single_Arrow 


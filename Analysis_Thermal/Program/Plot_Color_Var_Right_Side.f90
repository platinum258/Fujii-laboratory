

subroutine Plot_Color_Var_Right_Side & 
           ( File_Number, RGB, Number_Color_Level, Position_Maximum_Color_Var, Minimum_Value_Plotted, Maximum_Value_Plotted )

     !use omp_lib
     implicit none

     integer, intent(in) :: File_Number, Number_Color_Level 
     double precision, intent(in) :: RGB( 3, Number_Color_Level )
     double precision, intent(in) :: Position_Maximum_Color_Var( 2, 2, 2 )

     double precision, intent(in) :: Minimum_Value_Plotted, Maximum_Value_Plotted

     integer :: i, j, k 

     double precision, allocatable, dimension(:,:,:,:) :: Position_Color_Var_Level
     double precision, allocatable, dimension(:,:) :: Position_Boundary_Level
     double precision, allocatable, dimension(:) :: Width_Var

     allocate( Width_Var( 2 ) )

     do i= 1, 2
          Width_Var( i )= Position_Maximum_Color_Var( i, 2, 1 ) -Position_Maximum_Color_Var( i, 1, 1 )
     end do


     allocate( Position_Boundary_Level( 2, Number_Color_Level ) )

     do i= 1, Number_Color_Level
          Position_Boundary_Level( 1, i ) & 
          = ( Position_Maximum_Color_Var( 2, 2, 1 ) &
             -Position_Maximum_Color_Var( 2, 1, 1 ) )/dble( Number_Color_Level ) *dble( i-1 ) &
             +Position_Maximum_Color_Var( 2, 1, 1 )  
          Position_Boundary_Level( 2, i ) & 
          = ( Position_Maximum_Color_Var( 2, 2, 1 ) &
             -Position_Maximum_Color_Var( 2, 1, 1 ) )/dble( Number_Color_Level ) *dble( i ) &
             +Position_Maximum_Color_Var( 2, 1, 1 )  
     end do

 
     allocate( Position_Color_Var_Level( 2, 2, 2, Number_Color_Level ) )

     do i= 1, Number_Color_Level 
          do j= 1, 2
               do k= 1, 2
                    Position_Color_Var_Level( 1, k, j, i )= Position_Maximum_Color_Var( 1, 1, j )  
                    Position_Color_Var_Level( 2, k, j, i )= Position_Boundary_Level( k, i )
               end do
          end do
     end do

     deallocate( Position_Boundary_Level )

     !write( File_Number, * )'# Color Var --------------------------------------------'
     do i= 1, Number_Color_Level
          write( File_Number, * )'gsave'
          write( File_Number, * ) Position_Color_Var_Level( 1, 1, 1, i ), Position_Color_Var_Level( 2, 1, 1, i ), 'moveto'
          write( File_Number, * ) Position_Color_Var_Level( 1, 1, 2, i ), Position_Color_Var_Level( 2, 1, 2, i ), 'lineto'
          write( File_Number, * ) Position_Color_Var_Level( 1, 2, 2, i ), Position_Color_Var_Level( 2, 2, 2, i ), 'lineto'
          write( File_Number, * ) Position_Color_Var_Level( 1, 2, 1, i ), Position_Color_Var_Level( 2, 2, 1, i ), 'lineto'
          write( File_Number, * ) Position_Color_Var_Level( 1, 1, 1, i ), Position_Color_Var_Level( 2, 1, 1, i ), 'lineto'
          write( File_Number, * )'closepath'
          write( File_Number, * ) RGB( 1, i ), RGB( 2, i ), RGB( 3, i ), 'setrgbcolor' 
          write( File_Number, * )'fill'
          write( File_Number, * )'grestore'
          write( File_Number, * )'newpath'
     end do
     !write( File_Number, * )'# Color Var --------------------------------------------'

     deallocate( Position_Color_Var_Level )

     write( File_Number, * ) '/Helvetica findfont 18 scalefont setfont'
     write( File_Number, * ) Position_Maximum_Color_Var( 1, 1, 1 ) -60d0, &!*Width_Var( 1 ), &
                             Position_Maximum_Color_Var( 2, 1, 1 ) -Width_Var( 2 )/10d0, 'moveto'
     write( File_Number, 73 )'(', Minimum_Value_Plotted, ') show'
     !write( File_Number, * ) Minimum_Value_Plotted 

     write( File_Number, * ) Position_Maximum_Color_Var( 1, 2, 1 ) -60d0, &!*Width_Var( 1 ), &
                             Position_Maximum_Color_Var( 2, 2, 1 ) +Width_Var( 2 )/20d0, 'moveto'
     write( File_Number, 73 )'(', Maximum_Value_Plotted, ') show'
     !write( File_Number, * ) Maximum_Value_Plotted

     73 format(A1,es15.2,A6) 


     deallocate( Width_Var )

     return
end subroutine Plot_Color_Var_Right_Side



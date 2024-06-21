
subroutine Plot_Mesh_Configuration &
       ( Position_Node, Number_Node, &
         Index_Element_2_Node, Number_Element, &
         Class_Element, & 
         ID_Element_Material, ID_Element_OpenRegion, ID_Element_FixedDomain, ID_Element_Air, &
         ID_Element_Material_Exterior, ID_Element_Air_Exterior, &
         ID_Element_PML_Xp, ID_Element_PML_Xm, ID_Element_PML_Yp, ID_Element_PML_Ym, &
         ID_Element_PML_XpYp,  ID_Element_PML_XmYp, ID_Element_PML_XpYm,  ID_Element_PML_XmYm )

    !$use omp_lib
    use Parameters
    implicit none

    integer, intent(in) :: Number_Node, Number_Element
    double precision, intent(in) :: Position_Node( 2, Number_Node )
    integer, intent(in) :: Index_Element_2_Node( 4, Number_Element )
    integer, intent(in) :: Class_Element( Number_Element )
    integer, intent(in) :: ID_Element_Material, ID_Element_OpenRegion, ID_Element_FixedDomain, ID_Element_Air 
    integer, intent(in) :: ID_Element_Material_Exterior, ID_Element_Air_Exterior  
    integer, intent(in) :: ID_Element_PML_Xp, ID_Element_PML_Xm, ID_Element_PML_Yp, ID_Element_PML_Ym 
    integer, intent(in) :: ID_Element_PML_XpYp,  ID_Element_PML_XmYp, ID_Element_PML_XpYm,  ID_Element_PML_XmYm 

    integer :: i, j

    double precision, allocatable, dimension(:) :: Position_Minimum 
    double precision, allocatable, dimension(:) :: Position_Maximum 

    double precision, allocatable, dimension(:,:) :: Position_Node_Plot_Mesh 

    !=======================================================
    write(*,*)'      call Plot_Mesh_Configuration '
    !=======================================================
    !================================================
    !write(*,*)'         Initial Settings'
    !================================================

    allocate( Position_Minimum( 2 ) )
    allocate( Position_Maximum( 2 ) )
 
    !Minimum
    Position_Minimum( 1 )= Integer_Initialization_Plus 
    Position_Minimum( 2 )= Integer_Initialization_Plus
    !Maximum
    Position_Maximum( 1 )= Integer_Initialization_Minus
    Position_Maximum( 2 )= Integer_Initialization_Minus

    do j= 1, Number_Element
       do i= 1, 4
          if( Position_Minimum( 1 ) > Position_Node( 1, Index_Element_2_Node( i, j ) ) )then
             Position_Minimum( 1 ) = Position_Node( 1, Index_Element_2_Node( i, j ) )
          end if

          if( Position_Maximum( 1 ) < Position_Node( 1, Index_Element_2_Node( i, j ) ) )then
             Position_Maximum( 1 ) = Position_Node( 1, Index_Element_2_Node( i, j ) )
          end if

          if( Position_Minimum( 2 ) > Position_Node( 2, Index_Element_2_Node( i, j ) ) )then
             Position_Minimum( 2 ) = Position_Node( 2, Index_Element_2_Node( i, j ) )
          end if

          if( Position_Maximum( 2 ) < Position_Node( 2, Index_Element_2_Node( i, j ) ) )then
             Position_Maximum( 2 ) = Position_Node( 2, Index_Element_2_Node( i, j ) )
          end if
       end do
    end do

    allocate( Position_Node_Plot_Mesh( 2, Number_Node ) )
 
    do j= 1, Number_Node
       Position_Node_Plot_Mesh( 1, j )= ( Position_Node( 1, j ) -Position_Minimum( 1 ) )  &
                         /( Position_Maximum( 1 ) -Position_Minimum( 1 ) )*dble( WidthPS_X ) 
       Position_Node_Plot_Mesh( 2, j )= ( Position_Node( 2, j ) -Position_Minimum( 2 ) )  &
                         /( Position_Maximum( 2 ) -Position_Minimum( 2 ) )*dble( WidthPS_Y ) 
    end do

    do j= 1, Number_Node
       Position_Node_Plot_Mesh( 1, j )= Position_Node_Plot_Mesh( 1, j ) +Translation_X 
       Position_Node_Plot_Mesh( 2, j )= Position_Node_Plot_Mesh( 2, j ) +Translation_Y
    end do

    deallocate( Position_Minimum )
    deallocate( Position_Maximum )

    !================================================
    !write(*,*)'         Plot Configuration & Mesh'
    !================================================
    open( 200, file='Configuration_Mesh.ps', status='replace' )
     write(200,*)'%!PS-Adobe-2.0'
     write(200,*)'%BoundingBox: 0 0', WidthPS_X, WidthPS_Y
     write(200,*)'%EndHeader'
     write(200,*)'stroke'
     write(200,*)'0.02 setlinewidth'
  
       write(200,*)dble( Translation_X ), dble( Translation_Y ), 'moveto'
       write(200,*)dble( Translation_X ) +WidthPS_X, dble( Translation_Y ) +0D0, 'lineto'
       write(200,*)dble( Translation_X ) +WidthPS_X, dble( Translation_Y ) +WidthPS_Y, 'lineto'
       write(200,*)dble( Translation_X ) +0D0, dble( Translation_Y ) +WidthPS_Y, 'lineto'
       write(200,*)dble( Translation_X ), dble( Translation_Y ), 'lineto'
       write(200,*)'stroke'

       write(200,*)'0.3 setgray'

       do j= 1, Number_Element
          if( Class_Element(j) == ID_Element_Material )then

             write(200,*)'gsave'
             write(200,*) Position_Node_Plot_Mesh( 1, Index_Element_2_Node( 1, j ) ), &
                    Position_Node_Plot_Mesh( 2, Index_Element_2_Node( 1, j ) ), &
                    'moveto'
 
             if( Index_Element_2_Node( 3, j )==Index_Element_2_Node( 4, j ) )then
                do i= 2, 3
                   write(200,*) Position_Node_Plot_Mesh( 1, Index_Element_2_Node( i, j ) ), &
                          Position_Node_Plot_Mesh( 2, Index_Element_2_Node( i, j ) ), &
                          'lineto'
                end do
  
             else
                do i= 2, 4
                   write(200,*) Position_Node_Plot_Mesh( 1, Index_Element_2_Node( i, j ) ), &
                          Position_Node_Plot_Mesh( 2, Index_Element_2_Node( i, j ) ), &
                          'lineto'
                end do
 
             end if
 
             write(200,*) Position_Node_Plot_Mesh( 1, Index_Element_2_Node( 1, j ) ), &
                    Position_Node_Plot_Mesh( 2, Index_Element_2_Node( 1, j ) ), &
                    'lineto'

             write(200,*)'closepath'
             !write(200,*)'0 255 255 setrgbcolor' ! cyan
             write(200,*)'0 0 255 setrgbcolor' ! blue
             !write(200,*)'1 0 0 setrgbcolor' ! red 
             write(200,*)'fill'
             write(200,*)'grestore'
             write(200,*)'newpath'

          else if( Class_Element(j) == ID_Element_Air )then

             write(200,*)'gsave'
             write(200,*) Position_Node_Plot_Mesh( 1, Index_Element_2_Node( 1, j ) ), &
                    Position_Node_Plot_Mesh( 2, Index_Element_2_Node( 1, j ) ), &
                    'moveto'
 
             if( Index_Element_2_Node( 3, j )==Index_Element_2_Node( 4, j ) )then
                do i= 2, 3
                   write(200,*) Position_Node_Plot_Mesh( 1, Index_Element_2_Node( i, j ) ), &
                          Position_Node_Plot_Mesh( 2, Index_Element_2_Node( i, j ) ), &
                          'lineto'
                end do
  
             else
                do i= 2, 4
                   write(200,*) Position_Node_Plot_Mesh( 1, Index_Element_2_Node( i, j ) ), &
                          Position_Node_Plot_Mesh( 2, Index_Element_2_Node( i, j ) ), &
                          'lineto'
                end do
 
             end if
 
             write(200,*) Position_Node_Plot_Mesh( 1, Index_Element_2_Node( 1, j ) ), &
                    Position_Node_Plot_Mesh( 2, Index_Element_2_Node( 1, j ) ), &
                    'lineto'

             write(200,*)'closepath'
             write(200,*)'0 255 255 setrgbcolor' ! cyan
             !write(200,*)'0 0 255 setrgbcolor' ! blue
             !write(200,*)'1 0 0 setrgbcolor' ! red 
             write(200,*)'fill'
             write(200,*)'grestore'
             write(200,*)'newpath'

          else if( Class_Element(j) == ID_Element_FixedDomain )then

             write(200,*)'gsave'
             write(200,*) Position_Node_Plot_Mesh( 1, Index_Element_2_Node( 1, j ) ), &
                    Position_Node_Plot_Mesh( 2, Index_Element_2_Node( 1, j ) ), &
                    'moveto'
 
             if( Index_Element_2_Node( 3, j )==Index_Element_2_Node( 4, j ) )then
                do i= 2, 3
                   write(200,*) Position_Node_Plot_Mesh( 1, Index_Element_2_Node( i, j ) ), &
                          Position_Node_Plot_Mesh( 2, Index_Element_2_Node( i, j ) ), &
                          'lineto'
                end do
  
             else
                do i= 2, 4
                   write(200,*) Position_Node_Plot_Mesh( 1, Index_Element_2_Node( i, j ) ), &
                          Position_Node_Plot_Mesh( 2, Index_Element_2_Node( i, j ) ), &
                          'lineto'
                end do
 
             end if
  
             write(200,*) Position_Node_Plot_Mesh( 1, Index_Element_2_Node( 1, j ) ), &
                    Position_Node_Plot_Mesh( 2, Index_Element_2_Node( 1, j ) ), &
                    'lineto'

             write(200,*)'closepath'
             !write(200,*)'0 255 255 setrgbcolor' ! cyan
             !write(200,*)'0 0 255 setrgbcolor' ! blue
             write(200,*)'1 0 0 setrgbcolor' ! red 
             write(200,*)'fill'
             write(200,*)'grestore'
             write(200,*)'newpath'

          else if( Class_Element(j) == ID_Element_OpenRegion )then

             write(200,*)'gsave'
             write(200,*) Position_Node_Plot_Mesh( 1, Index_Element_2_Node( 1, j ) ), &
                    Position_Node_Plot_Mesh( 2, Index_Element_2_Node( 1, j ) ), &
                    'moveto'
 
             if( Index_Element_2_Node( 3, j )==Index_Element_2_Node( 4, j ) )then
                do i= 2, 3
                   write(200,*) Position_Node_Plot_Mesh( 1, Index_Element_2_Node( i, j ) ), &
                          Position_Node_Plot_Mesh( 2, Index_Element_2_Node( i, j ) ), &
                          'lineto'
                end do
             else
                do i= 2, 4
                   write(200,*) Position_Node_Plot_Mesh( 1, Index_Element_2_Node( i, j ) ), &
                          Position_Node_Plot_Mesh( 2, Index_Element_2_Node( i, j ) ), &
                          'lineto'
                end do
             end if
  
             write(200,*) Position_Node_Plot_Mesh( 1, Index_Element_2_Node( 1, j ) ), &
                    Position_Node_Plot_Mesh( 2, Index_Element_2_Node( 1, j ) ), &
                    'lineto'

             write(200,*)'closepath'
             !write(200,*)'0 255 255 setrgbcolor' ! cyan
             !write(200,*)'0 0 255 setrgbcolor' ! blue
             !write(200,*)'1 0 0 setrgbcolor' ! red 
             write(200,*)'0.1 0.1 0.1 setrgbcolor' ! yellow 
             write(200,*)'fill'
             write(200,*)'grestore'
             write(200,*)'newpath'

          else if( Class_Element(j) == ID_Element_PML_Xp .or. Class_Element(j) == ID_Element_PML_Xm )then

             write(200,*)'gsave'
             write(200,*) Position_Node_Plot_Mesh( 1, Index_Element_2_Node( 1, j ) ), &
                    Position_Node_Plot_Mesh( 2, Index_Element_2_Node( 1, j ) ), &
                    'moveto'
 
             if( Index_Element_2_Node( 3, j )==Index_Element_2_Node( 4, j ) )then
                do i= 2, 3
                   write(200,*) Position_Node_Plot_Mesh( 1, Index_Element_2_Node( i, j ) ), &
                          Position_Node_Plot_Mesh( 2, Index_Element_2_Node( i, j ) ), &
                          'lineto'
                end do
             else
                do i= 2, 4
                   write(200,*) Position_Node_Plot_Mesh( 1, Index_Element_2_Node( i, j ) ), &
                          Position_Node_Plot_Mesh( 2, Index_Element_2_Node( i, j ) ), &
                          'lineto'
                end do
             end if
  
             write(200,*) Position_Node_Plot_Mesh( 1, Index_Element_2_Node( 1, j ) ), &
                    Position_Node_Plot_Mesh( 2, Index_Element_2_Node( 1, j ) ), &
                    'lineto'

             write(200,*)'closepath'
             !write(200,*)'0 255 0 setrgbcolor' ! cyan
             !write(200,*)'0 255 255 setrgbcolor' ! cyan
             !write(200,*)'0 0 255 setrgbcolor' ! blue
             write(200,*)'255 0 255 setrgbcolor' ! red 
             !write(200,*)'1 1 0 setrgbcolor' ! yellow 
             write(200,*)'fill'
             write(200,*)'grestore'
             write(200,*)'newpath'

          else if( Class_Element(j) == ID_Element_PML_Yp .or. Class_Element(j) == ID_Element_PML_Ym )then

             write(200,*)'gsave'
             write(200,*) Position_Node_Plot_Mesh( 1, Index_Element_2_Node( 1, j ) ), &
                    Position_Node_Plot_Mesh( 2, Index_Element_2_Node( 1, j ) ), &
                    'moveto'
 
             if( Index_Element_2_Node( 3, j )==Index_Element_2_Node( 4, j ) )then
                do i= 2, 3
                   write(200,*) Position_Node_Plot_Mesh( 1, Index_Element_2_Node( i, j ) ), &
                          Position_Node_Plot_Mesh( 2, Index_Element_2_Node( i, j ) ), &
                          'lineto'
                end do
             else
                do i= 2, 4
                   write(200,*) Position_Node_Plot_Mesh( 1, Index_Element_2_Node( i, j ) ), &
                          Position_Node_Plot_Mesh( 2, Index_Element_2_Node( i, j ) ), &
                          'lineto'
                end do
             end if
  
             write(200,*) Position_Node_Plot_Mesh( 1, Index_Element_2_Node( 1, j ) ), &
                    Position_Node_Plot_Mesh( 2, Index_Element_2_Node( 1, j ) ), &
                    'lineto'

             write(200,*)'closepath'
             !write(200,*)'0 255 0 setrgbcolor' ! cyan
             !write(200,*)'0 255 255 setrgbcolor' ! cyan
             !write(200,*)'0 0 255 setrgbcolor' ! blue
             write(200,*)'255 1 255 setrgbcolor' ! red 
             !write(200,*)'1 1 0 setrgbcolor' ! yellow 
             write(200,*)'fill'
             write(200,*)'grestore'
             write(200,*)'newpath'

          else if( Class_Element(j) == ID_Element_PML_XpYp .or. &
               Class_Element(j) == ID_Element_PML_XpYm .or. &
               Class_Element(j) == ID_Element_PML_XmYp .or. &
               Class_Element(j) == ID_Element_PML_XmYm )then

             write(200,*)'gsave'
             write(200,*) Position_Node_Plot_Mesh( 1, Index_Element_2_Node( 1, j ) ), &
                    Position_Node_Plot_Mesh( 2, Index_Element_2_Node( 1, j ) ), &
                    'moveto'
 
             if( Index_Element_2_Node( 3, j )==Index_Element_2_Node( 4, j ) )then
                do i= 2, 3
                   write(200,*) Position_Node_Plot_Mesh( 1, Index_Element_2_Node( i, j ) ), &
                          Position_Node_Plot_Mesh( 2, Index_Element_2_Node( i, j ) ), &
                          'lineto'
                end do
             else
                do i= 2, 4
                   write(200,*) Position_Node_Plot_Mesh( 1, Index_Element_2_Node( i, j ) ), &
                          Position_Node_Plot_Mesh( 2, Index_Element_2_Node( i, j ) ), &
                          'lineto'
                end do
             end if
  
             write(200,*) Position_Node_Plot_Mesh( 1, Index_Element_2_Node( 1, j ) ), &
                    Position_Node_Plot_Mesh( 2, Index_Element_2_Node( 1, j ) ), &
                    'lineto'

             write(200,*)'closepath'
             !write(200,*)'0 255 0 setrgbcolor' ! green 
             !write(200,*)'0 255 255 setrgbcolor' ! cyan
             !write(200,*)'0 0 255 setrgbcolor' ! blue
             write(200,*)'255 0 0 setrgbcolor' ! red 
             !write(200,*)'1 1 0 setrgbcolor' ! yellow 
             write(200,*)'fill'
             write(200,*)'grestore'
             write(200,*)'newpath'

          else 

             write(200,*)'gsave'
             write(200,*) Position_Node_Plot_Mesh( 1, Index_Element_2_Node( 1, j ) ), &
                    Position_Node_Plot_Mesh( 2, Index_Element_2_Node( 1, j ) ), &
                    'moveto'
 
             if( Index_Element_2_Node( 3, j )==Index_Element_2_Node( 4, j ) )then
                do i= 2, 3
                   write(200,*) Position_Node_Plot_Mesh( 1, Index_Element_2_Node( i, j ) ), &
                          Position_Node_Plot_Mesh( 2, Index_Element_2_Node( i, j ) ), &
                          'lineto'
                end do
             else
                do i= 2, 4
                   write(200,*) Position_Node_Plot_Mesh( 1, Index_Element_2_Node( i, j ) ), &
                          Position_Node_Plot_Mesh( 2, Index_Element_2_Node( i, j ) ), &
                          'lineto'
                end do
             end if
  
             write(200,*) Position_Node_Plot_Mesh( 1, Index_Element_2_Node( 1, j ) ), &
                    Position_Node_Plot_Mesh( 2, Index_Element_2_Node( 1, j ) ), &
                    'lineto'

             write(200,*)'closepath'
             write(200,*)'0 0 0 setrgbcolor' ! black 
             !write(200,*)'0 255 255 setrgbcolor' ! cyan
             !write(200,*)'0 0 255 setrgbcolor' ! blue
            ! write(200,*)'255 0 0 setrgbcolor' ! red 
             !write(200,*)'1 1 0 setrgbcolor' ! yellow 
             write(200,*)'fill'
             write(200,*)'grestore'
             write(200,*)'newpath'


          end if
       end do

       do j= 1, Number_Element
          write(200,*) Position_Node_Plot_Mesh( 1, Index_Element_2_Node( 1, j ) ), &
                 Position_Node_Plot_Mesh( 2, Index_Element_2_Node( 1, j ) ), &
                 'moveto'

          if( Index_Element_2_Node( 3, j )==Index_Element_2_Node( 4, j ) )then
             do i= 2, 3
                write(200,*) Position_Node_Plot_Mesh( 1, Index_Element_2_Node( i, j ) ), &
                       Position_Node_Plot_Mesh( 2, Index_Element_2_Node( i, j ) ), &
                       'lineto'
             end do
          else
             do i= 2, 4
                write(200,*) Position_Node_Plot_Mesh( 1, Index_Element_2_Node( i, j ) ), &
                       Position_Node_Plot_Mesh( 2, Index_Element_2_Node( i, j ) ), &
                       'lineto'
             end do
          end if

          write(200,*) Position_Node_Plot_Mesh( 1, Index_Element_2_Node( 1, j ) ), &
                 Position_Node_Plot_Mesh( 2, Index_Element_2_Node( 1, j ) ), &
                 'lineto'
          write(200,*)'stroke'
       end do


       write(*,*)'         ==> Configuration_Mesh.ps'

       close( 200 )

       deallocate( Position_Node_Plot_Mesh )

    return
end subroutine Plot_Mesh_Configuration 



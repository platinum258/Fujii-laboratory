
subroutine Plot_Configuration &
       ( Position_Node, Number_Node, &
         Index_Element_2_Node, Number_Element, &
         Class_Element, & 
         ID_Element_Material, ID_Element_OpenRegion, ID_Element_FixedDomain, ID_Element_Base_Material, &
         ID_Element_Material_Exterior, ID_Element_Base_Material_Exterior )   

   !$use omp_lib
   use Parameters
   implicit none

   integer, intent(in) :: Number_Node, Number_Element
   double precision, intent(in) :: Position_Node( 2, Number_Node )
   integer, intent(in) :: Index_Element_2_Node( 4, Number_Element )
   integer, intent(in) :: Class_Element( Number_Element )
   integer, intent(in) :: ID_Element_Material, ID_Element_OpenRegion, ID_Element_FixedDomain, ID_Element_Base_Material 
   integer, intent(in) :: ID_Element_Material_Exterior, ID_Element_Base_Material_Exterior  
   
   integer :: i, j
   integer :: Counter_Percent
   
   double precision, allocatable, dimension(:) :: Position_Minimum 
   double precision, allocatable, dimension(:) :: Position_Maximum 
   
   double precision, allocatable, dimension(:,:) :: Position_Node_Plot_Mesh 

   !================================================
   write(*,*)'===================================================='
   write(*,*)'      call Plot_Configuration'
   write(*,*)'===================================================='
   !================================================
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
    !write(*,*)'         Plot Configuration'
    !================================================
   
    open( 201, file='Configuration.ps', status='replace' )
     write(201,*)'%!PS-Adobe-2.0'
     write(201,*)'%BoundingBox: 0 0', WidthPS_X, WidthPS_Y
     write(201,*)'%EndHeader'
     write(201,*)'stroke'
     write(201,*)'0.1 setlinewidth'
   
     Counter_Percent= 1 

     do j= 1, Number_Element

      if( j== Number_Element/100*Counter_Percent )then
         !write(*,'(a,i3,a)', advance='no')'[', Counter_Percent, '%]'
         write(*,'(a,i3,a,$)')'[', Counter_Percent, '%]'
         Counter_Percent= Counter_Percent +1
      end if 

      !====================================================================================
      if( Class_Element(j) == ID_Element_Material )then
      !====================================================================================
   
       write(201,*)'gsave'
       write(201,*) Position_Node_Plot_Mesh( 1, Index_Element_2_Node( 1, j ) ), &
              Position_Node_Plot_Mesh( 2, Index_Element_2_Node( 1, j ) ), &
              'moveto'
    
           if( Index_Element_2_Node( 3, j )==Index_Element_2_Node( 4, j ) )then
              do i= 2, 3
                 write(201,*) Position_Node_Plot_Mesh( 1, Index_Element_2_Node( i, j ) ), &
                        Position_Node_Plot_Mesh( 2, Index_Element_2_Node( i, j ) ), &
                        'lineto'
              end do
           else
              do i= 2, 4
                 write(201,*) Position_Node_Plot_Mesh( 1, Index_Element_2_Node( i, j ) ), &
                        Position_Node_Plot_Mesh( 2, Index_Element_2_Node( i, j ) ), &
                        'lineto'
              end do
           end if
   
   
       write(201,*) Position_Node_Plot_Mesh( 1, Index_Element_2_Node( 1, j ) ), &
              Position_Node_Plot_Mesh( 2, Index_Element_2_Node( 1, j ) ), &
              'lineto'
   
       write(201,*)'closepath'
       if( Flag_Thermal_Device==0 ) write(201,*)'0.765625  0.4375   0.0 setrgbcolor'
       if( Flag_Thermal_Device==10 ) write(201,*)'0.2 0.2 0.2 setrgbcolor'
       if( Flag_Thermal_Device==20 ) write(201,*)'1 1 1 setrgbcolor'
       if( Flag_Thermal_Device==21 ) write(201,*)'1 1 1 setrgbcolor'
       if( Flag_Thermal_Device==30 ) write(201,*)'1 1 1 setrgbcolor'
       write(201,*)'fill'
       write(201,*)'grestore'
       write(201,*)'newpath'
   
      !====================================================================================
      else if( Class_Element(j) == ID_Element_Base_Material )then
      !====================================================================================
   
       write(201,*)'gsave'
       write(201,*) Position_Node_Plot_Mesh( 1, Index_Element_2_Node( 1, j ) ), &
              Position_Node_Plot_Mesh( 2, Index_Element_2_Node( 1, j ) ), &
              'moveto'
    
           if( Index_Element_2_Node( 3, j )==Index_Element_2_Node( 4, j ) )then
              do i= 2, 3
                 write(201,*) Position_Node_Plot_Mesh( 1, Index_Element_2_Node( i, j ) ), &
                        Position_Node_Plot_Mesh( 2, Index_Element_2_Node( i, j ) ), &
                        'lineto'
              end do
           else
              do i= 2, 4
                 write(201,*) Position_Node_Plot_Mesh( 1, Index_Element_2_Node( i, j ) ), &
                        Position_Node_Plot_Mesh( 2, Index_Element_2_Node( i, j ) ), &
                        'lineto'
              end do
           end if
   
    
       write(201,*) Position_Node_Plot_Mesh( 1, Index_Element_2_Node( 1, j ) ), &
              Position_Node_Plot_Mesh( 2, Index_Element_2_Node( 1, j ) ), &
              'lineto'
   
       write(201,*)'closepath'
       if( Flag_Thermal_Device==0 ) write(201,*)'0.2 0.2 0.2 setrgbcolor'
       if( Flag_Thermal_Device==10 ) write(201,*)'0.765625  0.4375   0.0 setrgbcolor' ! yellow
       if( Flag_Thermal_Device==20 ) write(201,*)'0.97265625 0.7109375 0.0390625 setrgbcolor' ! yellow
       if( Flag_Thermal_Device==21 ) write(201,*)'0.97265625 0.7109375 0.0390625 setrgbcolor' ! yellow
       if( Flag_Thermal_Device==30 ) write(201,*)'0.97265625 0.7109375 0.0390625 setrgbcolor' ! yellow
       !write(201,*)'0 0 1 setrgbcolor' ! blue
       !write(201,*)'1 0 0 setrgbcolor' ! red 
       write(201,*)'fill'
       write(201,*)'grestore'
       write(201,*)'newpath'
   
      !====================================================================================
      else if( Class_Element(j) == ID_Element_FixedDomain )then
      !====================================================================================
   
       write(201,*)'gsave'
       write(201,*) Position_Node_Plot_Mesh( 1, Index_Element_2_Node( 1, j ) ), &
              Position_Node_Plot_Mesh( 2, Index_Element_2_Node( 1, j ) ), &
              'moveto'
    
           if( Index_Element_2_Node( 3, j )==Index_Element_2_Node( 4, j ) )then
              do i= 2, 3
                 write(201,*) Position_Node_Plot_Mesh( 1, Index_Element_2_Node( i, j ) ), &
                        Position_Node_Plot_Mesh( 2, Index_Element_2_Node( i, j ) ), &
                        'lineto'
              end do
           else
              do i= 2, 4
                 write(201,*) Position_Node_Plot_Mesh( 1, Index_Element_2_Node( i, j ) ), &
                        Position_Node_Plot_Mesh( 2, Index_Element_2_Node( i, j ) ), &
                        'lineto'
              end do
           end if
   
    
       write(201,*) Position_Node_Plot_Mesh( 1, Index_Element_2_Node( 1, j ) ), &
              Position_Node_Plot_Mesh( 2, Index_Element_2_Node( 1, j ) ), &
              'lineto'
   
       write(201,*)'closepath'
       !write(201,*)'0 1 1 setrgbcolor' ! cyan
       if( Flag_Thermal_Device==220 )then
          write(201,*)'0 0 1 setrgbcolor' ! Blue

       else if( Flag_Thermal_Device==100 )then

          write(201,*)'0 1 0 setrgbcolor' ! Yellow-Green

       else if( Flag_Thermal_Device==0  )then
          !write(201,*)'0 0 0 setrgbcolor' ! Black 
          !write(201,*)'0.57 0.57 0.57 setrgbcolor' ! Gray
          write(201,*)'1 1 1 setrgbcolor' ! Black 
       else if( Flag_Thermal_Device==10  )then
          write(201,*)'1 1 1 setrgbcolor' ! Black 
       else if( Flag_Thermal_Device==51 .or. Flag_Thermal_Device==50 )then
          write(201,*)'1 1 1 setrgbcolor' ! Black 
       end if

       if( Flag_Thermal_Device==0 ) write(201,*)'1 1 1 setrgbcolor'
       if( Flag_Thermal_Device==10 ) write(201,*)'1 1 1 setrgbcolor'
       if( Flag_Thermal_Device==20 ) write(201,*)'0.2 0.2 0.2 setrgbcolor' ! yellow
       if( Flag_Thermal_Device==21 ) write(201,*)'0.2 0.2 0.2 setrgbcolor' ! yellow
       if( Flag_Thermal_Device==30 ) write(201,*)'0.2 0.2 0.2 setrgbcolor' ! yellow

       write(201,*)'fill'
       write(201,*)'grestore'
       write(201,*)'newpath'
   
      !====================================================================================
      else if( Class_Element(j) == ID_Element_OpenRegion )then
      !====================================================================================
   
       write(201,*)'gsave'
       write(201,*) Position_Node_Plot_Mesh( 1, Index_Element_2_Node( 1, j ) ), &
              Position_Node_Plot_Mesh( 2, Index_Element_2_Node( 1, j ) ), &
              'moveto'
    
           if( Index_Element_2_Node( 3, j )==Index_Element_2_Node( 4, j ) )then
              do i= 2, 3
                 write(201,*) Position_Node_Plot_Mesh( 1, Index_Element_2_Node( i, j ) ), &
                        Position_Node_Plot_Mesh( 2, Index_Element_2_Node( i, j ) ), &
                        'lineto'
              end do
           else
              do i= 2, 4
                 write(201,*) Position_Node_Plot_Mesh( 1, Index_Element_2_Node( i, j ) ), &
                        Position_Node_Plot_Mesh( 2, Index_Element_2_Node( i, j ) ), &
                        'lineto'
              end do
           end if
   
    
       write(201,*) Position_Node_Plot_Mesh( 1, Index_Element_2_Node( 1, j ) ), &
              Position_Node_Plot_Mesh( 2, Index_Element_2_Node( 1, j ) ), &
              'lineto'
   
       write(201,*)'closepath'
       if( Flag_Thermal_Device==0 ) write(201,*)'0.2 0.2 0.2 setrgbcolor'
       if( Flag_Thermal_Device==10 ) write(201,*)'0.765625  0.4375   0.0 setrgbcolor'
       !if( Flag_Thermal_Device==20 ) write(201,*)'0.2 0.2 0.2 setrgbcolor' ! yellow
       if( Flag_Thermal_Device==20 ) write(201,*)'1 1 1 setrgbcolor' ! yellow
       if( Flag_Thermal_Device==21 ) write(201,*)'1 1 1 setrgbcolor' ! yellow
       if( Flag_Thermal_Device==30 ) write(201,*)'1 1 1 setrgbcolor' ! yellow
       write(201,*)'fill'
       write(201,*)'grestore'
       write(201,*)'newpath'
   
      !====================================================================================
      else if( Class_Element(j) == ID_Element_Material_Exterior )then
      !====================================================================================
   
       write(201,*)'gsave'
       write(201,*) Position_Node_Plot_Mesh( 1, Index_Element_2_Node( 1, j ) ), &
              Position_Node_Plot_Mesh( 2, Index_Element_2_Node( 1, j ) ), &
              'moveto'
    
           if( Index_Element_2_Node( 3, j )==Index_Element_2_Node( 4, j ) )then
              do i= 2, 3
                 write(201,*) Position_Node_Plot_Mesh( 1, Index_Element_2_Node( i, j ) ), &
                        Position_Node_Plot_Mesh( 2, Index_Element_2_Node( i, j ) ), &
                        'lineto'
              end do
           else
              do i= 2, 4
                 write(201,*) Position_Node_Plot_Mesh( 1, Index_Element_2_Node( i, j ) ), &
                        Position_Node_Plot_Mesh( 2, Index_Element_2_Node( i, j ) ), &
                        'lineto'
              end do
           end if
   
    
       write(201,*) Position_Node_Plot_Mesh( 1, Index_Element_2_Node( 1, j ) ), &
              Position_Node_Plot_Mesh( 2, Index_Element_2_Node( 1, j ) ), &
              'lineto'
   
       write(201,*)'closepath'
       if( Flag_Thermal_Device==0 ) write(201,*)'0.2 0.2 0.2 setrgbcolor'
       if( Flag_Thermal_Device==10 ) write(201,*)'0.765625  0.4375   0.0 setrgbcolor'
       !if( Flag_Thermal_Device==20 ) write(201,*)'0.2 0.2 0.2 setrgbcolor' ! yellow
       if( Flag_Thermal_Device==20 ) write(201,*)'1 1 1 setrgbcolor' ! yellow
       if( Flag_Thermal_Device==21 ) write(201,*)'1 1 1 setrgbcolor' ! yellow
       if( Flag_Thermal_Device==30 ) write(201,*)'1 1 1 setrgbcolor' ! yellow
       write(201,*)'fill'
       write(201,*)'grestore'
       write(201,*)'newpath'
    
      !====================================================================================
      else if( Class_Element(j) == ID_Element_Base_Material_Exterior )then
      !====================================================================================
   
       write(201,*)'gsave'
       write(201,*) Position_Node_Plot_Mesh( 1, Index_Element_2_Node( 1, j ) ), &
              Position_Node_Plot_Mesh( 2, Index_Element_2_Node( 1, j ) ), &
              'moveto'
    
           if( Index_Element_2_Node( 3, j )==Index_Element_2_Node( 4, j ) )then
              do i= 2, 3
                 write(201,*) Position_Node_Plot_Mesh( 1, Index_Element_2_Node( i, j ) ), &
                        Position_Node_Plot_Mesh( 2, Index_Element_2_Node( i, j ) ), &
                        'lineto'
              end do
           else
              do i= 2, 4
                 write(201,*) Position_Node_Plot_Mesh( 1, Index_Element_2_Node( i, j ) ), &
                        Position_Node_Plot_Mesh( 2, Index_Element_2_Node( i, j ) ), &
                        'lineto'
              end do
           end if
   
    
       write(201,*) Position_Node_Plot_Mesh( 1, Index_Element_2_Node( 1, j ) ), &
              Position_Node_Plot_Mesh( 2, Index_Element_2_Node( 1, j ) ), &
              'lineto'
   
       write(201,*)'closepath'
       if( Flag_Thermal_Device==0 ) write(201,*)'0.2 0.2 0.2 setrgbcolor'
       if( Flag_Thermal_Device==10 ) write(201,*)'0.765625  0.4375   0.0 setrgbcolor'
       !if( Flag_Thermal_Device==20 ) write(201,*)'0.2 0.2 0.2 setrgbcolor' ! yellow
       if( Flag_Thermal_Device==20 ) write(201,*)'1 1 1 setrgbcolor' ! yellow
       if( Flag_Thermal_Device==21 ) write(201,*)'1 1 1 setrgbcolor' ! yellow
       if( Flag_Thermal_Device==30 ) write(201,*)'1 1 1 setrgbcolor' ! yellow
       write(201,*)'fill'
       write(201,*)'grestore'
       write(201,*)'newpath'
   
      !====================================================================================
      end if
      !====================================================================================

     end do

     write(201,*)dble( Translation_X ), dble( Translation_Y ), 'moveto'
     write(201,*)dble( Translation_X ) +WidthPS_X, dble( Translation_Y ) +0D0, 'lineto'
     write(201,*)dble( Translation_X ) +WidthPS_X, dble( Translation_Y ) +WidthPS_Y, 'lineto'
     write(201,*)dble( Translation_X ) +0D0, dble( Translation_Y ) +WidthPS_Y, 'lineto'
     write(201,*)dble( Translation_X ), dble( Translation_Y ), 'lineto'
     write(201,*)'stroke'


      !==================================================================================================
      ! Plot Design Domain
      !==================================================================================================

      if( Flag_DesignDomain==1 )then
         if( Flag_Thermal_Device==50 )then
            !write( 201, * )'1 1 1 setrgbcolor 0.3 setlinewidth'
            write( 201, * )'0 0 0 setrgbcolor 1 setlinewidth'
            write( 201, * ) 'newpath'
            write( 201, * ) Position_Center_DesignDomain_X_PS, &
                            Position_Center_DesignDomain_Y_PS -Radius_DesignDomain_Plot_PS/2d0, &
                            Radius_DesignDomain_Plot_PS,' 0 180 arc'
            write( 201, * )'closepath'
            write( 201, * )'stroke'
         else
            !write( 201, * )'1 1 1 setrgbcolor 0.3 setlinewidth'
            write( 201, * )'0 0 0 setrgbcolor 1 setlinewidth'
            write( 201, * ) 'newpath'
            write( 201, * ) Position_Center_DesignDomain_X_PS, Position_Center_DesignDomain_Y_PS, &
                            Radius_DesignDomain_Plot_PS,' 0 360 arc'
            write( 201, * )'closepath'
            write( 201, * )'stroke'
         end if

      else if( Flag_DesignDomain==2 )then
         write( 201, * )'1 1 1 setrgbcolor 0.3 setlinewidth'
         write( 201, * ) Position_Corner_DesignDomain_X1_PS, Position_Corner_DesignDomain_Y1_PS, 'moveto' 
         write( 201, * ) Position_Corner_DesignDomain_X2_PS, Position_Corner_DesignDomain_Y1_PS, 'lineto' 
         write( 201, * ) Position_Corner_DesignDomain_X2_PS, Position_Corner_DesignDomain_Y2_PS, 'lineto' 
         write( 201, * ) Position_Corner_DesignDomain_X1_PS, Position_Corner_DesignDomain_Y2_PS, 'lineto' 
         write( 201, * ) Position_Corner_DesignDomain_X1_PS, Position_Corner_DesignDomain_Y1_PS, 'lineto' 
         write( 201, * )'stroke'
      end if

   
     write(201,*)'0.3 setgray'


    close( 201 )
   
    deallocate( Position_Node_Plot_Mesh )
   
   return
end subroutine Plot_Configuration



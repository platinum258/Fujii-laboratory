
subroutine Detect_Element_inwhich_Point_Exists &
   ( Number_Element, Position_Node_Plot_Element, Position_Line, &
     Element_Number_include_Center )

   implicit none

   integer, intent(in) :: Number_Element
   double precision, intent(in) :: Position_Line( 2 )
   double precision, intent(in) :: Position_Node_Plot_Element( 2, 3, Number_Element ) 
   integer, intent(out) :: Element_Number_include_Center

   integer e, i, j   
   double precision :: Vector1( 2, 3 ), Vector2( 2, 3 ) 
   integer :: Flag_Vector( 3 )
   double precision :: Distance_tmp
   double precision, allocatable, dimension(:) :: Distance_Center_2_Vector 
   double precision, allocatable, dimension(:,:) :: Position_Center_Element
   double precision :: Min_X, Max_X, Min_Y, Max_Y

   !==================================================================================================
   !write(*,*)'Detect Element in Which Point Exists'
   !==================================================================================================
   Element_Number_include_Center=923653459

   Min_X= minval( Position_Node_Plot_Element( 1, :, : ) )
   Max_X= maxval( Position_Node_Plot_Element( 1, :, : ) )
   Min_Y= minval( Position_Node_Plot_Element( 2, :, : ) )
   Max_Y= maxval( Position_Node_Plot_Element( 2, :, : ) )

   if( Position_Line( 1 ) < Min_X .or. Max_X < Position_Line( 1 ) )then
      write(*,*)'======================================================'
      write(*,*)'Error at X, out of range'
      write(*,*)'Position_Line( 1 )=', Position_Line( 1 ) 
      write(*,*)'Min_X=', Min_X
      write(*,*)'Max_X=', Max_X
      write(*,*)'Detect_Element_inwhich_Point_Exists.f90 33'
      write(*,*)'======================================================'
   end if
   if( Position_Line( 2 ) < Min_Y .or. Max_Y < Position_Line( 2 ) )then
      write(*,*)'======================================================'
      write(*,*)'Error at Y, out of range'
      write(*,*)'Position_Line( 2 )=', Position_Line( 2 ) 
      write(*,*)'Min_Y=', Min_Y
      write(*,*)'Max_Y=', Max_Y
      write(*,*)'Detect_Element_inwhich_Point_Exists.f90 38'
      write(*,*)'======================================================'
   end if

   do e= 1, Number_Element 
      do i= 1, 3 
         do j= 1, 2 
            if( i==3 )then
               Vector1( j, i )=Position_Node_Plot_Element( j, 3, e ) -Position_Line( j ) 
               Vector2( j, i )=Position_Node_Plot_Element( j, 1, e ) -Position_Line( j )
            else
               Vector1( j, i )=Position_Node_Plot_Element( j, i, e ) -Position_Line( j )
               Vector2( j, i )=Position_Node_Plot_Element( j, i+1, e ) -Position_Line( j )
            end if
         end do
      end do
      
      call Judge_Cross_Product( 3, Vector1, Vector2, 0, Flag_Vector )
      
      if( Flag_Vector( 1 )==1 .and. Flag_Vector( 2 )==1 .and. Flag_Vector( 3 )==1 )then
         Element_Number_include_Center= e
         go to 596
      else if( Flag_Vector( 1 )==-1 .and. Flag_Vector( 2 )==-1 .and. Flag_Vector( 3 )==-1 )then
         Element_Number_include_Center= e
         go to 596
      end if
   end do
      
   596  continue
   if( Element_Number_include_Center==923653459 )then
      do e= 1, Number_Element 
         do i= 1, 3 
            do j= 1, 2 
               if( i==3 )then
                  Vector1( j, i )=Position_Node_Plot_Element( j, 3, e ) -Position_Line( j ) 
                  Vector2( j, i )=Position_Node_Plot_Element( j, 1, e ) -Position_Line( j )
               else
                  Vector1( j, i )=Position_Node_Plot_Element( j, i, e ) -Position_Line( j )
                  Vector2( j, i )=Position_Node_Plot_Element( j, i+1, e ) -Position_Line( j )
               end if
            end do
         end do
         
         call Judge_Cross_Product( 3, Vector1, Vector2, 1, Flag_Vector )
         
         if( Flag_Vector( 1 )==1 .and. Flag_Vector( 2 )==1 .and. Flag_Vector( 3 )==1 )then
            Element_Number_include_Center= e
            go to 983
         else if( Flag_Vector( 1 )==-1 .and. Flag_Vector( 2 )==-1 .and. Flag_Vector( 3 )==-1 )then
            Element_Number_include_Center= e
            go to 983 
         end if
      end do
   end if

   983  continue
   if( Element_Number_include_Center==923653459 )then
      allocate( Position_Center_Element( 2, Number_Element ) )
   
      Position_Center_Element= 0.0d0
      
      do e= 1, Number_Element 
         do i= 1, 3 
            do j= 1, 2 
               Position_Center_Element( j, e )&
               = Position_Center_Element( j, e ) +Position_Node_Plot_Element( j, i, e ) 
            end do
         end do
      end do
      
      do e= 1, Number_Element 
         do i= 1, 2 
            Position_Center_Element( i, e )= Position_Center_Element( i, e )/3.0d0
         end do
      end do
      
      allocate( Distance_Center_2_Vector( Number_Element ) )
      
      do e= 1, Number_Element 
         Distance_Center_2_Vector( e )&
         = ( Position_Center_Element( 1, e ) -Position_Line( 1 ) )**2 & 
          +( Position_Center_Element( 2, e ) -Position_Line( 2 ) )**2 
      end do
      deallocate( Position_Center_Element )
   
      Distance_tmp= 1d10
      do e= 1, Number_Element 
         if( Distance_tmp >  Distance_Center_2_Vector( e ) )then
            Distance_tmp = Distance_Center_2_Vector( e )
            Element_Number_include_Center= e 
         end if
      end do
    
      deallocate( Distance_Center_2_Vector )
   end if

   !write(*,*) ' ' 
   !write(*,*) '     end Detect_Element_inwhich_Point_Exists' 
   return
end subroutine Detect_Element_inwhich_Point_Exists 



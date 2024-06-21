
subroutine Plot_Stream_Line &
     ( FN, Number_Element, Number_Line_X, Number_Line_Y, Ratio_Default, & 
       Position_Minimum_Plot, Position_Maximum_Plot, &
       Position_Node_Plot_Element, Gradient_in_Element, & 
       Optimization_Step )

  !use omp_lib
  use Parameters
  implicit none

  integer, intent(in) :: FN, Number_Element 
  integer, intent(in) :: Number_Line_X, Number_Line_Y
  double precision, intent(in) :: Ratio_Default
  double precision, intent(in) :: Position_Minimum_Plot( 2 ), Position_Maximum_Plot( 2 )
  double precision, intent(in) :: Gradient_in_Element( 2, Number_Element ) 
  double precision, intent(in) :: Position_Node_Plot_Element( 2, 3, Number_Element )
  integer, intent(in) :: Optimization_Step

  integer :: Loop_Line_Y, Loop_Line_X
  integer :: e, i, j

  integer :: Counter_Percent

  double precision, allocatable, dimension(:,:,:) :: Position_Line
  double precision, allocatable, dimension(:) :: Number_Line, Counter_Line
  double precision :: Line_Length_X, Line_Length_X_Default, Vector_Length


  double precision, allocatable, dimension(:,:) :: Position_Center_Element
  double precision, allocatable, dimension(:) :: Distance_Center_2_Vector 
  double precision, allocatable, dimension(:,:,:) :: Vector_tmp, Vector_Unit_tmp 
  !integer :: Number_Arrow, Number_Line_X, Number_Line_Y, Loop_Line_Y, Loop_Line_X
  !double precision :: Stemthick, StemLength, Headlength, Headangle, Headrad
  double precision :: Distance_tmp, Real_tmp
  double precision :: Vector1( 2, 3 ), Vector2( 2, 3 ) 
  integer :: Flag_Vector( 3 )
  integer, allocatable, dimension(:,:) :: Element_Number_include_Center
  !integer :: minloc
  character(len=256) :: Filename_Result

  character(len=Length_Character_Optimization_Step) :: Optimization_Step_Character
  character(len=8) :: Format_Filenumber
  
  write( Format_Filenumber, '(A2,I1,A1,I1,A1)' ) & 
  '(i', Length_Character_Optimization_Step,'.',Length_Character_Optimization_Step,')'
  
  write( Optimization_Step_Character, Format_Filenumber ) Optimization_Step
 
  Filename_Result= trim( "Heat_Flux_OC_"//trim(Optimization_Step_Character)//".ps" )

  !==================================================================================================
  write(*,*)'     call Plot_Stream_Line' 
  !==================================================================================================

  !Number_Arrow= Number_Line_X*Number_Line_Y
  !Line_Length_X=( Position_Maximum_Plot( 1 ) -Position_Minimum_Plot( 1 ) )/dble( Number_Line_X )
  Line_Length_X_Default=( Position_Maximum_Plot( 1 ) -Position_Minimum_Plot( 1 ) )/dble( Number_Line_X )
  Line_Length_X=Line_Length_X_Default/Ratio_Default

  !allocate( Position_Line( 2, Number_Line_Y, Number_Line_X+1 ) )
  allocate( Position_Line( 2, Number_Line_Y, Number_Line_X*2*int(Ratio_Default) ) )
  allocate( Number_Line( Number_Line_Y ) )
  allocate( Counter_Line( Number_Line_Y ) )
  Position_Line= 0.0d0
  Number_Line(:)= Number_Line_X+1
  Counter_Line(:)= 1 

  do j= 1, Number_Line_Y
     Position_Line( 1, j, 1 ) = Position_Maximum_Plot( 1 ) 
     Position_Line( 2, j, 1 ) &
     = Position_Minimum_Plot( 2 ) +( Position_Maximum_Plot( 2 ) -Position_Minimum_Plot( 2 ) )*( j-0.5d0 )/dble( Number_Line_Y )
  end do

  allocate( Element_Number_include_Center( Number_Line_Y, Number_Line_X+1 ) )
  allocate( Vector_tmp( 2, Number_Line_Y, Number_Line_X+1 ) )
  allocate( Vector_Unit_tmp( 2, Number_Line_Y, Number_Line_X+1 ) )

  Element_Number_include_Center= 0
  Counter_Percent= 1
  
  do Loop_Line_X= 1, ( Number_Line_X+1 )*int(Ratio_Default)

     if( dble( Loop_Line_X ) >= dble( ( Number_Line_X+1 )*int(Ratio_Default) )*Counter_Percent/100d0 )then
        !write(*,'(a,i3,a)', advance='no')'[', Counter_Percent, '%]'
        write(*,'(a,i3,a,$)')'[', Counter_Percent, '%]'
        Counter_Percent= Counter_Percent +1
     end if 

     do Loop_Line_Y= 1, Number_Line_Y

        590 continue

        !==================================================================================================
        !write(*,*)'Detect Element in Which Point Exists'
        !==================================================================================================
        do e= 1, Number_Element 
           do i= 1, 3 
              do j= 1, 2 
                 !if( i==3 )then
                 !   Vector1( j, i )=Position_Node_Plot_Element( j, 3, e ) -Position_Line( j, Loop_Line_Y, Loop_Line_X )
                 !   Vector2( j, i )=Position_Node_Plot_Element( j, 1, e ) -Position_Line( j, Loop_Line_Y, Loop_Line_X )
                 !else
                 !   Vector1( j, i )=Position_Node_Plot_Element( j, i, e ) -Position_Line( j, Loop_Line_Y, Loop_Line_X )
                 !   Vector2( j, i )=Position_Node_Plot_Element( j, i+1, e ) -Position_Line( j, Loop_Line_Y, Loop_Line_X )
                 !end if
!
                 if( i==3 )then
                    Vector1( j, i )=Position_Node_Plot_Element( j, 3, e ) -Position_Line( j, Loop_Line_Y, Counter_Line( Loop_Line_Y ) ) 
                    Vector2( j, i )=Position_Node_Plot_Element( j, 1, e ) -Position_Line( j, Loop_Line_Y, Counter_Line( Loop_Line_Y ) )
                 else
                    Vector1( j, i )=Position_Node_Plot_Element( j, i, e ) -Position_Line( j, Loop_Line_Y, Counter_Line( Loop_Line_Y ) )
                    Vector2( j, i )=Position_Node_Plot_Element( j, i+1, e ) -Position_Line( j, Loop_Line_Y, Counter_Line( Loop_Line_Y ) )
                 end if
              end do
           end do
   
           call Judge_Cross_Product( 3, Vector1, Vector2, 0, Flag_Vector )
   
           if( Flag_Vector( 1 )==1 .and. Flag_Vector( 2 )==1 .and. Flag_Vector( 3 )==1 )then
              Element_Number_include_Center( Loop_Line_Y, Loop_Line_X )= e
              go to 596
           else if( Flag_Vector( 1 )==-1 .and. Flag_Vector( 2 )==-1 .and. Flag_Vector( 3 )==-1 )then
              Element_Number_include_Center( Loop_Line_Y, Loop_Line_X )= e
              go to 596
           end if
        end do
   
        596  continue
        if( Element_Number_include_Center( Loop_Line_Y, Loop_Line_X )==0 )then
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
              = ( Position_Center_Element( 1, e ) -Position_Line( 1, Loop_Line_Y, Counter_Line( Loop_Line_Y ) ) )**2 & 
               +( Position_Center_Element( 2, e ) -Position_Line( 2, Loop_Line_Y, Counter_Line( Loop_Line_Y ) ) )**2 
           end do
           deallocate( Position_Center_Element )

           Distance_tmp= 1d10
           do e= 1, Number_Element 
              if( Distance_tmp >  Distance_Center_2_Vector( e ) )then
                 Distance_tmp = Distance_Center_2_Vector( e )
                 Element_Number_include_Center( Loop_Line_Y, Loop_Line_X )= e 
              end if
           end do
 
           !Element_Number_include_Center( Loop_Line_Y, Loop_Line_X )=minloc( Distance_Center_2_Vector )
  
           deallocate( Distance_Center_2_Vector )
        end if

        !==================================================================================================
        !write(*,*)'Compute Vector for Gradient Direction'
        !==================================================================================================
        do i= 1, 2 
           Vector_tmp( i, Loop_Line_Y, Loop_Line_X )= -Gradient_in_Element( i, Element_Number_include_Center( Loop_Line_Y, Loop_Line_X ) )
        end do
        Real_tmp= sqrt( Vector_tmp( 1, Loop_Line_Y, Loop_Line_X )*Vector_tmp( 1, Loop_Line_Y, Loop_Line_X ) &
                       +Vector_tmp( 2, Loop_Line_Y, Loop_Line_X )*Vector_tmp( 2, Loop_Line_Y, Loop_Line_X ) )
        do i= 1, 2 
           Vector_Unit_tmp( i, Loop_Line_Y, Loop_Line_X )= Vector_tmp( i, Loop_Line_Y, Loop_Line_X )/Real_tmp
        end do

        !==================================================================================================
        !write(*,*)'Update Position_Line'
        !==================================================================================================

        Vector_Length= ( Vector_Unit_tmp( 1, Loop_Line_Y, Loop_Line_X )*Line_Length_X &
                         /abs( Vector_Unit_tmp( 1, Loop_Line_Y, Loop_Line_X ) ) )**2 &
                      +( Vector_Unit_tmp( 2, Loop_Line_Y, Loop_Line_X )*Line_Length_X &
                         /abs( Vector_Unit_tmp( 1, Loop_Line_Y, Loop_Line_X ) ) )**2 

        !if( Line_Length_X*sqrt(2.0d0) <= sqrt( Vector_Length ) )then
        !if( Line_Length_X*1.3d0 <= sqrt( Vector_Length ) )then
        if( Line_Length_X*2.0d0/sqrt(3.0d0) <= sqrt( Vector_Length ) )then

           Number_Line( Loop_Line_Y )= Number_Line( Loop_Line_Y ) +1
           Counter_Line( Loop_Line_Y )= Counter_Line( Loop_Line_Y ) +1

           !Position_Line( 1, Loop_Line_Y, Counter_Line( Loop_Line_Y )+1 )&
           != Position_Line( 1, Loop_Line_Y, Counter_Line( Loop_Line_Y ) ) & 
           !  +Vector_Unit_tmp( 1, Loop_Line_Y, Loop_Line_X )*Line_Length_X 

           !Position_Line( 2, Loop_Line_Y, Counter_Line( Loop_Line_Y )+1 )&
           != Position_Line( 2, Loop_Line_Y, Counter_Line( Loop_Line_Y ) ) &
           !  +Vector_Unit_tmp( 2, Loop_Line_Y, Loop_Line_X )*Line_Length_X 

           Position_Line( 1, Loop_Line_Y, Counter_Line( Loop_Line_Y ) )&
           = Position_Line( 1, Loop_Line_Y, Counter_Line( Loop_Line_Y )-1 ) & 
             +Vector_Unit_tmp( 1, Loop_Line_Y, Loop_Line_X )*Line_Length_X 

           Position_Line( 2, Loop_Line_Y, Counter_Line( Loop_Line_Y ) )&
           = Position_Line( 2, Loop_Line_Y, Counter_Line( Loop_Line_Y )-1 ) &
             +Vector_Unit_tmp( 2, Loop_Line_Y, Loop_Line_X )*Line_Length_X

           goto 590 

        else if( Loop_Line_X <= Number_Line_X )then
        !if( Loop_Line_X <= Number_Line_X )then

           Counter_Line( Loop_Line_Y )= Counter_Line( Loop_Line_Y ) +1

           Position_Line( 1, Loop_Line_Y, Counter_Line( Loop_Line_Y ) )&
           = Position_Line( 1, Loop_Line_Y, Counter_Line( Loop_Line_Y )-1 ) & 
             +Vector_Unit_tmp( 1, Loop_Line_Y, Loop_Line_X )*Line_Length_X &
             /abs( Vector_Unit_tmp( 1, Loop_Line_Y, Loop_Line_X ) )  

           Position_Line( 2, Loop_Line_Y, Counter_Line( Loop_Line_Y ) )&
           = Position_Line( 2, Loop_Line_Y, Counter_Line( Loop_Line_Y )-1 ) &
             +Vector_Unit_tmp( 2, Loop_Line_Y, Loop_Line_X )*Line_Length_X &
             /abs( Vector_Unit_tmp( 1, Loop_Line_Y, Loop_Line_X ) ) 

           !Position_Line( 1, Loop_Line_Y, Counter_Line( Loop_Line_Y )+1 )&
           != Position_Line( 1, Loop_Line_Y, Counter_Line( Loop_Line_Y ) ) & 
           !  +Vector_Unit_tmp( 1, Loop_Line_Y, Loop_Line_X )*Line_Length_X &
           !  /abs( Vector_Unit_tmp( 1, Loop_Line_Y, Loop_Line_X ) )  

           !Position_Line( 2, Loop_Line_Y, Counter_Line( Loop_Line_Y )+1 )&
           != Position_Line( 2, Loop_Line_Y, Counter_Line( Loop_Line_Y ) ) &
           !  +Vector_Unit_tmp( 2, Loop_Line_Y, Loop_Line_X )*Line_Length_X &
           !  /abs( Vector_Unit_tmp( 1, Loop_Line_Y, Loop_Line_X ) ) 

           !Position_Line( 1, Loop_Line_Y, Loop_Line_X+1 )&
           != Position_Line( 1, Loop_Line_Y, Loop_Line_X ) & 
           !  +Vector_Unit_tmp( 1, Loop_Line_Y, Loop_Line_X )*Line_Length_X &
           !  /abs( Vector_Unit_tmp( 1, Loop_Line_Y, Loop_Line_X ) ) 

           !Position_Line( 2, Loop_Line_Y, Loop_Line_X+1 )&
           != Position_Line( 2, Loop_Line_Y, Loop_Line_X ) &
           !  +Vector_Unit_tmp( 2, Loop_Line_Y, Loop_Line_X )*Line_Length_X &
           !  /abs( Vector_Unit_tmp( 1, Loop_Line_Y, Loop_Line_X ) ) 
        end if
     end do
  end do

  deallocate( Vector_tmp )
  deallocate( Vector_Unit_tmp )
  deallocate( Element_Number_include_Center )


  !do Loop_Line_Y= 1, Number_Line_Y
  do Loop_Line_Y= 1, Number_Line_Y
     write( FN, * )'1 1 1 setrgbcolor'
     !write(FN,*)'0 0 0 setrgbcolor'
     write(FN,*) '1 setlinewidth'
     write(FN,*)'newpath'
     write(FN,*) Position_Line( 1, Loop_Line_Y, 1 ), Position_Line( 2, Loop_Line_Y, 1 ), 'moveto'
     do Loop_Line_X= 1, Number_Line( Loop_Line_Y )!Number_Line_X+1
        write(FN,*) Position_Line( 1, Loop_Line_Y, Loop_Line_X ), Position_Line( 2, Loop_Line_Y, Loop_Line_X ), 'lineto'
     end do
     write(FN,*)'stroke'
  end do

  deallocate( Position_Line )
  deallocate( Number_Line )
  deallocate( Counter_Line )

  2272 continue

  write(*,*) ' ' 
  write(*,*) '     end subroutine Plot_Stream_Line' 
  return
end subroutine Plot_Stream_Line




subroutine Compute_Position_and_Element_PV &
           ( Number_Point_PV_X, Number_Point_PV_Y, &
             Position_Node, Number_Node,  &
             Index_Element_2_Node, Number_Element, &
             !====================================================================================================
             Position_PV, Element_Number_PV, Number_PV )

   !$ use omp_lib
   implicit none

   integer, intent(in) :: Number_Point_PV_X, Number_Point_PV_Y  
   integer, intent(in) :: Number_Node, Number_Element
   double precision, intent(in) :: Position_Node( 2, Number_Node )
   integer, intent(in) :: Index_Element_2_Node( 4, Number_Element ) 

   double precision, intent(out) :: Position_PV( 2, Number_Point_PV_X*Number_Point_PV_Y )
   integer, intent(out) :: Element_Number_PV( Number_Point_PV_X*Number_Point_PV_Y )
   integer, intent(out) :: Number_PV 

   !double precision, allocatable, dimension(:,:,:) :: Position_PV_tmp( 2, Number_Point_PV_X, Number_Point_PV_Y )
   double precision, allocatable, dimension(:,:,:) :: Position_PV_tmp
   integer :: e, i, j, k, l, m, n, Loop_X, Loop_Y
   integer( kind=8 ) :: Integer_tmp
   !double precision :: Maximum_Position( 2 ), Minimum_Position( 2 )
   double precision, allocatable, dimension(:) :: Maximum_Position, Minimum_Position
   !double precision :: Position_Node_in_Element( 2, 3, Number_Element )
   double precision,  allocatable, dimension(:,:,:) :: Position_Node_in_Element
   !integer :: Number_Devision( 2 )
   !integer :: Flag_Direction_Vector( 3 )
   integer, allocatable, dimension(:) :: Number_Devision, Flag_Direction_Vector
   double precision :: Ratio_Number_Devision
   !double precision :: Vector_1( 2, 3 ), Vector_2( 2, 3 )
   double precision, allocatable, dimension(:,:)  :: Vector_1, Vector_2
   !double precision :: Position_Center_G_Element( 2, Number_Element ), Width( 2 )
   double precision, allocatable, dimension(:,:) :: Position_Center_G_Element
   double precision, allocatable, dimension(:) :: Width
   double precision, allocatable, dimension(:,:,:) :: Position_Division_Low, Position_Division_High
   integer, allocatable, dimension(:,:,:) :: Element_Number_in_Department, Department_Number 
   integer, allocatable, dimension(:,:,:) :: Element_Number_PV_tmp_1 
   integer, allocatable, dimension(:,:) :: Element_Number_PV_tmp_2
   integer, allocatable, dimension(:,:) :: Counter_Department 
   integer :: Size_Department, Maximum_Department
   !integer :: Department_Number_Previous( 2 )
   integer, allocatable, dimension(:) :: Department_Number_Previous
   double precision :: Minimum_Distance_G_2_PV( Number_Point_PV_X, Number_Point_PV_Y ), Distance_G_2_PV 

   integer :: Counter_PV 
   integer, parameter :: Flag_PV=0 

   !====================================================================================================
   write(*,*)'            call Compute_Position_and_Element_PV' 
   !====================================================================================================

   allocate( Position_Node_in_Element( 2, 3, Number_Element ) )

   do e= 1, Number_Element
      do i= 1, 3 
         do j= 1, 2 
            Position_Node_in_Element( j, i, e )= Position_Node( j, Index_Element_2_Node( i, e ) )
         end do
      end do
   end do

   allocate( Position_Center_G_Element( 2, Number_Element ) )

   do e= 1, Number_Element
      do i= 1, 2 
         Position_Center_G_Element( i, e )&
         = ( Position_Node_in_Element( i, 1, e ) +Position_Node_in_Element( i, 2, e ) +Position_Node_in_Element( i, 3, e ) )/3d0 
      end do
   end do

   Ratio_Number_Devision= Number_Point_PV_X/Number_Point_PV_Y

   write(*,*) Number_Point_PV_X, Number_Point_PV_Y, Number_Element
   Integer_tmp= Number_Point_PV_X*Number_Point_PV_Y*Number_Element
   write(*,*) Integer_tmp, Number_Point_PV_X*Number_Point_PV_Y +Number_Element

   allocate( Number_Devision( 2 ) )

   Number_Devision( 2 )= int( sqrt( 1d0/Ratio_Number_Devision &
                         *sqrt( dble( Integer_tmp )&
                         /dble( Number_Point_PV_X*Number_Point_PV_Y +Number_Element ) ) ) )

   Number_Devision( 1 )= Number_Devision( 2 ) *Ratio_Number_Devision

   write(*,*)'          =================================================='
   write(*,*)'            Number_Devision( 1 )=', Number_Devision( 1 )
   write(*,*)'            Number_Devision( 2 )=', Number_Devision( 2 )
   write(*,*)'          =================================================='

   !====================================================================================================
   write(*,*)'            Compute Maximum and Minimum Value of Position'
   !====================================================================================================

   allocate( Maximum_Position( 2 ) )
   allocate( Minimum_Position( 2 ) )

   Minimum_Position( 1 )= 1d10 
   Minimum_Position( 2 )= 1d10 

   Maximum_Position( 1 )= -1d10 
   Maximum_Position( 2 )= -1d10 

   do e= 1, Number_Element
      do i= 1, 3
         do j= 1, 2
            if( Maximum_Position( j ) < Position_Node_in_Element( j, i, e ) )then
               Maximum_Position( j )= Position_Node_in_Element( j, i, e ) 
            end if
            if( Minimum_Position( j ) > Position_Node_in_Element( j, i, e ) )then
               Minimum_Position( j )= Position_Node_in_Element( j, i, e ) 
            end if
         end do
      end do
   end do

   write(*,*)'Minimum_Position=', Minimum_Position( 1 ), Minimum_Position( 2 )
   write(*,*)'Maximum_Position=', Maximum_Position( 1 ), Maximum_Position( 2 )


   !====================================================================================================
   write(*,*)'            Compute Maximum and Minimum Value of Position'
   !====================================================================================================

   allocate( Position_Division_Low( 2, Number_Devision( 1 ), Number_Devision( 2 ) ) )
   allocate( Position_Division_High( 2, Number_Devision( 1 ), Number_Devision( 2 ) ) )
   allocate( Width( 2 ) )

   do i= 1, 2
      Width( i )=( Maximum_Position( i )-Minimum_Position( i ) )/Number_Devision( i )
   end do

   do i= 1, Number_Devision( 2 )
      do j= 1, Number_Devision( 1 ) 
         Position_Division_Low( 1, j, i )= Minimum_Position( 1 ) +Width( 1 )*( j-1 ) !-Width( 1 )/100d0
         Position_Division_Low( 2, j, i )= Minimum_Position( 2 ) +Width( 2 )*( i-1 ) !-Width( 2 )/100d0
         Position_Division_High( 1, j, i )= Minimum_Position( 1 ) +Width( 1 )*( j ) !+Width( 1 )/100d0
         Position_Division_High( 2, j, i )= Minimum_Position( 2 ) +Width( 2 )*( i ) !+Width( 2 )/100d0
      end do
   end do

   deallocate( Width )

   !====================================================================================================
   write(*,*)'            Check Position_Division'
   !====================================================================================================

!   ! OK
!   write(100,*) '#, Number_Devision=', Number_Devision( 1 ), Number_Devision( 2 ) 
!   write(101,*) '#, Number_Devision=', Number_Devision( 1 ), Number_Devision( 2 ) 
!   do i= 1, Number_Devision( 2 ) 
!      do j= 1, Number_Devision( 1 )
!         write(100,*) Position_Division_Low( 1, j, i ), Position_Division_Low( 2, j, i )
!         write(101,*) Position_Division_High( 1, j, i ), Position_Division_High( 2, j, i )
!      end do
!   end do

   !====================================================================================================
   write(*,*)'            Investigate Department of Element'
   !====================================================================================================

   allocate( Counter_Department( Number_Devision( 1 ), Number_Devision( 2 ) ) )

   do i= 1, Number_Devision( 2 ) 
      do j= 1, Number_Devision( 1 ) 
         Counter_Department( j, i )= 0
      end do
   end do

   Size_Department= Number_Element/( Number_Devision( 1 )*Number_Devision( 2 ) )*3d0

   allocate( Element_Number_in_Department( Size_Department, Number_Devision( 1 ), Number_Devision( 2 ) ) )

   Element_Number_in_Department= 0

   do e= 1, Number_Element
      do i= 1, Number_Devision( 2 ) 
         do j= 1, Number_Devision( 1 ) 
            if( Position_Division_Low( 2, j, i ) <= Position_Center_G_Element( 2, e ) .and. &
                Position_Center_G_Element( 2, e ) <= Position_Division_High( 2, j, i ) .and. &
                Position_Division_Low( 1, j, i ) <= Position_Center_G_Element( 1, e ) .and. &
                Position_Center_G_Element( 1, e ) <= Position_Division_High( 1, j, i )  )then
                
                Counter_Department( j, i )= Counter_Department( j, i ) +1
                Element_Number_in_Department( Counter_Department( j, i ), j, i )= e
            end if
         end do
      end do
   end do

   Maximum_Department= -1d8
   do i= 1, Number_Devision( 2 ) 
      do j= 1, Number_Devision( 1 ) 
         if( Maximum_Department < Counter_Department( j, i ) )then
            Maximum_Department= Counter_Department( j, i )
         end if
      end do
   end do

   write(*,*)'          =================================================='
   write(*,*)'            Maximum_Department=', Maximum_Department
   write(*,*)'          =================================================='

   !====================================================================================================
   write(*,*)'            Check Element_Number_in_Department'
   !====================================================================================================

   ! OK
   do i= 1, Number_Devision( 2 ) 
      do j= 1, Number_Devision( 1 )
         do k= 1, Counter_Department( j, i )
            if( Element_Number_in_Department( k, j, i )==0 )then
               call Output_Error( 'Compute_Position_and_Element_PV', 156 ) 
            end if

   !         l= 200 +10* i +j
   !         write(l,*) Position_Node_in_Element( 1, 1, Element_Number_in_Department( k, j, i ) ), & 
   !                Position_Node_in_Element( 2, 1, Element_Number_in_Department( k, j, i ) )  
         end do
      end do
   end do

   !====================================================================================================
   write(*,*)'            Determine the Position of PV Computation'
   !====================================================================================================

   allocate( Position_PV_tmp( 2, Number_Point_PV_X, Number_Point_PV_Y ) )

   do i= 1, Number_Point_PV_Y
      do j= 1, Number_Point_PV_X
         Position_PV_tmp( 1, j, i )= Minimum_Position( 1 ) +( Maximum_Position( 1 ) -Minimum_Position( 1 ) )/Number_Point_PV_X*( j-0.5d0 )
         Position_PV_tmp( 2, j, i )= Minimum_Position( 2 ) +( Maximum_Position( 2 ) -Minimum_Position( 2 ) )/Number_Point_PV_Y*( i-0.5d0 )
      end do
   end do

   deallocate( Maximum_Position )
   deallocate( Minimum_Position )

   !====================================================================================================
   write(*,*)'            Check Position_PV_tmp'
   !====================================================================================================

   !! OK
   !do i= 1, Number_Point_PV_Y
   !   do j= 1, Number_Point_PV_X
   !      write(300,*) Position_PV_tmp( 1, j, i ), Position_PV_tmp( 2, j, i ) 
   !   end do
   !end do

   !====================================================================================================
   write(*,*)'            Investigate Department of Position'
   !====================================================================================================

   allocate( Department_Number( 2, Number_Point_PV_X, Number_Point_PV_Y ) )

   Department_Number= 0

   do Loop_Y= 1, Number_Point_PV_Y
      do Loop_X= 1, Number_Point_PV_X
         do i= 1, Number_Devision( 2 ) 
            do j= 1, Number_Devision( 1 )
               if( Position_Division_Low( 2, j, i ) <= Position_PV_tmp( 2, Loop_X, Loop_Y ) .and. &
                   Position_PV_tmp( 2, Loop_X, Loop_Y ) <= Position_Division_High( 2, j, i ) .and. &
                   Position_Division_Low( 1, j, i ) <= Position_PV_tmp( 1, Loop_X, Loop_Y ) .and. &
                   Position_PV_tmp( 1, Loop_X, Loop_Y ) <= Position_Division_High( 1, j, i ) )then
                   
                   Department_Number( 1, Loop_X, Loop_Y )= j
                   Department_Number( 2, Loop_X, Loop_Y )= i
                   go to 212
               end if
            end do
         end do
         212 continue
      end do
   end do

   deallocate( Number_Devision )

   !====================================================================================================
   write(*,*)'            Check Department_Number '
   !====================================================================================================

   !! OK
   !do Loop_Y= 1, Number_Point_PV_Y
   !   do Loop_X= 1, Number_Point_PV_X
   !      write(400,*) Position_PV_tmp( 1, Loop_X, Loop_Y ), Position_PV_tmp( 2, Loop_X, Loop_Y ), & 
   !               Department_Number( 1, Loop_X, Loop_Y ), Department_Number( 2, Loop_X, Loop_Y )

   !      do i= 1, 2
   !          if( Department_Number( i, Loop_X, Loop_Y )==0 )then
   !             call Output_Error( 'Compute_Position_and_Element_PV', 190 )
   !          end if
   !      end do
   !   end do
   !end do

   deallocate( Position_Division_Low )
   deallocate( Position_Division_High )


   !====================================================================================================
   write(*,*)'            Check Department_Number '
   !====================================================================================================

   allocate( Element_Number_PV_tmp_1( Maximum_Department, Number_Point_PV_X, Number_Point_PV_Y ) )

   do Loop_Y= 1, Number_Point_PV_Y
      do Loop_X= 1, Number_Point_PV_X
         do k= 1, Counter_Department( Department_Number( 1, Loop_X, Loop_Y ), Department_Number( 2, Loop_X, Loop_Y ) ) 
            Element_Number_PV_tmp_1( k, Loop_X, Loop_Y )&
            = Element_Number_in_Department( k, Department_Number( 1, Loop_X, Loop_Y ), Department_Number( 2, Loop_X, Loop_Y ) ) 
         end do
      end do
   end do

   allocate( Department_Number_Previous( 2 ) )

   Department_Number_Previous( 1 )= 0
   Department_Number_Previous( 2 )= 0

   do Loop_Y= 1, Number_Point_PV_Y
      do Loop_X= 1, Number_Point_PV_X
!         if( Department_Number_Previous( 1 ) /= Department_Number( 1, Loop_X, Loop_Y ) .or. &
!             Department_Number_Previous( 2 ) /= Department_Number( 2, Loop_X, Loop_Y )  )then
!
!            do k= 1, Counter_Department( Department_Number( 1, Loop_X, Loop_Y ), Department_Number( 2, Loop_X, Loop_Y ) ) 
!               m= 500 +10*Department_Number( 2, Loop_X, Loop_Y ) +Department_Number( 1, Loop_X, Loop_Y ) 
!               write(m,*) Position_Node_in_Element( 1, 1, Element_Number_PV_tmp_1( k, Loop_X, Loop_Y ) ), & 
!                      Position_Node_in_Element( 2, 1, Element_Number_PV_tmp_1( k, Loop_X, Loop_Y ) )  
!            end do
!         end if
         Department_Number_Previous( 1 )= Department_Number( 1, Loop_X, Loop_Y )
         Department_Number_Previous( 2 )= Department_Number( 2, Loop_X, Loop_Y )
      end do
   end do

   deallocate( Department_Number_Previous )

   !====================================================================================================
   write(*,*)'            Detect Element including PV Points'
   !====================================================================================================

   Element_Number_PV= 0
   Minimum_Distance_G_2_PV= 1d8

   allocate( Element_Number_PV_tmp_2( Number_Point_PV_X, Number_Point_PV_Y ) )
   allocate( Flag_Direction_Vector( 3 ) )
   allocate( Vector_1( 2, 3 ) )
   allocate( Vector_2( 2, 3 ) )

   Element_Number_PV_tmp_2= 0

   do Loop_Y= 1, Number_Point_PV_Y
      do Loop_X= 1, Number_Point_PV_X
         do k= 1, Counter_Department( Department_Number( 1, Loop_X, Loop_Y ), Department_Number( 2, Loop_X, Loop_Y ) ) 

            ! OUT
            !m= 700 +10* Loop_Y +Loop_X
            !write(m,*) Position_Node_in_Element( 1, 1, Element_Number_PV_tmp_1( k, Loop_X, Loop_Y ) ), & 
            !       Position_Node_in_Element( 2, 1, Element_Number_PV_tmp_1( k, Loop_X, Loop_Y ) )  

            if( Flag_PV==0 )then

               Distance_G_2_PV &
               = ( Position_PV_tmp( 1, Loop_X, Loop_Y )-Position_Center_G_Element( 1, Element_Number_PV_tmp_1( k, Loop_X, Loop_Y ) ) )**2 &
                +( Position_PV_tmp( 2, Loop_X, Loop_Y )-Position_Center_G_Element( 2, Element_Number_PV_tmp_1( k, Loop_X, Loop_Y ) ) )**2

               if( Minimum_Distance_G_2_PV( Loop_X, Loop_Y ) > Distance_G_2_PV )then
                  Minimum_Distance_G_2_PV( Loop_X, Loop_Y ) = Distance_G_2_PV
                  Element_Number_PV_tmp_2( Loop_X, Loop_Y )= Element_Number_PV_tmp_1( k, Loop_X, Loop_Y ) 
               end if

            else if( Flag_PV==1 )then
               do l= 1, 2 
                  Vector_1( l, 1 )= Position_Node_in_Element( l, 2, Element_Number_PV_tmp_1( k, Loop_X, Loop_Y ) ) &
                             -Position_Node_in_Element( l, 1, Element_Number_PV_tmp_1( k, Loop_X, Loop_Y ) ) 
                  Vector_1( l, 2 )= Position_Node_in_Element( l, 3, Element_Number_PV_tmp_1( k, Loop_X, Loop_Y ) ) &
                             -Position_Node_in_Element( l, 2, Element_Number_PV_tmp_1( k, Loop_X, Loop_Y ) ) 
                  Vector_1( l, 3 )= Position_Node_in_Element( l, 1, Element_Number_PV_tmp_1( k, Loop_X, Loop_Y ) ) &
                             -Position_Node_in_Element( l, 3, Element_Number_PV_tmp_1( k, Loop_X, Loop_Y ) ) 
                  Vector_2( l, 1 )= Position_PV_tmp( l, Loop_X, Loop_Y ) &
                             -Position_Node_in_Element( l, 1, Element_Number_PV_tmp_1( k, Loop_X, Loop_Y ) ) 
                  Vector_2( l, 2 )= Position_PV_tmp( l, Loop_X, Loop_Y ) &
                             -Position_Node_in_Element( l, 2, Element_Number_PV_tmp_1( k, Loop_X, Loop_Y ) ) 
                  Vector_2( l, 3 )= Position_PV_tmp( l, Loop_X, Loop_Y ) &
                             -Position_Node_in_Element( l, 3, Element_Number_PV_tmp_1( k, Loop_X, Loop_Y ) ) 
               end do
 
               call Judge_Cross_Product &
                  ( 3, Vector_1, Vector_2, 0, Flag_Direction_Vector )

               if( Flag_Direction_Vector( 2 )==1 .and. Flag_Direction_Vector( 3 )==1  )then
                  Element_Number_PV_tmp_2( Loop_X, Loop_Y )= Element_Number_PV_tmp_1( k, Loop_X, Loop_Y ) 
                  !go to 114
               end if
            end if

         end do
         114 continue
      end do
   end do

   deallocate( Vector_1 )
   deallocate( Vector_2 )
   deallocate( Flag_Direction_Vector )
   deallocate( Position_Node_in_Element )
   deallocate( Element_Number_in_Department )
   deallocate( Element_Number_PV_tmp_1 )
   deallocate( Department_Number )
   deallocate( Counter_Department )

   Counter_PV= 0

   if( Flag_PV==0 )then
      do i= 1, Number_Point_PV_Y
         do j= 1, Number_Point_PV_X
            if( Minimum_Distance_G_2_PV( j, i ) <= 1d0/200d0 .and. Element_Number_PV_tmp_2( j, i )/=0 )then
               Counter_PV= Counter_PV +1
               Element_Number_PV( Counter_PV )= Element_Number_PV_tmp_2( j, i )

               do k= 1, 2
                  Position_PV( k, Counter_PV )= Position_Center_G_Element( k, Element_Number_PV_tmp_2( j, i ) ) 
               end do
            end if
         end do
      end do
   else if( Flag_PV==1 )then
      do i= 1, Number_Point_PV_Y
         do j= 1, Number_Point_PV_X
            if( Element_Number_PV_tmp_2( j, i )==0 )then
               call Output_Error( 'Compute_Position_and_Element_PV', 238 ) 
            else
               Counter_PV= Counter_PV +1
               Element_Number_PV( Counter_PV )= Element_Number_PV_tmp_2( j, i )

               do k= 1, 2
                  Position_PV( k, Counter_PV )=Position_PV_tmp( k, j, i ) 
               end do
            end if
         end do
      end do
   end if

   Number_PV= Counter_PV

   deallocate( Position_Center_G_Element )
   deallocate( Position_PV_tmp )
   deallocate( Element_Number_PV_tmp_2 )

   !====================================================================================================
   write(*,*)'            end Compute_Position_and_Element_PV'
   !====================================================================================================

   return
end subroutine Compute_Position_and_Element_PV




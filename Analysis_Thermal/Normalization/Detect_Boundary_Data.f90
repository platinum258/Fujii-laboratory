
subroutine Detect_Boundary_Data &
       ( ID_Element_Boundary, &
         Number_Node, Position_Node, &
         Number_Element, Shape_Element, Index_Element_2_Node, Class_Element, &
         Flag_Renumber_Local_Node, &
       !====================================================================================================
         Number_Edge_on_Boundary, Element_Number_Boundary, Local_Node_Number )

   !$ use omp_lib
   implicit none

   integer, intent(in) :: ID_Element_Boundary 
   integer, intent(in) :: Number_Node, Number_Element, Shape_Element
   double precision, intent(in) :: Position_Node( 2, Number_Node )
   integer, intent(in) :: Index_Element_2_Node( Shape_Element, Number_Element ), Class_Element( Number_Element )
   integer, intent(in) :: Flag_Renumber_Local_Node

   integer, intent(out) :: Number_Edge_on_Boundary
   integer, intent(out) :: Element_Number_Boundary( Number_Element )
   integer, intent(out) :: Local_Node_Number( 2, Number_Element )

   integer :: e, i, j, k, l
   integer, allocatable, dimension(:) :: Counter_Node_in_Element_1
   integer, allocatable, dimension(:) :: Flag_Node_In, Flag_Node_Out, Flag_Node
   integer, allocatable, dimension(:) :: Flag_Node_In_Plus, Flag_Node_Out_Plus, Flag_Node_Plus
   integer, allocatable, dimension(:,:) :: Local_Node_Number_tmp 
   integer, allocatable, dimension(:,:) :: Index_Node_2_Element
   integer, allocatable, dimension(:) :: Counter_Element_having_Node 
   integer, allocatable, dimension(:) :: Element_Number_Boundary_Common 

   integer :: Counter_Edge_1, Counter_Element_2, Counter_Node_Flag, Counter_Percent
 
   double precision, allocatable, dimension(:,:,:) :: Position_LocalNode_in_Elements 
   double precision, allocatable, dimension(:) :: Cross_Product_1 
   integer, allocatable, dimension(:) :: Integer_tmp
         
   integer, parameter :: Flag_Off= 0
   integer, parameter :: Flag_On= 1

   ! Flag_Renumber_Local_Node
   integer, parameter :: Flag_Position_X= 1
   integer, parameter :: Flag_Position_Y= 2
   integer, parameter :: Flag_Position_Clock_Rotation= 3

   double precision, parameter :: Position_Origin_X= 0d0
   double precision, parameter :: Position_Origin_Y= 0d0

   integer :: Number_Element_Common, Max_Number_Element_Common


   double precision, allocatable, dimension(:) :: Area_Element 
   double precision, allocatable, dimension(:,:,:) :: Position_Node_Element_Triangle 
   double precision, allocatable, dimension(:,:,:,:) :: Difference_Position_Node_Element_Triangle

   double precision, parameter :: Limit_Area = 0d0 !( 1d0/200d0 )**2 *2.4d-4 !1d-4 < 2d-4 < 2.25d-4 < 2.3d-4 < x =2.4d-4 < 2.5d-4 < 3d-4 < 5d-4 < 1d-3
   !double precision, parameter :: Limit_Length_Edge = 1d0/200d0*1d-5 !   1d-5 < 5d-5 < 7.5d-5 < x < 1d-4
   double precision, parameter :: Limit_Length_Edge = 1d0/200d0*1d-4 !   1d-5 < 5d-5 < 7.5d-5 < x < 1d-4
   double precision, allocatable, dimension(:,:,:) :: Length_Edge_Element 
 
   !====================================================================================================
   write(*,*)'            call Detect_Boundary_Data', ID_Element_Boundary 
   !====================================================================================================

   allocate( Flag_Node_In( Number_Node ) )
   allocate( Flag_Node_Out( Number_Node ) )
   allocate( Flag_Node_In_Plus( Number_Node ) )
   allocate( Flag_Node_Out_Plus( Number_Node ) )
   allocate( Counter_Element_having_Node( Number_Node ) )

   do i= 1, Number_Node
      Flag_Node_In( i )= 0 
      Flag_Node_Out( i )= 0 
      Flag_Node_In_Plus( i )= 0 
      Flag_Node_Out_Plus( i )= 0 
      Counter_Element_having_Node( i )= 0
   end do

   Number_Element_Common= 15
   Max_Number_Element_Common= 0

   allocate( Index_Node_2_Element( Number_Element_Common, Number_Node ) )

   do i= 1, Number_Node
      do j= 1, Number_Element_Common 
         Index_Node_2_Element( j, i )= 0
      end do
   end do

   do e= 1, Number_Element
      do i= 1, Shape_Element
         Counter_Element_having_Node( Index_Element_2_Node( i, e ) )&
         = Counter_Element_having_Node( Index_Element_2_Node( i, e ) ) +1

         Index_Node_2_Element( Counter_Element_having_Node( Index_Element_2_Node( i, e ) ), Index_Element_2_Node( i, e ) ) &
         = e

         if( Max_Number_Element_Common < Counter_Element_having_Node( Index_Element_2_Node( i, e ) ) )then
            Max_Number_Element_Common= Counter_Element_having_Node( Index_Element_2_Node( i, e ) )
         end if
      end do
   end do

   !====================================================================================================
   write(*,*)'               Detect Elements on Boundary'
   !====================================================================================================

   do e= 1, Number_Element
      if( Class_Element( e )==ID_Element_Boundary )then
         do i= 1, Shape_Element
            Flag_Node_In( Index_Element_2_Node( i, e ) )= 1
         end do
         do i= 1, Shape_Element
            Flag_Node_In_Plus( Index_Element_2_Node( i, e ) )= Flag_Node_In_Plus( Index_Element_2_Node( i, e ) ) +1
         end do
      else if( Class_Element( e )/=ID_Element_Boundary )then
         do i= 1, Shape_Element
            Flag_Node_Out( Index_Element_2_Node( i, e ) )= 1
         end do
         do i= 1, Shape_Element
            Flag_Node_Out_Plus( Index_Element_2_Node( i, e ) )= Flag_Node_Out_Plus( Index_Element_2_Node( i, e ) ) +10
         end do
      end if
   end do

   allocate( Flag_Node( Number_Node ) )
   allocate( Flag_Node_Plus( Number_Node ) )

   do i= 1, Number_Node
      Flag_Node( i )= Flag_Node_In( i ) +Flag_Node_Out( i )
      Flag_Node_Plus( i )= Flag_Node_In_Plus( i ) +Flag_Node_Out_Plus( i )
   end do

   deallocate( Flag_Node_In )
   deallocate( Flag_Node_Out )

   !====================================================================================================
   write(*,*)'               Check the Number of Detected Nodes on Boundary'
   !====================================================================================================
   Counter_Node_Flag= 0

   do i= 1, Number_Node
      if( Flag_Node( i )==2 )then
         Counter_Node_Flag= Counter_Node_Flag +1
      end if
   end do

   if( Counter_Node_Flag==0 )then
      write(*,*)'Counter_Node_Flag=', Counter_Node_Flag
      !call Output_Error( 'Detect_Boundary_Data', 93 )
   end if

   do i= 1, Number_Node
      if( Flag_Node( i ) < 0 .or. 2 < Flag_Node( i ) )then
         call Output_Error( 'Detect_Boundary_Data', 55 )
      end if
   end do

   !====================================================================================================
   write(*,*)'               Count Nodes Overlapping Boundary in an Element'
   !====================================================================================================

   allocate( Counter_Node_in_Element_1( Number_Element ) )
   allocate( Local_Node_Number_tmp( 3, Number_Element ) )

   do e= 1, Number_Element
      Counter_Node_in_Element_1( e )= 0
      do i= 1, 2
         Local_Node_Number( i, e )= 0
         Local_Node_Number_tmp( i, e )= 0
      end do
   end do

   do e= 1, Number_Element
      if( Class_Element( e )==ID_Element_Boundary )then 
         do i= 1, Shape_Element
            if( Flag_Node( Index_Element_2_Node( i, e ) )==2 )then
               Counter_Node_in_Element_1( e )= Counter_Node_in_Element_1( e ) +1
               Local_Node_Number_tmp( Counter_Node_in_Element_1( e ), e )= i
            end if
         end do
      end if
   end do

   deallocate( Flag_Node )

   !====================================================================================================
   write(*,*)'               Compute Area_Element'
   !====================================================================================================
   allocate( Position_Node_Element_Triangle( 2, 3, Number_Element ) )
    
   !$omp parallel do default( none ) &
   !$omp private( e, i, j ) & 
   !$omp shared( Number_Element, Position_Node, Position_Node_Element_Triangle ) & 
   !$omp shared( Index_Element_2_Node ) 
   do e= 1, Number_Element
      do i= 1, 3
         do j= 1, 2
            Position_Node_Element_Triangle( j, i, e )= Position_Node( j, Index_Element_2_Node( i, e ) )
         end do
      end do
   end do

   allocate( Difference_Position_Node_Element_Triangle( 2, 3, 3, Number_Element ) )

   do e= 1, Number_Element
      do i= 1, 3
         do j= 1, 3
            do k= 1, 2
               Difference_Position_Node_Element_Triangle( k, j, i, e ) &
               =Position_Node_Element_Triangle( k, j, e ) -Position_Node_Element_Triangle( k, i, e )
            end do
         end do
      end do
   end do

   allocate( Area_Element( Number_Element ) )

   !$omp parallel do default( none ) &
   !$omp private( e ) & 
   !$omp shared( Number_Element, Area_Element, Difference_Position_Node_Element_Triangle ) 
   do e= 1, Number_Element
      Area_Element( e ) &
      =abs( Difference_Position_Node_Element_Triangle( 1, 1, 3, e ) & 
         *Difference_Position_Node_Element_Triangle( 2, 2, 3, e ) &
         -Difference_Position_Node_Element_Triangle( 1, 2, 3, e ) & 
         *Difference_Position_Node_Element_Triangle( 2, 1, 3, e ) )/2d0
   end do 
    
   do e= 1, Number_Element
      if( Area_Element( e ) <= Limit_Area .and. Class_Element( e )==ID_Element_Boundary )then
         do i= 1, 3
            write( 990, * ) Position_Node( 1, Index_Element_2_Node( i, e ) ), Position_Node( 2, Index_Element_2_Node( i, e ) )
         end do
      end if
   end do 

   !====================================================================================================
   write(*,*)'               Compute Length of Edge in Element on Boundary'
   !====================================================================================================

   allocate( Length_Edge_Element( 3, 3, Number_Element ) )

   do e= 1, Number_Element
      do i= 1, 3
         do j= 1, 3
            Length_Edge_Element( j, i, e )&
            = sqrt( Difference_Position_Node_Element_Triangle( 1, j, i, e )**2 +Difference_Position_Node_Element_Triangle( 1, j, i, e )**2 )
         end do
      end do
   end do 

   deallocate( Position_Node_Element_Triangle )
   deallocate( Difference_Position_Node_Element_Triangle )

   !====================================================================================================
   write(*,*)'               Detect Element on Boundary'
   !====================================================================================================

   allocate( Element_Number_Boundary_Common( Number_Element ) )

   Counter_Edge_1= 0
   Counter_Element_2= 0

   Counter_Percent= 1
   do e= 1, Number_Element
      if( e== Number_Element/100*Counter_Percent )then
         write(*,'(a,i3,a,$)')'[', Counter_Percent, '%]'
         Counter_Percent= Counter_Percent +1
      end if 

      if( Counter_Node_in_Element_1( e ) >= 2 .and. Class_Element( e )==ID_Element_Boundary .and. Area_Element( e ) > Limit_Area )then

         do i= 1, Shape_Element-1 
            do j= 1, Max_Number_Element_Common 
               if( Index_Node_2_Element( j, Index_Element_2_Node( i, e ) ) /= 0 )then 
                  if( Class_Element( Index_Node_2_Element( j, Index_Element_2_Node( i, e ) ) ) /= ID_Element_Boundary )then

                     do k= i+1, Shape_Element 
                        if( Length_Edge_Element( k, i, e ) >= Limit_Length_Edge )then
                           do l= 1, Max_Number_Element_Common
                              if( Index_Node_2_Element( l, Index_Element_2_Node( k, e ) ) /= 0 )then 
                                 if( Class_Element( Index_Node_2_Element( l, Index_Element_2_Node( k, e ) ) ) &
                                     /= ID_Element_Boundary )then
                                    if(   Index_Node_2_Element( j, Index_Element_2_Node( i, e ) )&
                                        ==Index_Node_2_Element( l, Index_Element_2_Node( k, e ) ) )then

                                       Counter_Edge_1= Counter_Edge_1 +1
                                       Element_Number_Boundary( Counter_Edge_1 )= e
                                       Element_Number_Boundary_Common( Counter_Edge_1 )&
                                       = Index_Node_2_Element( j, Index_Element_2_Node( i, e ) ) 
                                       Local_Node_Number( 1, Counter_Edge_1 )= i 
                                       Local_Node_Number( 2, Counter_Edge_1 )= k
                                    end if
                                 end if
                              end if
                           end do
                        end if 
                     end do

                  end if 
               end if 
            end do
         end do

      end if 
   end do

   deallocate( Length_Edge_Element )
   write(*,*)' '
   !====================================================================================================
   write(*,*)'               Check Boundary Data'
   !====================================================================================================

   do e= 1, Counter_Edge_1 
      if( Element_Number_Boundary( e ) == Element_Number_Boundary_Common( e ) )then
         write(*,*)'Element_Number_Boundary( e )=', Element_Number_Boundary( e )
         write(*,*)'Element_Number_Boundary_Common( e )=', Element_Number_Boundary_Common( e )
         call Output_Error( 'Detect_Boundary_Data', 252 )
      end if
      if( Local_Node_Number( 1, e )==0 .or. Local_Node_Number( 2, e )==0 )then
         write(*,*)' '
         write(*,*)'Local_Node_Number( 1, e )=', Local_Node_Number( 1, e )
         write(*,*)'Local_Node_Number( 2, e )=', Local_Node_Number( 2, e )
         call Output_Error( 'Detect_Boundary_Data', 258 )
      end if
      if( Class_Element( Element_Number_Boundary( e ) ) /= ID_Element_Boundary )then
         write(*,*)'Class_Element( Element_Number_Boundary( e ) )=', Class_Element( Element_Number_Boundary( e ) )
         write(*,*)'Element_Number_Boundary( e )=', Element_Number_Boundary( e )
         call Output_Error( 'Detect_Boundary_Data', 263 )
      end if
      if( Class_Element( Element_Number_Boundary_Common( e ) ) == ID_Element_Boundary )then
         write(*,*)'Class_Element( Element_Number_Boundary_Common( e ) )=', Class_Element( Element_Number_Boundary_Common( e ) )
         write(*,*)'Element_Number_Boundary_Common( e )=', Element_Number_Boundary_Common( e )
         call Output_Error( 'Detect_Boundary_Data', 268 )
      end if
   end do

   deallocate( Counter_Node_in_Element_1 )
   deallocate( Local_Node_Number_tmp )

   deallocate( Flag_Node_In_Plus )
   deallocate( Flag_Node_Out_Plus )
   deallocate( Flag_Node_Plus )

   deallocate( Element_Number_Boundary_Common )

   Number_Edge_on_Boundary= Counter_Edge_1

   !====================================================================================================
   write(*,*)'               Renumber Local Node Number'
   !====================================================================================================

   allocate( Position_LocalNode_in_Elements( 2, 2, Number_Edge_on_Boundary ) )

   do e= 1, Number_Edge_on_Boundary
      do i= 1, 2
         do j= 1, 2
            Position_LocalNode_in_Elements( j, i, e ) & 
            = Position_Node( j, Index_Element_2_Node( Local_Node_Number( i, e ), Element_Number_Boundary( ( e ) ) ) )
         end do
      end do
   end do

   !====================================================================================================
   if( Flag_Renumber_Local_Node==Flag_Position_X )then
   write(*,*)'               Based on X Position'
   !====================================================================================================

      allocate( Integer_tmp( Number_Edge_on_Boundary ) )

      do e= 1, Number_Edge_on_Boundary
         if( Position_LocalNode_in_Elements( 1, 1, e ) > Position_LocalNode_in_Elements( 1, 2, e ) )then
            Integer_tmp( e )= Local_Node_Number( 1, e ) 
            Local_Node_Number( 1, e )= Local_Node_Number( 2, e )
            Local_Node_Number( 2, e )= Integer_tmp( e )
         else if( Position_LocalNode_in_Elements( 1, 1, e ) == Position_LocalNode_in_Elements( 1, 2, e ) )then
            if( Position_LocalNode_in_Elements( 2, 1, e ) > Position_LocalNode_in_Elements( 2, 2, e ) )then
               Integer_tmp( e )= Local_Node_Number( 1, e ) 
               Local_Node_Number( 1, e )= Local_Node_Number( 2, e )
               Local_Node_Number( 2, e )= Integer_tmp( e )
            end if
         end if
      end do

      deallocate( Integer_tmp )

   !====================================================================================================
   else if( Flag_Renumber_Local_Node==Flag_Position_Y )then
   write(*,*)'               Based on Y Position'
   !====================================================================================================

      allocate( Integer_tmp( Number_Edge_on_Boundary ) )

      do e= 1, Number_Edge_on_Boundary
         if( Position_LocalNode_in_Elements( 2, 1, e ) > Position_LocalNode_in_Elements( 2, 2, e ) )then
            Integer_tmp( e )= Local_Node_Number( 1, e ) 
            Local_Node_Number( 1, e )= Local_Node_Number( 2, e )
            Local_Node_Number( 2, e )= Integer_tmp( e )
         else if( Position_LocalNode_in_Elements( 2, 1, e ) == Position_LocalNode_in_Elements( 2, 2, e ) )then
            if( Position_LocalNode_in_Elements( 1, 1, e ) > Position_LocalNode_in_Elements( 1, 2, e ) )then
               Integer_tmp( e )= Local_Node_Number( 1, e ) 
               Local_Node_Number( 1, e )= Local_Node_Number( 2, e )
               Local_Node_Number( 2, e )= Integer_tmp( e )
            end if
         end if
      end do

      deallocate( Integer_tmp )
   !====================================================================================================
   else if( Flag_Renumber_Local_Node==Flag_Position_Clock_Rotation )then
   write(*,*)'               Based on Clock Rotation'
   !====================================================================================================

      do e= 1, Number_Edge_on_Boundary
         do i= 1, 2
            Position_LocalNode_in_Elements( 1, i, e )=  Position_LocalNode_in_Elements( 1, i, e ) -Position_Origin_X
            Position_LocalNode_in_Elements( 2, i, e )=  Position_LocalNode_in_Elements( 2, i, e ) -Position_Origin_Y
         end do
      end do
 
      !====================================================================================================
      write(*,*)'               Compute Cross Products'
      !====================================================================================================

      allocate( Cross_Product_1( Number_Edge_on_Boundary ) )

      do e= 1, Number_Edge_on_Boundary
         Cross_Product_1( e ) &
         = Position_LocalNode_in_Elements( 1, 1, e ) *Position_LocalNode_in_Elements( 2, 2, e ) &
          -Position_LocalNode_in_Elements( 2, 1, e ) *Position_LocalNode_in_Elements( 1, 2, e ) 
      end do
      allocate( Integer_tmp( Number_Edge_on_Boundary ) )

      do e= 1, Number_Edge_on_Boundary
         if( Cross_Product_1( e ) < 0d0  )then
            Integer_tmp( e )= Local_Node_Number( 1, e ) 
            Local_Node_Number( 1, e )= Local_Node_Number( 2, e )
            Local_Node_Number( 2, e )= Integer_tmp( e )
         !else if( Cross_Product_1( e )==0d0 )then
         !   call Output_Error( 'Detect_Boundary_Data', 339 )
         end if
      end do

      deallocate( Cross_Product_1 )
      deallocate( Integer_tmp )
   !====================================================================================================
   end if
   !====================================================================================================


   deallocate( Counter_Element_having_Node )
   deallocate( Index_Node_2_Element )
   deallocate( Position_LocalNode_in_Elements )

   !====================================================================================================
   write(*,*)'            end Detect_Boundary_Data'
   !====================================================================================================

   return
end subroutine Detect_Boundary_Data




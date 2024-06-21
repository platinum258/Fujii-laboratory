
subroutine Detect_Line_on_Waveguide &
       ( ID_Element, Range_Position_Node, &
         Number_Node, Position_Node, &
         Number_Element, Shape_Element, Index_Element_2_Node, Class_Element, &
         Flag_Renumber_Local_Node, &
       !====================================================================================================
         Number_Element_on_Boundary_1, Element_Number_Boundary_1, Local_Node_Number_1, &
         Number_Element_on_Boundary_2, Element_Number_Boundary_2, Local_Node_Number_2 )

   !$ use omp_lib
   implicit none

   integer, intent(in) :: ID_Element 
   integer, intent(in) :: Number_Node, Number_Element, Shape_Element
   double precision, intent(in) :: Range_Position_Node( 2, 2 ) ! ( x or y, min or max)
   double precision, intent(in) :: Position_Node( 2, Number_Node )
   integer, intent(in) :: Index_Element_2_Node( Shape_Element, Number_Element ), Class_Element( Number_Element )
   integer, intent(in) :: Flag_Renumber_Local_Node

   integer, intent(out) :: Number_Element_on_Boundary_1, Number_Element_on_Boundary_2
   integer, intent(out) :: Element_Number_Boundary_1( Number_Element ), Element_Number_Boundary_2( Number_Element )
   integer, intent(out) :: Local_Node_Number_1( 2, Number_Element ), Local_Node_Number_2( 2, Number_Element )

   integer :: e, i, j
   integer, allocatable, dimension(:) :: Counter_Node_in_Element_1, Counter_Node_in_Element_2
   integer, allocatable, dimension(:) :: Flag_Node
   integer, allocatable, dimension(:,:) :: Local_Node_Number_1_tmp, Local_Node_Number_2_tmp

   integer :: Counter_Element_1, Counter_Element_2, Counter_Node_Flag
 
   double precision, allocatable, dimension(:,:,:) :: Position_LocalNode_in_Elements_1, Position_LocalNode_in_Elements_2
   double precision, allocatable, dimension(:) :: Cross_Product_1, Cross_Product_2 
   integer, allocatable, dimension(:) :: Integer_tmp
         
   integer, parameter :: Flag_Off= 0
   integer, parameter :: Flag_On= 1

   ! Flag_Renumber_Local_Node
   integer, parameter :: Flag_Position_X= 1
   integer, parameter :: Flag_Position_Y= 2
   integer, parameter :: Flag_Position_Clock_Rotation= 3

   double precision, parameter :: Position_Origin_X= 0d0
   double precision, parameter :: Position_Origin_Y= 0d0

   !====================================================================================================
   write(*,*)'   call Detect_Line_on_Waveguide'
   !====================================================================================================

   allocate( Flag_Node( Number_Node ) )

   do i= 1, Number_Node
      Flag_Node( i )= 0 
   end do

   !====================================================================================================
   write(*,*)'      Detect Elements on Boundary'
   !====================================================================================================

   do e= 1, Number_Element
      if( Class_Element( e )==ID_Element )then
         do i= 1, Shape_Element
            Flag_Node( Index_Element_2_Node( i, e ) )= 1
         end do
      end if
   end do


   !====================================================================================================
   write(*,*)'      Detect Elements on Boundary'
   !====================================================================================================

   do i= 1, Number_Node 
      if( Range_Position_Node( 1, 1 ) <= Position_Node( 1, i ) .and. & 
          Position_Node( 1, i ) <= Range_Position_Node( 1, 2 ) )then

         if( Range_Position_Node( 2, 1 ) <= Position_Node( 2, i ) .and. & 
             Position_Node( 2, i ) <= Range_Position_Node( 2, 2 ) )then

            Flag_Node( i )= Flag_Node( i ) +1

         end if
      end if
   end do

   !====================================================================================================
   write(*,*)'      Check the Number of Detected Nodes on Boundary'
   !====================================================================================================
   Counter_Node_Flag= 0

   do i= 1, Number_Node
      if( Flag_Node( i )==2 )then
         Counter_Node_Flag= Counter_Node_Flag +1
      end if
   end do

   if( Counter_Node_Flag==0 )then
      write(*,*)'Counter_Node_Flag=', Counter_Node_Flag
      call Output_Error( 'Detect_Line_on_Waveguide', 93 )
   end if

   do i= 1, Number_Node
      if( Flag_Node( i ) <= -1 .or. Flag_Node( i ) >= 3 )then
         call Output_Error( 'Detect_Line_on_Waveguide', 55 )
      end if
   end do

   !====================================================================================================
   write(*,*)'      Count Nodes Overlapping Boundary in an Element'
   !====================================================================================================

   allocate( Counter_Node_in_Element_1( Number_Element ) )
   allocate( Counter_Node_in_Element_2( Number_Element ) )
   allocate( Local_Node_Number_1_tmp( 2, Number_Element ) )
   allocate( Local_Node_Number_2_tmp( 2, Number_Element ) )

   do e= 1, Number_Element
      Counter_Node_in_Element_1( e )= 0
      Counter_Node_in_Element_2( e )= 0
      do i= 1, 2
         Local_Node_Number_1( i, e )= 0
         Local_Node_Number_2( i, e )= 0
         Local_Node_Number_1_tmp( i, e )= 0
         Local_Node_Number_2_tmp( i, e )= 0
      end do
   end do

   do e= 1, Number_Element
      if( Class_Element( e )==ID_Element )then
         do i= 1, Shape_Element
            if( Flag_Node( Index_Element_2_Node( i, e ) )==2 )then
               Counter_Node_in_Element_1( e )= Counter_Node_in_Element_1( e ) +1

               Local_Node_Number_1_tmp( Counter_Node_in_Element_1( e ), e )= i
            end if
         end do
      end if
   end do

   deallocate( Flag_Node )

   !====================================================================================================
   write(*,*)'      Detect Element on Boundary'
   !====================================================================================================

   Counter_Element_1= 0
   Counter_Element_2= 0

   do e= 1, Number_Element
      if( Counter_Node_in_Element_1( e )==2 )then
         Counter_Element_1= Counter_Element_1 +1

         Element_Number_Boundary_1( Counter_Element_1 )= e
         do i= 1, 2
            Local_Node_Number_1( i, Counter_Element_1 )= Local_Node_Number_1_tmp( i, e )
         end do
      else if( Counter_Node_in_Element_2( e )==2 )then
         Counter_Element_2= Counter_Element_2 +1

         Element_Number_Boundary_2( Counter_Element_2 )= e
         do i= 1, 2
            Local_Node_Number_2( i, Counter_Element_2 )= Local_Node_Number_2_tmp( i, e )
         end do
      end if 
   end do

   deallocate( Counter_Node_in_Element_1 )
   deallocate( Counter_Node_in_Element_2 )
   deallocate( Local_Node_Number_1_tmp )
   deallocate( Local_Node_Number_2_tmp )

   Number_Element_on_Boundary_1= Counter_Element_1
   Number_Element_on_Boundary_2= Counter_Element_2

   !if( Number_Element_on_Boundary_1/=Number_Element_on_Boundary_2 )then
   !   call Output_Error( 'Detect_Line_on_Waveguide', 140 )
   !end if

   !====================================================================================================
   write(*,*)'      Renumber Local Node Number'
   !====================================================================================================

   allocate( Position_LocalNode_in_Elements_1( 2, 2, Number_Element_on_Boundary_1 ) )
   allocate( Position_LocalNode_in_Elements_2( 2, 2, Number_Element_on_Boundary_2 ) )

   do e= 1, Number_Element_on_Boundary_1
      do i= 1, 2
         do j= 1, 2
            Position_LocalNode_in_Elements_1( j, i, e ) & 
            = Position_Node( j, Index_Element_2_Node( Local_Node_Number_1( i, e ), Element_Number_Boundary_1( ( e ) ) ) )
         end do
      end do
   end do

   do e= 1, Number_Element_on_Boundary_2
      do i= 1, 2
         do j= 1, 2
            Position_LocalNode_in_Elements_2( j, i, e ) &
            = Position_Node( j, Index_Element_2_Node( Local_Node_Number_2( i, e ), Element_Number_Boundary_2( ( e ) ) ) )
         end do
      end do
   end do

   !====================================================================================================
   if( Flag_Renumber_Local_Node==Flag_Position_X )then
   write(*,*)'         Based on X Position'
   !====================================================================================================

      allocate( Integer_tmp( Number_Element_on_Boundary_1 ) )

      do e= 1, Number_Element_on_Boundary_1
         if( Position_LocalNode_in_Elements_1( 1, 1, e ) > Position_LocalNode_in_Elements_1( 1, 2, e ) )then
            Integer_tmp( e )= Local_Node_Number_1( 1, e ) 
            Local_Node_Number_1( 1, e )= Local_Node_Number_1( 2, e )
            Local_Node_Number_1( 2, e )= Integer_tmp( e )
         else if( Position_LocalNode_in_Elements_1( 1, 1, e ) == Position_LocalNode_in_Elements_1( 1, 2, e ) )then
            if( Position_LocalNode_in_Elements_1( 2, 1, e ) > Position_LocalNode_in_Elements_1( 2, 2, e ) )then
               Integer_tmp( e )= Local_Node_Number_1( 1, e ) 
               Local_Node_Number_1( 1, e )= Local_Node_Number_1( 2, e )
               Local_Node_Number_1( 2, e )= Integer_tmp( e )
            end if
         end if
      end do

      deallocate( Integer_tmp )
      allocate( Integer_tmp( Number_Element_on_Boundary_2 ) )

      do e= 1, Number_Element_on_Boundary_2
         if( Position_LocalNode_in_Elements_2( 1, 1, e ) > Position_LocalNode_in_Elements_2( 1, 2, e ) )then
            Integer_tmp( e )= Local_Node_Number_2( 1, e ) 
            Local_Node_Number_2( 1, e )= Local_Node_Number_2( 2, e )
            Local_Node_Number_2( 2, e )= Integer_tmp( e )
         else if( Position_LocalNode_in_Elements_2( 1, 1, e ) == Position_LocalNode_in_Elements_2( 1, 2, e ) )then
            if( Position_LocalNode_in_Elements_2( 2, 1, e ) > Position_LocalNode_in_Elements_2( 2, 2, e ) )then
               Integer_tmp( e )= Local_Node_Number_2( 1, e ) 
               Local_Node_Number_2( 1, e )= Local_Node_Number_2( 2, e )
               Local_Node_Number_2( 2, e )= Integer_tmp( e )
            end if
         end if
      end do

      deallocate( Integer_tmp )

   !====================================================================================================
   else if( Flag_Renumber_Local_Node==Flag_Position_Y )then
   write(*,*)'         Based on Y Position'
   !====================================================================================================

      allocate( Integer_tmp( Number_Element_on_Boundary_1 ) )

      do e= 1, Number_Element_on_Boundary_1
         if( Position_LocalNode_in_Elements_1( 2, 1, e ) > Position_LocalNode_in_Elements_1( 2, 2, e ) )then
            Integer_tmp( e )= Local_Node_Number_1( 1, e ) 
            Local_Node_Number_1( 1, e )= Local_Node_Number_1( 2, e )
            Local_Node_Number_1( 2, e )= Integer_tmp( e )
         else if( Position_LocalNode_in_Elements_1( 2, 1, e ) == Position_LocalNode_in_Elements_1( 2, 2, e ) )then
            if( Position_LocalNode_in_Elements_1( 1, 1, e ) > Position_LocalNode_in_Elements_1( 1, 2, e ) )then
               Integer_tmp( e )= Local_Node_Number_1( 1, e ) 
               Local_Node_Number_1( 1, e )= Local_Node_Number_1( 2, e )
               Local_Node_Number_1( 2, e )= Integer_tmp( e )
            end if
         end if
      end do

      deallocate( Integer_tmp )
      allocate( Integer_tmp( Number_Element_on_Boundary_2 ) )

      do e= 1, Number_Element_on_Boundary_2
         if( Position_LocalNode_in_Elements_2( 2, 1, e ) > Position_LocalNode_in_Elements_2( 2, 2, e ) )then
            Integer_tmp( e )= Local_Node_Number_2( 1, e ) 
            Local_Node_Number_2( 1, e )= Local_Node_Number_2( 2, e )
            Local_Node_Number_2( 2, e )= Integer_tmp( e )
         else if( Position_LocalNode_in_Elements_2( 2, 1, e ) == Position_LocalNode_in_Elements_2( 2, 2, e ) )then
            if( Position_LocalNode_in_Elements_2( 1, 1, e ) > Position_LocalNode_in_Elements_2( 1, 2, e ) )then
               Integer_tmp( e )= Local_Node_Number_2( 1, e ) 
               Local_Node_Number_2( 1, e )= Local_Node_Number_2( 2, e )
               Local_Node_Number_2( 2, e )= Integer_tmp( e )
            end if
         end if
      end do

      deallocate( Integer_tmp )

   !====================================================================================================
   else if( Flag_Renumber_Local_Node==Flag_Position_Clock_Rotation )then
   write(*,*)'         Based on Clock Rotation'
   !====================================================================================================

      do e= 1, Number_Element_on_Boundary_1
         do i= 1, 2
            Position_LocalNode_in_Elements_1( 1, i, e )=  Position_LocalNode_in_Elements_1( 1, i, e ) -Position_Origin_X
            Position_LocalNode_in_Elements_1( 2, i, e )=  Position_LocalNode_in_Elements_1( 2, i, e ) -Position_Origin_Y
         end do
      end do
 
      do e= 1, Number_Element_on_Boundary_2
         do i= 1, 2
            Position_LocalNode_in_Elements_2( 1, i, e )=  Position_LocalNode_in_Elements_2( 1, i, e ) -Position_Origin_X
            Position_LocalNode_in_Elements_2( 2, i, e )=  Position_LocalNode_in_Elements_2( 2, i, e ) -Position_Origin_Y
         end do
      end do

      !====================================================================================================
      write(*,*)'            Compute Cross Products'
      !====================================================================================================

      allocate( Cross_Product_1( Number_Element_on_Boundary_1 ) )
      allocate( Cross_Product_2( Number_Element_on_Boundary_2 ) )

      do e= 1, Number_Element_on_Boundary_1
         Cross_Product_1( e ) &
         = Position_LocalNode_in_Elements_1( 1, 1, e ) *Position_LocalNode_in_Elements_1( 2, 2, e ) &
          -Position_LocalNode_in_Elements_1( 2, 1, e ) *Position_LocalNode_in_Elements_1( 1, 2, e ) 
      end do

      do e= 1, Number_Element_on_Boundary_2
         Cross_Product_2( e ) &
         = Position_LocalNode_in_Elements_2( 1, 1, e ) *Position_LocalNode_in_Elements_2( 2, 2, e ) &
          -Position_LocalNode_in_Elements_2( 2, 1, e ) *Position_LocalNode_in_Elements_2( 1, 2, e ) 
      end do

      allocate( Integer_tmp( Number_Element_on_Boundary_1 ) )

      do e= 1, Number_Element_on_Boundary_1
         if( Cross_Product_1( e ) < 0d0  )then
            Integer_tmp( e )= Local_Node_Number_1( 1, e ) 
            Local_Node_Number_1( 1, e )= Local_Node_Number_1( 2, e )
            Local_Node_Number_1( 2, e )= Integer_tmp( e )
         else if( Cross_Product_1( e )==0d0 )then
            call Output_Error( 'Detect_Line_on_Waveguide', 339 )
         end if
      end do

      deallocate( Integer_tmp )
      allocate( Integer_tmp( Number_Element_on_Boundary_2 ) )

      do e= 1, Number_Element_on_Boundary_2
         if( Cross_Product_2( e ) < 0d0  )then
            Integer_tmp( e )= Local_Node_Number_2( 1, e ) 
            Local_Node_Number_2( 1, e )= Local_Node_Number_2( 2, e )
            Local_Node_Number_2( 2, e )= Integer_tmp( e )
         else if( Cross_Product_2( e )==0d0 )then
            call Output_Error( 'Detect_Line_on_Waveguide', 352 )
         end if
      end do

      deallocate( Integer_tmp )
      deallocate( Cross_Product_1 )
      deallocate( Cross_Product_2 )

   !====================================================================================================
   end if
   !====================================================================================================

   deallocate( Position_LocalNode_in_Elements_1 )
   deallocate( Position_LocalNode_in_Elements_2 )

   return
end subroutine Detect_Line_on_Waveguide 



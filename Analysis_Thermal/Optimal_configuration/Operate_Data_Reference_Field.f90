
subroutine Operate_Data_Reference_Field &
       ( Flag_Structure, File_Pass_Reference_Field, &
         Number_Node, Number_Node_Reference, Number_Node_tmp, Temperature_Solution, &
         Circular_Frequency_Normalized, Incident_Angle_Degree, &
         Temperature_Reference, Number_Node_Reference_Read )

   !$use omp_lib
   use Parameters
   implicit none

   integer, intent(in) :: Flag_Structure, Number_Node, Number_Node_Reference, Number_Node_tmp 
   character(len=*), intent(in) :: File_Pass_Reference_Field 
   double precision, intent(in) :: Temperature_Solution( Number_Node )
   double precision, intent(out) :: Temperature_Reference( Number_Node_tmp ) 
   integer, intent(out) :: Number_Node_Reference_Read
   double precision, intent(in) :: Circular_Frequency_Normalized
   double precision, intent(in) :: Incident_Angle_Degree 

   integer :: i, Integer_tmp
   integer :: Integer_part_Circular_Frequency_Normalized 
   integer :: Integer_part_Incident_Angle
   double precision :: Dicimal_part_Circular_Frequency_Normalized 
   double precision :: Dicimal_part_Incident_Angle 

   !character(len=256) :: Filename_Reference_Field
   character(:), allocatable :: Filename_Reference_Field

   integer :: Number_Power_Decimal_Part_Frequency
   integer :: Number_Power_Decimal_Part_Incident_Angle

   integer :: Length_Character_Integer_input, Length_Character_Frequency
   integer :: Length_Character_Incident_Angle, Length_Filename_Reference_Field

   !character, allocatable(:) :: Character_Number_Grid_Scale_DesignDomain
   character(:), allocatable :: Character_Number_Grid_Scale_DesignDomain
   character(:), allocatable :: Character_Frequency_Normalized
   character(:), allocatable :: Character_Incident_Angle_Degree

   !====================================================================================
   write(*,*)'    call Operate_Data_Reference_Field'
   !====================================================================================

   !====================================================================================
   write(*,*)'    Number_Grid_Scale_DesignDomain --> Character_Number_Grid_Scale_DesignDomain '
   !====================================================================================

   Length_Character_Integer_input= log10( dble( Number_Grid_Scale_DesignDomain ) ) +1
    
   allocate( character(Length_Character_Integer_input)::Character_Number_Grid_Scale_DesignDomain )
 
   call Transform_Integer_2_Character &
       ( Number_Grid_Scale_DesignDomain, Length_Character_Integer_input, &
         Character_Number_Grid_Scale_DesignDomain(:) )

   !write(*,*) Character_Number_Grid_Scale_DesignDomain, Length_Character_Integer_input

   !====================================================================================
   write(*,*)'    Circular_Frequency_Normalized --> Character_Frequency_Normalized'
   !====================================================================================

   Integer_part_Circular_Frequency_Normalized= int( Circular_Frequency_Normalized +1d-6 )

   Dicimal_part_Circular_Frequency_Normalized= Circular_Frequency_Normalized -Integer_part_Circular_Frequency_Normalized

   if( Dicimal_part_Circular_Frequency_Normalized < 1d-8 )then
      if( Circular_Frequency_Normalized==0d0 )then
         Length_Character_Frequency= 1
      else
         Length_Character_Frequency= log10( dble( Circular_Frequency_Normalized ) ) +1
      end if

      allocate( character(Length_Character_Frequency)::Character_Frequency_Normalized )

      call Transform_Integer_2_Character &
          ( Integer_part_Circular_Frequency_Normalized, Length_Character_Frequency, &
          Character_Frequency_Normalized(:) )
   else

      call Obtain_Number_Power_Decimal_Part &
          ( Circular_Frequency_Normalized, Number_Power_Decimal_Part_Frequency )
   
      Length_Character_Frequency= log10( dble( Circular_Frequency_Normalized ) ) +2 +Number_Power_Decimal_Part_Frequency
   
      allocate( character(Length_Character_Frequency)::Character_Frequency_Normalized )
   
      call Transform_Real_2_Character &
          ( Circular_Frequency_Normalized, Length_Character_Frequency, &
            Character_Frequency_Normalized(:) )
   
   end if

   !write(*,*) Character_Frequency_Normalized, Length_Character_Frequency

   !====================================================================================
   write(*,*)'   Incident_Angle_Degree --> Character_Incident_Angle_Degree'
   !====================================================================================

   if( Incident_Angle_Degree >= 0d0 )then
      Integer_part_Incident_Angle= int( Incident_Angle_Degree +1d-6 )
   else 
      Integer_part_Incident_Angle= int( Incident_Angle_Degree -1d-6 )
   end if 
   Dicimal_part_Incident_Angle= Incident_Angle_Degree -dble( Integer_part_Incident_Angle )

   write(*,*)'   Incident_Angle_Degree=', Incident_Angle_Degree
   write(*,*)'   Integer_part_Incident_Angle=', Integer_part_Incident_Angle
   write(*,*)'   Dicimal_part_Incident_Angle=', Dicimal_part_Incident_Angle

   if( abs( Incident_Angle_Degree ) < 1d-8 )then

      if( Incident_Angle_Degree==0d0 )then
         Length_Character_Incident_Angle= 1
      else
         Length_Character_Incident_Angle= log10( dble( Circular_Frequency_Normalized ) ) +1
      end if

      allocate( character(Length_Character_Incident_Angle)::Character_Incident_Angle_Degree )

      call Transform_Integer_2_Character &
          ( Integer_part_Incident_Angle, Length_Character_Incident_Angle, &
          Character_Incident_Angle_Degree(:) )

   else if( abs( Dicimal_part_Incident_Angle ) < 1d-8 .and. Incident_Angle_Degree >= 0d0  )then

      if( Integer_part_Incident_Angle==0 )then
         Length_Character_Incident_Angle= 1
      else
         Length_Character_Incident_Angle= log10( dble( Integer_part_Incident_Angle ) ) +1
      end if

      allocate( character(Length_Character_Incident_Angle)::Character_Incident_Angle_Degree )
      call Transform_Integer_2_Character &
          ( Integer_part_Incident_Angle, Length_Character_Incident_Angle, &
          Character_Incident_Angle_Degree(:) )

   else if( abs( Dicimal_part_Incident_Angle ) < 1d-8 .and. Incident_Angle_Degree < 0d0  )then

      if( Integer_part_Incident_Angle==0 )then
         Length_Character_Incident_Angle= 2
      else
         Length_Character_Incident_Angle= log10( dble( abs( Integer_part_Incident_Angle ) ) ) +2
      end if
 


      allocate( character(Length_Character_Incident_Angle)::Character_Incident_Angle_Degree )

      call Transform_Integer_2_Character &
          ( Integer_part_Incident_Angle, Length_Character_Incident_Angle, &
          Character_Incident_Angle_Degree(:) )

   else if( abs( Dicimal_part_Incident_Angle ) >= 1d-8 .and. Incident_Angle_Degree >= 0d0  )then
      call Obtain_Number_Power_Decimal_Part &
         ( Incident_Angle_Degree, Number_Power_Decimal_Part_Incident_Angle )

      if( Integer_part_Incident_Angle==0 )then
         Length_Character_Incident_Angle &
         = 2 &
          +Number_Power_Decimal_Part_Incident_Angle
      else
         Length_Character_Incident_Angle &
         = log10( dble( abs( Integer_part_Incident_Angle ) ) ) +2 &
          +Number_Power_Decimal_Part_Incident_Angle
      end if
   
      allocate( character(Length_Character_Incident_Angle)::Character_Incident_Angle_Degree )
   
      call Transform_Real_2_Character &
          ( Incident_Angle_Degree, Length_Character_Incident_Angle, &
            Character_Incident_Angle_Degree(:) )

   else if( abs( Dicimal_part_Incident_Angle ) >= 1d-8 .and. Incident_Angle_Degree < 0d0  )then

      call Obtain_Number_Power_Decimal_Part &
         ( Incident_Angle_Degree, Number_Power_Decimal_Part_Incident_Angle )

write(*,*)'Integer_part_Incident_Angle=', Integer_part_Incident_Angle

      if( Integer_part_Incident_Angle==0 )then
         Length_Character_Incident_Angle &
         = 3 &
          +Number_Power_Decimal_Part_Incident_Angle
      else
         Length_Character_Incident_Angle &
         = log10( dble( abs( Integer_part_Incident_Angle ) ) ) +3 &
          +Number_Power_Decimal_Part_Incident_Angle
      end if

      allocate( character(Length_Character_Incident_Angle)::Character_Incident_Angle_Degree )
   
      call Transform_Real_2_Character &
          ( Incident_Angle_Degree, Length_Character_Incident_Angle, &
            Character_Incident_Angle_Degree(:) )

   end if

   !====================================================================================
   write(*,*)'    Set Filename_Reference_Field '
   !====================================================================================

   Length_Filename_Reference_Field= 19 +Length_Character_Integer_input +11 +Length_Character_Frequency +4
    
   allocate( character(Length_Filename_Reference_Field)::Filename_Reference_Field )

   Filename_Reference_Field= &
   trim("Reference_Field_OF_"//trim(Character_Number_Grid_Scale_DesignDomain)//&
      "_Frequency_"//trim(Character_Frequency_Normalized)//&
      "_IncidentAngle_"//trim(Character_Incident_Angle_Degree)//".ref")

write(*,*)'Flag_Structure=', Flag_Structure
   if( Flag_Structure==2 )then
   !====================================================================================
   write(*,*)'   Output ', Filename_Reference_Field
   !====================================================================================
      open( 653, file=Filename_Reference_Field, status='replace' )

      write(653,*) Number_Node_tmp 
      do i= 1, Number_Node_tmp 
          !write(653,*) i, Field_Electric_Incident( i ) 
          write(653,*) i, Temperature_Solution( i ) 
          !write(653,*) i, Field_Electric_All( i )
          !658 format(I10,1X,es15.8,1X,es15.8,1X,es15.8)
      end do
      close( 653 )

      do i= 1, Number_Node_tmp
          Temperature_Reference( i )= Temperature_Solution( i )
      end do

   else 
   !====================================================================================
   write(*,*)'   Read ', Filename_Reference_Field
   !====================================================================================
      !open( 653, file=Filename_Reference_Field, action='read' )
      open( 653, file=trim(File_Pass_Reference_Field)//trim(Filename_Reference_Field), action='read' ) 

      read(653,*) Number_Node_Reference_Read 

          if( Number_Node_Reference_Read/=Number_Node_Reference )then
             write(*,*)'Number_Node_Reference=', Number_Node_Reference
             write(*,*)'Number_Node_Reference_Read=', Number_Node_Reference_Read
             call Output_Error( 'Operate_Data_Reference_Field', 682 )
          end if

          do i= 1, Number_Node_tmp
             read(653,*) Integer_tmp, Temperature_Reference( i ) 
          end do
          close( 653 )
   end if

   deallocate( Character_Number_Grid_Scale_DesignDomain )
   deallocate( Character_Frequency_Normalized )
   deallocate( Character_Incident_Angle_Degree )
   deallocate( Filename_Reference_Field )

   return
end subroutine Operate_Data_Reference_Field


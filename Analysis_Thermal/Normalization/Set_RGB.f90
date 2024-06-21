
subroutine Set_RGB( Flag_RGB_Color, Number_Level, RGB, RGB_Lowest, RGB_Highest ) 

   !use omp_lib
   implicit none

   integer, intent(in) :: Flag_RGB_Color, Number_Level 
   integer :: i 
   double precision, intent(out) :: RGB( 3, Number_Level )
   double precision, intent(out) :: RGB_Lowest( 3 ), RGB_Highest( 3 ) 
   integer, allocatable, dimension(:) :: Level_Boundary 
   
   !======================================================================================================================
   ! call Set_RGB
   !======================================================================================================================
   ! Flag_RGB_Color 
   ! 0:Blue-Green_Red
   ! 1:Blue-Green_Red
   ! 2:darkBlue-darkGreen_darkRed
   ! 3:darkBlue-darkGreen_darkRed
   ! 5:DarkBlue-Blue-Green_Red-DarkRed
   !======================================================================================================================

   
   !======================================================================================================================
   if( Flag_RGB_Color==0 )then
   !======================================================================================================================

      allocate( Level_Boundary( 5 ) )
   
      do i= 1, 5
         Level_Boundary( i )= 1.0d0 +( dble( Number_Level ) -1.0d0 )*( i-1 )/4.0d0
      end do
   
      !$omp parallel do default( none ) &
      !$omp private( i ) &
      !$omp shared( Number_Level, Level_Boundary, RGB ) 
      do i= 1, Number_Level
         ! Blue
         if( Level_Boundary( 1 ) <= i .and. i <= Level_Boundary( 2 )+1 )then
            RGB( 3, i )= 1.0d0 
         else if( Level_Boundary( 2 )+1 <= i .and. i <= Level_Boundary( 3 )-1 )then
            RGB( 3, i )= -1.0d0/( Level_Boundary( 3 )-1-Level_Boundary( 2 )-1 )*( i-Level_Boundary( 3 )+1 ) 
         else
            RGB( 3, i )= 0.0d0 
         end if
   
         ! Green
         if( Level_Boundary( 1 )+1 <= i .and. i <= Level_Boundary( 2 ) )then
            RGB( 2, i )= 1.0d0/( Level_Boundary( 2 )-Level_Boundary( 1 )-1 )*( i-Level_Boundary( 1 )-1 ) 
         else if( Level_Boundary( 2 )+1 <= i .and. i <= Level_Boundary( 4 )-1 )then
            RGB( 2, i )= 1.0d0 
         else if( Level_Boundary( 4 ) <= i .and. i <= Level_Boundary( 5 )-1 )then
            RGB( 2, i )= -1.0d0/( Level_Boundary( 5 )-1-Level_Boundary( 4 ) )*( i-Level_Boundary( 5 )+1 ) 
         else
            RGB( 2, i )= 0.0d0 
         end if
    
         ! Red 
         if( Level_Boundary( 3 )+1 <= i .and. i <= Level_Boundary( 4 )-1 )then
            RGB( 1, i )= 1.0d0/( Level_Boundary( 4 )-1-Level_Boundary( 3 )-1 )*( i-Level_Boundary( 3 )-1 ) 
         else if( Level_Boundary( 4 )-1 <= i .and. i <= Level_Boundary( 5 ) )then
            RGB( 1, i )= 1.0d0 
         else
            RGB( 1, i )= 0.0d0 
         end if
      end do
   
      deallocate( Level_Boundary )

      RGB_Lowest( 1 )= 0.0d0 
      RGB_Lowest( 2 )= 0.0d0
      RGB_Lowest( 3 )= 1.0d0
      RGB_Highest( 1 )= 1.0d0
      RGB_Highest( 2 )= 0.0d0
      RGB_Highest( 3 )= 0.0d0
 
   else if( Flag_RGB_Color==1 )then

      allocate( Level_Boundary( 5 ) )

      do i= 1, 5
         Level_Boundary( i )= 1.0d0 +( dble( Number_Level ) -1.0d0 )*( i-1 )/4.0d0
      end do

      !$omp parallel do default( none ) &
      !$omp private( i ) &
      !$omp shared( Number_Level, Level_Boundary, RGB ) 
      do i= 1, Number_Level
         ! Blue
         if( Level_Boundary( 1 ) <= i .and. i <= Level_Boundary( 3 )+1 )then
            RGB( 3, i )= 1.0d0
         else if( Level_Boundary( 3 )+1 <= i .and. i <= Level_Boundary( 4 )-1 )then
            RGB( 3, i )= -1.0d0/( Level_Boundary( 4 )-1-Level_Boundary( 3 )-1 )*( i-Level_Boundary( 4 )+1 )
         else
            RGB( 3, i )= 0.0d0
         end if

         ! Green
         if( Level_Boundary( 1 )+1 <= i .and. i <= Level_Boundary( 2 ) )then
            RGB( 2, i )= 1.0d0/( Level_Boundary( 2 )-Level_Boundary( 1 )-1 )*( i-Level_Boundary( 1 )-1 )
         else if( Level_Boundary( 2 )+1 <= i .and. i <= Level_Boundary( 4 )-1 )then
            RGB( 2, i )= 1.0d0
         else if( Level_Boundary( 4 ) <= i .and. i <= Level_Boundary( 5 )-1 )then
            RGB( 2, i )= -1.0d0/( Level_Boundary( 5 )-1-Level_Boundary( 4 ) )*( i-Level_Boundary( 5 )+1 )
         else
            RGB( 2, i )= 0.0d0
         end if

         ! Red 
         if( Level_Boundary( 2 )+1 <= i .and. i <= Level_Boundary( 3 )-1 )then
            RGB( 1, i )= 1.0d0/( Level_Boundary( 3 )-1-Level_Boundary( 2 )-1 )*( i-Level_Boundary( 2 )-1 )
         else if( Level_Boundary( 3 )-1 <= i .and. i <= Level_Boundary( 5 ) )then
            RGB( 1, i )= 1.0d0
         else
            RGB( 1, i )= 0.0d0
         end if
      end do

      deallocate( Level_Boundary )

      RGB_Lowest( 1 )= 0.0d0 
      RGB_Lowest( 2 )= 0.0d0
      RGB_Lowest( 3 )= 0.5d0
      RGB_Highest( 1 )= 0.5d0
      RGB_Highest( 2 )= 0.0d0
      RGB_Highest( 3 )= 0.0d0

   else if( Flag_RGB_Color==2 )then

      !$omp parallel do default( none ) &
      !$omp private( i ) &
      !$omp shared( Number_Level, RGB ) 
      do i= 1, Number_Level
         RGB( 1, i )= 1.0d0/( Number_Level-1 )*( i-1 )
         RGB( 3, i )= -1.0d0/( Number_Level-1 )*( i-Number_Level ) 
      
         if( dble( i ) <= Number_Level/2.0d0 )then
            RGB( 2, i )= 2.0d0/( Number_Level-2.0d0 )*( i-1 ) 
         else
            RGB( 2, i )= -2.0d0/( Number_Level-2.0d0 )*( i -Number_Level ) 
         end if 
      end do
   
   else if( Flag_RGB_Color==3 )then

      allocate( Level_Boundary( 7 ) )
   
      do i= 1, 7
         Level_Boundary( i )= 1.0d0 +( dble( Number_Level ) -1.0d0 )*( i-1 )/6.0d0
      end do
   
      !$omp parallel do default( none ) &
      !$omp private( i ) &
      !$omp shared( Number_Level, Level_Boundary, RGB ) 
      do i= 1, Number_Level
         if( Level_Boundary( 1 ) <= i .and. i <= Level_Boundary( 4 ) )then
            RGB( 3, i )= -1.0d0/( Level_Boundary( 4 )-Level_Boundary( 1 ) )*( i-Level_Boundary( 4 ) ) 
         else
            RGB( 3, i )= 0.0d0 
         end if
   
         if( Level_Boundary( 1 ) <= i .and. i <= Level_Boundary( 4 ) )then
            RGB( 2, i )= 1.0d0/( Level_Boundary( 4 )-Level_Boundary( 1 ) )*( i-Level_Boundary( 1 ) ) 
         else if( Level_Boundary( 4 ) <= i .and. i <= Level_Boundary( 7 ) )then
            RGB( 2, i )= -1.0d0/( Level_Boundary( 7 )-Level_Boundary( 4 ) )*( i-Level_Boundary( 7 ) ) 
         else
            RGB( 2, i )= 0.0d0 
         end if
   
         if( Level_Boundary( 4 ) <= i .and. i <= Level_Boundary( 7 ) )then
            RGB( 1, i )= 1.0d0/( Level_Boundary( 7 )-Level_Boundary( 4 ) )*( i-Level_Boundary( 4 ) ) 
         else
            RGB( 1, i )= 0.0d0 
         end if
      end do
      
      deallocate( Level_Boundary )
   
      RGB_Lowest( 1 )= 0.0d0 
      RGB_Lowest( 2 )= 0.0d0
      RGB_Lowest( 3 )= 1.0d0
      RGB_Highest( 1 )= 1.0d0
      RGB_Highest( 2 )= 0.0d0
      RGB_Highest( 3 )= 0.0d0

    else if( Flag_RGB_Color==4 )then

      allocate( Level_Boundary( 7 ) )

      do i= 1, 7
         Level_Boundary( i )= 1.0d0 +( dble( Number_Level ) -1.0d0 )*( i-1 )/6.0d0
      end do

      !$omp parallel do default( none ) &
      !$omp private( i ) &
      !$omp shared( Number_Level, Level_Boundary, RGB ) 
      do i= 1, Number_Level
         if( Level_Boundary( 4 ) <= i .and. i <= Level_Boundary( 5 ) )then
            RGB( 1, i )= 1.0d0/( Level_Boundary( 5 )-Level_Boundary( 4 ) )*( i-Level_Boundary( 4 ) )
         else if( Level_Boundary( 5 ) <= i )then
            RGB( 1, i )= 1.0d0
         else
            RGB( 1, i )= 0.0d0
         end if

         if( Level_Boundary( 1 ) <= i .and. i <= Level_Boundary( 2 ) )then
            RGB( 2, i )= 0.0d0
         else if( Level_Boundary( 2 ) <= i .and. i <= Level_Boundary( 3 ) )then
            RGB( 2, i )= 1.0d0/( Level_Boundary( 3 )-Level_Boundary( 2 ) )*( i-Level_Boundary( 2 ) )
         else if( Level_Boundary( 3 ) <= i .and. i <= Level_Boundary( 5 ) )then
            RGB( 2, i )= 1.0d0
         else if( Level_Boundary( 5 ) <= i .and. i <= Level_Boundary( 6 ) )then
            RGB( 2, i )= 1d0 -1.0d0/( Level_Boundary( 6 )-Level_Boundary( 5 ) )*( i-Level_Boundary( 5 ) )
         else if( Level_Boundary( 6 ) <= i .and. i <= Level_Boundary( 7 ) )then
            RGB( 2, i )= 1.0d0/( Level_Boundary( 7 )-Level_Boundary( 6 ) )*( i-Level_Boundary( 6 ) )
         end if

         if( Level_Boundary( 1 ) <= i .and. i <= Level_Boundary( 2 ) )then
            RGB( 3, i )= 1.0d0/( Level_Boundary( 2 )-Level_Boundary( 1 ) )*( i-Level_Boundary( 1 ) )
         else if( Level_Boundary( 2 ) <= i .and. i <= Level_Boundary( 3 ) )then
            RGB( 3, i )= 1.0d0
         else if( Level_Boundary( 3 ) <= i .and. i <= Level_Boundary( 4 ) )then
            RGB( 3, i )= 1d0 -1.0d0/( Level_Boundary( 4 )-Level_Boundary( 3 ) )*( i-Level_Boundary( 3 ) )
         else if( Level_Boundary( 4 ) <= i .and. i <= Level_Boundary( 6 ) )then
            RGB( 3, i )= 0.0d0
         else if( Level_Boundary( 6 ) <= i .and. i <= Level_Boundary( 7 ) )then
            RGB( 3, i )= 1.0d0/( Level_Boundary( 7 )-Level_Boundary( 6 ) )*( i-Level_Boundary( 6 ) )
         end if

      end do

      deallocate( Level_Boundary )

      RGB_Lowest( 1 )= 0.0d0 
      RGB_Lowest( 2 )= 0.0d0
      RGB_Lowest( 3 )= 1.0d0
      RGB_Highest( 1 )= 1.0d0
      RGB_Highest( 2 )= 0.0d0
      RGB_Highest( 3 )= 0.0d0

    else if( Flag_RGB_Color==5 )then

      allocate( Level_Boundary( 7 ) )

      do i= 1, 7
         Level_Boundary( i )= 1.0d0 +( dble( Number_Level ) -1.0d0 )*( i-1 )/6.0d0
      end do

      !$omp parallel do default( none ) &
      !$omp private( i ) &
      !$omp shared( Number_Level, Level_Boundary, RGB ) 
      do i= 1, Number_Level
         if( Level_Boundary( 4 ) <= i .and. i <= Level_Boundary( 5 ) )then
            RGB( 1, i )= 1.0d0/( Level_Boundary( 5 )-Level_Boundary( 4 ) )*( i-Level_Boundary( 4 ) )
         else if( Level_Boundary( 5 ) <= i .and. i <= Level_Boundary( 6 ) )then
            RGB( 1, i )= 1.0d0
         else if( Level_Boundary( 6 ) <= i )then
            RGB( 1, i )= ( 0.5d0 -1.0d0  )/( Level_Boundary( 7 )-Level_Boundary( 6 ) )*( i-Level_Boundary( 6 ) ) +1.0d0
         else
            RGB( 1, i )= 0.0d0
         end if

         if( Level_Boundary( 1 ) <= i .and. i <= Level_Boundary( 2 ) )then
            RGB( 2, i )= 0.0d0
         else if( Level_Boundary( 2 ) <= i .and. i <= Level_Boundary( 3 ) )then
            RGB( 2, i )= 1.0d0/( Level_Boundary( 3 )-Level_Boundary( 2 ) )*( i-Level_Boundary( 2 ) )
         else if( Level_Boundary( 3 ) <= i .and. i <= Level_Boundary( 5 ) )then
            RGB( 2, i )= 1.0d0
         else if( Level_Boundary( 5 ) <= i .and. i <= Level_Boundary( 6 ) )then
            RGB( 2, i )= 1d0 -1.0d0/( Level_Boundary( 6 )-Level_Boundary( 5 ) )*( i-Level_Boundary( 5 ) )
         else if( Level_Boundary( 6 ) <= i .and. i <= Level_Boundary( 7 ) )then
            !RGB( 2, i )= 1.0d0/( Level_Boundary( 7 )-Level_Boundary( 6 ) )*(
            !i-Level_Boundary( 6 ) )
            RGB( 2, i )= 0.0d0
         end if

         if( Level_Boundary( 1 ) <= i .and. i <= Level_Boundary( 2 ) )then
            !RGB( 3, i )= 1.0d0/( Level_Boundary( 2 )-Level_Boundary( 1 ) )*(
            !i-Level_Boundary( 1 ) )
            RGB( 3, i )= ( 1.0d0 -0.5d0 )/( Level_Boundary( 2 )-Level_Boundary( 1 ) )*( i-Level_Boundary( 1 ) ) +0.5d0
         else if( Level_Boundary( 2 ) <= i .and. i <= Level_Boundary( 3 ) )then
            RGB( 3, i )= 1.0d0
         else if( Level_Boundary( 3 ) <= i .and. i <= Level_Boundary( 4 ) )then
            RGB( 3, i )= 1d0 -1.0d0/( Level_Boundary( 4 )-Level_Boundary( 3 ) )*( i-Level_Boundary( 3 ) )
         else if( Level_Boundary( 4 ) <= i .and. i <= Level_Boundary( 6 ) )then
            RGB( 3, i )= 0.0d0
         else if( Level_Boundary( 6 ) <= i .and. i <= Level_Boundary( 7 ) )then
            !RGB( 3, i )= 1.0d0/( Level_Boundary( 7 )-Level_Boundary( 6 ) )*(
            !i-Level_Boundary( 6 ) )
            RGB( 3, i )= 0.0d0
         end if

      end do

      deallocate( Level_Boundary )

      RGB_Lowest( 1 )= 0.0d0 
      RGB_Lowest( 2 )= 0.0d0
      RGB_Lowest( 3 )= 0.5d0
      RGB_Highest( 1 )= 0.5d0
      RGB_Highest( 2 )= 0.0d0
      RGB_Highest( 3 )= 0.0d0

    else if( Flag_RGB_Color==6 )then

      allocate( Level_Boundary( 7 ) )

      do i= 1, 7
         Level_Boundary( i )= 1.0d0 +( dble( Number_Level ) -1.0d0 )*( i-1 )/6.0d0
      end do

      !$omp parallel do default( none ) &
      !$omp private( i ) &
      !$omp shared( Number_Level, Level_Boundary, RGB ) 
      do i= 1, Number_Level
         if( Level_Boundary( 1 ) <= i .and. i <= Level_Boundary( 3 ) )then
            RGB( 1, i )= 1d0 -1.0d0/( Level_Boundary( 3 )-Level_Boundary( 1 ) )*( i-Level_Boundary( 1 ) )
         else if( Level_Boundary( 3 ) <= i .and. i <= Level_Boundary( 5 ) )then
            RGB( 1, i )= 0.0d0
         else if( Level_Boundary( 5 ) <= i .and. i <= Level_Boundary( 7 ) )then
            RGB( 1, i )= 1.0d0/( Level_Boundary( 7 )-Level_Boundary( 5 ) )*( i-Level_Boundary( 5 ) )
         end if

         if( Level_Boundary( 1 ) <= i .and. i <= Level_Boundary( 3 ) )then
            RGB( 2, i )= 1d0 -1.0d0/( Level_Boundary( 3 )-Level_Boundary( 1 ) )*( i-Level_Boundary( 1 ) )
         else if( Level_Boundary( 3 ) <= i .and. i <= Level_Boundary( 5 ) )then
            RGB( 2, i )= 1.0d0/( Level_Boundary( 5 )-Level_Boundary( 3 ) )*( i-Level_Boundary( 3 ) )
         elseif( Level_Boundary( 5 ) <= i .and. i <= Level_Boundary( 7 ) )then
            RGB( 2, i )= 1d0 -1.0d0/( Level_Boundary( 7 )-Level_Boundary( 5 ) )*( i-Level_Boundary( 5 ) )
         end if

         if( i <= Level_Boundary( 3 ) )then
            RGB( 3, i )= 1.0d0
         else if( Level_Boundary( 3 ) <= i .and. i <= Level_Boundary( 5 ) )then
            RGB( 3, i )= 1d0 -1.0d0/( Level_Boundary( 5 )-Level_Boundary( 3 ) )*( i-Level_Boundary( 3 ) )
         else if( Level_Boundary( 5 ) <= i )then
            RGB( 3, i )= 0.0d0
         end if

      end do

      deallocate( Level_Boundary )

      RGB_Lowest( 1 )= 1.0d0
      RGB_Lowest( 2 )= 1.0d0
      RGB_Lowest( 3 )= 1.0d0
      RGB_Highest( 1 )= 0.5d0
      RGB_Highest( 2 )= 0.0d0
      RGB_Highest( 3 )= 0.0d0

    else if( Flag_RGB_Color==7 )then

      allocate( Level_Boundary( 7 ) )

      do i= 1, 7
         Level_Boundary( i )= 1.0d0 +( dble( Number_Level ) -1.0d0 )*( i-1 )/6.0d0
      end do

      !$omp parallel do default( none ) &
      !$omp private( i ) &
      !$omp shared( Number_Level, Level_Boundary, RGB ) 
      do i= 1, Number_Level
         if( Level_Boundary( 4 ) <= i .and. i <= Level_Boundary( 5 ) )then
            RGB( 1, i )= 1.0d0/( Level_Boundary( 5 )-Level_Boundary( 4 ) )*( i-Level_Boundary( 4 ) )
         else if( Level_Boundary( 5 ) <= i .and. i <= Level_Boundary( 6 ) )then
            RGB( 1, i )= 1.0d0
         else if( Level_Boundary( 6 ) <= i )then
            RGB( 1, i )= ( 0.5d0 -1.0d0  )/( Level_Boundary( 7 )-Level_Boundary( 6 ) )*( i-Level_Boundary( 6 ) ) +1.0d0
         else
            RGB( 1, i )= 0.0d0
         end if

         if( Level_Boundary( 1 ) <= i .and. i <= Level_Boundary( 2 ) )then
            RGB( 2, i )= 0.0d0
         else if( Level_Boundary( 2 ) <= i .and. i <= Level_Boundary( 3 ) )then
            RGB( 2, i )= 1.0d0/( Level_Boundary( 3 )-Level_Boundary( 2 ) )*( i-Level_Boundary( 2 ) )
         else if( Level_Boundary( 3 ) <= i .and. i <= Level_Boundary( 5 ) )then
            RGB( 2, i )= 1.0d0
         else if( Level_Boundary( 5 ) <= i .and. i <= Level_Boundary( 6 ) )then
            RGB( 2, i )= 1d0 -1.0d0/( Level_Boundary( 6 )-Level_Boundary( 5 ) )*( i-Level_Boundary( 5 ) )
         else if( Level_Boundary( 6 ) <= i .and. i <= Level_Boundary( 7 ) )then
            !RGB( 2, i )= 1.0d0/( Level_Boundary( 7 )-Level_Boundary( 6 ) )*(
            !i-Level_Boundary( 6 ) )
            RGB( 2, i )= 0.0d0
         end if

         if( Level_Boundary( 1 ) <= i .and. i <= Level_Boundary( 2 ) )then
            !RGB( 3, i )= 1.0d0/( Level_Boundary( 2 )-Level_Boundary( 1 ) )*(
            !i-Level_Boundary( 1 ) )
            RGB( 3, i )= ( 1.0d0 -0.5d0 )/( Level_Boundary( 2 )-Level_Boundary( 1 ) )*( i-Level_Boundary( 1 ) ) +0.5d0
         else if( Level_Boundary( 2 ) <= i .and. i <= Level_Boundary( 3 ) )then
            RGB( 3, i )= 1.0d0
         else if( Level_Boundary( 3 ) <= i .and. i <= Level_Boundary( 4 ) )then
            RGB( 3, i )= 1d0 -1.0d0/( Level_Boundary( 4 )-Level_Boundary( 3 ) )*( i-Level_Boundary( 3 ) )
         else if( Level_Boundary( 4 ) <= i .and. i <= Level_Boundary( 6 ) )then
            RGB( 3, i )= 0.0d0
         else if( Level_Boundary( 6 ) <= i .and. i <= Level_Boundary( 7 ) )then
            !RGB( 3, i )= 1.0d0/( Level_Boundary( 7 )-Level_Boundary( 6 ) )*(
            !i-Level_Boundary( 6 ) )
            RGB( 3, i )= 0.0d0
         end if

      end do

      deallocate( Level_Boundary )

      RGB_Lowest( 1 )= 0.0d0 
      RGB_Lowest( 2 )= 0.0d0
      RGB_Lowest( 3 )= 0.0d0
      RGB_Highest( 1 )= 1.0d0
      RGB_Highest( 2 )= 1.0d0
      RGB_Highest( 3 )= 1.0d0

    else if( Flag_RGB_Color==8 )then

      allocate( Level_Boundary( 4 ) )

      do i= 1, 4
         Level_Boundary( i )= 1.0d0 +( dble( Number_Level ) -1.0d0 )*( i-1 )/3.0d0
      end do

      !$omp parallel do default( none ) &
      !$omp private( i ) &
      !$omp shared( Number_Level, Level_Boundary, RGB ) 
      do i= 1, Number_Level
         if( Level_Boundary( 1 ) <= i .and. i <= Level_Boundary( 2 ) )then
            !RGB( 1, i )= 1.0d0/( Level_Boundary( 2 )-Level_Boundary( 1 ) )*( i-Level_Boundary( 1 ) )
            RGB( 1, i )= 0.95d0/( Level_Boundary( 2 )-Level_Boundary( 1 ) )*( i-Level_Boundary( 1 ) )
         else if( Level_Boundary( 2 ) <= i .and. i <= Level_Boundary( 4 ) )then
            !RGB( 1, i )= 1.0d0
            RGB( 1, i )= 0.95d0
         end if

         if( Level_Boundary( 1 ) <= i .and. i <= Level_Boundary( 2 ) )then
            RGB( 2, i )= 0.0d0
         else if( Level_Boundary( 2 ) <= i .and. i <= Level_Boundary( 3 ) )then
            !RGB( 2, i )= 1.0d0/( Level_Boundary( 3 )-Level_Boundary( 2 ) )*( i-Level_Boundary( 2 ) )
            RGB( 2, i )= 0.95d0/( Level_Boundary( 3 )-Level_Boundary( 2 ) )*( i-Level_Boundary( 2 ) )
         else if( Level_Boundary( 3 ) <= i .and. i <= Level_Boundary( 4 ) )then
            !RGB( 2, i )= 1.0d0
            RGB( 2, i )= 0.95d0
         end if

         if( Level_Boundary( 3 ) <= i .and. i <= Level_Boundary( 4 ) )then
            !RGB( 3, i )= 1.0d0/( Level_Boundary( 4 )-Level_Boundary( 3 ) )*( i-Level_Boundary( 3 ) ) 
            RGB( 3, i )= 0.95d0/( Level_Boundary( 4 )-Level_Boundary( 3 ) )*( i-Level_Boundary( 3 ) ) 
         else if( Level_Boundary( 1 ) <= i .and. i <= Level_Boundary( 3 ) )then
            RGB( 3, i )= 0.0d0
         end if

      end do

      deallocate( Level_Boundary )

      RGB_Lowest( 1 )= 0.0d0 
      RGB_Lowest( 2 )= 0.0d0
      RGB_Lowest( 3 )= 0.0d0
      RGB_Highest( 1 )= 1.0d0
      RGB_Highest( 2 )= 1.0d0
      RGB_Highest( 3 )= 1.0d0

    else if( Flag_RGB_Color==9 )then

      allocate( Level_Boundary( 7 ) )
   
      do i= 1, 7
         Level_Boundary( i )= 1.0d0 +( dble( Number_Level ) -1.0d0 )*( i-1 )/6.0d0
      end do
   
      !$omp parallel do default( none ) &
      !$omp private( i ) &
      !$omp shared( Number_Level, Level_Boundary, RGB ) 
      do i= 1, Number_Level
         if( Level_Boundary( 4 ) <= i .and. i <= Level_Boundary( 5 ) )then
            RGB( 1, i )= 1.0d0/( Level_Boundary( 5 )-Level_Boundary( 4 ) )*( i-Level_Boundary( 4 ) ) 
         else if( Level_Boundary( 5 ) <= i )then
            RGB( 1, i )= 1.0d0 
         else
            RGB( 1, i )= 0.0d0 
         end if
 
         if( Level_Boundary( 1 ) <= i .and. i <= Level_Boundary( 2 ) )then
            RGB( 2, i )= 0.0d0 
         else if( Level_Boundary( 2 ) <= i .and. i <= Level_Boundary( 3 ) )then
            RGB( 2, i )= 1.0d0/( Level_Boundary( 3 )-Level_Boundary( 2 ) )*( i-Level_Boundary( 2 ) ) 
         else if( Level_Boundary( 3 ) <= i .and. i <= Level_Boundary( 5 ) )then
            RGB( 2, i )= 1.0d0 
         else if( Level_Boundary( 5 ) <= i .and. i <= Level_Boundary( 6 ) )then
            RGB( 2, i )= 1d0 -1.0d0/( Level_Boundary( 6 )-Level_Boundary( 5 ) )*( i-Level_Boundary( 5 ) ) 
         else if( Level_Boundary( 6 ) <= i .and. i <= Level_Boundary( 7 ) )then
            RGB( 2, i )= 1.0d0/( Level_Boundary( 7 )-Level_Boundary( 6 ) )*( i-Level_Boundary( 6 ) ) 
         end if
   
         if( Level_Boundary( 1 ) <= i .and. i <= Level_Boundary( 2 ) )then
            RGB( 3, i )= 1.0d0/( Level_Boundary( 2 )-Level_Boundary( 1 ) )*( i-Level_Boundary( 1 ) ) 
         else if( Level_Boundary( 2 ) <= i .and. i <= Level_Boundary( 3 ) )then
            RGB( 3, i )= 1.0d0 
         else if( Level_Boundary( 3 ) <= i .and. i <= Level_Boundary( 4 ) )then
            RGB( 3, i )= 1d0 -1.0d0/( Level_Boundary( 4 )-Level_Boundary( 3 ) )*( i-Level_Boundary( 3 ) ) 
         else if( Level_Boundary( 4 ) <= i .and. i <= Level_Boundary( 6 ) )then
            RGB( 3, i )= 0.0d0 
         else if( Level_Boundary( 6 ) <= i .and. i <= Level_Boundary( 7 ) )then
            RGB( 3, i )= 1.0d0/( Level_Boundary( 7 )-Level_Boundary( 6 ) )*( i-Level_Boundary( 6 ) ) 
         end if
   
      end do
      
      deallocate( Level_Boundary )

      RGB_Lowest( 1 )= 0.0d0 
      RGB_Lowest( 2 )= 0.0d0
      RGB_Lowest( 3 )= 0.0d0
      !RGB_Highest( 1 )= 1.0d0
      !RGB_Highest( 2 )= 1.0d0
      !RGB_Highest( 3 )= 1.0d0
      RGB_Highest( 1 )= 0.95d0
      RGB_Highest( 2 )= 0.95d0
      RGB_Highest( 3 )= 0.95d0


   else
      write(*,*)'Flag_RGB_Color=', Flag_RGB_Color
      call Output_Error( 'Set_RGB', 259 ) 
   end if 

   return
end subroutine Set_RGB




subroutine Format_Global_aa_DMUMPS &
           ( GlobalMatrix, J_GlobalMatrix, & 
             Number_Node, Number_Node_InOneLow, Number_NonZero, &
             Flag_Symmetric, Flag_Matrix_Upper, &
             !====================================================
             aa_zmumps, ia_zmumps, ja_zmumps )

   !$ use omp_lib
   use Parameters
   implicit none  
  
   integer :: i, j 
   integer :: Counter_NonZero
   
   integer, intent(in) :: Number_Node, Number_Node_InOneLow
   integer, intent(in) :: Flag_Symmetric, Flag_Matrix_Upper 
   integer, intent(in) :: Number_NonZero
   double precision, intent(in) :: GlobalMatrix( Number_Node, Number_Node_InOneLow )
   integer, intent(inout) :: J_GlobalMatrix( Number_Node, Number_Node_InOneLow )
  
   double precision, intent(out) :: aa_zmumps( Number_NonZero )
   integer, intent(out) :: ia_zmumps( Number_NonZero ), ja_zmumps( Number_NonZero )

   !====================================================================================
   write(*,*)'      call Format_Global_aa_DMUMPS'
   !====================================================================================

   !do i= 1, Number_Node 
   !   do j= 1, Number_Node_InOneLow
   !      if( 0 < J_GlobalMatrix( i, j ) .and. J_GlobalMatrix( i, j ) <= Number_Node )then
   !      if( i > J_GlobalMatrix( i, j ) )then
   !         if( GlobalMatrix( i, j )==Zero )then
   !            write(*,*)'GlobalMatrix( i, j )=', GlobalMatrix( i, j )
   !            write(*,*)'i=', i, 'j=', j
   !            write(*,*)'J_GlobalMatrix( i, j )=', J_GlobalMatrix( i, j )
   !            write(*,*)'Format_Global_aa_DMUMPS.f90 31 '
   !            stop
   !         end if
   !         end if
   !      end if
   !   end do
   !end do

   !====================================================================================
   if( Flag_Symmetric==Flag_On .and. Flag_Matrix_Upper==Flag_Off )then
   !====================================================================================
  
      Counter_NonZero= 0 
    
      do i= 1, Number_NonZero 
         aa_zmumps( i )= Zero
         ia_zmumps( i )= 0
         ja_zmumps( i )= 0
      end do
   
      do i= 1, Number_Node 
         do j= 1, Number_Node_InOneLow
            if( 0 < J_GlobalMatrix( i, j ) .and. J_GlobalMatrix( i, j ) <= Number_Node )then
               if( i >= J_GlobalMatrix( i, j ) )then ! Symmetric
                  Counter_NonZero= Counter_NonZero +1
                  if( Counter_NonZero > Number_NonZero )then
                     write(*,*)'ERROR :Counter_NonZero=' ,Counter_NonZero, '>', & 
                           'Number_NonZero', Number_NonZero
                     write(*,*)'Format_Global_aa_DMUMPS.f90 65'
                     stop
                  end if
               
                  aa_zmumps( Counter_NonZero )= GlobalMatrix( i, j )
                  ia_zmumps( Counter_NonZero )= i
                  ja_zmumps( Counter_NonZero )= J_GlobalMatrix( i, j )
               end if
            end if
         end do
      end do
     
      if( Counter_NonZero /= Number_NonZero)then
         write(*,*)'Counter_NonZero /= Number_NonZero'
         write(*,*)'Format_Global_aa_DMUMPS.f90 59 '
         stop
      end if
     
      do i= 1, Number_NonZero
         if( ja_zmumps( i ) <= 0 )then
     
            write(*,*)' ja_zmumps( i ) <= 0 '
            write(*,*)'Format_Global_aa_DMUMPS.f90 67 '
            stop
         end if       
         if( Number_Node < ja_zmumps( i )   )then  
    
            write(*,*)' Number_Node < ja_zmumps( i )  '
            write(*,*)'Format_Global_aa_DMUMPS.f90 73 '
            stop
   
         end if
      end do

   !====================================================================================
   else if( Flag_Symmetric==Flag_On .and. Flag_Matrix_Upper==Flag_On )then
   !====================================================================================
  
      Counter_NonZero= 0 
    
      do i= 1, Number_NonZero 
         aa_zmumps( i )= Zero
         ia_zmumps( i )= 0
         ja_zmumps( i )= 0
      end do
   
      do i= 1, Number_Node 
         do j= 1, Number_Node_InOneLow
            if( 0 < J_GlobalMatrix( i, j ) .and. J_GlobalMatrix( i, j ) <= Number_Node )then
               if( i <= J_GlobalMatrix( i, j ) )then ! Symmetric
                  Counter_NonZero= Counter_NonZero +1
                  if( Counter_NonZero > Number_NonZero )then
                     write(*,*)'ERROR :Counter_NonZero=' ,Counter_NonZero, '>', & 
                           'Number_NonZero', Number_NonZero
                     write(*,*)'Format_Global_aa_DMUMPS.f90 65'
                     stop
                  end if
               
                  aa_zmumps( Counter_NonZero )= GlobalMatrix( i, j )
                  ia_zmumps( Counter_NonZero )= i
                  ja_zmumps( Counter_NonZero )= J_GlobalMatrix( i, j )
               end if
            end if
         end do
      end do
     
      if( Counter_NonZero /= Number_NonZero)then
         write(*,*)'Counter_NonZero /= Number_NonZero'
         write(*,*)'Format_Global_aa_DMUMPS.f90 59 '
         stop
      end if
     
      do i= 1, Number_NonZero
         if( ja_zmumps( i ) <= 0 )then
     
            write(*,*)' ja_zmumps( i ) <= 0 '
            write(*,*)'Format_Global_aa_DMUMPS.f90 67 '
            stop
         end if       
         if( Number_Node < ja_zmumps( i )   )then  
    
            write(*,*)' Number_Node < ja_zmumps( i )  '
            write(*,*)'Format_Global_aa_DMUMPS.f90 73 '
            stop
   
         end if
      end do

   !====================================================================================
   else if( Flag_Symmetric==Flag_Off )then
   !====================================================================================

      Counter_NonZero= 0 
    
      do i= 1, Number_NonZero 
         aa_zmumps( i )= Zero
         ia_zmumps( i )= 0
         ja_zmumps( i )= 0
      end do
   
      do i= 1, Number_Node 
         do j= 1, Number_Node_InOneLow
            if( 0 < J_GlobalMatrix( i, j ) .and. J_GlobalMatrix( i, j ) <= Number_Node )then
               !if( i >= J_GlobalMatrix( i, j ) )then ! Symmetric
                  Counter_NonZero= Counter_NonZero +1
                  if( Counter_NonZero > Number_NonZero )then
                     write(*,*)'ERROR :Counter_NonZero=' ,Counter_NonZero, '>', & 
                           'Number_NonZero', Number_NonZero
                     write(*,*)'Format_Global_aa_DMUMPS.f90 119'
                     stop
                  end if
               
                  aa_zmumps( Counter_NonZero )= GlobalMatrix( i, j )
                  ia_zmumps( Counter_NonZero )= i
                  ja_zmumps( Counter_NonZero )= J_GlobalMatrix( i, j )
               !end if
            end if
         end do
      end do
     
      if( Counter_NonZero /= Number_NonZero)then
         write(*,*)'Counter_NonZero /= Number_NonZero'
         write(*,*)'Counter_NonZero=', Counter_NonZero
         write(*,*)'Number_NonZero=', Number_NonZero
         write(*,*)'Format_Global_aa_DMUMPS.f90 113 '
         stop
      end if

      do i= 1, Number_NonZero
         if( ja_zmumps( i ) <= 0 )then
     
            write(*,*)' ja_zmumps( i ) <= 0 '
            write(*,*)' ja_zmumps( i )=', ja_zmumps( i )
            write(*,*)'Format_Global_aa_DMUMPS.f90 121 '
            stop
         end if       
         if( Number_Node < ja_zmumps( i )   )then  
    
            write(*,*)' Number_Node < ja_zmumps( i )  '
            write(*,*)' Number_Node=', Number_Node
            write(*,*)' ja_zmumps( i )=', ja_zmumps( i )
            write(*,*)' Number_NonZero=', Number_NonZero
            write(*,*)' i=', i 
            write(*,*)'Format_Global_aa_DMUMPS.f90 127 '
            stop
   
         end if
      end do

   !====================================================================================
   else
   !====================================================================================
      write(*,*)'Flag_Symmetric=', Flag_Symmetric
      write(*,*)'Format_Global_aa_DMUMPS.f90 82 '
      stop
   end if

   !====================================================================================
   write(*,*)'      end Format_Global_aa_DMUMPS'
   !====================================================================================

   return
end subroutine Format_Global_aa_DMUMPS  


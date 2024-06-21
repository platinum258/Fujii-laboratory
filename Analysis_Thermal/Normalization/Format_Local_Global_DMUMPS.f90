
subroutine Format_Local_Global_DMUMPS &
       ( LocalMatrix, Node_Element, & 
         Shape_Element, Number_Node, Number_Element, & 
         Flag_Initialization, Number_Node_InOneLow, Flag_Symmetric,  &
         !Class_Element, ID_Element_Removed, Flag_Remove_Element, &
         Class_Element, Flag_Element_Removed, Flag_Remove_Element, &
         !=====================================================
         GlobalMatrix, J_GlobalMatrix, Number_NonZero )

   !$ use omp_lib
   implicit none  
  
   integer :: i, j, k, e
   integer, intent(in) :: Number_Node, Number_Node_InOneLow, Number_Element
   integer, intent(in) :: Shape_Element ! 3 or 4
   integer, intent(in) :: Flag_Initialization ! 0:Off, 1:On  
   integer, intent(in) :: Node_Element( Shape_Element, Number_Element )
   integer, intent(in) :: Flag_Symmetric 
   integer, intent(in) :: Class_Element( Number_Element ) 
   !integer, intent(in) :: ID_Element_Removed, Flag_Remove_Element 
   integer, intent(in) :: Flag_Element_Removed( Number_Element ) 
   integer, intent(in) :: Flag_Remove_Element 
   double precision, intent(in) :: LocalMatrix( Shape_Element, Shape_Element, Number_Element )
  
   double precision, intent(out) :: GlobalMatrix( Number_Node, Number_Node_InOneLow )
   integer, intent(out) :: J_GlobalMatrix( Number_Node, Number_Node_InOneLow )
   integer, intent(out) :: Number_NonZero
  
   integer :: integer_tmp, Counter_J, Counter_J_MAX
   double precision :: complex_tmp

   double precision, parameter :: Zero= 0d0 !dcmplx( 0.0d0, 0.0d0 )
  
   integer, parameter :: Flag_On= 1
   integer, parameter :: Flag_Off= 0
   
   !==================================================================
   write(*,*)'         call Format_Local_Global_DMUMPS'
   !==================================================================
  
   if( Flag_Initialization==1 )then
      !!$omp parallel do default( none ) &
      !!$omp private( i, j ) &
      !!$omp shared( Number_Node_InOneLow, Number_Node ) & 
      !!$omp shared( GlobalMatrix, J_GlobalMatrix ) 
      do j= 1, Number_Node_InOneLow
         do i= 1, Number_Node 
            GlobalMatrix( i, j )= Zero 
            J_GlobalMatrix( i, j )= 0
         end do
      end do
   end if

   !==============================================================================================
   if( Flag_Remove_Element==Flag_Off )then
   !==============================================================================================
     
      Counter_J_MAX= 1
     
      do e= 1, Number_Element 
         do i= 1, Shape_Element 
            do j= 1, Shape_Element
               Counter_J= 1
    
             591   if( J_GlobalMatrix( Node_Element( i, e ), Counter_J )==Node_Element( j, e ) .and. & 
                   J_GlobalMatrix( Node_Element( i, e ), Counter_J )/=0 )then
             
                  GlobalMatrix( Node_Element( i, e ), Counter_J ) &
                  = GlobalMatrix( Node_Element( i, e ), Counter_J )+ LocalMatrix( j, i, e )
             
               else if( J_GlobalMatrix( Node_Element( i, e ), Counter_J )==0 )then
                  GlobalMatrix( Node_Element( i, e ), Counter_J )= LocalMatrix( j, i, e )
                  J_GlobalMatrix( Node_Element( i, e ), Counter_J )= Node_Element( j, e )
             
               else
                  Counter_J= Counter_J+1
                  if( Counter_J > Number_Node_InOneLow )then
                     write(*,*)'Counter_J=', Counter_J,'> Number_Node_InOneLow=', Number_Node_InOneLow, 'triangle'
                     write(*,*)'Class_Element =', Class_Element( e ) 
                     write(*,*)'Row Number =', Node_Element( i, e ) 
                     do k= 1, Number_Node_InOneLow
                        write(*,*)'J_GlobalMatrix( Node_Element( i, e ), k )=', J_GlobalMatrix( Node_Element( i, e ), k )
                     end do
                     call Output_Error( 'Format_Local_Global_DMUMPS', 70 )
                  end if
                  if( Counter_J > Counter_J_MAX )then
                     Counter_J_MAX= Counter_J
                  end if
                  go to 591
               end if
   
            end do
         end do
      end do
     
      Number_NonZero=0
   
      if( Flag_Symmetric==Flag_On )then  
         do j= 1, Number_Node_InOneLow
            do i= 1, Number_Node 
               if( J_GlobalMatrix( i, j )/=0 )then
                  if( i >= J_GlobalMatrix( i, j ) )then 
                     Number_NonZero= Number_NonZero +1
                  end if
               end if
            end do
         end do
      else if( Flag_Symmetric==Flag_Off )then  
         do j= 1, Number_Node_InOneLow
            do i= 1, Number_Node 
               if( J_GlobalMatrix( i, j )/=0 )then
                  !if( i >= J_GlobalMatrix( i, j ) )then 
                     Number_NonZero= Number_NonZero +1
                  !end if
               end if
            end do
         end do
      else
         call Output_Error( 'Format_Local_Global_DMUMPS', 106 )
      end if
   
      do i= 1, Number_Node 
         do j= 1, Number_Node_InOneLow-1
            do k= J+1, Number_Node_InOneLow
               if( J_GlobalMatrix( i, j )/=0 .and. J_GlobalMatrix( i, k )/=0 )then
                  if( J_GlobalMatrix( i, j ) > J_GlobalMatrix( i, k ))then
               
                     complex_tmp= GlobalMatrix( i, j )
                     GlobalMatrix( i, j )= GlobalMatrix( i, k )
                     GlobalMatrix( i, k )= complex_tmp
                 
                     integer_tmp= J_GlobalMatrix( i, j )
                     J_GlobalMatrix( i, j )= J_GlobalMatrix( i, k )
                     J_GlobalMatrix( i, k )= integer_tmp
               
                  end if
               end if
            end do
         end do
      end do
     
      !$omp parallel do default( none ) &
      !$omp private( i, j ) &
      !$omp shared( Number_Node_InOneLow, Number_Node, J_GlobalMatrix ) 
      do i= 1, Number_Node 
         do j= 1, Number_Node_InOneLow -1
            if( J_GlobalMatrix( i, j )/=0 .and. J_GlobalMatrix( i, j+1 )/=0 )then
               if( J_GlobalMatrix( i, j ) > J_GlobalMatrix( i, j+1 ) )then
                  write(*,*)'ERROR:  J_GlobalMatrix( i, j ) > J_GlobalMatrix( i, j+1 )'
                  call Output_Error( 'Format_Local_Global_DMUMPS', 113)
               end if
            end if
         end do
      end do

   !==============================================================================================
   else if( Flag_Remove_Element==Flag_On )then
   !==============================================================================================
     
      Counter_J_MAX= 1
     
      do e= 1, Number_Element 
         !if( Class_Element( e ) /= ID_Element_Removed )then
         if( Flag_Element_Removed( e ) == Flag_Off )then
            do i= 1, Shape_Element 
               do j= 1, Shape_Element
                  Counter_J= 1
       
                691   if( J_GlobalMatrix( Node_Element( i, e ), Counter_J )==Node_Element( j, e ) .and. & 
                      J_GlobalMatrix( Node_Element( i, e ), Counter_J )/=0 )then
                
                     GlobalMatrix( Node_Element( i, e ), Counter_J ) &
                     = GlobalMatrix( Node_Element( i, e ), Counter_J )+ LocalMatrix( j, i, e )
                
                  else if( J_GlobalMatrix( Node_Element( i, e ), Counter_J )==0 )then
                     GlobalMatrix( Node_Element( i, e ), Counter_J )= LocalMatrix( j, i, e )
                     J_GlobalMatrix( Node_Element( i, e ), Counter_J )= Node_Element( j, e )
                
                  else
                     Counter_J= Counter_J+1
                     if( Counter_J > Number_Node_InOneLow )then
                        write(*,*) 'Counter_J=', Counter_J, '> Number_Node_InOneLow=', Number_Node_InOneLow, 'triangle'
                        write(*,*)'Row Number =', Node_Element( i, e ) 
                        do k= 1, Number_Node_InOneLow
                           write(*,*)'J_GlobalMatrix( Node_Element( i, e ), k )=', J_GlobalMatrix( Node_Element( i, e ), k )
                        end do
                        call Output_Error( 'Format_Local_Global_DMUMPS', 179 )
                     end if
                     if( Counter_J > Counter_J_MAX )then
                        Counter_J_MAX= Counter_J
                     end if
                     go to 691
                  end if
      
               end do
            end do
         end if
      end do
     
      Number_NonZero=0
   
      if( Flag_Symmetric==Flag_On )then  
         do j= 1, Number_Node_InOneLow
            do i= 1, Number_Node 
               if( J_GlobalMatrix( i, j )/=0 )then
                  if( i >= J_GlobalMatrix( i, j ) )then 
                     Number_NonZero= Number_NonZero +1
                  end if
               end if
            end do
         end do
      else if( Flag_Symmetric==Flag_Off )then  
         do j= 1, Number_Node_InOneLow
            do i= 1, Number_Node 
               if( J_GlobalMatrix( i, j )/=0 )then
                  !if( i >= J_GlobalMatrix( i, j ) )then 
                     Number_NonZero= Number_NonZero +1
                  !end if
               end if
            end do
         end do
      else
         call Output_Error( 'Format_Local_Global_DMUMPS', 215 )
      end if
    
      do i= 1, Number_Node 
         do j= 1, Number_Node_InOneLow-1
            do k= J+1, Number_Node_InOneLow
               if( J_GlobalMatrix( i, j )/=0 .and. J_GlobalMatrix( i, k )/=0 )then
                  if( J_GlobalMatrix( i, j ) > J_GlobalMatrix( i, k ))then
               
                     complex_tmp= GlobalMatrix( i, j )
                     GlobalMatrix( i, j )= GlobalMatrix( i, k )
                     GlobalMatrix( i, k )= complex_tmp
                 
                     integer_tmp= J_GlobalMatrix( i, j )
                     J_GlobalMatrix( i, j )= J_GlobalMatrix( i, k )
                     J_GlobalMatrix( i, k )= integer_tmp
               
                  end if
               end if
            end do
         end do
      end do
     
      !$omp parallel do default( none ) &
      !$omp private( i, j ) &
      !$omp shared( Number_Node_InOneLow, Number_Node, J_GlobalMatrix ) 
      do i= 1, Number_Node 
         do j= 1, Number_Node_InOneLow -1
            if( J_GlobalMatrix( i, j )/=0 .and. J_GlobalMatrix( i, j+1 )/=0 )then
               if( J_GlobalMatrix( i, j ) > J_GlobalMatrix( i, j+1 ) )then
                  write(*,*)'ERROR:  J_GlobalMatrix( i, j ) > J_GlobalMatrix( i, j+1 )'
                  call Output_Error( 'Format_Local_Global_DMUMPS', 246 )
               end if
            end if
         end do
      end do
   !==============================================================================================
   else
   !==============================================================================================
      write(*,*)'Flag_Remove_Element=', Flag_Remove_Element
      call Output_Error( 'Format_Local_Global_DMUMPS', 253 )
   end if


  
   return
end subroutine Format_Local_Global_DMUMPS 


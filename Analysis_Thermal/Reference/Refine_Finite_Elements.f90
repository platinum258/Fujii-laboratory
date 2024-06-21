
subroutine Refine_Finite_Elements&
           ( Class_Element_Refine, &
             Number_Node, Position_Node, &
             Number_Element_Triangle, Index_Element_2_Node_Triangle, Class_Element_Triangle, &
             Width_Matrix_LHS )
  
   implicit none
   
   integer, intent(in) :: Class_Element_Refine
   integer, intent(inout) :: Number_Node, Number_Element_Triangle, Width_Matrix_LHS
   integer, intent(inout) :: Index_Element_2_Node_Triangle( 3, Number_Element_Triangle*3 )
   integer, intent(inout) :: Class_Element_Triangle( Number_Element_Triangle*3 )
   double precision, intent(inout) :: Position_Node( 2, Number_Node+Number_Element_Triangle )

   integer e, i, j, Counter, Number_Element_Refine
   double precision, allocatable, dimension(:,:,:) :: Position_Node_Element
   double precision, allocatable, dimension(:,:) :: Position_Center_Element
   integer, allocatable, dimension(:,:) :: Index_Element_2_Node_Refine
   integer, allocatable, dimension(:) :: Node_Number_Center_Element, Element_Number_Original

   !===========================================
   write(*,*)'   call Refine_Finite_Elements'
   !===========================================

   if( Class_Element_Refine==999 )then
      Counter= Number_Element_Triangle
   else
      Counter= 0
      do e= 1, Number_Element_Triangle 
         if( Class_Element_Triangle( e )==Class_Element_Refine )then
            Counter= Counter +1
         end if
      end do 
   end if

   if( Counter==0 ) goto 94 

   Number_Element_Refine= Counter

   allocate( Index_Element_2_Node_Refine( 3, Number_Element_Refine ) )
   allocate( Element_Number_Original( Number_Element_Refine ) )

   Counter= 0
   do e= 1, Number_Element_Triangle 
      if( Class_Element_Triangle( e )==Class_Element_Refine .or. Class_Element_Refine==999 )then
         Counter= Counter +1
         Element_Number_Original( Counter )= e
         do i= 1, 3
            Index_Element_2_Node_Refine( i, Counter )= Index_Element_2_Node_Triangle( i, e )
         end do 
      end if
   end do 

   allocate( Position_Node_Element( 2, 3, Number_Element_Refine ) )

   do e= 1, Number_Element_Refine
      do i= 1, 3
         do j= 1, 2
            Position_Node_Element( j, i, e )= Position_Node( j, Index_Element_2_Node_Refine( i, e ) )
         end do
      end do
   end do

   allocate( Position_Center_Element( 2, Number_Element_Refine ) )
   allocate( Node_Number_Center_Element( Number_Element_Refine ) )

   do e= 1, Number_Element_Refine
      Node_Number_Center_Element( e )= Number_Node +e 
      do j= 1, 2
         Position_Center_Element( j, e ) &
         = ( Position_Node_Element( j, 1, e ) +Position_Node_Element( j, 2, e ) +Position_Node_Element( j, 3, e ) )/3.0d0  
      end do
   end do

   do e= 1, Number_Element_Refine
      if( Class_Element_Refine==999 )then
         Class_Element_Triangle( 2*e +Number_Element_Triangle -1 )= Class_Element_Triangle( Element_Number_Original( e ) ) 
         Class_Element_Triangle( 2*e +Number_Element_Triangle )= Class_Element_Triangle( Element_Number_Original( e ) )
      else
         Class_Element_Triangle( 2*e +Number_Element_Triangle -1 )= Class_Element_Refine
         Class_Element_Triangle( 2*e +Number_Element_Triangle )= Class_Element_Refine
      end if
      Index_Element_2_Node_Triangle( 3, Element_Number_Original( e ) )= Node_Number_Center_Element( e ) 

      Index_Element_2_Node_Triangle( 1, 2*e +Number_Element_Triangle -1 )&
      = Index_Element_2_Node_Refine( 2, e ) 
      Index_Element_2_Node_Triangle( 2, 2*e +Number_Element_Triangle -1 )&
      = Index_Element_2_Node_Refine( 3, e ) 
      Index_Element_2_Node_Triangle( 3, 2*e +Number_Element_Triangle -1 )&
      = Node_Number_Center_Element( e ) 

      Index_Element_2_Node_Triangle( 1, 2*e +Number_Element_Triangle )&
      = Index_Element_2_Node_Refine( 3, e ) 
      Index_Element_2_Node_Triangle( 2, 2*e +Number_Element_Triangle )&
      = Index_Element_2_Node_Refine( 1, e ) 
      Index_Element_2_Node_Triangle( 3, 2*e +Number_Element_Triangle )&
      = Node_Number_Center_Element( e ) 
   end do

   do e= 1, Number_Element_Refine 
      do j= 1, 2
         Position_Node( j, e +Number_Node )= Position_Center_Element( j, e )
      end do
   end do

   Number_Node= Number_Node +Number_Element_Refine
   Number_Element_Triangle= Number_Element_Triangle +Number_Element_Refine*2

   Width_Matrix_LHS= Width_Matrix_LHS*2

   deallocate( Element_Number_Original )

   deallocate( Position_Center_Element )
   deallocate( Node_Number_Center_Element )

   deallocate( Position_Node_Element )
   deallocate( Index_Element_2_Node_Refine )

   !===========================================
   write(*,*)'   end Refine_Finite_Elements'
   !===========================================

   !stop

   94 continue
   return
end subroutine Refine_Finite_Elements


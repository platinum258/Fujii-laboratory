
program thermalcloak_ITRfree
   !$ use omp_lib
   implicit none
   
!============================================================================================   
     
   integer :: num_threads, SizeA, Node_Omega_out, Node_Omega_in
     
   double precision,allocatable :: rko( :,: ), rhso( : ), ipmkl( : )
   integer,allocatable :: rko_id( : ), T_rhs_ID( : ), Omega_out_ID( : ), Omega_in_ID( : )
   
   integer,allocatable :: Element_Design_Area_T( : ), Element_Design_Area_S( : )
   double precision, allocatable :: x_beta( : ), density_updated( : ), diff_density( : )
   double precision, allocatable :: Triangle_density_vec( : ), Best_Triangle_density_vec( : )
   double precision, allocatable :: Square_density_vec( : ) 
	
!============================================================================================   
   
   character filename*128
   integer,parameter :: nelx = 300   !全体領域のx方向の分割数
   integer,parameter :: nely = 200   !全体領域のy方向の分割数
   
   integer,parameter :: nd = 3
   double precision, parameter :: vol_constraint = 0.4d0 !初期値
   integer,parameter :: penal = 3
 
   integer,parameter :: Number_Node = ( nely+1 )*( nelx+1 )
   integer,parameter :: Square_Number_Element   = nely*nelx              
   integer,parameter :: Triangle_Number_Element = nelx*nely*2

   integer,parameter :: g_limit = 1000
   integer g_min
   
   double precision,parameter :: move_limit = 1.0d-3 !設計変数の増減値
 
   double precision,parameter :: high = 1.0d0       !クロークの高さ
   double precision,parameter :: thickness = 1.0d0  !平板の厚さ
   double precision,parameter :: t_boundary_higher_side=1.0d0 
   double precision,parameter :: t_boundary_lower_side=0.0d0          
   
   integer :: counter( Number_Node ), counter1( Number_Node )
   double precision :: xcoord( Number_Node ), ycoord( Number_Node ), nesc( Square_Number_Element,2 )
   double precision :: Position_Node( 2, Number_Node )
   double precision :: xmin = - nelx/2.0d0 , xmax = nelx/2.0d0 
   double precision :: ymin = - nely/2.0d0 , ymax = nely/2.0d0
   double precision :: Radius_2Area = (nely/8.0d0)*1.5d0 + 0.5d0  !断熱部分の半径
   double precision :: Radius_4Area = (nely/8.0d0)*3.0d0 + 0.5d0  !設計領域の半径

   integer,parameter :: Type_Element_Insulator = 2
   integer,parameter :: Type_Element_DesignDomain = 4
   integer,parameter :: Type_Element_Outer = 5
   integer,parameter :: Width_Matrix_LHS=7
   
   integer,parameter :: Loop_Solution = 1
   integer,parameter :: Loop_MatParam_T = 1
   integer,parameter :: Loop_Obj_Func = 1
   double precision,parameter :: Perimeter_Structure = 0.0d0
!============================================================================================  
   
   double precision diff_OF( 2:g_limit )

!============================================================================================  
   integer :: Number_EDesignArea  !断熱部分を除いた三角形要素数
   integer :: Number_Eins
   integer :: Number_SEDesignArea
   integer :: Number_SEins
   integer :: nodex_s( Square_Number_Element,4 ), nodex( Triangle_Number_Element,3 ), eid( Triangle_Number_Element )
   integer :: Index_Element_2_Node( 3, Triangle_Number_Element ) 
   double precision :: dis
             
!============================================================================================
   integer :: g, Best_generation
   integer :: i,j,k,l,m,n,e,t,ie,e1,e2,info,node,int_tmp
   integer :: position_in_element( Triangle_Number_Element, 3, 2 )
   
   double precision :: l1,l2,lmid
   double precision :: change
   double precision :: PD_Objective_function( Triangle_Number_Element )
   double precision :: beta( Triangle_Number_Element )
   double precision :: Objective_function( g_limit )
   double precision :: x( nd ), y( nd ), A( nd ), B( nd )
   double precision :: Area
   double precision :: Local_matrix( nd,nd,Triangle_Number_Element )
   double precision :: shape_function( Triangle_Number_Element, 3 )

   double precision :: psi_n,psi,volume,diff_T
   double precision :: nec( Triangle_Number_Element,2 )
   double precision :: kappa( Triangle_Number_Element )       !設計変数から等価変換した熱伝導率
   double precision :: T_ref( Number_Node ), T_bare( Number_Node )
   double precision :: T_rhs( Number_Node ), Global_matrix( Number_Node,Number_Node )
   double precision :: epsi_n( Number_Node ), epsi( Number_Node )
   double precision :: T_e( Triangle_Number_Element ), T_mat( g_limit, Number_Node )   
   double precision :: T_n( Number_Node ), Td( Number_Node )
   double precision :: E_str( Triangle_Number_Element )
  
  !============================================================================================

   double precision :: Le( Triangle_Number_Element, nd )
   double precision :: Lamda( Triangle_Number_Element, nd )
   double precision :: PD_Lamda( Triangle_Number_Element, 2 )
   double precision :: PD_T( Triangle_Number_Element, 2 )
      
   double precision :: lm( nd ), gm( Number_Node )
   double precision :: w_n( Number_Node ), Q( Number_Node ), Q1( Number_Node )
   double precision :: PD_W( Triangle_Number_Element,2 )
      
   !===========================================================================================
   !Set number threads  
   !===========================================================================================         
   !$call omp_set_num_threads(28)
   !call mkl_set_num_threads(2)
      
   !=====================================================
   write(*,*)'**Create Finite Element Data**'
   !=====================================================
   
   do i = 1, nelx + 1
      do j = 1, nely + 1
         k = j + ( nely + 1 )*( i - 1 )
         xcoord( k ) = xmin + dble( i - 1 )
      end do
   end do

   do i = 1, nelx + 1
      do j = 1, nely + 1
         k = j + ( nely + 1 )*( i - 1 )
         ycoord( k ) = ymin + dble( j - 1 )
      end do
   end do

   do i = 1, Number_Node
      Position_Node( 1,i ) = xcoord( i )
	  Position_Node( 2,i ) = ycoord( i )
   enddo	  

   do i = 0, nelx - 1
      do j = 1, nely
         m = j + i *nely
         nodex_s( m,1 ) = j + i*( nely + 1 )
         nodex_s( m,2 ) = j + i*( nely + 1 ) + ( nely + 1 )
         nodex_s( m,3 ) = j + i*( nely + 1 ) + ( nely + 2 )
         nodex_s( m,4 ) = j + i*( nely + 1 ) + 1
      end do
   end do

   do i = 1, Square_Number_Element
     nodex( i*2,1 ) = nodex_s( i,4 )
     nodex( i*2,2 ) = nodex_s( i,2 )
     nodex( i*2,3 ) = nodex_s( i,3 )
     nodex( i*2-1,1 ) = nodex_s( i,4 )
     nodex( i*2-1,2 ) = nodex_s( i,1 )
     nodex( i*2-1,3 ) = nodex_s( i,2 )
   end do

   do i = 1, nelx 
      do j = 1, nely
         n = j + nely *( i - 1 )
         nesc( n,1 ) = i + ( xmin - 0.5d0 )
         nesc( n,2 ) = j + ( ymin - 0.5d0 )
      end do
   end do

   do i = 1, Triangle_Number_Element
      eid( i ) = Type_Element_Outer
   end do

   do i = 1, Square_Number_Element
      dis = sqrt( nesc( i,1 )*nesc( i,1 ) + nesc( i,2 )*nesc( i,2 ) )
      if( dis <= Radius_2Area ) then
         eid( i*2 ) = Type_Element_Insulator
         eid( i*2 - 1 ) = Type_Element_Insulator
      else if( dis <= Radius_4Area .and.  dis > Radius_2Area ) then
         eid( i*2 ) = Type_Element_DesignDomain
         eid( i*2 - 1 ) = Type_Element_DesignDomain
      end if
   end do

   Number_EDesignArea = 0
   do i = 1, Triangle_Number_Element
      if( eid( i ) == Type_Element_DesignDomain ) then
         Number_EDesignArea = Number_EDesignArea + 1
      end if 
   end do 
   Number_SEDesignArea = Number_EDesignArea/2

   Number_Eins = 0
   do i = 1, Triangle_Number_Element
      if( eid( i ) == Type_Element_Insulator ) then
         Number_Eins = Number_Eins + 1
      end if
   end do
   Number_SEins = Number_Eins/2

   allocate( Element_Design_Area_T( Number_EDesignArea ))
   allocate( Element_Design_Area_S( Number_EDesignArea/2 ))
   allocate( x_beta( Number_EDesignArea ))
   allocate( density_updated( Number_EDesignArea ), diff_density( Number_EDesignArea ))
   allocate( Triangle_density_vec( Number_EDesignArea ))
   allocate( Best_Triangle_density_vec( Number_EDesignArea ))
   allocate( Square_density_vec( Number_EDesignArea/2 ))
   
   !=====================================================
   write(*,*)'Initialization'
   !=====================================================
   
   !Reset Temperature
   T_e( : ) = 0.0d0
   
   !Thermal conductivity value
   kappa( : ) = 1.0d0
   
   SizeA = Number_Node

   Index_Element_2_Node = TRANSPOSE( nodex )
   
   allocate( rko( SizeA,SizeA ), rhso( SizeA ), ipmkl( SizeA ))
   allocate( rko_id( SizeA ) , T_rhs_ID( SizeA ))
   allocate( Omega_out_ID( SizeA ), Omega_in_ID( SizeA ))
   
   !========================================================================================
   write(*,*)'Set any point in the Triangle element'
   !========================================================================================     

   do i = 1, nelx 
     do j = 1, nely
       e2 = ( j + nely *( i - 1 ) )*2 
       e1 = e2 - 1
       nec( e2,1 ) = i + ( xmin - 0.25d0 ) 
       nec( e2,2 ) = j + ( ymin - 0.25d0 )
       nec( e1,1 ) = i + ( xmin - 0.75d0 )
       nec( e1,2 ) = j + ( ymin - 0.75d0 )
     end do
   end do
   
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Reading Temperature without cloak
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   
		!======================================================== 
		write(*,*)'Making Local Matrix'
		!========================================================
		call mkGSM_0 (nd,Triangle_Number_Element,Number_Node,x,y,A,B,nodex,xcoord,ycoord,Local_matrix )
		!======================================================== 
		!write(*,*)'Making Global Matrix and set boundary conditions'
		!========================================================
		!call setBC_0( Triangle_Number_Element,Number_Node,nd,SizeA,xcoord,ycoord,Global_matrix,&
		!				Local_matrix,nodex,T_rhs,rko_id,rko,rhso,xmin,xmax,ymin,ymax )
		!===================================================================================
		!write(*,*)'Solving Linear Equation' 
		!===================================================================================
		!call dgesv( SizeA, 1, rko, SizeA, ipmkl, rhso, SizeA, info )
      !===================================================================================
       call analyze_thermal_distribution(Local_matrix, Index_Element_2_Node, eid, &
              0, Number_Node, Triangle_Number_Element, Width_Matrix_LHS, &
              xcoord, ycoord, t_boundary_higher_side, t_boundary_lower_side, &
              xmax,ymax,xmin,ymin, &
              !===========================================================================
              T_ref)
   
		
	!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~	
	 
	! do i = 1, SizeA
   !     T_rhs( rko_id( i ) ) = rhso( i )
   !  end do
!
   !  if( info/=0 )then
   !     write(*,*)'info of dgesv:', info, 'at 154'
   !     stop
   !  end if
   !  !======================================================== 
   !  write(*,*)'Converts to temperature for the element'
   !  !========================================================
!
   !do i = 1, SizeA      
   !   T_ref( i ) = T_rhs( i )
   !enddo  
   
   open( 200, file = 'T_ref.txt' )
      do i = 1, SizeA
         write( 200,* ) xcoord( i ), ycoord( i ),  T_ref( i )
      end do
   close( 200 )
   
   !initialization
   !=======================================
   !Compute SizeA
   !=======================================
 
   counter( : ) = 1 
   do i = 1, Triangle_Number_Element
      if( eid( i ) == Type_Element_Insulator ) then
         do j = 1, nd
            counter( nodex( i,j ) ) = 0
         enddo
      endif 
   enddo
   do i = 1, Triangle_Number_Element   
      if( eid( i ) /= Type_Element_Insulator ) then
         do j = 1, nd
            counter( nodex( i,j ) ) = 1
         enddo
      endif
   enddo 
   
   j = 1
   do i = 1, Number_Node
     if( counter( i ) == 1 ) then
        T_rhs_ID( j ) = i
        j = j + 1
     end if
   end do

   SizeA = 0
   do i = 1, Number_Node
      if ( counter( i ) == 1 ) then
         SizeA = SizeA + 1
      endif
   enddo
   
   !=======================================
   !Compute Node_Omega_out
   !=======================================

   counter1( : ) = 1 
   do i = 1, Triangle_Number_Element
      if( eid( i ) /= Type_Element_Outer ) then
         do j = 1, nd
            counter1( nodex( i,j ) ) = 0
         enddo
      endif 
   enddo
   do i = 1, Triangle_Number_Element   
      if( eid( i ) == Type_Element_Outer ) then
         do j = 1, nd
            counter1( nodex( i,j ) ) = 1
         enddo
      endif
   enddo 
   
   j = 1
   do i = 1, Number_Node
     if( counter1( i ) == 1 ) then
        Omega_out_ID( j ) = i
        j = j + 1
     end if
   end do
   
   j = 1
   do i = 1, Number_Node
     if( counter1( i ) == 0 ) then
        Omega_in_ID( j ) = i
        j = j + 1
     end if
   end do

   Node_Omega_out = 0
   do i = 1, Number_Node
      if ( counter1( i ) == 1 ) then
         Node_Omega_out = Node_Omega_out + 1
      endif
   enddo
   
   Node_Omega_in = 0
   do i = 1, Number_Node
      if ( counter1( i ) == 0 ) then
         Node_Omega_in = Node_Omega_in + 1
      endif
   enddo

   write(*,*)'Node_Omega_out=',Node_Omega_out
   write(*,*)'Node_Omega_in=', Node_Omega_in   
   
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Reading temperature before cloaking
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      Local_matrix=0.0d0
	  
		call mkGSM(nd,Type_Element_Insulator,Triangle_Number_Element,&
						Number_Node,kappa,eid,x,y,A,B,nodex,xcoord,ycoord,Local_matrix)
		!call setBC( Triangle_Number_Element,Type_Element_Insulator,Number_Node,nd,SizeA,xcoord,ycoord,nodex,&
		!				Local_matrix,Global_matrix,eid,T_rhs,rko_id,rko,rhso,xmin,xmax,ymin,ymax )
		!call dgesv( SizeA, 1, rko, SizeA, ipmkl, rhso, SizeA, info )
      call analyze_thermal_distribution(Local_matrix, Index_Element_2_Node, eid, &
         1, Number_Node, Triangle_Number_Element, Width_Matrix_LHS, &
         xcoord, ycoord, t_boundary_higher_side, t_boundary_lower_side, &
         xmax,ymax,xmin,ymin, &
         !===========================================================================
         T_bare)      
   
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   
	! do i = 1, SizeA
   !     T_rhs( rko_id( i ) ) = rhso( i )
   !  end do
   !
   !  if( info/=0 )then
   !     write(*,*)'info of dgesv:', info, 'at 185'
   !     stop
   !  end if

   open( 205, file = 'T_bare.txt' )
      do i = 1, Number_Node
         write( 205,* ) xcoord( i ), ycoord( i ),  T_bare( i )
      end do
   close( 205 )

   !=====================================================
   write(*,*)'**Read Complete**'
   !=====================================================

   !========================================================================================
   !Set element design Area
   !========================================================================================     
  
   j = 1
   do i = 1, Triangle_Number_Element
     if( eid( i ) == Type_Element_DesignDomain ) then
        Element_Design_Area_T( j ) = i
        j = j + 1
     end if
   end do
   
   j = 1
   do i = 1, Triangle_Number_Element, 2
     if( eid( i ) == Type_Element_DesignDomain ) then
        Element_Design_Area_S( j ) = i
        j = j + 1
     end if
   end do   
   
   !========================================================================================
   !Calculate psi_n
   !========================================================================================       
   
   !T_bare( : ) = T_rhs( : )
   
   do i = 1, Number_Node
      if ( counter1( i ) == 1 ) then
         epsi_n( i ) = T_bare( i ) - T_ref( i )
      else
	     epsi_n( i ) = 0.0d0
      end if   	  
   end do	  
  
   call Compute_Objective_Function &
      ( Type_Element_Outer, &
        epsi_n, Number_Node, & 
        Position_Node, Number_Node, &
        Index_Element_2_Node, eid, Triangle_Number_Element, &
        g, Loop_Solution, &
        Perimeter_Structure, &
   !  =================================================================
        psi_n )

   write(*,*) 'psi_n=',psi_n

   !===========================================================================================
   write(*,*)'Initialization'
   !===========================================================================================
  
   !density value
   Triangle_density_vec( : ) = vol_constraint
   Square_density_vec( : ) = vol_constraint
   
   !psi value
   epsi( : ) = 0.0d0
   
   !parameter
   change = 1.0d0 
   
   !===========================================================================================
   write(*,*)'Optimization Computation'
   !===========================================================================================

   g= 0
   g_min= 1
   do while( change > 1.0d-4 )
      g = g + 1

      write(*,*)'=============================================================='
      write(*,*)'g=', g, 'change=', change
      write(*,*)'=============================================================='

      !========================================================
      write(*,*) 'Set kappa'
      !========================================================     

      do i = 1, Number_EDesignArea
         j = Element_Design_Area_T( i )
         kappa( j ) = 1.0d0 + ( high/thickness ) *Triangle_density_vec( i )**penal 
      end do
    
      !========================================================
      write(*,*) 'Finite Element Analysis'
      !========================================================     
	  
	  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ! Reading temperature after cloaking
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~	  
   
			call mkGSM(nd,Type_Element_Insulator,Triangle_Number_Element,&
						Number_Node,kappa,eid,x,y,A,B,nodex,xcoord,ycoord,Local_matrix) 
			!call setBC( Triangle_Number_Element,Type_Element_Insulator,Number_Node,nd,SizeA,xcoord,ycoord,nodex,&
			!			Local_matrix,Global_matrix,eid,T_rhs,rko_id,rko,rhso,xmin,xmax,ymin,ymax )
			!call dgesv( SizeA, 1, rko, SizeA, ipmkl, rhso, SizeA, info )
         call analyze_thermal_distribution(Local_matrix, Index_Element_2_Node, eid, &
            1, Number_Node, Triangle_Number_Element, Width_Matrix_LHS, &
            xcoord, ycoord, t_boundary_higher_side, t_boundary_lower_side, &
            xmax,ymax,xmin,ymin, &
            !===========================================================================
            T_rhs)      
      
	  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	  
!	  if( info/=0 )then
!         write(*,*)'info of dgesv:', info, 'at 295'
!         stop
!      end if
!	  
!   	 do i = 1, SizeA
!        T_rhs( rko_id( i ) ) = rhso( i )
!     end do 
  
      !=========================================================================================
      write(*,*)'Compute Objective Function and Volume'
      !=========================================================================================

	  T_n( : ) = T_rhs( : )

	  do i = 1, Number_Node
         if( counter( i ) == 1 ) then
            T_mat( g,i ) = T_n( i )
         end if
      end do
	 	  		  
	  do i = 1, Number_Node
         if ( counter1( i ) == 1 ) then
            epsi( i ) = T_n( i ) - T_ref( i )
         else
	        epsi( i ) = 0.0d0
         end if   	  
      end do	   
	  	  
      call Compute_Objective_Function &
         ( Type_Element_Outer, &
           epsi, Number_Node, & 
		   Position_Node, Number_Node, &
		   Index_Element_2_Node, eid, Triangle_Number_Element, &
		   g, Loop_Solution, &
           Perimeter_Structure, &
   !  =================================================================
           psi )
		   
      write(*,*) 'psi=',psi		   

      Objective_function( g ) = psi / psi_n
      volume = sum( Square_density_vec ) / dble( Number_EDesignArea/2 )

      write(*,*)'===================================================='
      write(*,*)'   Objective_Function=', Objective_function( g )
      write(*,*)'   Volume=', volume
      write(*,*)'===================================================='

      write(*,*) g, Objective_function( g ), volume
      open( 500, file='ObjectiveFunction.dat', position='append' )
         write(500,*) g, Objective_function( g ), volume
      close(500)
	  
	  !=========================================================================================
      write(*,*)'Generation of temperature distribution per element'
      !=========================================================================================
	  
	  do i = 1, SizeA
	     T_rhs( rko_id( i ) ) = rhso( i )
      end do   
	  	  
	  call etemp( Triangle_Number_Element,Number_Node,Type_Element_Insulator,nd,SizeA,Area,eid,&
                      xcoord,ycoord,nec,x,y,rhso,rko_id,nodex,shape_function,T_rhs,T_e )   

      open(400,file='etemperature.txt',status='replace')
         do i = 1, Triangle_Number_Element
            if( eid( i ) /= Type_Element_Insulator ) then
               write(400,*) nec( i,1 ), nec( i,2 ), T_e( i )
            end if
         end do
      close(400)

      !Preparing for optimization

      !========================================================================================
      write(*,*)'Calculation of the accompanying field'
	  !========================================================================================      
    
  !    call mkGSM(nd,Type_Element_Insulator,Triangle_Number_Element,&
!						Number_Node,kappa,eid,x,y,A,B,nodex,xcoord,ycoord,Local_matrix)
			!call setBC_sv( Triangle_Number_Element,Number_Node,Type_Element_Insulator,&
			!				nd,SizeA,xcoord,ycoord,Local_matrix,Global_matrix,T_rhs,&
			!				rko_id,rko,rhso,xmin,xmax,ymin,ymax,ie,eid,counter,x,y,&
			!				nodex,Area,a,b,psi_n,T_n,T_ref,lm,gm,Q,Q1 )
       call analyze_thermal_distribution(Local_matrix, Index_Element_2_Node, eid, &
          1, Number_Node, Triangle_Number_Element, Width_Matrix_LHS, &
          xcoord, ycoord, t_boundary_higher_side, t_boundary_lower_side, &
          xmax,ymax,xmin,ymin, &
          !===========================================================================
          w_n)      
	
      !========================================================================================
		
	  !if( info/=0 )then
     !    write(*,*)'info of dgesv:', info, 'at 999'
     !    stop
     ! end if
	  
!	  w_n( : ) = 0.0d0	  
!   	  do i = 1, SizeA
!         w_n( rko_id( i ) ) = rhso( i )
!      end do 

      open( 500, file = 'participating_venue.txt' )
         do i = 1, Number_Node
            write( 500,* ) xcoord( i ), ycoord( i ), W_n( i )
         end do
      close( 500 )	  
	  
      !========================================================================================
      write(*,*)'Calculating Sensitivity Functions with Objective Functions'
	  !========================================================================================  	  

      do i = 1, Triangle_Number_Element
	     if( eid( i ) /= Type_Element_Insulator ) then 
		 
		    do j = 1, nd
               x( j ) = xcoord( nodex( i,j ) )
               y( j ) = ycoord( nodex( i,j ) )
            end do
          
            B(1) = y(2) - y(3)
            B(2) = y(3) - y(1)
            B(3) = y(1) - y(2)
            A(1) = x(3) - x(2)
            A(2) = x(1) - x(3)
            A(3) = x(2) - x(1)          
		 
			PD_W( i,1 ) = B(1)*w_n( nodex( i,1 )) + B(2)*w_n( nodex( i,2 )) + B(3)*w_n( nodex( i,3 ))
			PD_W( i,2 ) = A(1)*w_n( nodex( i,1 )) + A(2)*w_n( nodex( i,2 )) + A(3)*w_n( nodex( i,3 ))
			PD_T( i,1 ) = B(1)*T_n( nodex( i,1 )) + B(2)*T_n( nodex( i,2 )) + B(3)*T_n( nodex( i,3 ))
			PD_T( i,2 ) = A(1)*T_n( nodex( i,1 )) + A(2)*T_n( nodex( i,2 )) + A(3)*T_n( nodex( i,3 ))
 
            PD_Objective_function( i ) = high*( PD_W( i,1 )*PD_T( i,1 ) + PD_W( i,2 )*PD_T( i,2 ) ) / thickness
         end if
      end do	
	  
      !=========================================================================================
      write(*,*)'Update Density Field'
      !=========================================================================================

      l1 = 0.0d0
      l2 = 1.0d9
 
      do while( ( l2-l1 )/( l1+l2 ) > 1.0d-2 )

        lmid= 0.5d0*( l2+l1 )

        do i = 1, Triangle_Number_Element, 2
           if( eid( i ) /= Type_Element_Insulator ) then
              beta( i ) = ( ( PD_Objective_function( i ) ) / lmid )**7.5d-1 
           end if 
        end do    		

        do i = 1, Number_EDesignArea/2
           j = Element_Design_Area_S( i )
           x_beta( i ) = Square_density_vec( i ) *beta( j )
        end do 
 
        do i = 1, Number_EDesignArea/2
           if(      x_beta( i )  <= max( 0.0d0, Square_density_vec( i ) -move_limit ) )then
              density_updated( i )= max( 0.0d0, Square_density_vec( i ) -move_limit )
           else if( x_beta( i )  >= min( 1.0d0, Square_density_vec( i ) +move_limit ) )then
              density_updated( i )= min( 1.0d0, Square_density_vec( i ) +move_limit )
           else
              density_updated( i )= x_beta( i )
           end if
        end do

        volume = sum( density_updated )/dble( Number_EDesignArea/2 )
		
        if( volume > vol_constraint )then
           l1= lmid
        else
           l2= lmid
        end if
    
      end do !do while((l2-l1)/(l1+l2)>1.0d-3)

      !=========================================================================================
      write(*,*)'Check Convergence' 
      !=========================================================================================

      diff_density( : ) = abs( density_updated( : ) - Square_density_vec( : ) )
      change = maxval( diff_density )
      write(*,*)'   change=', change
	  
      !=========================================================================================
      !Make .ps
      !=========================================================================================

      j = 1 
      do i = 1 , Triangle_Number_Element
         if ( eid( i ) == Type_Element_DesignDomain ) then
            E_str( i ) = Triangle_density_vec( j )
            j = j + 1    
         end if
      end do
	        
      do e = 1, Triangle_Number_Element
         do i = 1, 3
            position_in_element( e,i,1 ) = xcoord( nodex( e,i ) )
            position_in_element( e,i,2 ) = ycoord( nodex( e,i ) )
         end do
      end do

      write(filename,'(i5.5,"_Structure.ps")') g
      !=========================================================================================
      write(*,*)'Output .ps >>>>> ', filename
      !=========================================================================================  
	  
      open(999,file=filename,status='replace')
         write(999,*)'%!PS-Adobe-3.0 EPSF-3.0'
         write(999,*)'%%BoundingBox:',0,0,600,400
         write(999,*)'0.0 setgray'
         do e = 1, Triangle_Number_Element
           if ( eid( e ) == Type_Element_DesignDomain ) then
              write(999,*)'newpath'
              write(999,*) 5.0d0*position_in_element(e,1,1) + 300.0d0, 5.0d0*position_in_element(e,1,2) + 200.0d0, 'moveto'
              write(999,*) 5.0d0*position_in_element(e,2,1) + 300.0d0, 5.0d0*position_in_element(e,2,2) + 200.0d0, 'lineto'
              write(999,*) 5.0d0*position_in_element(e,3,1) + 300.0d0, 5.0d0*position_in_element(e,3,2) + 200.0d0, 'lineto'
              write(999,*) 'closepath'
              write(999,*) 1.00 - E_str(e),'setgray'
              write(999,*) 'fill'
           else if ( eid( e ) == Type_Element_Insulator ) then 
              write(999,*)'newpath'
              write(999,*) 5.0d0*position_in_element(e,1,1) + 300.0d0, 5.0d0*position_in_element(e,1,2) + 200.0d0, 'moveto'
              write(999,*) 5.0d0*position_in_element(e,2,1) + 300.0d0, 5.0d0*position_in_element(e,2,2) + 200.0d0, 'lineto'
              write(999,*) 5.0d0*position_in_element(e,3,1) + 300.0d0, 5.0d0*position_in_element(e,3,2) + 200.0d0, 'lineto'
              write(999,*) 'closepath'
              write(999,*) '1.0 0.0 0.0 setrgbcolor'
              write(999,*) 'fill'
            end if
         end do    
      close(999)
	  
      !=========================================================================================
      !convergence judgment
      !=========================================================================================

      if( g >= 2 )then
         diff_OF( g )= abs( Objective_function( g ) - Objective_function( g-1 ) )/Objective_function( g )
         write(*,*)'diff_OF( g )=', diff_OF( g ) 

         if( Objective_function( g ) < Objective_function( g_min ) ) then
			g_min = g
			do i = 1, Number_EDesignArea
		       Best_Triangle_density_vec( i ) = Triangle_density_vec( i )
		    end do
		 end if

         open( 574, file='Convergence.dat', position='append' )
            write( 574, * ) g, diff_OF( g ), Objective_function( g )
         close( 574 )

         if( diff_OF( g ) <= 1.0d-4 )then
            write(*,*)'Convergence : diff_OF( g ) < 0.0001'
            write(*,*)'g_min=', g_min
            call Output_the_best_value( g_limit,g_min,Number_EDesignArea,Best_Triangle_density_vec,&
										Number_Node,counter,xcoord,ycoord,T_mat )
			call Compute_difference_in_temperature( g_limit,g_min,Number_Node,diff_T,counter1,&
													xcoord,ycoord,T_ref,T_mat,Td)
            call Make_Analysis_date( i,j,nd,Triangle_Number_Element,Number_Node,Number_EDesignArea&
						,nodex,eid,xmin,xmax,ymin,ymax,xcoord,ycoord&
						,Best_Triangle_density_vec,high,thickness,penal,g_min)
			change = 0.0d0 !stop
         else if( g > g_min +10 )then
            write(*,*)'Convergence : g > g_min +10'
            write(*,*)'g_min=', g_min
            call Output_the_best_value( g_limit,g_min,Number_EDesignArea,Best_Triangle_density_vec,&
										Number_Node,counter,xcoord,ycoord,T_mat )
			call Compute_difference_in_temperature( g_limit,g_min,Number_Node,diff_T,counter1,&
													xcoord,ycoord,T_ref,T_mat,Td)										
            call Make_Analysis_date( i,j,nd,Triangle_Number_Element,Number_Node,Number_EDesignArea&
						,nodex,eid,xmin,xmax,ymin,ymax,xcoord,ycoord&
						,Best_Triangle_density_vec,high,thickness,penal,g_min)
			change = 0.0d0 !stop
         end if 
      end if 

      !=========================================================================================
      write(*,*)'Update Design variable'
      !=========================================================================================	  
 
      Square_density_vec( : ) = density_updated( : )
	  
	  do j = 1, Number_EDesignArea/2
	     i = j*2
         Triangle_density_vec( i )   = density_updated( j ) 
		 Triangle_density_vec( i-1 ) = density_updated( j )
      end do		 
 
 
      if( g+1 > g_limit ) change = 0.0d0 
   enddo !do while( change > 1.0d-2 ) 

   !stop 

   deallocate( rko )
   deallocate( rhso )
   deallocate( ipmkl )
   deallocate( rko_id )
   deallocate( T_rhs_ID )
   deallocate( Omega_out_ID )
   deallocate( Omega_in_ID )
   deallocate( Element_Design_Area_T )
   deallocate( Element_Design_Area_S )
   deallocate( x_beta )
   deallocate( density_updated )
   deallocate( diff_density )
   deallocate( Triangle_density_vec )
   deallocate( Best_Triangle_density_vec )
   deallocate( Square_density_vec )
   
end program 

  subroutine mkGSM_0 (nd,Triangle_Number_Element,Number_Node,x,y,A,B,nodex,xcoord,ycoord,Local_matrix )
     implicit none
     integer i,j,k,l,ie
     integer nd,Number_Node,Triangle_Number_Element
     integer :: nodex( Triangle_Number_Element,nd )
     double precision :: Area
     double precision :: x( nd ), y( nd ), A( nd ), B( nd ), Local_matrix( nd,nd,Triangle_Number_Element )
     double precision :: xcoord( Number_Node ), ycoord( Number_Node )
 
     do ie = 1, Triangle_Number_Element

         do i = 1, nd
            x( i ) = xcoord( nodex( ie,i ) )
            y( i ) = ycoord( nodex( ie,i ) )
         end do
           
         B(1) = y(2) - y(3)
         B(2) = y(3) - y(1)
         B(3) = y(1) - y(2)
         A(1) = x(3) - x(2)
         A(2) = x(1) - x(3)
         A(3) = x(2) - x(1)
        
         Area = ( A(3)*B(2) - B(3)*A(2) )/ 2.0d0

         do i = 1, nd
            B(i) = B(i) / ( Area*2.0d0 )
            A(i) = A(i) / ( Area*2.0d0 )
         end do
     
         do i = 1, nd
            do j = 1, nd
               Local_matrix( i,j,ie ) = ( B(i)*B(j) + A(i)*A(j) )*Area 
            end do
         end do
                   
      end do

      return
   end subroutine mkGSM_0
 
   subroutine setBC_0( Triangle_Number_Element,Number_Node,nd,SizeA,xcoord,ycoord,Global_matrix,&
						Local_matrix,nodex,T_rhs,rko_id,rko,rhso,xmin,xmax,ymin,ymax )
      implicit none
      integer :: i,j,k,l,ie
      integer :: nd, Triangle_Number_Element, Number_Node, SizeA
      integer :: rko_id( SizeA ), nodex( Triangle_Number_Element,nd )
      double precision :: xmin, xmax, ymin, ymax 
      double precision :: xcoord( Number_Node ), ycoord( Number_Node ), T_rhs( Number_Node )
      double precision :: shape_function( Triangle_Number_Element,nd ), T_e( Triangle_Number_Element )
      double precision :: rko( SizeA,SizeA ), rhso( SizeA )
      double precision :: Local_matrix( nd,nd,Triangle_Number_Element ), Global_matrix( Number_Node,Number_Node )
  
      !======================================================== 
      !Making Global Matrix
      !========================================================  
  
      Global_matrix( :,: ) = 0.0d0
  
      do ie = 1, Triangle_Number_Element
	     do k = 1, nd
            i = nodex( ie,k ) 
            do l = 1, nd
               j = nodex( ie,l )  
               Global_matrix( i,j ) = Global_matrix( i,j ) + Local_matrix( k,l,ie )
            end do
         end do   
      end do
	  
      !Reset
      T_rhs( : ) = 0.0d0
  
      !========================================================	
      !Set Dirichlet boundary
      !========================================================
 
      do i = 1, Number_Node
        if( xcoord( i ) == xmin ) then
    
           do j = 1, Number_Node
              Global_matrix( i,j ) = 0.0d0
           end do
 
              Global_matrix( i,i ) = 1.0d0
              T_rhs( i ) = 0.0d0 
  
        else if( xcoord( i ) == xmax ) then
        
           do j = 1, Number_Node
              Global_matrix( i,j ) = 0.0d0
           end do
   
              Global_matrix( i,i ) = 1.0d0
              T_rhs( i ) = 1.0d0
   
        end if
      end do
   
      j = 1
      do i = 1, Number_Node
        if( Global_matrix( i,i ) /= 0.0d0 ) then
           rko_id( i ) = i
           j = j + 1
        end if
      end do

	  do i = 1, SizeA
         do j = 1, SizeA
            rko( i,j ) = Global_matrix( rko_id( i ), rko_id( j ) )
         end do
         rhso( i ) = T_rhs( rko_id( i ) )
      end do   
	             
      return 
   end subroutine setBC_0   
   

   subroutine mkGSM (nd,Type_Element_Insulator,Triangle_Number_Element,&
						Number_Node,kappa,eid,x,y,A,B,nodex,xcoord,ycoord,Local_matrix)
      implicit none
      integer i,j,k,l,ie
      integer nd,Number_Node,Triangle_Number_Element,Type_Element_Insulator
      integer :: nodex( Triangle_Number_Element,nd ),eid( Triangle_Number_Element )
      double precision :: Area
      double precision :: x( nd ), y( nd ), A( nd ), B( nd ), Local_matrix( nd,nd,Triangle_Number_Element )
      double precision :: xcoord( Number_Node ), ycoord( Number_Node )
      double precision :: kappa( Triangle_Number_Element )
  
      do ie = 1, Triangle_Number_Element
         if( eid( ie ) /= Type_Element_Insulator ) then
  
            do i = 1, nd
               x( i ) = xcoord( nodex( ie,i ) )
               y( i ) = ycoord( nodex( ie,i ) )
            end do
          
            B(1) = y(2) - y(3)
            B(2) = y(3) - y(1)
            B(3) = y(1) - y(2)
            A(1) = x(3) - x(2)
            A(2) = x(1) - x(3)
            A(3) = x(2) - x(1)
     
            Area = ( A(3)*B(2) - B(3)*A(2) ) / 2.0d0
            do i = 1, nd
               B(i) = B(i) / ( Area*2.0d0 ) 
               A(i) = A(i) / ( Area*2.0d0 )
            end do
   
            do i = 1, nd
               do j = 1, nd
                  Local_matrix( i,j,ie ) = ( B(i)*B(j) + A(i)*A(j) )*( kappa(ie) )*Area
               end do
            end do
            
         end if
      end do
 
      return
   end subroutine mkGSM
      
   subroutine setBC( Triangle_Number_Element,Type_Element_Insulator,Number_Node,nd,SizeA,xcoord,ycoord,nodex,&
						Local_matrix,Global_matrix,eid,T_rhs,rko_id,rko,rhso,xmin,xmax,ymin,ymax )
      implicit none
      integer :: i,j,k,l,ie
      integer :: nd, Triangle_Number_Element, Number_Node, SizeA, Type_Element_Insulator
      integer :: rko_id( SizeA ), nodex( Triangle_Number_Element,nd )
	  integer :: eid( Triangle_Number_Element )
      double precision :: xmin, xmax, ymin, ymax 
      double precision :: xcoord( Number_Node ), ycoord( Number_Node ), T_rhs( Number_Node )
      double precision :: shape_function( Triangle_Number_Element,nd ), T_e( Triangle_Number_Element )
      double precision :: rko( SizeA,SizeA ), rhso( SizeA )
      double precision :: Local_matrix( nd,nd,Triangle_Number_Element ), Global_matrix( Number_Node,Number_Node )
  
      !======================================================== 
      !Making Global Matrix
      !========================================================
	  
      Global_matrix( :,: ) = 0.0d0

      do ie = 1, Triangle_Number_Element
	     if ( eid( ie ) /= Type_Element_Insulator ) then
	     do k = 1, nd
            i = nodex( ie,k ) 
            do l = 1, nd
               j = nodex( ie,l )  
               Global_matrix( i,j ) = Global_matrix( i,j ) + Local_matrix( k,l,ie )
            end do
         end do   
         end if			
      end do

      !Reset
      T_rhs( : ) = 0.0d0
  
      !======================================================== 
      !Set Norman boundary
      !========================================================
 
      do i = 1, Number_Node
         if( ycoord( i ) == ymin .or. ycoord( i ) == ymax ) then
            T_rhs( i ) = 0.0d0
         end if
      end do
  
      !========================================================	
      !Set Dirichlet boundary
      !========================================================
 
      do i = 1, Number_Node
        if( xcoord( i ) == xmin ) then
    
           do j = 1, Number_Node
              Global_matrix( i,j ) = 0.0d0
           end do
 
              Global_matrix( i,i ) = 1.0d0
              T_rhs( i ) = 0.0d0 
  
        else if( xcoord( i ) == xmax ) then
        
           do j = 1, Number_Node
              Global_matrix( i,j ) = 0.0d0
           end do
   
              Global_matrix( i,i ) = 1.0d0
              T_rhs( i ) = 1.0d0
   
        end if
      end do
     
      j = 1
      do i = 1, Number_Node
        if( Global_matrix( i,i ) /= 0.0d0 ) then
           rko_id( j ) = i
           j = j + 1
        end if
      end do
 
      do i = 1, SizeA
         do j = 1, SizeA
            rko( i,j ) = Global_matrix( rko_id( i ), rko_id( j ) )
         end do
         rhso( i ) = T_rhs( rko_id( i ) )
      end do   
	             
      return 
   end subroutine setBC
   
   subroutine setBC_sv( Triangle_Number_Element,Number_Node,Type_Element_Insulator,&
							nd,SizeA,xcoord,ycoord,Local_matrix,Global_matrix,T_rhs,&
							rko_id,rko,rhso,xmin,xmax,ymin,ymax,ie,eid,counter,x,y,&
							nodex,Area,a,b,psi_n,T_n,T_ref,lm,gm,Q,Q1 )
      implicit none
      integer :: i,j,k,l,ie
      integer :: nd, Triangle_Number_Element, Number_Node, SizeA,Type_Element_Insulator, Type_Element_Outer
      integer :: eid( Triangle_Number_Element ), rko_id( SizeA )
	  integer :: nodex( Triangle_Number_Element,nd ),counter( Number_Node )
      double precision :: xmin, xmax, ymin, ymax, Area, psi_n 
	  double precision :: x( nd ), y( nd )
	  double precision :: a( nd ), b( nd ), c( nd )
      double precision :: xcoord( Number_Node ), ycoord( Number_Node ), T_rhs( Number_Node )
      double precision :: T_n( Number_Node ), T_ref( Number_Node )
      double precision :: rko( SizeA,SizeA ), rhso( SizeA )
	  double precision :: Local_matrix( nd,nd,Triangle_Number_Element ),Global_matrix( Number_Node,Number_Node )
	  double precision :: lm( nd ), gm( Number_Node ), Q( Number_Node ), Q1( Number_Node )
	  
	  gm( : ) = 0.0d0
	  lm( : ) = 0.0d0
	  
	  do ie = 1, Triangle_Number_Element
         if( eid( ie ) /= Type_Element_Insulator ) then

            do i = 1, nd
               x( i ) = xcoord( nodex( ie,i ) )
               y( i ) = ycoord( nodex( ie,i ) )
            end do
           
		    Area = ( x(1)*( y(2)-y(3) ) + x(2)*( y(3)-y(1) ) + x(3)*( y(1)-y(2) )) / 2.0d0
		   
		    a(1) = ( x(2)*y(3) - x(3)*y(2) ) / ( 2.0d0 * Area )
		    b(1) = ( y(2) - y(3) ) / ( 2.0d0 * Area )
			c(1) = ( x(3) - x(2) ) / ( 2.0d0 * Area )
			
		    a(2) = ( x(3)*y(1) - x(1)*y(3) ) / ( 2.0d0 * Area )
		    b(2) = ( y(3) - y(1) ) / ( 2.0d0 * Area )
			c(2) = ( x(1) - x(3) ) / ( 2.0d0 * Area )			
		  
		    a(3) = ( x(1)*y(2) - x(2)*y(1) ) / ( 2.0d0 * Area )
		    b(3) = ( y(1) - y(2) ) / ( 2.0d0 * Area )
			c(3) = ( x(2) - x(1) ) / ( 2.0d0 * Area )		  
			
			if( eid( ie ) == 5 ) then		     
			   do i = 1, nd 
			      lm( i ) =  ( 2.0d0/psi_n )&
					  			    *abs( 	a(i)*( ( T_n( nodex( ie,i )) - T_ref( nodex( ie,i ))))&
			   					    + b(i)*x(i)*( ( T_n( nodex( ie,i )) - T_ref( nodex( ie,i ))))&
								    + c(i)*y(i)*( ( T_n( nodex( ie,i )) - T_ref( nodex( ie,i )))))
			   end do	
		   
               do j = 1, nd
                  i = nodex( ie,j ) 
                  gm( i ) = gm( i ) + lm( j )
               end do
			
			end if	
         end if 
      end do

      !========================================================	
      !Making Global Matrix
      !========================================================	  
	  
      Global_matrix( :,: ) = 0.0d0

      do ie = 1, Triangle_Number_Element
	     if ( eid( ie ) /= Type_Element_Insulator ) then
	     do k = 1, nd
            i = nodex( ie,k ) 
            do l = 1, nd
               j = nodex( ie,l )  
               Global_matrix( i,j ) = Global_matrix( i,j ) + Local_matrix( k,l,ie )
            end do
         end do   
         end if			
      end do	  

      !========================================================	
      !Set Dirichlet boundary
      !========================================================

      !Reset
      T_rhs( : ) = 0.0d0
	  Q( : ) = 0.0d0	  	  	  
 
      do i = 1, Number_Node
        if( xcoord( i ) == xmin ) then
    
           do j = 1, Number_Node
              Global_matrix( i,j ) = 0.0d0
           end do
		   
		   do j = 1, Number_Node	
              Q( j ) = Q( j ) - Global_matrix( j,i )*0.0d0
           end do
		   
		   do j = 1, Number_Node
		      Global_matrix( j,i ) = 0.0d0
		   end do	  
 
              Global_matrix( i,i ) = 1.0d0
              Q( i ) = 0.0d0 
			    
        else if( xcoord( i ) == xmax ) then
        
           do j = 1, Number_Node
              Global_matrix( i,j ) = 0.0d0
           end do
		   
           do j = 1, Number_Node	
              Q( j ) = Q( j ) - Global_matrix( j,i )*1.0d0
           end do				   
		   
		   do j = 1, Number_Node
		      Global_matrix( j,i ) = 0.0d0
		   end do	 
		   
              Global_matrix( i,i ) = 1.0d0
              Q( i ) = 1.0d0  
   
        end if
      end do
   
      j = 1
      do i = 1, Number_Node
        if( Global_matrix( i,i ) /= 0.0d0 ) then
           rko_id( j ) = i
           j = j + 1
        end if
      end do
 
      do i = 1, SizeA
         do j = 1, SizeA
            rko( i,j ) = Global_matrix( rko_id( i ), rko_id( j ) )
         end do
         rhso( i ) = Q( rko_id( i ) )
      end do   		  

      return 
   end subroutine setBC_sv   
 
   subroutine etemp( Triangle_Number_Element,Number_Node,Type_Element_Insulator,nd,SizeA,Area,eid,&
                      xcoord,ycoord,nec,x,y,rhso,rko_id,nodex,shape_function,T_rhs,T_e )  
      implicit none
      integer :: i,j
      integer :: Triangle_Number_Element, Number_Node, nd, SizeA, Type_Element_Insulator
      integer :: nodex( Triangle_Number_Element,nd ), eid( Triangle_Number_Element )
      integer :: rko_id( SizeA )
      double precision :: Area
      double precision :: x( nd ), y( nd ), rhso( SizeA )
      double precision :: xcoord( Number_Node ), ycoord( Number_Node ), nec( Triangle_Number_Element,2 ) ,T_rhs( Number_Node )
      double precision :: shape_function( Triangle_Number_Element,nd ), T_e( Triangle_Number_Element )

      do i = 1, Triangle_Number_Element
         if( eid( i ) /= Type_Element_Insulator ) then
            do j = 1, nd
               x( j ) = xcoord( nodex( i,j ) )
               y( j ) = ycoord( nodex( i,j ) )
            end do     
                   
            Area = (( x( 2 ) - x( 1 ) )*( y( 3 ) - y( 1 ) )) - (( y( 1 ) - y( 2 ))*( x( 1 ) - x( 3 ) ))
            shape_function( i,1 ) = ( x( 2 )*y( 3 ) - x( 3 )*y( 2 ) ) /Area&
                                  + ( nec( i,1 )*( y( 2 ) - y( 3 ) )) /Area&
                                  + ( nec( i,2 )*( x( 3 ) - x( 2 ) )) /Area 
                                                                     
            shape_function( i,2 ) = ( x( 3 )*y( 1 ) - x( 1 )*y( 3 ) ) /Area&
                                  + ( nec( i,1 )*( y( 3 ) - y( 1 ) )) /Area&
                                  + ( nec( i,2 )*( x( 1 ) - x( 3 ) )) /Area
                                                                     
            shape_function( i,3 ) = ( x( 1 )*y( 2 ) - x( 2 )*y( 1 ) ) /Area&
                                  + ( nec( i,1 )*( y( 1 ) - y( 2 ) )) /Area&
                                  + ( nec( i,2 )*( x( 2 ) - x( 1 ) )) /Area 
                
            T_e( i ) = shape_function( i,1 ) *T_rhs( nodex( i,1 ) )&
                     + shape_function( i,2 ) *T_rhs( nodex( i,2 ) )&
                     + shape_function( i,3 ) *T_rhs( nodex( i,3 ) )
            
         end if
      end do 

      return 
   end subroutine etemp 
 
subroutine Compute_Objective_Function &
       ( ID_Element, &
         Value_Obj_Func, Number_Node_Value_Obj_Func, &
         Position_Node, Number_Node, &
         Index_Element_2_Node, Class_Element, Number_Element, &
         Optimization_Step, Loop_Solution, &
         Perimeter_Structure, &
         !=================================================================
         Obj_Func )

   !$use omp_lib
 !  use Parameters
   implicit none

   integer, intent(in) :: ID_Element 
   integer, intent(in) :: Number_Node_Value_Obj_Func 
   integer, intent(in) :: Number_Node, Number_Element

   !complex( kind( 0d0 ) ), intent(in) :: Value_Obj_Func( Number_Node_Value_Obj_Func )
   double precision, intent(in) :: Value_Obj_Func( Number_Node_Value_Obj_Func )

   double precision, intent(in) :: Position_Node( 2, Number_Node )
   integer, intent(in) :: Index_Element_2_Node( 3, Number_Element ) 
   integer, intent(in) :: Class_Element( Number_Element ) 

   integer, intent(in) :: Optimization_Step, Loop_Solution

   double precision, intent(in) :: Perimeter_Structure

   double precision, intent(out) :: Obj_Func

   integer :: e, i, j, k
   double precision, allocatable, dimension(:) :: Obj_Func_Element
   double precision, allocatable, dimension(:) :: Obj_Func_Element_tmp
   double precision, allocatable, dimension(:,:,:) :: Position_Node_Element_Triangle
   double precision, allocatable, dimension(:,:,:,:) :: Difference_Position_Node_Element_Triangle
   double precision, allocatable, dimension(:) :: Area_Element  
   double precision, allocatable, dimension(:,:) :: BasisFunction_A, BasisFunction_B, BasisFunction_C
   double precision, allocatable, dimension(:,:,:) :: BasisFunctionBC_Times_Difference_Position_Node_TE
   double precision, allocatable, dimension(:,:) :: Basis_Function_Times_Position_Node_TE
   !complex( kind( 0d0 ) ), allocatable, dimension(:,:,:) :: Coefficient_C, Coefficient_xi, Coefficient_eta 
   !complex( kind( 0d0 ) ), allocatable, dimension(:,:,:) :: Coefficient_xi_eta, Coefficient_xi_xi, Coefficient_eta_eta
   !complex( kind( 0d0 ) ), allocatable, dimension(:,:,:) :: Integral_Element_Triangle
   double precision, allocatable, dimension(:,:,:) :: Coefficient_C, Coefficient_xi, Coefficient_eta 
   double precision, allocatable, dimension(:,:,:) :: Coefficient_xi_eta, Coefficient_xi_xi, Coefficient_eta_eta
   double precision, allocatable, dimension(:,:,:) :: Integral_Element_Triangle

   integer :: Counter_Element_OnjectiveFunction, Number_Element_ObjectiveFunction
   integer, allocatable, dimension(:,:) :: Index_Element_2_Node_ObjectiveFunction, Index_Element_2_Node_ObjectiveFunction_tmp

   double precision, allocatable, dimension(:,:,:) :: Intensity_Value_Obj_Func

   integer, allocatable, dimension( : ) :: Counter_Element_OnjectiveFunction_ID

   !=======================================================================================================
   write(*,*)'   call Compute_Objective_Function ' 
   !=======================================================================================================

   allocate( Index_Element_2_Node_ObjectiveFunction_tmp( 3, Number_Element ) ) 
      
   Counter_Element_OnjectiveFunction= 0
   ! No Parallel, Vectorize

   do e= 1, Number_Element
      if( Class_Element( e )==ID_Element )then

         Counter_Element_OnjectiveFunction= Counter_Element_OnjectiveFunction +1

         do i= 1, 3
            Index_Element_2_Node_ObjectiveFunction_tmp( i, Counter_Element_OnjectiveFunction ) & 
            = Index_Element_2_Node( i, e ) 
         end do
      end if
   end do

   Number_Element_ObjectiveFunction= Counter_Element_OnjectiveFunction

   allocate( Index_Element_2_Node_ObjectiveFunction( 3, Number_Element_ObjectiveFunction ) ) 
   Index_Element_2_Node_ObjectiveFunction( 1:3, 1:Number_Element_ObjectiveFunction )&
   = Index_Element_2_Node_ObjectiveFunction_tmp( 1:3, 1:Number_Element_ObjectiveFunction )
   

   deallocate( Index_Element_2_Node_ObjectiveFunction_tmp ) 

   !===================================================================================
   write(*,*)'      Position_Node --> Position_Node_Element_Triangle '
   !===================================================================================
   allocate( Position_Node_Element_Triangle( 2, 3, Number_Element_ObjectiveFunction ) )
   
   !$omp parallel do default( none )   &
   !$omp private( e, i, j ) &
   !$omp shared( Number_Element_ObjectiveFunction, Position_Node_Element_Triangle, Position_Node, Index_Element_2_Node_ObjectiveFunction ) 
   do e= 1, Number_Element_ObjectiveFunction
      do i= 1, 3
         do j= 1, 2
            Position_Node_Element_Triangle( j, i, e )= Position_Node( j, Index_Element_2_Node_ObjectiveFunction( i, e ) )
         end do
      end do
   end do

   allocate( Difference_Position_Node_Element_Triangle( 2, 3, 3, Number_Element_ObjectiveFunction ) )

   do e= 1, Number_Element_ObjectiveFunction
      do i= 1, 3
         do j= 1, 3
            do k= 1, 2
               Difference_Position_Node_Element_Triangle( k, j, i, e ) &
               =Position_Node_Element_Triangle( k, j, e ) -Position_Node_Element_Triangle( k, i, e )
            end do
         end do
      end do
   end do

   allocate( Area_Element( Number_Element_ObjectiveFunction ) )
   Area_Element( : ) &
   =abs( Difference_Position_Node_Element_Triangle( 1, 1, 3, : ) & 
      *Difference_Position_Node_Element_Triangle( 2, 2, 3, : ) &
      -Difference_Position_Node_Element_Triangle( 1, 2, 3, : ) & 
      *Difference_Position_Node_Element_Triangle( 2, 1, 3, : ) )/2.0D0

   !$omp parallel do default( none ) &
   !$omp private( e ) & 
   !$omp shared( Number_Element_ObjectiveFunction, Area_Element, Position_Node_Element_Triangle ) 
   do e= 1, Number_Element_ObjectiveFunction
      if( Area_Element( e )==0.0d0 )then
         write(*,*)'Area_Element( e )==0.0d0'
         write(*,*)'Element Number=', e
         write(*,*)'Number_Element_ObjectiveFunction=', Number_Element_ObjectiveFunction
        ! call Output_Error( 'Compute_Objective_Function', 166 )
		stop
      end if
   end do 
    
   allocate( BasisFunction_A( 3, Number_Element_ObjectiveFunction ) )
   allocate( BasisFunction_B( 3, Number_Element_ObjectiveFunction ) )
   allocate( BasisFunction_C( 3, Number_Element_ObjectiveFunction ) )

   BasisFunction_A( 1, : ) &
   = ( Position_Node_Element_Triangle( 1, 2, : ) *Position_Node_Element_Triangle( 2, 3, : ) &
      -Position_Node_Element_Triangle( 1, 3, : ) *Position_Node_Element_Triangle( 2, 2, : ) ) &
    /( 2.0D0 *Area_Element( : ) )

   BasisFunction_A( 2, : ) &
   = ( Position_Node_Element_Triangle( 1, 3, : ) *Position_Node_Element_Triangle( 2, 1, : ) & 
      -Position_Node_Element_Triangle( 1, 1, : ) *Position_Node_Element_Triangle( 2, 3, : ) ) &
    /( 2.0D0 *Area_Element( : ) )

   BasisFunction_A( 3, : ) & 
   = ( Position_Node_Element_Triangle( 1, 1, : ) *Position_Node_Element_Triangle( 2, 2, : ) & 
      -Position_Node_Element_Triangle( 1, 2, : ) *Position_Node_Element_Triangle( 2, 1, : ) ) & 
    /( 2.0D0 *Area_Element( : ) )

   BasisFunction_B( 1, : ) &
   = Difference_Position_Node_Element_Triangle( 2, 2, 3, : )/( 2.0D0 *Area_Element( : ) )

   BasisFunction_B( 2, : ) & 
   = Difference_Position_Node_Element_Triangle( 2, 3, 1, : )/( 2.0D0 *Area_Element( : ) )

   BasisFunction_B( 3, : ) & 
   = Difference_Position_Node_Element_Triangle( 2, 1, 2, : )/( 2.0D0 *Area_Element( : ) )
      
   BasisFunction_C( 1, : ) & 
   = Difference_Position_Node_Element_Triangle( 1, 3, 2, : )/( 2.0D0 *Area_Element( : ) )

   BasisFunction_C( 2, : ) & 
   = Difference_Position_Node_Element_Triangle( 1, 1, 3, : )/( 2.0D0 *Area_Element( : ) )

   BasisFunction_C( 3, : ) & 
   = Difference_Position_Node_Element_Triangle( 1, 2, 1, : )/( 2.0D0 *Area_Element( : ) )

   allocate( Basis_Function_Times_Position_Node_TE( 3, Number_Element_ObjectiveFunction ) ) 

   !$omp parallel do default( none )   &
   !$omp private( e, i ) &
   !$omp shared( Number_Element_ObjectiveFunction, Basis_Function_Times_Position_Node_TE ) &
   !$omp shared( BasisFunction_A, BasisFunction_B, BasisFunction_C, Position_Node_Element_Triangle )
   do e= 1, Number_Element_ObjectiveFunction
      do i= 1, 3
         Basis_Function_Times_Position_Node_TE( i, e ) & 
         = BasisFunction_A( i, e ) & 
          +BasisFunction_B( i, e ) *Position_Node_Element_Triangle( 1, 3, e ) & 
          +BasisFunction_C( i, e ) *Position_Node_Element_Triangle( 2, 3, e ) 
      end do
   end do

   allocate( BasisFunctionBC_Times_Difference_Position_Node_TE( 2, 3, Number_Element_ObjectiveFunction ) )

   !$omp parallel do default( none )   &
   !$omp private( e, i, j ) &
   !$omp shared( Number_Element_ObjectiveFunction, BasisFunctionBC_Times_Difference_Position_Node_TE ) &
   !$omp shared( BasisFunction_A, BasisFunction_B, BasisFunction_C, Difference_Position_Node_Element_Triangle )
   do e= 1, Number_Element_ObjectiveFunction
      do i= 1, 3
         do j= 1, 2
            BasisFunctionBC_Times_Difference_Position_Node_TE( j, i, e ) &
            = BasisFunction_B( i, e ) *Difference_Position_Node_Element_Triangle( 1, j, 3, e ) &
             +BasisFunction_C( i, e ) *Difference_Position_Node_Element_Triangle( 2, j, 3, e )
         end do
      end do
   end do
 
   deallocate( Difference_Position_Node_Element_Triangle )
 
   allocate( Coefficient_C( 3, 3, Number_Element_ObjectiveFunction ) )
   allocate( Coefficient_xi( 3, 3, Number_Element_ObjectiveFunction ) )
   allocate( Coefficient_eta( 3, 3, Number_Element_ObjectiveFunction ) )
   allocate( Coefficient_xi_eta( 3, 3, Number_Element_ObjectiveFunction ) )
   allocate( Coefficient_xi_xi( 3, 3, Number_Element_ObjectiveFunction ) )
   allocate( Coefficient_eta_eta( 3, 3, Number_Element_ObjectiveFunction ) )

 
   !$omp parallel do default( none )   &
   !$omp private( e, i, j ) &
   !$omp shared( Number_Element_ObjectiveFunction ) &
   !$omp shared( Basis_Function_Times_Position_Node_TE, Coefficient_C ) 
   do e= 1, Number_Element_ObjectiveFunction
      do i= 1, 3
         do j= 1, 3
            Coefficient_C( j, i, e ) & 
            = Basis_Function_Times_Position_Node_TE( i, e )*Basis_Function_Times_Position_Node_TE( j, e ) 
         end do
      end do
   end do

   !$omp parallel do default( none )   &
   !$omp private( e, i, j ) &
   !$omp shared( Number_Element_ObjectiveFunction, Coefficient_Xi ) & 
   !$omp shared( Basis_Function_Times_Position_Node_TE ) &
   !$omp shared( BasisFunctionBC_Times_Difference_Position_Node_TE ) 
   do e= 1, Number_Element_ObjectiveFunction
      do i= 1, 3
         do j= 1, 3
            Coefficient_Xi( j, i, e ) & 
            = Basis_Function_Times_Position_Node_TE( j, e ) &
             *BasisFunctionBC_Times_Difference_Position_Node_TE( 1, i, e ) &
             +Basis_Function_Times_Position_Node_TE( i, e ) &
             *BasisFunctionBC_Times_Difference_Position_Node_TE( 1, j, e ) 
         end do
      end do
   end do

   !$omp parallel do default( none )   &
   !$omp private( e, i, j ) &
   !$omp shared( Number_Element_ObjectiveFunction, Coefficient_Eta ) & 
   !$omp shared( Basis_Function_Times_Position_Node_TE ) &
   !$omp shared( BasisFunctionBC_Times_Difference_Position_Node_TE ) 
   do e= 1, Number_Element_ObjectiveFunction
      do i= 1, 3
         do j= 1, 3
            Coefficient_Eta( j, i, e ) &
            = Basis_Function_Times_Position_Node_TE( j, e ) &
             *BasisFunctionBC_Times_Difference_Position_Node_TE( 2, i, e ) &
             +Basis_Function_Times_Position_Node_TE( i, e ) &
             *BasisFunctionBC_Times_Difference_Position_Node_TE( 2, j, e ) 
         end do
      end do
   end do

   deallocate( Basis_Function_Times_Position_Node_TE ) 

   !$omp parallel do default( none )   &
   !$omp private( e, i, j ) &
   !$omp shared( Number_Element_ObjectiveFunction, Coefficient_Xi_Eta ) & 
   !$omp shared( BasisFunctionBC_Times_Difference_Position_Node_TE ) 
   do e= 1, Number_Element_ObjectiveFunction
      do i= 1, 3
         do j= 1, 3
            Coefficient_Xi_Eta( j, i, e ) & 
            = BasisFunctionBC_Times_Difference_Position_Node_TE( 2, i, e ) &
             *BasisFunctionBC_Times_Difference_Position_Node_TE( 1, j, e ) &
             +BasisFunctionBC_Times_Difference_Position_Node_TE( 1, i, e ) &
             *BasisFunctionBC_Times_Difference_Position_Node_TE( 2, j, e ) 
         end do
      end do
   end do

   !$omp parallel do default( none )   &
   !$omp private( e, i, j ) &
   !$omp shared( Number_Element_ObjectiveFunction, Coefficient_Xi_Xi ) & 
   !$omp shared( BasisFunctionBC_Times_Difference_Position_Node_TE ) 
   do e= 1, Number_Element_ObjectiveFunction
      do i= 1, 3
         do j= 1, 3
            Coefficient_Xi_Xi( j, i, e ) & 
            = BasisFunctionBC_Times_Difference_Position_Node_TE( 1, i, e ) &
             *BasisFunctionBC_Times_Difference_Position_Node_TE( 1, j, e ) 
         end do
      end do
   end do

   !$omp parallel do default( none )   &
   !$omp private( e, i, j ) &
   !$omp shared( Number_Element_ObjectiveFunction, Coefficient_Eta_Eta ) & 
   !$omp shared( BasisFunctionBC_Times_Difference_Position_Node_TE ) 
   do e= 1, Number_Element_ObjectiveFunction
      do i= 1, 3
         do j= 1, 3
            Coefficient_Eta_Eta( j, i, e ) & 
            = BasisFunctionBC_Times_Difference_Position_Node_TE( 2, i, e ) &
             *BasisFunctionBC_Times_Difference_Position_Node_TE( 2, j, e ) 
         end do
      end do
   end do
   
   deallocate( BasisFunctionBC_Times_Difference_Position_Node_TE )
   deallocate( BasisFunction_A )
   deallocate( BasisFunction_B )
   deallocate( BasisFunction_C )
  
   allocate( Integral_Element_Triangle( 3, 3, Number_Element_ObjectiveFunction ) )
 
   !$omp parallel do default( none )   &
   !$omp private( e, i, j ) &
   !$omp shared( Number_Element_ObjectiveFunction, Integral_Element_Triangle ) & 
   !$omp shared( Coefficient_C, Coefficient_xi, Coefficient_eta, Coefficient_xi_eta, Coefficient_xi_xi, Coefficient_eta_eta ) 
   do e= 1, Number_Element_ObjectiveFunction
      do i= 1, 3
         do j= 1, 3
            Integral_Element_Triangle( j, i, e ) & 
            = Coefficient_C( j, i, e ) +Coefficient_xi( j, i, e )/3.0d0 +Coefficient_eta( j, i, e )/3.0d0 &
             +Coefficient_xi_xi( j, i, e )/6.0d0 +Coefficient_eta_eta( j, i, e )/6.0d0 +Coefficient_xi_eta( j, i, e )/12.0d0 
         end do
      end do
   end do

   deallocate( Coefficient_C )
   deallocate( Coefficient_xi )
   deallocate( Coefficient_eta )
   deallocate( Coefficient_xi_eta )
   deallocate( Coefficient_xi_xi )
   deallocate( Coefficient_eta_eta ) 
  
   allocate( Obj_Func_Element_tmp( Number_Element_ObjectiveFunction ) )
   allocate( Intensity_Value_Obj_Func( 3, 3, Number_Element_ObjectiveFunction ) )

   !$omp parallel do default( none )   &
   !$omp private( e, i, j ) &
   !$omp shared( Number_Element_ObjectiveFunction, Intensity_Value_Obj_Func, Value_Obj_Func, Index_Element_2_Node_ObjectiveFunction ) 
   do e= 1, Number_Element_ObjectiveFunction
      do i= 1, 3
         do j= 1, 3
            Intensity_Value_Obj_Func( j, i, e ) &
		     = Value_Obj_Func( Index_Element_2_Node_ObjectiveFunction( i, e ) ) & 
			 *( Value_Obj_Func( Index_Element_2_Node_ObjectiveFunction( j, e ) ) ) 
         end do
 	  end do
   end do

   Obj_Func_Element_tmp= 0.0d0 
 
   !$omp parallel do default( none )   &
   !$omp private( e, i, j ) &
   !$omp shared( Number_Element_ObjectiveFunction, Obj_Func_Element_tmp ) &
   !$omp shared( Intensity_Value_Obj_Func, Integral_Element_Triangle ) 
   do e= 1, Number_Element_ObjectiveFunction
      do i= 1, 3
         do j= 1, 3
            Obj_Func_Element_tmp( e ) &
            = Obj_Func_Element_tmp( e ) &
              +Intensity_Value_Obj_Func( j, i, e ) *Integral_Element_Triangle( j, i, e ) 
         end do
      end do
   end do

   deallocate( Intensity_Value_Obj_Func )
   deallocate( Integral_Element_Triangle )
 
   allocate( Obj_Func_Element( Number_Element_ObjectiveFunction ) )

   !$omp parallel do default( none )   &
   !$omp private( e ) &
   !$omp shared( Number_Element_ObjectiveFunction, Obj_Func_Element, Area_Element, Obj_Func_Element_tmp ) 
   do e= 1, Number_Element_ObjectiveFunction
      Obj_Func_Element( e )= Area_Element( e ) *Obj_Func_Element_tmp( e )
   end do

   deallocate( Obj_Func_Element_tmp )

   Obj_Func= 0.0d0

   do e= 1, Number_Element_ObjectiveFunction
      Obj_Func= Obj_Func +Obj_Func_Element( e )
   end do

   deallocate( Position_Node_Element_Triangle )
   deallocate( Area_Element )
   deallocate( Obj_Func_Element )

   deallocate( Index_Element_2_Node_ObjectiveFunction ) 

   return
end subroutine Compute_Objective_Function   

subroutine Output_the_best_value( g_limit,g_min,Number_EDesignArea,Best_Triangle_density_vec,Number_Node,&
									counter,xcoord,ycoord,T_mat )
	implicit none
		integer :: i,g_limit,g_min
		integer :: Number_EDesignArea,Number_Node
		integer :: counter( Number_Node )
		
		double precision :: Best_Triangle_density_vec( Number_EDesignArea )
		double precision :: xcoord( Number_Node ), ycoord( Number_Node )
		double precision :: T_mat( g_limit,Number_Node )
	
			open( 500, file='ObjectiveFunction.dat', position='append' )
				write(500,*)'g_min=', g_min
			close(500)	
			
			open( 600,file='density.txt',status='replace' )
               do i = 1, Number_EDesignArea
                  write(600,*) Best_Triangle_density_vec( i )
               end do
            close(600)
            
			open( 700,file='temperature.txt',status='replace' )
               do i = 1, Number_Node
                  if( counter( i ) == 1  ) then
                     write(700,*)  xcoord( i ), ycoord( i ), T_mat( g_min, i )
                  end if
               end do
            close(700)


end subroutine Output_the_best_value

subroutine Compute_difference_in_temperature( g_limit,g_min,Number_Node,diff_T,counter1,xcoord,ycoord,T_ref,T_mat,Td)
	implicit none
	
	integer :: i,g_limit,g_min,Number_Node
	double precision :: diff_T

	integer :: counter1( Number_Node )	
	double precision :: xcoord( Number_Node ), ycoord( Number_Node )
	double precision :: T_ref( Number_Node ), T_mat( g_limit,Number_Node ), Td( Number_Node )
		
	!==============================================================================

	do i = 1, Number_Node
		if ( counter1( i ) == 1 ) then
			Td( i ) = abs( T_mat( g_min, i ) - T_ref( i ) )  
		end if
	end do

    diff_T = sum( Td )

	write(*,*)'++++++++++++++++++++++++++++++++++++++++++++++'
	write(*,*)' Total temperature difference at each element = ',diff_T
	write(*,*)'++++++++++++++++++++++++++++++++++++++++++++++'
	
	open( 444, file='diff_T.dat' )
		do i = 1, Number_Node
			if( counter1( i ) == 1 ) then 
				write(444,*) xcoord( i ), ycoord( i ), Td( i )
			endif
		end do
	close(444)	

	write(*,*)'========= diff_T Finish =========='

end subroutine Compute_difference_in_temperature


subroutine Make_Analysis_date( i,j,nd,Triangle_Number_Element,Number_Node,Number_EDesignArea&
						,nodex,eid,xmin,xmax,ymin,ymax,xcoord,ycoord&
						,Best_Triangle_density_vec,high,thickness,penal,g_min)
	implicit none
	integer :: i,j,nd,penal,g_min
	integer :: Triangle_Number_Element, Number_Node, Number_EDesignArea
    integer :: nodex( Triangle_Number_Element,nd ), eid( Triangle_Number_Element )
	double precision :: high, thickness
	double precision :: xmin, xmax, ymin, ymax
    double precision :: xcoord( Number_Node ), ycoord( Number_Node )
	double precision :: Best_Triangle_density_vec( Number_EDesignArea ) 

	open( 888, file='FEM_Data_Interval_001000' )
		write( 888,* ) 1,'# Flag_Electrical_Insulation_BC '
		write( 888,* ) Number_Node,'# The Number of Nodes'
		write( 888,* ) Triangle_Number_Element,'# The Number of Triangle Element'
		write( 888,* ) 6,'# Maximum Number of Elements sharing one Nodes '
		write( 888,* ) 0,'# The Number of Nodes on PEC Boundary'
		write( 888,* ) 0,'# The Number of Nodes in PEC'
		write( 888,* ) 0,'# The Number of Element on PEC Boundary'
		write( 888,* ) 5,'# ID_Element_Obj_Func_1'
		write( 888,* ) 1,5,2,4,'# ID_Element'
		write( 888,* ) 6,5,'# ID_Element_Exterior'
		write( 888,* ) high,'# Cloakroom Height '
		write( 888,* ) thickness,'# Plate thickness '
		write( 888,* ) penal,'# Penalty for density method'
		write( 888,* ) g_min,'# Convergent Generations'		
		write( 888,* ) xmin,xmax,ymin,ymax

		
		do i = 1, Number_Node
			write( 888,* ) i, xcoord( i ), ycoord( i )
		end do	
		
		j = 1
		do i = 1, Triangle_Number_Element
			if ( eid( i ) == 4 ) then
				write( 888,* ) i, eid( i ), Best_Triangle_density_vec( j ),&
								nodex( i,1 ), nodex( i,2 ), nodex( i,3 )
				j = j + 1
			else 
				write( 888,* ) i, eid( i ), 0.0d0, nodex( i,1 ), nodex( i,2 ), nodex( i,3 )		
			end if
		end do
	close( 888 )

end subroutine Make_Analysis_date

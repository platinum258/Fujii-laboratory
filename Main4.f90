
program thermalcloak_ITRfree
   !$ use omp_lib
   implicit none
   
!============================================================================================   
     
   integer :: num_threads, SizeA
     
   double precision,allocatable :: rko( :,: ), rhso( : ), ipmkl( : )
   integer,allocatable :: rko_id( : )
   
   integer,allocatable :: Element_Design_Area( : )
   double precision,allocatable :: x_beta( : ), density_updated( : ), diff_density( : )
   double precision,allocatable :: Triangle_density_vec( : ), Best_Triangle_density_vec( : )
   
!============================================================================================   
   
   character filename*128
   integer,parameter :: nelx = 120  !全体領域のx方向の分割数
   integer,parameter :: nely = 80   !全体領域のy方向の分割数
   
   integer,parameter :: nd = 3
   double precision, parameter :: vol_constraint = 0.5d0 !0.5d0
   integer,parameter :: penal = 3, rmin = 3
 
   integer,parameter :: Number_Node = ( nely+1 )*( nelx+1 )
   integer,parameter :: Square_Number_Element   = nely*nelx             !9600 
   integer,parameter :: Triangle_Number_Element = nelx*nely*2

 !  integer,parameter :: SizeA = 9508 !9801-293     !全体の接点数ー断熱部分の接点数
   integer,parameter :: g_limit = 300
   integer,parameter :: e_Ip = 293                  !断熱部分の接点数
 
   double precision,parameter :: high = 2.0d0       !クロークの高さ
   double precision,parameter :: thickness = 1.0d0  !平板の厚さ
   
   integer :: counter( Number_Node )
   double precision :: xcoord( Number_Node ), ycoord( Number_Node ), nesc( Square_Number_Element,2 )
   double precision,parameter :: mnx = - 60.0d0, mxx = 60.0d0, mny = - 40.0d0, mxy = 40.0d0
!============================================================================================  
   
   integer :: Number_EDesignArea  !断熱部分を除いた三角形要素数
   integer :: Number_Eins
   integer :: Number_SEDesignArea
   integer :: Number_SEins
   integer :: nodex_s( Square_Number_Element,4 ), nodex( Triangle_Number_Element,3 ), eid( Triangle_Number_Element )
   double precision :: dis
             
!============================================================================================
   integer :: g, Best_generation
   integer :: i,j,k,l,m,n,e,t,ie,e1,e2,info,node,int_tmp
   integer :: position_in_element( Triangle_Number_Element, 3, 2 )
   
   double precision :: l1,l2,lmid
   double precision :: move_limit, change
   double precision :: PD_Objective_function( Triangle_Number_Element )
   double precision :: beta( Triangle_Number_Element )
   double precision :: Objective_function( g_limit )
   double precision :: x( nd ), y( nd ), A( nd ), B( nd )
   double precision :: det
   double precision :: sk( nd,nd )
   double precision :: shape_function( Triangle_Number_Element, 3 )

   double precision :: psi_n,psi,volume  
   double precision :: nec( Triangle_Number_Element,2 )
   double precision :: kappa( Triangle_Number_Element )       !設計変数から等価変換した熱伝導率
   double precision :: T_ref( Triangle_Number_Element ), T_bare( Triangle_Number_Element )
   double precision :: T_rhs( Number_Node ), rk( Number_Node,Number_Node )
   double precision :: epsi_n( Triangle_Number_Element ), epsi( Triangle_Number_Element )
   double precision :: T_e( Triangle_Number_Element ), T1( g_limit, Triangle_Number_Element )   
   double precision :: E_str( Triangle_Number_Element )
  
  !============================================================================================

   double precision :: Le( Triangle_Number_Element, nd )
   double precision :: Lamda( Triangle_Number_Element, nd )
   double precision :: PD_Lamda( Triangle_Number_Element, 2 )
   double precision :: PD_T( Triangle_Number_Element, 2 )
   
   !===========================================================================================
   !Set number threads  
   !===========================================================================================         
   !$call omp_set_num_threads(28)
 !  call mkl_set_num_threads(28)
      
   !=====================================================
   write(*,*)'**Reading each date**'
   !=====================================================
   
   do i = 1, nelx + 1
      do j = 1, nely + 1
         k = j + ( nely + 1 )*( i - 1 )
         xcoord( k ) = mnx + dble( i - 1 )
      end do
   end do

   do i = 1, nelx + 1
      do j = 1, nely + 1
         k = j + ( nely + 1 )*( i - 1 )
         ycoord( k ) = mny + dble( j - 1 )
      end do
   end do

   do i = 0, nelx - 1
      do j = 1, nely
         m = j + i *nely
         nodex_s( m,1 ) = j + i*( nely + 1 )
	     nodex_s( m,2 ) = j + i*( nely + 1 ) + ( nely + 1 )
	     nodex_s( m,3 ) = j + i*( nely + 1 ) + ( nely + 2 )
	     nodex_s( m,4 ) = j + i*( nely + 1 ) + 1
      end do
!   write(*,*) "loop" ,i
   end do

   do i = 1, Square_Number_Element
	     nodex( i*2,1 ) = nodex_s( i,4 )
	     nodex( i*2,2 ) = nodex_s( i,2 )
	     nodex( i*2,3 ) = nodex_s( i,3 )	
	     nodex( i*2-1,1 ) = nodex_s( i,4 )
	     nodex( i*2-1,2 ) = nodex_s( i,1 )
	     nodex( i*2-1,3 ) = nodex_s( i,2 )
 !        write(*,*) "loop2" ,i
   end do

   do i = 1, nelx 
      do j = 1, nely
         n = j + nely *( i - 1 )
         nesc( n,1 ) = i + ( mnx - 0.5d0 )
         nesc( n,2 ) = j + ( mny - 0.5d0 )
      end do
   end do

   do i = 1, Triangle_Number_Element
      eid( i ) = 5
   end do

   do i = 1, Square_Number_Element
      dis = sqrt( nesc( i,1 )*nesc( i,1 ) + nesc( i,2 )*nesc( i,2 ) )
      if( dis <= 10.5d0 ) then
         eid( i*2 ) = 2
         eid( i*2 - 1 ) = 2
      else if( dis <= 25.5d0 .and.  dis > 10.5d0 ) then
         eid( i*2 ) = 4
         eid( i*2 - 1 ) = 4
      end if
   end do

   Number_EDesignArea = 0
   do i = 1, Triangle_Number_Element
      if( eid( i ) == 4 ) then
         Number_EDesignArea = Number_EDesignArea + 1
      end if 
   end do 
   Number_SEDesignArea = Number_EDesignArea/2

   Number_Eins = 0
   do i = 1, Triangle_Number_Element
      if( eid( i ) == 2 ) then
         Number_Eins = Number_Eins + 1
      end if
   end do
   Number_SEins = Number_Eins/2

   allocate( Element_Design_Area( Number_EDesignArea ))
   allocate( x_beta( Number_EDesignArea ))
   allocate( density_updated( Number_EDesignArea ), diff_density( Number_EDesignArea ))
   allocate( Triangle_density_vec( Number_EDesignArea ))
   allocate( Best_Triangle_density_vec( Number_EDesignArea ))
   
   !=====================================================
   write(*,*)'Initializtion'
   !=====================================================
   
   !Reset Temperature
   T_e( : ) = 0.0d0
   
   !Thermal conductivity value
   kappa( : ) = 1.0d0
   
   SizeA = Number_Node
   
   allocate( rko( SizeA,SizeA ), rhso( SizeA ), ipmkl( SizeA ))
   allocate( rko_id( SizeA ) )
   
   !========================================================================================
   write(*,*)'Set element center'
   !========================================================================================     
   
   e1 = 0
   e2 = 0
   do i = 0, nelx-1
     do j = 1, nely
       e2 = i *160 + j *2
       e1 = e2 - 1
       nec( e2,1 ) = dble( i ) - 59.25d0
       nec( e2,2 ) = dble( j ) - 40.25d0
       nec( e1,1 ) = dble( i ) - 59.75d0
       nec( e1,2 ) = dble( j ) - 40.75d0
     end do
   end do
   
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Reading Temperature without cloak
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   
     !======================================================== 
     write(*,*)'Making Global Matrix'
     !========================================================
     call mkGSM_0 ( nd,Triangle_Number_Element,Number_Node,x,y,A,B,nodex,xcoord,ycoord,sk,rk )
     !======================================================== 
     write(*,*)'Set boundary conditions'
     !========================================================
     call setBC( Triangle_Number_Element,Number_Node,nd,SizeA,xcoord,ycoord,rk,T_rhs,rko_id,rko,rhso,mnx,mxx,mny,mxy )
     !===================================================================================
     write(*,*)'Solving Linear Equation' 
     !===================================================================================
      call dgesv( SizeA, 1, rko, SizeA, ipmkl, rhso, SizeA, info )

     if( info/=0 )then
        write(*,*)'info of dgesv:', info, 'at 154'
        stop
     end if
     !======================================================== 
     write(*,*)'Converts to temperature for the element'
     !========================================================
      call etemp_0( Triangle_Number_Element,Number_Node,nd,SizeA,det,eid, &
                    xcoord,ycoord,nec,x,y,rhso,rko_id,nodex,shape_function,T_rhs,T_e )  

   do i = 1, Triangle_Number_Element      
      T_ref( i ) = T_e( i )
   enddo  
        
   open( 200, file = 'Treference111.txt' )
      do i = 1, Triangle_Number_Element
         write( 200,* ) T_ref( i )
      end do
   close( 200 )

   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   ! Reading Temperature before cloaking
   !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   
   !initialization
!   SizeA = Number_Node - e_Ip

   counter( : ) = 1 
   do i = 1, Triangle_Number_Element
         if( eid( i ) == 2 ) then
            do j = 1, nd
               counter( nodex( i,j ) ) = 0
			enddo
	     endif 
   enddo	
   do i = 1, Triangle_Number_Element   
	     if( eid( i ) /= 2 ) then		
		    do j = 1, nd
			   counter( nodex( i,j ) ) = 1
            enddo
         endif		  
   enddo		 

   SizeA = 0
   do i = 1, Triangle_Number_Element
      if ( counter( i ) == 1 ) then
	     SizeA = SizeA + 1
      endif		 
   enddo

   write(*,*)'SizeA=',SizeA

   T_e( : ) = 0.0d0
  
   call mkGSM( nd,Triangle_Number_Element,Number_Node,kappa,eid,x,y,A,B,nodex,xcoord,ycoord,sk,rk )
   call setBC( Triangle_Number_Element,Number_Node,nd,SizeA,xcoord,ycoord,rk,T_rhs,rko_id,rko,rhso,mnx,mxx,mny,mxy )
   call dgesv( SizeA, 1, rko, SizeA, ipmkl, rhso, SizeA, info )
     if( info/=0 )then
        write(*,*)'info of dgesv:', info, 'at 185'
        stop
     end if
   call etemp( Triangle_Number_Element,Number_Node,nd,SizeA,det,eid, &
                     xcoord,ycoord,nec,x,y,rhso,rko_id,nodex,shape_function,T_rhs,T_e )  

   do i = 1, Triangle_Number_Element
      if ( eid( i ) /= 2 ) then
         T_bare( i ) = T_e( i ) 
      endif  
   enddo 

   !=====================================================
   write(*,*)'**Read Complete**'
   !=====================================================

   !========================================================================================
   !Set element design area
   !========================================================================================     
  
   j = 1
   do i = 1, Triangle_Number_Element
     if( eid( i ) == 4 ) then
        Element_Design_Area( j ) = i
        j = j + 1
     end if
   end do
   
   !========================================================================================
   !Calculate psi_n
   !========================================================================================     

   epsi_n( : ) = 0.0d0

   do i = 1, Triangle_Number_Element
      if( eid( i ) == 5 ) then 
         epsi_n( i ) = abs( T_bare( i ) - T_ref( i ) )**2.0d0
      end if
   end do   
   psi_n  = sum( epsi_n )
   write(*,*) 'psi_n=',psi_n

   !===========================================================================================
   write(*,*)'Initialization'
   !===========================================================================================
  
   !density value
   Triangle_density_vec( : ) = vol_constraint
   
   !psi value
   epsi( : ) = 0.0d0
   
   !parameter
   change = 1.0d0 
   
   !===========================================================================================
   write(*,*)'Optimization Computation'
   !===========================================================================================

   !write(*,*)'aho g=', g, 'change=', change
   g= 0
   do while( change > 1.0d-4 )
      g = g + 1

      write(*,*)'=============================================================='
      write(*,*)'g=', g, 'change=', change
      write(*,*)'=============================================================='

      !========================================================
      write(*,*)'Impose symmetry constraints'
      !========================================================     
 
      do i = 2, Number_EDesignArea, 2
         Triangle_density_vec( i-1 ) = ( Triangle_density_vec( i ) + Triangle_density_vec( i-1 ) )/2.0d0
         Triangle_density_vec( i )   = ( Triangle_density_vec( i ) + Triangle_density_vec( i-1 ) )/2.0d0
      enddo

      !========================================================
      write(*,*) 'Set kappa'
      !========================================================     

      do i = 1, Number_EDesignArea
         j = Element_Design_Area( i )
         !kappa( j ) = ( 1.0d0 + high/thickness ) *Triangle_density_vec( i ) 
         kappa( j ) = 1.0d0 + ( high/thickness ) *Triangle_density_vec( i )**penal 
      end do
  
      !========================================================
      write(*,*) 'Finite Element Analysis'
      !========================================================     

      call mkGSM( nd,Triangle_Number_Element,Number_Node,kappa,eid,x,y,A,B,nodex,xcoord,ycoord,sk,rk )
      call setBC( Triangle_Number_Element,Number_Node,nd,SizeA,xcoord,ycoord,rk,T_rhs,rko_id,rko,rhso,mnx,mxx,mny,mxy )
      call dgesv( SizeA, 1, rko, SizeA, ipmkl, rhso, SizeA, info )
      if( info/=0 )then
         write(*,*)'info of dgesv:', info, 'at 295'
         stop
      end if
      call etemp( Triangle_Number_Element,Number_Node,nd,SizeA,det,eid, &
                  xcoord,ycoord,nec,x,y,rhso,rko_id,nodex,shape_function,T_rhs,T_e )   

      open(400,file='etemperature.txt',status='replace')
         do i = 1, Triangle_Number_Element
            if( eid( i ) /= 2 ) then
               write(400,*) nec( i,1 ), nec( i,2 ), T_e( i )
            end if
         end do
      close(400)

      do i = 1, Triangle_Number_Element
         if( eid( i ) /=2  ) then
            T1( g,i ) = T_e( i )
         end if
      end do
 
      !=========================================================================================
      write(*,*)'Compute Objective Function and Volume'
      !=========================================================================================
     
      do i = 1, Triangle_Number_Element
         if( eid( i ) == 5 ) then
            epsi( i ) = abs( T_e( i ) - T_ref( i ) )**2.0d0
         end if
      end do

      Objective_function( g ) = sum( epsi ) / psi_n
      volume = sum( Triangle_density_vec ) / dble( Number_EDesignArea )

      write(*,*)'===================================================='
      write(*,*)'   Objective_Function=', Objective_function( g )
      write(*,*)'   Volume=', volume
      write(*,*)'===================================================='

      write(*,*) g, Objective_function( g ), volume
      open( 500, file='ObjectiveFunction.dat', position='append' )
         write(500,*) g, Objective_function( g ), volume
      close(500)

      !========================================================================================
      !Preparing for optimization
      !========================================================================================      

      do i = 1, Triangle_Number_Element
         if( eid( i ) /= 2 ) then
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
 
            Le( i,1 ) = sqrt( A( 2 )**2.0d0 + B( 2 )**2.0d0 )
            Le( i,2 ) = sqrt( A( 3 )**2.0d0 + B( 3 )**2.0d0 )
            Le( i,3 ) = sqrt( A( 1 )**2.0d0 + B( 1 )**2.0d0 )
              
            Lamda( i,1 ) = ( Le( i,1 ) + Le( i,3 ) ) / 2.0d0 
            Lamda( i,2 ) = ( Le( i,1 ) + Le( i,2 ) ) / 2.0d0
            Lamda( i,3 ) = ( Le( i,2 ) + Le( i,3 ) ) / 2.0d0 
              
            PD_Lamda( i,1 ) = B(1)*Lamda( i,1 ) + B(2)*Lamda( i,2 ) + B(3)*Lamda( i,3 )
            PD_Lamda( i,2 ) = A(1)*Lamda( i,1 ) + A(2)*Lamda( i,2 ) + A(3)*Lamda( i,3 )
            
            PD_T( i,1 ) = B(1)*T_rhs( nodex( i,1 )) + B(2)*T_rhs( nodex( i,2 )) + B(3)*T_rhs( nodex( i,3 ))
            PD_T( i,2 ) = A(1)*T_rhs( nodex( i,1 )) + A(2)*T_rhs( nodex( i,2 )) + A(3)*T_rhs( nodex( i,3 ))
            
            PD_Objective_function( i ) = PD_Lamda( i,1 )*PD_T( i,1 ) + PD_Lamda( i,2 )*PD_T( i,2 )

         end if 
      end do

      !=========================================================================================
      write(*,*)'Update Density Field'
      !=========================================================================================

      l1 = 0.0d0
      l2 = 1.0d9
      move_limit = 1.0d-2
 
      do while( ( l2-l1 )/( l1+l2 ) > 1.0d-3 )

        lmid= 0.5d0*( l2+l1 )
!write(*,*)'lmid=', lmid, ( l2-l1 )/( l1+l2 )

        do i = 1, Triangle_Number_Element
           if( eid( i ) /= 2 ) then
              beta( i ) = ( abs( PD_Objective_function( i ) ) / lmid )**7.5d-1 
           end if 
        end do           
                   
        do e = 1, Number_EDesignArea
           i = Element_Design_Area( e )
           x_beta( e ) = Triangle_density_vec( e ) *beta( i )
        enddo
 
        do e = 1, Number_EDesignArea
           if( x_beta( e ) <= max( 0.0d0, Triangle_density_vec( e ) -move_limit ) )then
              density_updated( e )= max( 0.0d0, Triangle_density_vec( e ) -move_limit )
           else if( x_beta( e )  >= min( 1.0d0, Triangle_density_vec( e ) -move_limit) )then
              density_updated( e )= min( 1.0d0, Triangle_density_vec( e ) +move_limit )
           else
              density_updated( e )= x_beta( e )
           endif
        enddo

        volume= sum( density_updated )/dble( Number_EDesignArea )
        if( volume > vol_constraint )then
           l1= lmid
        else
           l2= lmid
        endif
   
      enddo !do while((l2-l1)/(l1+l2)>1.0d-3)

      !=========================================================================================
      write(*,*)'Check Convergence' 
      !=========================================================================================

      diff_density( : ) = abs( density_updated( : ) - Triangle_density_vec( : ) )
      change = maxval( diff_density )
      write(*,*)'   change=', change
 
      Triangle_density_vec( : ) = density_updated( : )
 
      !=========================================================================================
      ! Output the best value
      !=========================================================================================
      
      if ( g >= 2 ) then
         if ( Objective_function( g ) < Objective_function( g-1 ) ) then
            Best_Triangle_density_vec( : ) = Triangle_density_vec( : )
            Best_generation = g
         endif
      else
         Best_generation = g
      endif
 
      open(600,file='density.txt',status='replace')
         do i = 1, Number_EDesignArea
            write(600,*) Best_Triangle_density_vec( i )
         end do
      close(600)

      open(700,file='temperature.txt',status='replace')
         do i = 1, Triangle_Number_Element
            if( eid( i ) /=2  ) then
               write(700,*)  T1( Best_generation, i )
            end if
         end do
      close(700)

      !=========================================================================================
      !Make .ps
      !=========================================================================================

      j = 1 
      do i = 1 , Triangle_Number_Element
         if ( eid( i ) == 4 ) then
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
           if ( eid( e ) == 4 ) then
              write(999,*)'newpath'
              write(999,*) 5.0d0*position_in_element(e,1,1) + 300.0d0, 5.0d0*position_in_element(e,1,2) + 200.0d0, 'moveto'
              write(999,*) 5.0d0*position_in_element(e,2,1) + 300.0d0, 5.0d0*position_in_element(e,2,2) + 200.0d0, 'lineto'
              write(999,*) 5.0d0*position_in_element(e,3,1) + 300.0d0, 5.0d0*position_in_element(e,3,2) + 200.0d0, 'lineto'
              write(999,*) 'closepath'
              write(999,*) 1.00 - E_str(e),'setgray'
              write(999,*) 'fill'
            end if
         end do    
      close(999)

      !=========================================================================================
      !convergence judgment
      !=========================================================================================

      !   if( Objective_function( g ) -  Objective_function( g-1 ) <= 1.0d-4 ) stop

      if( g+1 > g_limit ) change = 0.0d0 
   enddo !do while( change > 1.0d-2 ) 

   stop 
 !  deallocate( xcoord )
 !  deallocate( ycoord )
 !  deallocate( T_rhs )
 !  deallocate( rk )
   deallocate( rko )
   deallocate( rhso )
   deallocate( ipmkl )
   deallocate( rko_id )
   deallocate( Element_Design_Area )
   deallocate( x_beta )
   deallocate( density_updated )
   deallocate( diff_density )
   deallocate( Triangle_density_vec )
   deallocate( Best_Triangle_density_vec )
   
end program 

  subroutine mkGSM_0 (nd,Triangle_Number_Element,Number_Node,x,y,A,B,nodex,xcoord,ycoord,sk,rk)
     implicit none
     integer i,j,k,l,ie
     integer nd,Number_Node,Triangle_Number_Element
     integer :: nodex( Triangle_Number_Element,nd )
     double precision :: det
     double precision :: x( nd ), y( nd ), A( nd ), B( nd ), sk( nd,nd )
     double precision :: xcoord( Number_Node ), ycoord( Number_Node ), rk( Number_Node,Number_Node )
     double precision :: kappa( Triangle_Number_Element )

     rk( :,: ) = 0.0d0
 
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
        
         det = A(3)*B(2) - B(3)*A(2)

         do i = 1, nd
            B(i) = B(i) / det
            A(i) = A(i) / det
         end do
     
         do i = 1, nd
            do j = 1, nd
               sk( i,j ) = ( B(i)*B(j) + A(i)*A(j) )*det /2.0d0
            end do
         end do
                    
         do k = 1, nd
            i = nodex( ie,k ) 
            do l = 1, nd
               j = nodex( ie,l )  
               rk( i,j ) = rk( i,j ) + sk( k,l )
            end do
         end do

      end do

      return
   end subroutine mkGSM_0
   
   subroutine etemp_0( Triangle_Number_Element,Number_Node,nd,SizeA,det,eid, &
                      xcoord,ycoord,nec,x,y,rhso,rko_id,nodex,shape_function,T_rhs,T_e )  
      implicit none
      integer :: i,j
      integer :: Triangle_Number_Element, Number_Node, nd, SizeA
      integer :: nodex( Triangle_Number_Element,nd ), eid( Triangle_Number_Element )
      integer :: rko_id( SizeA )
      double precision :: det
      double precision :: x( nd ), y( nd ), rhso( SizeA )
      double precision :: xcoord( Number_Node ), ycoord( Number_Node ), nec( Triangle_Number_Element,2 ) ,T_rhs( Number_Node )
      double precision :: shape_function( Triangle_Number_Element,nd ), T_e( Triangle_Number_Element )
         
      do i = 1, SizeA
         T_rhs( rko_id( i ) ) = rhso( i )
      end do
   
      do i = 1, Triangle_Number_Element
            do j = 1, nd
               x( j ) = xcoord( nodex( i,j ) )
               y( j ) = ycoord( nodex( i,j ) )
            end do     
                
            det = (( x( 2 ) - x( 1 ) )*( y( 3 ) - y( 1 ) )) - (( y( 1 ) - y( 2 ))*( x( 1 ) - x( 3 ) ))
            shape_function( i,1 ) = ( x( 2 )*y( 3 ) - x( 3 )*y( 2 ) ) /det&
                                  + ( nec( i,1 )*( y( 2 ) - y( 3 ) )) /det&
                                  + ( nec( i,2 )*( x( 3 ) - x( 2 ) )) /det 
                                                                  
            shape_function( i,2 ) = ( x( 3 )*y( 1 ) - x( 1 )*y( 3 ) ) /det&
                                  + ( nec( i,1 )*( y( 3 ) - y( 1 ) )) /det&
                                  + ( nec( i,2 )*( x( 1 ) - x( 3 ) )) /det
                                                                  
            shape_function( i,3 ) = ( x( 1 )*y( 2 ) - x( 2 )*y( 1 ) ) /det&
                                  + ( nec( i,1 )*( y( 1 ) - y( 2 ) )) /det&
                                  + ( nec( i,2 )*( x( 2 ) - x( 1 ) )) /det 
             
            T_e( i ) = shape_function( i,1 ) *T_rhs( nodex( i,1 ) )&
                     + shape_function( i,2 ) *T_rhs( nodex( i,2 ) )&
                     + shape_function( i,3 ) *T_rhs( nodex( i,3 ) )
         
      end do             
 
      return
   end subroutine etemp_0

   subroutine mkGSM (nd,Triangle_Number_Element,Number_Node,kappa,eid,x,y,A,B,nodex,xcoord,ycoord,sk,rk)
      implicit none
      integer i,j,k,l,ie
      integer nd,Number_Node,Triangle_Number_Element
      integer :: nodex( Triangle_Number_Element,nd ),eid( Triangle_Number_Element )
      double precision :: det
      double precision :: x( nd ), y( nd ), A( nd ), B( nd ), sk( nd,nd )
      double precision :: xcoord( Number_Node ), ycoord( Number_Node ), rk( Number_Node,Number_Node )
      double precision :: kappa( Triangle_Number_Element )
 
      rk( :,: ) = 0.0d0
  
      do ie = 1, Triangle_Number_Element
         if( eid( ie ) /= 2 ) then
  
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
     
            det = A(3)*B(2) - B(3)*A(2)
            do i = 1, nd
               B(i) = B(i) / det
               A(i) = A(i) / det
            end do
   
            do i = 1, nd
               do j = 1, nd
                  sk( i,j ) = (( B(i)*B(j) + A(i)*A(j) )*( kappa(ie) )*det) /2.0d0
               end do
            end do
            
            do k = 1, nd
               i = nodex( ie,k ) 
               do l = 1, nd
                  j = nodex( ie,l )  
                  rk( i,j ) = rk( i,j ) + sk( k,l )
               end do
            end do
         end if
      end do

      return
   end subroutine mkGSM
   
   subroutine setBC( Triangle_Number_Element,Number_Node,nd,SizeA,xcoord,ycoord,rk,T_rhs,rko_id,rko,rhso,mnx,mxx,mny,mxy )
      implicit none
      integer :: i,j
      integer :: nd, Triangle_Number_Element, Number_Node, SizeA
      integer :: rko_id( SizeA )
      double precision :: mnx, mxx, mny, mxy 
      double precision :: xcoord( Number_Node ), ycoord( Number_Node ), T_rhs( Number_Node )
      double precision :: shape_function( Triangle_Number_Element,nd ), T_e( Triangle_Number_Element )
      double precision :: rko( SizeA,SizeA ), rhso( SizeA ), rk( Number_Node,Number_Node )
   
      !Reset
      T_rhs( : ) = 0.0d0
  
      !======================================================== 
      !Set Norman boundary
      !========================================================
 
      do i = 1, Number_Node
         if( ycoord( i ) == mny .or. ycoord( i ) == mxy ) then
            T_rhs( i ) = 0.0d0
         end if
      end do
  
      !========================================================	
      !Set Dirichlet boundary
      !========================================================
 
      do i = 1, Number_Node
        if( xcoord( i ) == mnx ) then
    
           do j = 1, Number_Node
              rk( i,j ) = 0.0d0
           end do
 
              rk( i,i ) = 1.0d0
              T_rhs( i ) = 0.0d0 
  
        else if( xcoord( i ) == mxx ) then
        
           do j = 1, Number_Node
              rk( i,j ) = 0.0d0
           end do
   
              rk( i,i ) = 1.0d0
              T_rhs( i ) = 1.0d0
   
        end if
      end do
     
      j = 1
      do i = 1, Number_Node
        if( rk( i,i ) /= 0.0d0 ) then
           rko_id( j ) = i
           j = j + 1
        end if
      end do
 
      do i = 1, SizeA
         do j = 1, SizeA
            rko( i,j ) = rk( rko_id( i ), rko_id( j ) )
         end do
         rhso( i ) = T_rhs( rko_id( i ) )
      end do   
 
      return 
   end subroutine setBC
 
   subroutine etemp( Triangle_Number_Element,Number_Node,nd,SizeA,det,eid, &
                      xcoord,ycoord,nec,x,y,rhso,rko_id,nodex,shape_function,T_rhs,T_e )  
      implicit none
      integer :: i,j
      integer :: Triangle_Number_Element, Number_Node, nd, SizeA
      integer :: nodex( Triangle_Number_Element,nd ), eid( Triangle_Number_Element )
      integer :: rko_id( SizeA )
      double precision :: det
      double precision :: x( nd ), y( nd ), rhso( SizeA )
      double precision :: xcoord( Number_Node ), ycoord( Number_Node ), nec( Triangle_Number_Element,2 ) ,T_rhs( Number_Node )
      double precision :: shape_function( Triangle_Number_Element,nd ), T_e( Triangle_Number_Element )
           
      do i = 1, SizeA
         T_rhs( rko_id( i ) ) = rhso( i )
      end do
      
      do i = 1, Triangle_Number_Element
         if( eid( i ) /= 2 ) then
            do j = 1, nd
               x( j ) = xcoord( nodex( i,j ) )
               y( j ) = ycoord( nodex( i,j ) )
            end do     
                   
            det = (( x( 2 ) - x( 1 ) )*( y( 3 ) - y( 1 ) )) - (( y( 1 ) - y( 2 ))*( x( 1 ) - x( 3 ) ))
            shape_function( i,1 ) = ( x( 2 )*y( 3 ) - x( 3 )*y( 2 ) ) /det&
                                  + ( nec( i,1 )*( y( 2 ) - y( 3 ) )) /det&
                                  + ( nec( i,2 )*( x( 3 ) - x( 2 ) )) /det 
                                                                     
            shape_function( i,2 ) = ( x( 3 )*y( 1 ) - x( 1 )*y( 3 ) ) /det&
                                  + ( nec( i,1 )*( y( 3 ) - y( 1 ) )) /det&
                                  + ( nec( i,2 )*( x( 1 ) - x( 3 ) )) /det
                                                                     
            shape_function( i,3 ) = ( x( 1 )*y( 2 ) - x( 2 )*y( 1 ) ) /det&
                                  + ( nec( i,1 )*( y( 1 ) - y( 2 ) )) /det&
                                  + ( nec( i,2 )*( x( 2 ) - x( 1 ) )) /det 
                
            T_e( i ) = shape_function( i,1 ) *T_rhs( nodex( i,1 ) )&
                     + shape_function( i,2 ) *T_rhs( nodex( i,2 ) )&
                     + shape_function( i,3 ) *T_rhs( nodex( i,3 ) )
            
         end if
      end do 

      return 
   end subroutine etemp


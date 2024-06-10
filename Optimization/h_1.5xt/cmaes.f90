
subroutine cmaes&
          ( x_box, ratio_mm, &
            flg_min, type_cma, g_opt, n_dv, popsize, min_val, max_val, fitness )

   use cmaes_subprogram 

   character(len=3), intent(in) :: flg_min
   ! min or max
   character(len=1), intent(in)  :: type_cma 
   ! N : Normal CMA, D : Separable CMA ( Diagonal ), S : Sparse CMA, 
   ! F : Sepctral LSM based CMA ( Fourier series expantion) 

   integer, intent(in) :: g_opt
   integer, intent(in) :: n_dv
   integer, intent(in) :: popsize
   double precision, intent(in) :: min_val, max_val
   double precision, intent(in) :: fitness( popsize )

   double precision, intent(inout) :: x_box( n_dv, popsize) ! -1< x <1
   double precision, intent(out) :: ratio_mm 

   double precision, allocatable :: f_box(:), f_min(:)
   double precision, allocatable :: m_tmp(:), diff_m(:)

   integer i, j, k
   integer size_wm, size_fh, size_work, n_tmp, ILAENV
   integer, parameter :: seed=1

   type(DistributionParameter) :: dp
   type(CMA_Parameter)  :: cma
   type(save_parameter), save :: sv

   !================================================================================
   write(*,*)'   call cmaes : g=', g_opt
   write(*,*)'      type_cma = ',  type_cma
   write(*,*)'      sigma=', cma%sigma 
   !================================================================================ 

   if( popsize==1 )then
      x_box = 0.0d0
      ratio_mm = 0.0d0
      go to 213
   end if

   size_fh=20+(3*n_dv)/max( popsize, 4+3*int( log(dble( n_save) ) ) )
   size_wm=max( n_dv, popsize )

   n_tmp = ILAENV( 1, 'DSYTRD', 'U', n_dv, -1, -1, -1 )
   size_work= ( n_tmp +2 )*n_dv
   !size_work= max( (NB+2)*cma%dim, 3*n_dv-1 )
   !write(*,*)'   3*n_save-1=', 3*n_save-1, 'size_work=', size_work 

   !sv%xmean( 1:n_dv, 2 )=sv%xmean( 1:n_dv, 1 )

   cma%type_cma = type_cma
   cma%lwork = size_work
   write(*,*)'   cma%type_cma = ',  cma%type_cma

   !================================================================================ 
   if( g_opt == 0 )then
   !================================================================================ 

      if( type_cma == 'S' )then
         open( 68, file='dep_dv.dat', status='old' )
            read(68,*) cma%num_dep_dv
            if( cma%num_dep_dv > 27 )then
               write(*,*)'==========================='
               write(*,*)'Error'
               write(*,*)'cma%num_dep_dv=', cma%num_dep_dv, '< 27'
               write(*,*)'cmaes.f90 82'
               write(*,*)'==========================='
               stop
            else
               allocate( dp%dep_dv( cma%num_dep_dv, n_dv ) ) 
               read(68,*) dp%dep_dv 
               !do i= 1, n_dv
               !   do j= 1, cma%num_dep_dv 
               !      read(68,*) dp%dep_dv ( j, i )
               !   end do
               !end do
            end if 
         close( 68 )
         cma%flag_eig = 0
      end if 

      call check_cma_parameters( n_save, n_dv, n_cs_save, popsize, sv%fhistory, size_fh )

      allocate( dp%x_work( n_dv, popsize ) )
      if( type_cma/='D' )then
         allocate( dp%matC( n_dv, n_dv ), dp%matB( n_dv, n_dv ), dp%matInvSqrtC( n_dv, n_dv ) )
         allocate( dp%workmat( n_dv, size_wm), dp%work( size_work ) )
      end if
      allocate( dp%vecD( n_dv ) )
      allocate( dp%weights( popsize/2 ), dp%diagSqrtC( n_dv ) )
      if( type_cma=='A' ) allocate( dp%neg_weights( popsize/2 ) )
      allocate( dp%xmean( n_dv ), dp%ps( n_dv ), dp%pc( n_dv ) )
      allocate( dp%ub( n_dv ), dp%lb( n_dv ), dp%penalty( n_dv ) )
      allocate( dp%fhistory( size_fh ) )

      call setseed( seed )
      call initialize_cmaes( cma, dp, n_dv, min_val, max_val, popsize )
      call generate_candidate_solutions( cma, dp, dp%x_work )
   
      if( min_val < max_val )then 
         call handle_box_boundary(cma, dp, dp%x_work, x_box)
      else
         x_box=dp%x_work
      end if

      call dpcma_2_sv( dp, sv, cma, type_cma, n_dv, popsize, size_fh )

      deallocate( dp%x_work )
      if( type_cma/='D' )then
         deallocate( dp%matC, dp%matB, dp%matInvSqrtC )
         deallocate( dp%workmat, dp%work )
      end if 
      deallocate( dp%weights, dp%diagSqrtC, dp%vecD )
      if( type_cma=='A' ) deallocate( dp%neg_weights )
      deallocate( dp%xmean, dp%ps, dp%pc )
      deallocate( dp%ub, dp%lb, dp%penalty )
      deallocate( dp%fhistory )

      ratio_mm=1.0d8

   !================================================================================
   else if( g_opt >= 1 )then
   !================================================================================

      ! save --> tmp
      if( type_cma/='D' )then
         allocate( dp%matC( n_dv, n_dv ), dp%matB( n_dv, n_dv ), dp%matInvSqrtC( n_dv, n_dv ) )
         allocate( dp%workmat( n_dv, size_wm), dp%work( size_work ) )
      end if

      allocate( dp%vecD( n_dv ) )
      allocate( dp%weights( popsize/2 ), dp%diagSqrtC( n_dv ) )
      if( type_cma=='A' ) allocate( dp%neg_weights( popsize/2 ) )
      allocate( dp%xmean( n_dv ), dp%ps( n_dv ), dp%pc( n_dv ) )
      allocate( dp%ub( n_dv ), dp%lb( n_dv ), dp%penalty( n_dv ) )
      allocate( dp%fhistory( size_work ) )

      allocate( dp%x_work( n_dv, popsize ) )

      call sv_2_dpcma( dp, sv, cma, type_cma, n_dv, popsize, size_fh )

      if( type_cma == 'S' ) allocate( dp%dep_dv( cma%num_dep_dv, n_dv ) ) 
      allocate( f_min(cma%lambda) )
   
      if ( flg_min=='min' ) then
         f_min(:) = fitness(:)
      else if ( flg_min=='max' ) then
         f_min(:) = - fitness(:)
      else 
         write(*,*)'==========================='
         write(*,*)'Error'
         write(*,*)'flg_min=', flg_min
         write(*,*)"flg_min must be 'min' or 'max'. "
         write(*,*)'cmaes.f90 156'
         write(*,*)'==========================='
         stop
      end if
      allocate( f_box(cma%lambda) )

      if( min_val < max_val )then 
         call compute_penalty_v2( cma, dp, x_box, f_min, dp%x_work, f_box )
      else 
         f_box=f_min
      end if

      deallocate( f_min )

      call update_parameters(cma, dp, dp%x_work, f_box)
   
      deallocate( f_box )
   
      call reset_randm( n_dv, popsize, g_opt, seed )
      call generate_candidate_solutions(cma, dp, dp%x_work )
   
      if( min_val < max_val )then 
         call handle_box_boundary(cma, dp, dp%x_work, x_box)
      else
         x_box=dp%x_work
      end if

      call dpcma_2_sv( dp, sv, cma, type_cma, n_dv, popsize, size_fh )

      deallocate( dp%x_work )

      if( type_cma/='D' )then
         deallocate( dp%matC, dp%matB, dp%matInvSqrtC )
         deallocate( dp%workmat, dp%work )
      end if

      deallocate( dp%weights, dp%diagSqrtC, dp%vecD )
      if( type_cma=='A' ) deallocate( dp%neg_weights )
      deallocate( dp%xmean, dp%ps, dp%pc )
      deallocate( dp%ub, dp%lb, dp%penalty )
      deallocate( dp%fhistory )

      allocate( m_tmp( n_dv ), diff_m( n_dv ) )
   
      diff_m( 1:n_dv )= sv%xmean( 1:n_dv, 1 )-sv%xmean( 1:n_dv, 2 )
      m_tmp( 1:n_dv )= sv%xmean( 1:n_dv, 1 )
   
      ratio_mm=sqrt( dot_product( diff_m, diff_m ) )/sqrt( dot_product( m_tmp, m_tmp ) )
   
      deallocate( m_tmp, diff_m )
   end if

   if( type_cma == 'S' ) deallocate( dp%dep_dv ) 

   213 continue
   !================================================================================
   write(*,*)'   end cmaes'
   !================================================================================

end subroutine cmaes
 
subroutine reset_randm( NumParam, popsize, g_opt, seed )
    use cmaes_subprogram

    implicit none

    integer, intent(in) :: NumParam, popsize, g_opt, seed
    integer, parameter :: Limit_Integer= 2147483647

    integer :: Reset_Rand_Generation

    !Reset_Rand_Generation= int( Limit_Integer/( NumParam*popsize*100d0 ) -1 )*100
    !Reset_Rand_Generation= int( Limit_Integer/( NumParam*popsize*1000d0 ) )*1000
    Reset_Rand_Generation= int( Limit_Integer/( NumParam*popsize*1000d0 )-1 )*1000

    write(*,*)'            Reset_Rand_Generation=', Reset_Rand_Generation

    if( g_opt==0 )then
       open( 366, file='./Generation_Reset_Rand.dat', position='append')
          write( 366,*) 'Reset_Rand_Generation=', Reset_Rand_Generation, ' cmaes.f90 : 366'
       close( 366 )
    else if( mod( g_opt,Reset_Rand_Generation )==0 )then
       open( 370, file='./Generation_Reset_Rand.dat', position='append')
          write( 370,*) g_opt, ' cmaes.f90 : 370'
       close( 370 )
       call setseed( seed )
    end if
 
end subroutine reset_randm

subroutine check_cma_parameters( n_save, n_dv, n_cs_save, popsize, fhistory, size_fh )
   implicit none

   integer, intent(in) :: n_save, n_dv, n_cs_save, popsize, size_fh
   double precision, intent(in) :: fhistory( 20+(3*n_save)/( 4+3*int( log(dble( n_save ) ) ) ) ) 

      if( n_save < n_dv )then
         write(*,*)'=========================================================================='
         write(*,*)'Error'
         write(*,*)'n_save=', n_save, ' < ', 'n_dv=', n_dv
         write(*,*)'Modify the vlaue of "n_save" in "cmaes_subprogram.f90"'
         write(*,*)'cmaes.f90 299'
         write(*,*)'=========================================================================='
         stop
      else if( n_cs_save < popsize )then
         write(*,*)'=========================================================================='
         write(*,*)'Error'
         write(*,*)'n_cs_save=', n_cs_save, ' < ', 'popsize=', popsize
         write(*,*)'Modify the vlaue of "n_cs_save" in "cmaes_subprogram.f90"'
         write(*,*)'cmaes.f90 307'
         write(*,*)'=========================================================================='
         stop
      else if( size( fhistory,1) < size_fh )then
         write(*,*)'=========================================================================='
         write(*,*)'Error'
         write(*,*)'size(fhistory,1)=', size( fhistory,1), ' < ', 'size_fh=', size_fh 
         write(*,*)'cmaes.f90 314'
         write(*,*)'=========================================================================='
         stop
      end if

      if ( popsize/=1 .and. popsize < 4 + int(3.0d0 * log( dble(n_dv) ) ) ) then
         write(*,*)'=========================================================================='
         write(*,*)'Error'
         write(*,*)'popsize=', popsize, '>', 4 + int(3.0d0 * log(dble(n_dv)))
         write(*,*)'cmaes.f90 314'
         write(*,*)'=========================================================================='
         stop
      end if
   return
end subroutine check_cma_parameters

subroutine dpcma_2_sv( dp, sv, cma, type_cma, n_dv, popsize, size_fh )

   use cmaes_subprogram 
   implicit none
   type(DistributionParameter), intent(in) :: dp
   type(CMA_Parameter), intent(in) :: cma
   type(save_parameter), intent(out) :: sv
   character(len=1), intent(in)  :: type_cma 
   integer, intent(in)  :: n_dv, popsize, size_fh 

      ! tmp --> save
      sv%x_work( 1:n_dv, 1:popsize )=dp%x_work( 1:n_dv, 1:popsize )
      sv%D( 1:n_dv )=dp%vecD( 1:n_dv )
      if( type_cma/='D' )then
         sv%C( 1:n_dv, 1:n_dv )=dp%matC( 1:n_dv, 1:n_dv )
         sv%B( 1:n_dv, 1:n_dv )=dp%matB( 1:n_dv, 1:n_dv )
         sv%InvSqrtC( 1:n_dv, 1:n_dv )=dp%matInvSqrtC( 1:n_dv, 1:n_dv )
         !sv%workmat( 1:n_dv, 1:size_wm )=dp%workmat( 1:n_dv, 1:size_wm )
         !sv%work( 1:size_work )=dp%work( 1:size_work )
      end if
      sv%weights( 1:popsize/2 )=dp%weights( 1:popsize/2 )
      if( type_cma=='A' ) sv%negweights( 1:popsize/2 )=dp%neg_weights( 1:popsize/2 )
      sv%diagSqrtC( 1:n_dv )=dp%diagSqrtC( 1:n_dv )
      sv%xmean( 1:n_dv, 1 )=dp%xmean( 1:n_dv )
      sv%ps( 1:n_dv )=dp%ps( 1:n_dv )
      sv%pc( 1:n_dv )=dp%pc( 1:n_dv )
      sv%ub( 1:n_dv )=dp%ub( 1:n_dv )
      sv%lb( 1:n_dv )=dp%lb( 1:n_dv )
      sv%penalty( 1:n_dv )=dp%penalty( 1:n_dv )
      sv%fhistory( 1:size_fh )=dp%fhistory( 1:size_fh )
      if( type_cma=='S' ) sv%dep_dv( 1:cma%num_dep_dv, 1:n_dv )=dp%dep_dv( 1:cma%num_dep_dv, 1:n_dv )

      if( type_cma=='S' ) sv%int_cma( 2 )= cma%num_dep_dv
      sv%int_cma( 3 )= cma%dim
      sv%int_cma( 4 )= cma%lambda
      sv%int_cma( 5 )= cma%mu
      sv%int_cma( 6 )= cma%fhistory_size
      sv%int_cma( 7 )= cma%eig_iter
      sv%int_cma( 8 )= cma%idx_fhistory
      sv%int_cma( 9 )= cma%flg_fhistory
      sv%int_cma( 10 )= cma%flg_winit
      sv%int_cma( 11 )= cma%cnt_winit
      sv%int_cma( 12 )= cma%flg_dsyev_err
      sv%int_cma( 13 )= cma%lwork

      sv%dp_cma( 1 )= cma%mueff 
      sv%dp_cma( 2 )= cma%cm
      sv%dp_cma( 3 )= cma%cc 
      sv%dp_cma( 4 )= cma%cone
      sv%dp_cma( 5 )= cma%cmu
      sv%dp_cma( 6 )= cma%cs
      sv%dp_cma( 7 )= cma%ds
      sv%dp_cma( 8 )= cma%dble_dim
      sv%dp_cma( 9 )= cma%chidim
      sv%dp_cma( 10 )= cma%sigma

   return
end subroutine dpcma_2_sv

subroutine sv_2_dpcma( dp, sv, cma, type_cma, n_dv, popsize, size_fh )
   use cmaes_subprogram 
   implicit none

   type(DistributionParameter), intent(inout) :: dp
   type(CMA_Parameter), intent(inout) :: cma
   type(save_parameter), intent(inout) :: sv
   character(len=1), intent(in)  :: type_cma 
   integer, intent(in)  :: n_dv, popsize, size_fh 

   sv%xmean( 1:n_dv, 2 )=sv%xmean( 1:n_dv, 1 )

      dp%x_work( 1:n_dv, 1:popsize ) = sv%x_work( 1:n_dv, 1:popsize )
      dp%vecD( 1:n_dv )=sv%D( 1:n_dv )

      if( type_cma/='D' )then
         dp%matC( 1:n_dv, 1:n_dv )=sv%C( 1:n_dv, 1:n_dv )
         dp%matB( 1:n_dv, 1:n_dv )=sv%B( 1:n_dv, 1:n_dv )
         dp%matInvSqrtC( 1:n_dv, 1:n_dv )=sv%InvSqrtC( 1:n_dv, 1:n_dv )
      end if
      dp%weights( 1:popsize/2 )=sv%weights( 1:popsize/2 )
      if( type_cma=='A' ) dp%neg_weights( 1:popsize/2 )=sv%negweights( 1:popsize/2 )
      dp%diagSqrtC( 1:n_dv )=sv%diagSqrtC( 1:n_dv )
      dp%xmean( 1:n_dv )=sv%xmean( 1:n_dv, 1 )
      dp%ps( 1:n_dv )=sv%ps( 1:n_dv )
      dp%pc( 1:n_dv )=sv%pc( 1:n_dv )
      dp%ub( 1:n_dv )=sv%ub( 1:n_dv )
      dp%lb( 1:n_dv )=sv%lb( 1:n_dv )
      dp%penalty( 1:n_dv )=sv%penalty( 1:n_dv )
      dp%fhistory( 1:size_fh )=sv%fhistory( 1:size_fh )

      if( type_cma=='S' ) cma%num_dep_dv=sv%int_cma( 2 ) 
      cma%dim=sv%int_cma( 3 )
      cma%lambda=sv%int_cma( 4 )
      cma%mu=sv%int_cma( 5 )
      cma%fhistory_size=sv%int_cma( 6 )
      cma%eig_iter=sv%int_cma( 7 )
      cma%idx_fhistory=sv%int_cma( 8 )
      cma%flg_fhistory=sv%int_cma( 9 )
      cma%flg_winit=sv%int_cma( 10 )
      cma%cnt_winit=sv%int_cma( 11 )
      cma%flg_dsyev_err=sv%int_cma( 12 )
      cma%lwork=sv%int_cma( 13 )

      cma%mueff=sv%dp_cma( 1 ) 
      cma%cm=sv%dp_cma( 2 )
      cma%cc=sv%dp_cma( 3 ) 
      cma%cone=sv%dp_cma( 4 )
      cma%cmu=sv%dp_cma( 5 )
      cma%cs=sv%dp_cma( 6 ) 
      cma%ds=sv%dp_cma( 7 )
      cma%dble_dim=sv%dp_cma( 8 )
      cma%chidim=sv%dp_cma( 9 )
      cma%sigma=sv%dp_cma( 10 ) 

   return
end subroutine sv_2_dpcma


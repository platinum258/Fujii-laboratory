! CMA-ES Module
! The Covariance Matrix Adaptation Evolution Strategy (CMA-ES) with
! - Separable CMA for large scale problems [Ros 2008]
! - (Improved) Box Constraint Handling [Hansen 2009]
! - (Improved) Active-CMA [Jastrebski 2006]
!
! The current version can solve box constrained optimization problem of the form
!   Minimize f(x), x = (x_1,...,x_n)
!   Subject to lb <= x_i <= ub for i = 1,...,n
!
! Written and maintained by Youhei Akimoto @ Shinshu University (yohe.aki@gmail.com)
!
! 
!
! BUG fixed (2016/12/26): effective_median_iqr.
! BUG fixed (2016/11/14): inout -> out. Thanks to Prof. Fujii and Mr. Takahashi.


module cmaes_subprogram
  use random, only : setseed, normal
  !use eigen, only : myddot, mydger, mydsyr, mydgemv, mydgemm, mydsyev 
  implicit none
  
  integer, parameter :: n_save=8911 !6452
  integer, parameter :: n_cs_save=280
  integer, parameter :: mu_save=n_cs_save/2

  type CMA_Parameter
     ! Static Parameters
     integer :: dim, lambda, mu, fhistory_size
     !!double precision :: flg_min     
     double precision :: mueff, cm, cc, cone, cmu, cs, ds
     double precision :: dble_dim, chidim
     ! Dynamic Parameters
     integer :: eig_iter
     double precision :: sigma
     ! Box Constraint Parameters
     integer :: idx_fhistory, flg_fhistory, flg_winit, cnt_winit
     ! Termination Condition
     integer :: flg_dsyev_err
     ! Work Space
     integer :: lwork
     
     character(len=1) :: type_cma 
     integer :: num_dep_dv, flag_eig
  end type CMA_Parameter
    
  type DistributionParameter
     double precision, allocatable :: matC(:,:)
     double precision, allocatable :: matB(:, :), matInvSqrtC(:, :)
     double precision, allocatable :: vecD(:)
     double precision, allocatable :: weights(:), neg_weights(:) 
     double precision, allocatable :: diagSqrtC(:)
     double precision, allocatable :: xmean(:), ps(:), pc(:) 
     double precision, allocatable :: ub(:), lb(:), penalty(:)  ! dim-D
     double precision, allocatable :: fhistory(:)
     double precision, allocatable :: work(:) 
     double precision, allocatable :: workmat(:, :), x_work(:,:) 
     integer, allocatable :: dep_dv(:,:)
  end type DistributionParameter

   type save_parameter
      double precision :: C( n_save, n_save ), B( n_save, n_save ), InvSqrtC( n_save, n_save ) 
      double precision :: x_work(n_save, n_cs_save) ! x 
      double precision :: D( n_save ), weights( mu_save ), negweights( mu_save ), diagSqrtC( n_save ) 
      double precision :: xmean( n_save, 2 ), ps( n_save ), pc( n_save ) 
      double precision :: ub( n_save ), lb( n_save ), penalty( n_save ) 
      double precision :: fhistory( 20+(3*n_save)/( 4+3*int( log(dble( n_save) ) ) ) ) 
      !double precision :: sv%workmat( n_save, max( n_save, n_cs_save ) ), sv%work( 3*n_save-1 ) 
      integer :: int_cma( 13 )
      double precision :: dp_cma( 10 ), dep_dv( 27, n_save ) 
   end type save_parameter 

contains

  subroutine initialize_cmaes( cma, dp, num_param, min_val, max_val, popsize )

    implicit none
    
    type(CMA_Parameter), intent(inout) :: cma    
    type(DistributionParameter), intent(inout) :: dp
    integer, intent(in) :: num_param
    double precision, intent(in) :: min_val
    double precision, intent(in) :: max_val
    integer, intent(in) :: popsize       ! default = 4 + |3 * ln(N)|

    integer :: i
    double precision :: neg_mueff, neg_amu

    ! Compute the default static parameters
    cma%dim = num_param
    cma%dble_dim = dble(cma%dim)

    if ( popsize >= 4 + int(3.0d0 * log(cma%dble_dim)) ) then
       cma%lambda = popsize
    else
       cma%lambda = 4 + int(3.0d0 * log(cma%dble_dim))
    endif    
    cma%mu = cma%lambda / 2
    cma%fhistory_size = 20+(3*cma%dim)/cma%lambda

    dp%lb(:) = min_val
    dp%ub(:) = max_val
    dp%penalty(:) = 0.0d0
    !dp%xmean(:) = (dp%ub(:) + dp%lb(:)) / 2.0d0
    dp%xmean(:) = 0.25d0
    !dp%vecD(:) = ((dp%ub(:) - dp%lb(:)) / 2.0d0) / 2.0d0
    dp%vecD(:) = 1d0 !((dp%ub(:) - dp%lb(:)) / 2.0d0) / 2.0d0
    dp%diagSqrtC(:) = dp%vecD
    if ( cma%type_cma == 'S' ) then
       dp%matB(:, :) = 0.0d0
       dp%matC(:, :) = 0.0d0       
       dp%matInvSqrtC(:, :) = 0.0d0
       do i = 1, cma%dim
          dp%matB(i, i) = 1.0d0
          dp%matC(1, i) = 1.0d0 
          dp%matInvSqrtC(i, i) = 1.0d0 
       enddo
    else if( cma%type_cma /= 'D' )then
       dp%matB(:, :) = 0.0d0
       dp%matC(:, :) = 0.0d0       
       dp%matInvSqrtC(:, :) = 0.0d0
       do i = 1, cma%dim
          dp%matB(i, i) = 1.0d0
          dp%matC(i, i) = dp%vecD(i) * dp%vecD(i)
          dp%matInvSqrtC(i, i) = 1.0d0 / dp%vecD(i)
       end do
    end if
    !cma%sigma = 1.0d0
    cma%sigma = 0.3d0*( max_val -min_val )
    dp%pc(:) = 0.0d0
    dp%ps(:) = 0.0d0

    ! Compute the positive weights
    do i = 1, cma%mu
       dp%weights(i) = log(dble(cma%lambda + 1) / 2.0d0) - log(dble(i))
    end do
    dp%weights(:) = dp%weights(:) / sum(dp%weights)
    cma%mueff = sum(dp%weights)**2 / sum(dp%weights * dp%weights)

    ! Learning Rates
    cma%cm = 1.0d0
    cma%cc = (4.0d0 + cma%mueff / cma%dble_dim) / (cma%dble_dim + 4.0d0 + 2.0d0 * cma%mueff / cma%dble_dim)
    cma%cs = (cma%mueff + 2.0d0) / (cma%mueff + cma%dble_dim + 5.0d0)
    cma%ds = 1.0d0 + cma%cs + &
         2.0d0 * max(0.0d0, sqrt((cma%mueff - 1.0d0) / (cma%dble_dim + 1.0d0)) - 1.0d0)

    if ( cma%type_cma == 'S' ) then
       write(*,*)'cma%cone and cma%cmu are set after conting DoF of C at cmaes_subprogram.f90 478 lines'
    else if( cma%type_cma == 'D' )then
       cma%cone = 2.0d0 / ((cma%dble_dim + 1.3d0)**2 + cma%mueff) * (cma%dble_dim + 2.0d0) / 3.0d0
       cma%cmu = min(1.0d0 - cma%cone, &
            2.0d0 * (cma%mueff - 2.0d0 + 1.0d0/cma%mueff) / ((cma%dble_dim + 2.0d0)**2 + cma%mueff) * &
            (cma%dble_dim + 2.0d0) / 3.0d0)
    else if( cma%type_cma == 'N' .or. cma%type_cma == 'A' )then
       cma%cone = 2.0d0 / ((cma%dble_dim + 1.3d0)**2 + cma%mueff)
       cma%cmu = min(1.0d0 - cma%cone, &
            2.0d0 * (cma%mueff - 2.0d0 + 1.0d0/cma%mueff) / ((cma%dble_dim + 2.0d0)**2 + cma%mueff))
    else
       write(*,*)'==========================='
       write(*,*)'Error'
       write(*,*)'cma%type_cma=', cma%type_cma
       write(*,*)'cma%type_cma must be N, A, D, and S'
       write(*,*)'cmaes_subprogram.f90 161'
       write(*,*)'==========================='
       stop
    end if

    ! Compute the negative weights
    if( cma%type_cma == 'A' )then
       do i = 1, cma%mu
          dp%neg_weights(i) = log(dble(cma%lambda + 1) / 2.0d0) - log(dble(cma%lambda + 1 - i))
       end do
       dp%neg_weights(:) = dp%neg_weights(:) / sum(abs(dp%neg_weights))
       neg_mueff = sum(dp%neg_weights)**2 / sum(dp%neg_weights * dp%neg_weights)
       neg_amu = 1.0d0 + min(cma%cone / cma%cmu, 2.0d0 * neg_mueff / (neg_mueff + 2.0d0))
       dp%neg_weights(:) = dp%neg_weights(:) * min(neg_amu, (1.0d0 - cma%cone - cma%cmu) / (cma%dble_dim * cma%cmu))
    endif
    
    ! Others
    cma%flg_winit = 0
    cma%cnt_winit = 3
    cma%idx_fhistory = 1
    cma%flg_fhistory = 0
    cma%eig_iter = 0
    cma%chidim = sqrt(cma%dble_dim) * &
         (1.0d0 - 0.25d0 / cma%dble_dim + 1.0d0 / (21.0d0 * cma%dble_dim * cma%dble_dim))
    cma%flg_dsyev_err = 0

  end subroutine initialize_cmaes

  subroutine randn( outmat )
    implicit none
    double precision, intent(out) :: outmat(:, :)

    integer :: i, j, m, n

    m = size(outmat, 1)
    n = size(outmat, 2)    
    do j = 1, n
       do i = 1, m
          outmat(i, j) = normal(0.0d0, 1.0d0)
       end do
    end do
  end subroutine randn

  subroutine generate_candidate_solutions( cma, dp, outmat )
    implicit none
    
    type(CMA_Parameter), intent(inout) :: cma
    type(DistributionParameter), intent(inout) :: dp
    double precision, intent(out) :: outmat(:, :)

    integer :: i, m, n

    m = size(outmat, 1)
    n = size(outmat, 2)    

    call randn(outmat)
    if( cma%type_cma /= 'D' )then
       
       do i = 1, n
          dp%workmat(:, i) = outmat(:, i) * dp%vecD(:)
          outmat(:, i) = dp%xmean(:)
       end do
       call dgemm( 'N','N',m,n,m,cma%sigma,dp%matB,m,dp%workmat,m,1.0d0,outmat,m)
    else if( cma%type_cma == 'D' )then
       do i = 1, n
          outmat(:, i) = dp%xmean + cma%sigma * dp%vecD(:) * outmat(:, i)
       end do
    end if    
  end subroutine generate_candidate_solutions

  subroutine update_parameters( cma, dp, arx, f_box )
    implicit none

    type(CMA_Parameter), intent(inout) :: cma    
    type(DistributionParameter), intent(inout) :: dp
    double precision, intent(in) :: arx(:, :)
    double precision, intent(in) :: f_box(:)
    
    integer :: i, j, k, l
    integer :: idx(cma%lambda)
    double precision :: hsig, normps, sum_weights, scaled_nw
    double precision :: ymean(cma%dim)
    double precision :: sary(cma%dim, cma%lambda)

    ! Sparse CMA
    integer :: num_nonzero, count_nonzero, nb, dofC
    integer, allocatable, dimension(:) :: irn_tmp, icn_tmp
    integer, allocatable, dimension(:) :: irn, icn
    double precision, allocatable, dimension(:) :: vecC, vecC_tmp
    double precision, allocatable, dimension(:,:) :: matyy, C_full 

    ! FEAST
    integer :: M0,M,info
    double precision :: Emin,Emax
    !!!!!!!!!!!!!!!!! Feast declaration variable
    integer,dimension(64) :: feastparam
    double precision :: epsout
    integer :: loop
    character(len=1) :: UPLO='F'
    double precision,dimension(:),allocatable :: sb
    integer,dimension(:),allocatable :: isa
    double precision,dimension(:),allocatable :: res

    idx = indexsort(f_box)
    do i = 1, cma%lambda
       sary(:, i) = (arx(:, idx(i)) - dp%xmean) / cma%sigma
    end do

    call dgemv('N',cma%dim,cma%mu,1.0d0,sary,cma%dim,dp%weights,1,0.0d0,ymean,1)
    ! ymean(:) = 0.0d0
    ! do i = 1, cma%mu
    !    ymean = ymean + dp%weights(i) * sary(:, i)
    ! end do

    if ( cma%type_cma == 'S' ) then
       dp%ps(:) &
       = (1.0d0 - cma%cs) * dp%ps(:) &
        +sqrt(cma%cs * (2.0d0 - cma%cs) * cma%mueff) * matmul(dp%matInvSqrtC, ymean)
    else if( cma%type_cma /= 'D' )then
       call dgemv('N', cma%dim, cma%dim, sqrt(cma%cs * (2.0d0 - cma%cs) * cma%mueff), &
            dp%matInvSqrtC, cma%dim, ymean, 1, (1.0d0 - cma%cs), dp%ps, 1)
       ! dp%ps = (1.0d0 - cma%cs) * dp%ps + &
       !      sqrt(cma%cs * (2.0d0 - cma%cs) * cma%mueff) * matmul(dp%matInvSqrtC, ymean)
    else if( cma%type_cma == 'D' )then
       dp%ps = (1.0d0 - cma%cs) * dp%ps + &
            sqrt(cma%cs * (2.0d0 - cma%cs) * cma%mueff) * (ymean / dp%vecD)
    end if
       
    normps = sqrt(dot_product(dp%ps, dp%ps))
    if (normps < (1.5d0 + 1.0 / (cma%dble_dim - 0.5)) * cma%chidim) then
       hsig = 1.0d0  !TODO:
    else
       hsig = 0.0d0
    end if
    !dp%pc = (1.0d0 - cma%cc) * dp%pc + &
    !     hsig * sqrt(cma%cc * (2.0d0 - cma%cc) * cma%mueff) * ymean
    dp%pc(:) = (1.0d0 - cma%cc) * dp%pc(:) + &
         hsig * sqrt(cma%cc * (2.0d0 - cma%cc) * cma%mueff) * ymean(:)

    dp%xmean(:) = dp%xmean(:) + (cma%cm * cma%sigma) * ymean(:)
    cma%sigma = cma%sigma * exp((cma%cs / cma%ds) * (normps / cma%chidim - 1.0d0))

    ! CMA
    if( cma%type_cma == 'A' )then
       sum_weights = sum(dp%weights) + sum(dp%neg_weights)
    else if( cma%type_cma /= 'A' )then
       sum_weights = sum(dp%weights)
    endif

    ! Full CMA
    if( cma%type_cma /= 'D' .and. cma%type_cma /= 'S' )then
       dp%matC(:, :) = (1.0d0 - hsig * cma%cone - cma%cmu * sum_weights) * dp%matC(:, :)
       call dsyr('U', cma%dim, (hsig * cma%cone), dp%pc, 1, dp%matC, cma%dim)
       ! do j = 1, cma%dim
       !    dp%matC(:, j) = dp%matC(:, j) + (hsig * cma%cone) * dp%pc(:) * dp%pc(j)
       ! end do
       do i = 1, cma%mu

          ! Negative Update          
          if( cma%type_cma == 'A' )then
             call dgemv('N',cma%dim,cma%dim,1.0d0,dp%matInvSqrtC,cma%dim,&
                  sary(:, cma%lambda + 1 - i),1,0.0d0,dp%work(1:cma%dim),1)
             scaled_nw = dp%neg_weights(i) * cma%dble_dim / dot_product(dp%work(1:cma%dim),dp%work(1:cma%dim))
             !scaled_nw = dp%neg_weights(i) * cma%dble_dim / sum(matmul(dp%matInvSqrtC, sary(:, cma%lambda + 1 - i)) ** 2)
             call dsyr('U', cma%dim, (cma%cmu * scaled_nw), sary(:, cma%lambda + 1 - i), 1, dp%matC, cma%dim)                    
          endif

          ! Positive Update
          call dsyr('U', cma%dim, (cma%cmu * dp%weights(i)), sary(:, i), 1, dp%matC, cma%dim)
          ! do j = 1, cma%dim
          !    dp%matC(:, j) = dp%matC(:, j) + cma%cmu * sary(:, i) * (sary(j, i) * dp%weights(i))
          !    dp%matC(:, j) = dp%matC(:, j) + cma%cmu * sary(:, cma%lambda + 1 - i) * (sary(j, cma%lambda + 1 - i) * scaled_nw)             
          ! end do
       end do

       ! Perform Eigendecomposition every after 1 / (cone + cmu) / (10 * N) iterations
       cma%eig_iter = cma%eig_iter + 1
       if (dble(cma%eig_iter * cma%dim * 10) * (cma%cone + cma%cmu) > 1.0d0) then
          cma%eig_iter = 0  ! initialize the counter

          ! Eigendecomposition of matC.
          ! Note that matC is replaced with its eigenvector matrix.
          call dsyev( 'V', 'U', cma%dim, dp%matC, cma%dim, dp%vecD, &
               dp%work, cma%lwork, cma%flg_dsyev_err)
          dp%matB(:, :) = dp%matC(:, :)
          !call eig(dp%matC, dp%vecD, dp%matB, size(dp%matC, 1), dp%work, cma%lwork)

          ! (Re)compute C = B * D * B'
          do i = 1, cma%dim
             dp%workmat(:, i) = dp%matB(:, i) * dp%vecD(i)
          end do
          call dgemm('N','T',cma%dim,cma%dim,cma%dim,1.0d0,dp%matB,cma%dim,&
               dp%workmat,cma%dim,0.0d0,dp%matC,cma%dim)

          ! D <- sqrt of eigenvalues
          dp%vecD(:) = sqrt(dp%vecD)

          ! Compute C^{-1/2}
          do i = 1, cma%dim
             dp%workmat(:, i) = dp%matB(:, i) / dp%vecD(i)
          end do
          call dgemm('N','T',cma%dim,cma%dim,cma%dim,1.0d0,dp%matB,cma%dim,&
               dp%workmat,cma%dim,0.0d0,dp%matInvSqrtC,cma%dim)
          ! do i = 1, cma%dim
          !    dp%matInvSqrtC(i, :) = dp%matB(:, i) / dp%vecD(i)
          ! end do
          ! dp%matInvSqrtC(:, :) = matmul(dp%matB, dp%matInvSqrtC)

          ! diag(C)^{1/2} for stopping criteria and visualization
          do i = 1, cma%dim
             dp%diagSqrtC(i) = sqrt(dp%matC(i, i))
          end do
       end if

    !=====================================================
    else if( cma%type_cma == 'D' )then ! Separable CMA
    !=====================================================
       dp%vecD(:) = (1.0d0 - hsig * cma%cone - cma%cmu * sum_weights) * dp%vecD(:) * dp%vecD(:)
       dp%vecD(:) = dp%vecD(:) + (hsig * cma%cone) * dp%pc(:) * dp%pc(:)
       do i = 1, cma%mu
          
          ! Negative Update
          if( cma%type_cma == 'A' )then
             scaled_nw = dp%neg_weights(i) * cma%dble_dim / sum((sary(:, cma%lambda + 1 - i) / dp%vecD(:)) ** 2)
             dp%vecD(:) = dp%vecD(:) + cma%cmu * scaled_nw * sary(:, cma%lambda + 1 - i) * sary(:, cma%lambda + 1 - i)
          endif
          ! Positive Update
          dp%vecD(:) = dp%vecD(:) + cma%cmu * dp%weights(i) * sary(:, i) * sary(:, i)
       end do

       ! D <- sqrt(D)
       dp%vecD(:) = sqrt(dp%vecD)
       dp%diagSqrtC(:) = dp%vecD(:)


    !=====================================================
    else if ( cma%type_cma == 'S' ) then ! Sparse CMA
    !=====================================================
       write(*,*)'      Update sparse CMA-ES'

       write(*,*)'         Check dp%dep_dv'
       count_nonzero = 0
       do i = 1, cma%dim
          do j = 1, cma%num_dep_dv
             if( dp%dep_dv( j, i )/=0 )then
                count_nonzero = count_nonzero +1
                k = dp%dep_dv( j, i )
                do l= 1, cma%num_dep_dv
                   if( dp%dep_dv( l, k ) == i )then 
                      !write(441,*) i, dp%dep_dv( j, i ), k, dp%dep_dv( l, k )
                      goto 445 
                   end if
                end do
                write(*,*)'=============================' 
                write(*,*)'error' 
                write(*,*)'=============================' 
                write(*,*)'i, j =', i, dp%dep_dv( j, i ) 
                do l= 1, cma%num_dep_dv
                   if( dp%dep_dv( j, i )/=0 ) write(*,*) dp%dep_dv( l, k ) 
                end do
                write(*,*)'cmaes_subprogram.f90 414' 
                stop

                445 continue
             end if
          end do
       end do

       dofC = ( count_nonzero - cma%dim )/2 +cma%dim 
       write(*,*)'         Reformulate Learning Parameters'
       write(*,*)'            dofC : ', dofC 
       !write(*,*)'            dofC : ', dofC, count_nonzero - cma%dim, count_nonzero, cma%dim*9

       write(*,*)'            cma%cone : ', cma%cone 
       write(*,*)'            cma%cmu : ', cma%cmu 
       write(*,*)'            cma%cc : ', cma%cc 
       cma%cone = min(1.0d0, dble(cma%lambda)/6.0d0) / ( dofC +2.0d0*sqrt(dble(dofC))+cma%mueff/dble(cma%dim) )
       cma%cmu = min(1.0d0 - cma%cone, &
             (0.3d0 +cma%mueff - 2.0d0 + 1.0d0/cma%mueff) / ( dofC+4.0d0**sqrt(dble(dofC))+cma%mueff/2.0d0 ))
       cma%cc = sqrt(cma%cone) 
       write(*,*)'            cma%cone : ', cma%cone 
       write(*,*)'            cma%cmu : ', cma%cmu 
       write(*,*)'            cma%cc : ', cma%cc 

       dp%matC(:, :) = ( 1.0d0 - hsig * cma%cone - cma%cmu * sum_weights ) * dp%matC(:, :)

       write(*,*)'         Rank 1 update'
       do i = 1, cma%dim
          do j = 1, cma%num_dep_dv
             if( dp%dep_dv( j, i )/=0 )then
                dp%matC( j, i ) = dp%matC( j, i ) + (hsig * cma%cone) * dp%pc( dp%dep_dv( j, i  ) )* dp%pc(i)
             end if
          end do
       end do

       write(*,*)'         Rank mu update'
       allocate( matyy( cma%dim, cma%dim ) )
       matyy= 0.0d0

       do i = 1, cma%dim
          do j = 1, cma%dim 
             do k = 1, cma%mu  
                matyy( j, i ) = matyy( j, i ) + sary(j, k)* sary(i, k)*dp%weights(k) 
             end do
          end do
       end do

       do i = 1, cma%dim
          do j = 1, cma%num_dep_dv
             if( dp%dep_dv( j, i )/=0 )then
                dp%matC(j, i) = dp%matC(j, i) + cma%cmu * matyy( dp%dep_dv( j, i ), i ) 
             end if 
          end do
       end do

       deallocate( matyy )

       cma%eig_iter = cma%eig_iter + 1
       if (dble(cma%eig_iter * cma%dim * 10) * (cma%cone + cma%cmu) > 1.0d0) then
          cma%eig_iter = 0  ! initialize the counter

          num_nonzero = cma%num_dep_dv *cma%dim
          write(*,*)'num_nonzero=', num_nonzero
          allocate( irn_tmp( num_nonzero ) )
          allocate( icn_tmp( num_nonzero ) )
          allocate( vecC_tmp( num_nonzero ) )

          !Compressed Sparse Row --> Vector 
          write(*,*)'         Vectorize C'
          irn_tmp = 0
          icn_tmp = 0
          vecC_tmp = 0.0d0
          count_nonzero = 0
          do i = 1, cma%dim
             do j = 1, cma%num_dep_dv
                if( dp%dep_dv( j, i )/=0 )then
                   count_nonzero = count_nonzero +1
 
                   irn_tmp( count_nonzero ) = i 
                   icn_tmp( count_nonzero ) = dp%dep_dv( j, i )

                   vecC_tmp( count_nonzero ) = dp%matC( j, i ) 
                end if
             end do
          end do

          num_nonzero = count_nonzero
          allocate( irn( num_nonzero ) )
          allocate( icn( num_nonzero ) )
          allocate( vecC( num_nonzero ) )

          write(*,*)'         irn_tmp --> irn & icn_tmp --> icn'
          irn( 1:num_nonzero ) = irn_tmp( 1:num_nonzero )
          icn( 1:num_nonzero ) = icn_tmp( 1:num_nonzero )
          vecC( 1:num_nonzero ) = vecC_tmp( 1:num_nonzero )

          deallocate( irn_tmp )
          deallocate( icn_tmp )
          deallocate( vecC_tmp )

          !write(31,*) cma%dim, num_nonzero 
          !write(31,*) '1.0d-2' 
          !write(31,*) num_nonzero
          !do i = 1, num_nonzero
          !   if( irn( i )==icn( i ) )then
          !      write(31,*) irn( i ), icn( i ), vecC( i ), 1.0d0 
          !   else
          !      write(31,*) irn( i ), icn( i ), vecC( i ), 0.0d0 
          !   end if 
          !end do
          !write(*,*)'cmaes_subprogram.f90 554'
          !stop

          ! Eigendecomposition of C.
          call Check_Matrix_EA22(  cma%dim, num_nonzero, irn, icn, vecC )

          write(*,*)'         Eigendecomposition of C'
          !=============================================================================
          if( cma%flag_eig==0 )then
          !=============================================================================
             allocate( C_full( cma%dim,cma%dim ) )
             C_full = 0.0d0
             do i= 1, num_nonzero 
                C_full( irn(i) , icn(i) ) = vecC( i )
             end do
             write(*,*)'   call dsyev'
     
             call dsyev( 'V', 'U', cma%dim, C_full, cma%dim, dp%vecD, dp%work, cma%lwork, cma%flg_dsyev_err )

             if( cma%flg_dsyev_err/=0 )then
                write(*,*)'cma%flg_dsyev_err=', cma%flg_dsyev_err
                write(*,*)'cmaes_subprogram.f90 575'
                stop
             end if
             write(*,*)'   end dsyev'
             dp%matB = C_full
             dp%vecD(:) = sqrt(dp%vecD(:))

             do i = 1, cma%dim
                dp%matInvSqrtC(i, :) = dp%matB(:, i) / dp%vecD(i)
             end do
             dp%matInvSqrtC = matmul(dp%matB, dp%matInvSqrtC)

             deallocate( C_full )

          !=============================================================================
          else if( cma%flag_eig==1 )then
          !=============================================================================
             call EA22( cma%dim, 1.0d-2, num_nonzero, irn, icn, vecC, dp%vecD, dp%matB )
             dp%vecD(:) = sqrt(dp%vecD(:))
          !=============================================================================
          else if( cma%flag_eig==2 )then
          !=============================================================================

             !!! search interval [Emin,Emax] including M eigenpairs
             Emin= 0.00d0
             Emax= 2.00d0 
             M0= cma%dim
           
             !!!!!!!!!!!!! ALLOCATE VARIABLE 
             allocate(res(1:M0))   ! Residual 
  allocate(isa(1:cma%dim+1))
  isa = 0
  isa( 1 )= 1
  do i= 1, num_nonzero
     isa( irn(i)+1)= isa( irn(i) +1 ) +1
  end do
  do i=2, cma%dim+1
     isa(i)=isa(i)+isa(i-1)
  enddo

  allocate(sb(num_nonzero))
  sb = 0.0d0
  ! diagonal element = 1.0d0
  do i= 1, cma%dim 
     sb( isa(i) ) = 1.0d0
  end do

             ! sa = vecC 
             ! jsa = icn

             call feastinit(feastparam)
             feastparam(1)=1
             call dfeast_scsrgv(UPLO,cma%dim,vecC,isa,icn,sb,isa,icn,feastparam,epsout,loop,Emin,Emax,M0,dp%vecD,dp%matB,M,res,info)

             deallocate(res)   ! Residual 
             deallocate(isa) 
             deallocate(sb) 

             dp%vecD(:) = sqrt(dp%vecD(:))

             do i = 1, cma%dim
                dp%matInvSqrtC(i, :) = dp%matB(:, i) / dp%vecD(i)
             end do
             dp%matInvSqrtC = matmul(dp%matB, dp%matInvSqrtC)

             !call FEASTINIT(fea)
          end if

          deallocate( vecC )
          deallocate( irn )
          deallocate( icn )

       end if

    end if

  end subroutine update_parameters

  subroutine handle_box_boundary( cma, dp, arx, x_box)
    implicit none
    !! argument declaration 
    type(CMA_Parameter), intent(in) :: cma    
    type(DistributionParameter), intent(inout) :: dp
    double precision, intent(in) :: arx(:, :)
    double precision, intent(out) :: x_box(:, :)
    !! variable declaration
    integer :: i, j 

    x_box(:, :) = arx(:, :)
    do j = 1, cma%lambda
       do i = 1, cma%dim
          if (arx(i, j) < dp%lb(i)) then
             x_box(i, j) = dp%lb(i)
          else if (arx(i, j) > dp%ub(i)) then
             x_box(i, j) = dp%ub(i)
          end if
       end do
    end do
  end subroutine handle_box_boundary
  
  subroutine compute_penalty( cma, dp, x_box, f_min, arx, f_box)
    !! Box Constraint Handling with Adaptive Penalty Function
    !! Author Version of "A Method for Handling Uncertainty ..." by Niko (2009)
    implicit none
    !! argument declaration 
    type(CMA_Parameter), intent(inout) :: cma
    type(DistributionParameter), intent(inout) :: dp
    double precision, intent(in) :: x_box(:, :)    
    double precision, intent(in) :: f_min(:)    
    double precision, intent(in) :: arx(:, :)
    double precision, intent(out) :: f_box(:)    
    !! variable declaration
    integer :: i, j, k, flg
    double precision :: delta, tmp, delth, dgam, delm
    integer :: fhidx(cma%fhistory_size), fidx(cma%lambda)
    double precision :: xi(cma%dim)

    ! Modified Penalty Box Constraint (based on the author's version of Hansen et al. 2009)
    cma%cnt_winit = max(cma%cnt_winit - 1, 0)
    
    ! Compute the trace of C
    tmp = trace( cma, dp )
        
    ! Update fhistory
    fidx = indexsort(f_min)
    !! CHECK POINT I: Normalized IQR
    dp%fhistory(cma%idx_fhistory) = &
         (f_min(fidx(cma%lambda*3/4)) - f_min(fidx(cma%lambda/4))) / &
         (cma%sigma * cma%sigma * tmp / cma%dble_dim)
    
    cma%idx_fhistory = cma%idx_fhistory + 1
    if (cma%idx_fhistory > cma%fhistory_size) then
       cma%flg_fhistory = 1
       cma%idx_fhistory = 1 !! BUG: fixed on Dec. 7, 2015, (thanks to Mr. Takahashi)
    end if

    ! Check the feasibility of xmean
    flg = 0
    do i = 1, cma%dim
       if (dp%xmean(i) < dp%lb(i)) then
          flg = 1
       else if (dp%xmean(i) > dp%ub(i)) then
          flg = 1
       end if
    end do

    ! Set the penalty coefficient            
    if (flg == 1) then
       if (cma%flg_fhistory == 1) then
          k = size(dp%fhistory)
       else
          k = cma%idx_fhistory - 1
       end if

       !! CHECK POINT II: delta computed only once or twice       
       if ((cma%flg_winit == 0) .or. (cma%cnt_winit > 0)) then
          fhidx(1:k) = indexsort(dp%fhistory(1:k))
          delta = dp%fhistory(fhidx(k/2))

          !! CHECK POINT I: Normalized IQR
          dp%penalty(:) = 2.0d0 * delta          
          cma%flg_winit = 1
       end if
    end if

    !! CHECK POINT III: different update
    ! Increase and decrease the penalty coefficients
    delth = 3.0d0 * max(1.0d0, sqrt(cma%dble_dim) / cma%mueff)
    dgam = min(1.0d0, cma%mueff / (1.0d1 * cma%dble_dim))
    do i = 1, cma%dim
       delm = max(dp%lb(i) - dp%xmean(i), 0.0d0) + max(dp%xmean(i) - dp%ub(i), 0.0d0)
       delm = delm / coordinate_length( cma, dp, i) ! BUG: fixed on Dec. 9, 2016
       dp%penalty(i) = dp%penalty(i) * exp((dgam / 2.0d0) * tanh(max(0.0d0, delm - delth) / 3.0d0))
       if (dp%penalty(i) > 5.0d0 * delta) then
          dp%penalty(i) = dp%penalty(i) * exp(- dgam / 3.0d0)
       end if
    end do

    do i = 1, cma%dim
       xi(i) = 2.0d0 * log(coordinate_length( cma, dp, i))
    end do
    xi(:) = exp(0.9d0 * (xi(:) - sum(xi(:))/cma%dble_dim))
       
    f_box(:) = 0.0d0
    do j = 1, cma%lambda
       do i = 1, cma%dim
          f_box(j) = f_box(j) + (x_box(i, j) - arx(i, j))**2 * dp%penalty(i) / xi(i)
       end do
    end do
    f_box(:) = f_min(:) + f_box(:) / cma%dble_dim
  end subroutine compute_penalty


  subroutine compute_penalty_v1( cma, dp, x_box, f_min, arx, f_box)
    !! Box Constraint Handling with Adaptive Penalty Function
    !! Based on "A Method for Handling Uncertainty ..." by Niko (2009),
    !! but the penalty coefficients are updated every iteration using IQR history.
    implicit none

    type(CMA_Parameter), intent(inout) :: cma
    type(DistributionParameter), intent(inout) :: dp
    double precision, intent(in) :: x_box(:, :)    
    double precision, intent(in) :: f_min(:)    
    double precision, intent(in) :: arx(:, :)
    double precision, intent(out) :: f_box(:)    

    integer :: i, j, k, flg
    double precision :: delta, tmp
    integer :: fhidx(size(dp%fhistory)), fidx(cma%lambda)
    double precision :: xi(cma%dim)

    ! Update fhistory
    fidx = indexsort(f_min)

    !! CHECK POINT I
    dp%fhistory(cma%idx_fhistory) = f_min(fidx(cma%lambda*3/4)) - f_min(fidx(cma%lambda/4))

    cma%idx_fhistory = cma%idx_fhistory + 1
    if (cma%idx_fhistory > size(dp%fhistory)) then
       cma%flg_fhistory = 1
       !cma%idx_fhistory = 0 !! BUG: fixed on Dec. 7, 2015, (thanks to Mr. Takahasi)
       cma%idx_fhistory = 1 
    end if

    flg = 0
    ! Set the penalty coefficient    
    do i = 1, cma%dim
       if (dp%xmean(i) < dp%lb(i)) then
          flg = 1
       else if (dp%xmean(i) > dp%ub(i)) then
          flg = 1
       end if
    end do
    
    if (flg == 1) then
       if (cma%flg_fhistory == 1) then
          k = size(dp%fhistory)
       else
          k = cma%idx_fhistory - 1
       end if

       !! CHECK POINT II
       fhidx(1:k) = indexsort(dp%fhistory(1:k))
       delta = dp%fhistory(fhidx(k/2))

       tmp = trace( cma, dp )

       !! CHECK POINT I
       dp%penalty(:) = 2.0d0 * delta / (cma%sigma * cma%sigma * tmp / cma%dble_dim)
    end if
    
    !! CHECK POINT III: different update
    ! Increase the penalty coefficient
    do i = 1, cma%dim
       if ((dp%xmean(i) < dp%lb(i)) .or. (dp%xmean(i) > dp%ub(i))) then
          if (max(dp%lb(i) - dp%xmean(i), dp%xmean(i) - dp%ub(i)) > &
               3.0d0 * coordinate_length( cma, dp, i) * &
               max(1.0d0, sqrt(cma%dble_dim)/cma%mueff)) then
             dp%penalty(i) = dp%penalty(i) * 1.1d0**min(1.0d0, cma%mueff/dble(10*cma%dim))
          end if
       end if
    end do    

    do i = 1, cma%dim
       xi(i) = 2.0d0 * log(coordinate_length( cma, dp, i))
    end do    
    xi(:) = exp(0.9d0 * (xi(:) - sum(xi(:))/cma%dble_dim))
       
    f_box(:) = 0.0d0
    do j = 1, cma%lambda
       do i = 1, cma%dim
          f_box(j) = f_box(j) + (x_box(i, j) - arx(i, j))**2 * dp%penalty(i) / xi(i)
       end do
    end do
    f_box(:) = f_min(:) + f_box(:) / cma%dble_dim
  end subroutine compute_penalty_v1

  subroutine compute_penalty_v2( cma, dp, x_box, f_min, arx, f_box)
    !! Box Constraint Handling with Adaptive Penalty Function
    !! Based on the author version of "A Method for Handling Uncertainty ..." by Niko (2009)
    !! Two main changes are:
    !! 1. median computation of the IQR history is replaced with the effective median
    !! 2. penalty coefficient is updated if their mean value is greater than three times
    !!    the effective median of the IQR.
    !! These modification improves the performance when the objective function is far from
    !! quadratic function such as an exponential function.
    implicit none
    !! argument declaration 
    type(CMA_Parameter), intent(inout) :: cma
    type(DistributionParameter), intent(inout) :: dp
    double precision, intent(in) :: x_box(:, :)    
    double precision, intent(in) :: f_min(:)    
    double precision, intent(in) :: arx(:, :)
    double precision, intent(out) :: f_box(:)    
    !! variable declaration
    integer :: i, j, k, flg
    double precision :: delta, tmp, delth, dgam, delm, iqr_f_min
    integer :: fhidx(cma%fhistory_size), fidx(cma%lambda)
    double precision :: xi(cma%dim)

    cma%cnt_winit = max(cma%cnt_winit - 1, 0)
    
    ! Compute the trace of C
    tmp = trace( cma, dp )
        
    ! Update fhistory
    fidx = indexsort(f_min)
    iqr_f_min = f_min(fidx(cma%lambda*3/4)) - f_min(fidx(cma%lambda/4))
    !print *, 'check A'
    call append_iqr_hist( cma, dp, iqr_f_min / (cma%sigma * cma%sigma * tmp / cma%dble_dim) )
    ! dp%fhistory(cma%idx_fhistory) = iqr_f_min / (cma%sigma * cma%sigma * tmp / cma%dble_dim)
    ! cma%idx_fhistory = cma%idx_fhistory + 1
    ! if (cma%idx_fhistory > size(dp%fhistory)) then
    !    cma%flg_fhistory = 1
    !    cma%idx_fhistory = 1 !! BUG: fixed on Dec. 7, 2015, (thanks to Mr. Takahashi)
    ! end if

    !print *, 'check B'
    ! Compute the effective median of IQR history
    delta = effective_median_iqr( cma, dp )

    !print *, 'check C'    
    ! Check the feasibility of xmean
    flg = 0
    do i = 1, cma%dim
       if (dp%xmean(i) < dp%lb(i)) then
          flg = 1
       else if (dp%xmean(i) > dp%ub(i)) then
          flg = 1
       end if
    end do

    !print *, 'check D'
    ! Set the penalty coefficient            
    if (flg == 1) then
       k = get_iqr_hist_size( cma )

       if ((cma%flg_winit == 0) .or. (cma%cnt_winit > 0)) then
          ! fhidx(1:k) = indexsort(dp%fhistory(1:k))
          ! delta = dp%fhistory(fhidx(k/2))
          dp%penalty(:) = 2.0d0 * delta          
          cma%flg_winit = 1
       end if
    end if

    !print *, 'check E'
    ! Increase the penalty coefficients
    delth = 3.0d0 * max(1.0d0, sqrt(cma%dble_dim) / cma%mueff)
    dgam = min(1.0d0, cma%mueff / (1.0d1 * cma%dble_dim))
    do i = 1, cma%dim
       delm = max(dp%lb(i) - dp%xmean(i), 0.0d0) + max(dp%xmean(i) - dp%ub(i), 0.0d0)
       delm = delm / coordinate_length( cma, dp, i) ! BUG: fixed on Dec. 9, 2016
       dp%penalty(i) = dp%penalty(i) * exp((dgam / 2.0d0) * tanh(max(0.0d0, delm - delth) / 3.0d0))
    end do

    !print *, 'check F'
    ! Decrease the penalty coefficients
    tmp = sum( dp%penalty(:) ) / cma%dble_dim
    if ( tmp > 3.0d0 * delta ) then
       dp%penalty(:) = dp%penalty(:) * (3.0d0 * delta / tmp)
    end if

    ! Compute the penalty factor
    do i = 1, cma%dim
       xi(i) = 2.0d0 * log(coordinate_length( cma, dp, i))
    end do
    xi(:) = exp(0.9d0 * (xi(:) - sum(xi(:)) / cma%dble_dim))
       
    f_box(:) = 0.0d0
    do j = 1, cma%lambda
       do i = 1, cma%dim
          f_box(j) = f_box(j) + (x_box(i, j) - arx(i, j))**2 * dp%penalty(i) / xi(i)
       end do
    end do
    f_box(:) = f_min(:) + f_box(:) / cma%dble_dim
  end subroutine compute_penalty_v2

  function effective_median_iqr( cma, dp )
    implicit none
    type(CMA_Parameter), intent(inout) :: cma
    type(DistributionParameter), intent(inout) :: dp
    double precision :: effective_median_iqr
    ! variable declaration
    integer :: i, t, k
    double precision :: med, tmp1, tmp2
    double precision :: iqrhist(cma%fhistory_size)
    integer :: idx(cma%fhistory_size)

    !print *, 'check I'
    ! step 0.
    k = get_iqr_hist_size( cma )    
    if ( k == 1 ) then
       effective_median_iqr = get_iqr_hist( cma, dp, 1 )
       return
    else if ( k == 2 ) then
       tmp1 = get_iqr_hist( cma, dp, 1 )
       tmp2 = get_iqr_hist( cma, dp, 2 )          
       effective_median_iqr = (tmp1 + tmp2) / 2.0d0
       return
    end if

    !print *, 'check II'    
    ! step 1. compute the median of latest three iterations
    iqrhist(1) = get_iqr_hist( cma, dp, 1 )
    iqrhist(2) = get_iqr_hist( cma, dp, 2 )
    iqrhist(3) = get_iqr_hist( cma, dp, 3 )    
    if ( iqrhist(1) < iqrhist(2) ) then
       if ( iqrhist(3) < iqrhist(1) ) then
          med = iqrhist(1)
       else if ( iqrhist(2) < iqrhist(3) ) then
          med = iqrhist(2)
       else
          med = iqrhist(3)          
       end if
    else if ( iqrhist(2) < iqrhist(1) ) then
       if ( iqrhist(1) < iqrhist(3) ) then
          med = iqrhist(1)
       else if ( iqrhist(3) < iqrhist(2) ) then
          med = iqrhist(2)
       else
          med = iqrhist(3)                    
       end if
    end if

    !print *, 'check III'
    ! step 2. compute effective iterations
    if ( k == 3 ) then
       effective_median_iqr = med
    else 
       tmp1 = log(med)
       t = k !! BUG: fixed on Dec. 26, 2016, (thanks to Prof. Fujii)
       do i = 4, k
          iqrhist(i) = get_iqr_hist( cma, dp, i )
          tmp2 = log(iqrhist(i))
          if ( abs(tmp2 - tmp1) > log(5.0d0) ) then
             t = i - 1
             exit
          end if
       end do

       ! step 3. comupte
       idx(1:t) = indexsort( iqrhist(1:t) )
       !print *, 'check V'
       !print *, size(idx), t
       if (mod(t, 2) == 1) then
          !print *, size(idx), t          
          effective_median_iqr = iqrhist( idx((t + 1) / 2) )
       else
          !print *, iqrhist( idx(1:t) ), t
          effective_median_iqr = iqrhist( idx(t / 2) ) + iqrhist( idx(t / 2 + 1) )
       end if
    end if
    !print *, 'check V'    
    return    
  end function effective_median_iqr

  subroutine append_iqr_hist( cma, dp, iqr )
    implicit none
    type(CMA_Parameter), intent(inout) :: cma
    type(DistributionParameter), intent(inout) :: dp
    double precision, intent(in) :: iqr

    dp%fhistory(cma%idx_fhistory) = iqr
    ! update the next index
    cma%idx_fhistory = cma%idx_fhistory + 1
    if (cma%idx_fhistory > cma%fhistory_size) then
       cma%flg_fhistory = 1
       cma%idx_fhistory = 1 !! BUG: fixed on Dec. 7, 2015, (thanks to Mr. Takahashi)
    end if
  end subroutine append_iqr_hist

  function get_iqr_hist_size( cma )
    implicit none
    type(CMA_Parameter), intent(in) :: cma
    integer :: get_iqr_hist_size
    
    if (cma%flg_fhistory == 1) then
       get_iqr_hist_size = cma%fhistory_size
    else
       get_iqr_hist_size = cma%idx_fhistory - 1
    end if
    return
  end function get_iqr_hist_size
    
  function get_iqr_hist( cma, dp, idx )
    ! idx >= 1: return the IQR of the (t + 1 - idx)th iteration
    implicit none
    type(CMA_Parameter), intent(in) :: cma
    type(DistributionParameter), intent(inout) :: dp
    integer, intent(in) :: idx
    double precision :: get_iqr_hist
    integer :: i, k

    k = get_iqr_hist_size( cma )
    i = cma%idx_fhistory - idx
    
    if ( idx < 1 .or. idx > k ) then
       ! should not occur
       ! TODO: error
       print *, 'error in get_iqr_hist', idx, k
       get_iqr_hist = -1.0d0
    else if ( i >= 1 ) then
       get_iqr_hist = dp%fhistory( i )
    else
       get_iqr_hist = dp%fhistory( k + i )
    end if

    return
  end function get_iqr_hist
    
  
  logical function check_convergence( cma, dp, tolx, tolf )
    implicit none

    type(CMA_Parameter), intent(inout) :: cma
    type(DistributionParameter), intent(inout) :: dp
    double precision, intent(in) :: tolx, tolf

    integer :: i, k
    double precision :: fdiff
    double precision :: dx(cma%dim)

    check_convergence = .false.
    if (cma%flg_dsyev_err .ne. 0) then
       check_convergence = .true.
       print *, 'dsyer_err_', cma%flg_dsyev_err
    end if
    ! check tolf
    if (cma%flg_fhistory == 1) then
       k = size(dp%fhistory)
    else
       k = cma%idx_fhistory - 1
    end if
    fdiff = maxval(dp%fhistory) * cma%sigma * cma%sigma * trace( cma, dp ) / cma%dble_dim
    if (fdiff <= tolf) then
       check_convergence = .true.
       print *, 'tolf'
    end if
    ! check tolx
    do i = 1, cma%dim
       dx(i) = coordinate_length( cma, dp, i)
    enddo
    if (maxval(dx(:)) <= tolx) then
       check_convergence = .true.
       print *, 'tolx'       
    end if
  end function check_convergence
  
  !----- indexsort: Return Sorted Indeces -----
  function indexsort( farray )
    implicit none

    double precision, intent(in), dimension(:) :: farray

    integer, dimension(size(farray)) :: indexsort
    integer, dimension(size(farray)) :: ranking
    integer i, j
    ranking = 1
    do i = 1, size(farray)-1
       do j = i+1, size(farray)
          if ( farray(i) .le. farray(j) ) then
             ranking(j) = ranking(j) + 1
          else
             ranking(i) = ranking(i) + 1
          end if
       end do
    end do
    do i = 1, size(farray)
       indexsort( ranking(i) ) = i
    end do
    return 
  end function indexsort

  !----- sort: Return Sorted Array -----
  function sort( inarr )
    implicit none
    double precision, intent(in), dimension(:) :: inarr

    integer :: i
    integer, dimension(size(inarr)) :: idx
    double precision, dimension(size(inarr)) :: sort

    idx = indexsort( inarr )
    do i = 1, size( inarr )
       sort(i) = inarr( idx(i) )
    end do
    return
  end function sort

  function trace( cma, dp )
    implicit none

    type(CMA_Parameter), intent(in) :: cma
    type(DistributionParameter), intent(in) :: dp

    integer :: k
    double precision :: trace

    trace = 0.0d0

    if( cma%type_cma == 'S' )then
      do k = 1, cma%dim
          trace = trace + dp%matC(1, k)
       end do
    else if( cma%type_cma /= 'D' )then
       do k = 1, cma%dim
          trace = trace + dp%matC(k, k)
       end do
    else if( cma%type_cma == 'D' )then
       do k = 1, cma%dim
          trace = trace + dp%vecD(k) * dp%vecD(k)
       end do
    end if
    return
  end function trace
      
  function coordinate_length( cma, dp, i )
    implicit none
    
    type(CMA_Parameter), intent(in) :: cma
    type(DistributionParameter), intent(inout) :: dp
    integer, intent(in) :: i
    double precision :: coordinate_length

    if( cma%type_cma == 'S' ) then            
       coordinate_length = cma%sigma * sqrt(dp%matC(1, i))
    else if( cma%type_cma /= 'D' )then
       coordinate_length = cma%sigma * sqrt(dp%matC(i, i))
    else if( cma%type_cma == 'D' )then
       coordinate_length = cma%sigma * dp%vecD(i)
    end if

    return
  end function coordinate_length
    
end module cmaes_subprogram 

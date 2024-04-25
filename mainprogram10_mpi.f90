program thermalcloak
      !$ use omp_lib
      !   use parameters
      !   use Objective_function_subprogram
   use mpi
   implicit none

!============================================================================================
   !解析領域分割サイズ設定
   integer,parameter :: nelx=200
   integer,parameter :: nely=150
   integer,parameter :: r_designarea=50
   integer,parameter :: r_ins=10
   real(8),parameter :: mnx=-(dble(nelx/2))
   real(8),parameter :: mxx=-mnx
   real(8),parameter :: mny=-(dble(nely/2))
   real(8),parameter :: mxy=-mny
!============================================================================================

!============================================================================================
   !fem関連の変数(要素、節点など)の定義
   integer,parameter :: nd=3
   integer,parameter :: Number_Node=(nelx+1)*(nely+1)
   integer Number_Node_for_fem!断熱領域の中にある節点を除いた値
   integer,parameter :: Number_Element=nelx*nely*2
   integer,parameter :: Number_SElement=nelx*nely
   integer Number_EDesignArea,Number_SEDesignArea,Number_Eins,Number_SEins
!============================================================================================

!============================================================================================
   !CMA-ES関連変数
   integer,parameter :: Number_sampling=100
   integer lambda
   integer Optimization_step
   integer Number_Design_Variable
   real(8) convergence_ratio,Convergence_Ratio_dowhile
   real(8),parameter :: Convergence_Error=1.0d-7
   real(8) sigma
   character(len=1), parameter :: Type_CMA = 'N'
!============================================================================================

!============================================================================================
   !最適化関連変数
   integer Best_Objective_function_Loc,Best_Objective_function_Loc_Temp
   integer Number_matrix_design_variable
   real(8),parameter :: t_flat=1.0d0
   real(8) psi_n,psi!目的関数計算用
   real(8) :: objective_function_average,objective_function_pre,Best_Objective_function
!============================================================================================

!============================================================================================
   !その他、プログラム処理用変数
   integer,parameter :: total_process=10!mpi関連
   integer,parameter :: size_flat_nodex=Number_Element*nd!mpi送受信用
   integer :: process_id, num_procs, namelen, rc, ierr, number_firsttask,tasks_per_process!mpi関連
   integer num_threads!スレッド指定用
   integer info!配列計算用
   integer i,i1,i2,j,k,e,e1,e2,d1,d2,d3,d4,lo,lo1,lo2,lo3,lo4,lo5
   integer sample
   integer ir,iw!入出力関連
   character filename*128
!============================================================================================

!============================================================================================
   !fem関連配列
   integer :: nodex(Number_Element,nd),nodex_s(Number_SElement,4),eid(Number_Element)
   integer :: SElement_coordinate(Number_SElement,Number_SElement),Reverse_SElement_coordinate(nelx,nely)
   integer :: Number_Node_for_fem_count(Number_Node)
   integer,allocatable :: g_matrix2_id(:)
   integer,allocatable :: SElement_Design_Area(:),SElement_Design_Area_quarter(:,:)
   real(8) :: xcoord(Number_Node),ycoord(Number_Node),T_rhs(Number_Node),ipmkl_0(Number_Node)
   real(8) :: l_matrix(nd,nd),g_matrix(Number_Node,Number_Node)
   real(8) :: nec(Number_Element,2),nesc(Number_SElement,2)
   real(8),allocatable :: g_matrix2(:,:),rhs2(:),ipmkl(:)
!============================================================================================

!============================================================================================
   !最適化関連配列
   real(8) :: lambda_va(Number_Element,Number_Sampling),Best_lambda_va(Number_Element)!設計変数から等価変換した熱伝導率
   real(8) :: Objective_function(Number_Sampling)
   real(8) :: T_ref(Number_Node),T_bare(Number_Node)
   real(8) :: npsi_n(Number_Node),npsi(Number_Node)
   real(8),allocatable :: design_variable(:,:),density_filtered_vec(:)
   real(8),allocatable :: Vec_D(:),mat_B(:,:),average(:)
   real(8),allocatable :: fitness(:)
!============================================================================================

!============================================================================================
   !mpi送受信用配列
   integer :: flat_nodex(size_flat_nodex)
   integer,allocatable :: flat_SElement_Design_Area_quarter(:)
   real(8),allocatable :: flat_design_variable(:)
!============================================================================================

!============================================================================================
   !結果出力用配列
   integer :: position_in_element(Number_Element,3,2)
   real(8) :: Density_plot(Number_Element)
!============================================================================================   

!============================================================================================
  !start mpi
  call MPI_INIT(ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD,process_id,ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,num_procs,ierr)

  if(process_id==0) then 
   !======================================================
   write(*,*) 'generate finite element'
   !======================================================
    
    call Generate_Finite_Element(nelx,nely,mnx,mny,mxx,mxy,Number_Element,Number_SElement,Number_Node,&
    xcoord,ycoord,nodex,nodex_s,nesc,eid,r_designarea,r_ins)
 
    Number_EDesignArea = 0
     do i = 1,Number_Element
        if( eid( i ) == 4 ) then
           Number_EDesignArea = Number_EDesignArea + 1
        end if 
     end do 
    Number_SEDesignArea = Number_EDesignArea/2
    Number_Design_Variable=Number_SEDesignArea/4
 
   !======================================================
    write(*,*) 'Design Variable=', Number_Design_Variable
    write(*,*) 'generate finite element complete'
   !======================================================
    flat_nodex=reshape(nodex,[size_flat_nodex])
  end if
   !==========================================================================================
   !broadcast data
    call MPI_BCAST(xcoord,Number_Node,MPI_DOUBLE,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(ycoord,Number_Node,MPI_DOUBLE,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(flat_nodex,size_flat_nodex,MPI_INT,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(eid,Number_Element,MPI_INT,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(Number_EDesignArea,1,MPI_INT,0,MPI_COMM_WORLD,ierr) 
    call MPI_BCAST(Number_SEDesignArea,1,MPI_INT,0,MPI_COMM_WORLD,ierr) 
    call MPI_BCAST(Number_Design_Variable,1,MPI_INT,0,MPI_COMM_WORLD,ierr) 
   !==========================================================================================

!============================================================================================= 
  
!============================================================================================= 
  if(process_id/=0) then
    nodex=reshape(flat_nodex,[Number_Element,nd])
  end if
!=============================================================================================

!=============================================================================================  
  call MPI_Barrier(MPI_COMM_WORLD,ierr)
!=============================================================================================

    if(Number_Sampling>=4+int(3.0d0*log(dble(Number_Element)))) then
      lambda=Number_Sampling
    else
      lambda=4+int(3.0d0*log(dble(Number_Element)))
    end if
 
    allocate(design_variable(Number_Design_Variable,lambda))
    allocate(fitness(lambda))
    allocate(SElement_Design_Area(Number_SEDesignArea),SElement_Design_Area_quarter(Number_SEDesignArea/4,4))
    allocate(flat_SElement_Design_Area_quarter(Number_SEDesignArea))   
    allocate(density_filtered_vec(Number_Design_Variable))
    allocate(average(Number_Design_Variable))
    allocate(Vec_D(Number_Design_Variable),mat_B(Number_Design_Variable,Number_Design_Variable))
 
   !============================================================================================
   !count Number_Node_for_fem
  if(process_id==0) then 
    do i=1,Number_Node
      Number_Node_for_fem_count(i)=0
    end do
 
    do j=1,nd
      do i=1,Number_Element
        if(eid(i)/=2)then
         Number_Node_for_fem_count(nodex(i,j))=1
        end if
      end do
    end do
 
    Number_Node_for_fem=0
    do i=1,Number_Node
      if(Number_Node_for_fem_count(i)==1)then 
        Number_Node_for_fem=Number_Node_for_fem+1
      end if
    end do
 
   !======================================================
    write(*,*)"Number_Node_for_fem=",Number_Node_for_fem
   !======================================================
  end if

!=============================================================================================
!broadcast data
call MPI_BCAST(Number_Node_for_fem,1,MPI_INT,0,MPI_COMM_WORLD,ierr) 
!=============================================================================================

!=============================================================================================
call MPI_Barrier(MPI_COMM_WORLD,ierr)
!============================================================================================

   allocate(g_matrix2(Number_Node_for_fem,Number_Node_for_fem),rhs2(Number_Node_for_fem),ipmkl(Number_Node_for_fem))
   allocate(g_matrix2_id(Number_Node_for_fem))

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Cauculate T-reference
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  if(process_id==0) then 
   !======================================================== 
    write(*,*)'Making Global Matrix'
   !========================================================
    call mkGSM_0 ( nd,Number_Element,Number_Node,nodex,xcoord,ycoord,g_matrix)
   !======================================================== 
    write(*,*)'Set boundary conditions'
   !========================================================
    call setBC_0(Number_Element,Number_Node,nd,xcoord,ycoord,g_matrix,T_ref,mnx,mxx,mny,mxy) 
   !===================================================================================
    write(*,*)'Solving Linear Equation' 
   !===================================================================================
    call dgesv(Number_Node,1,g_matrix,Number_Node,ipmkl_0,T_ref,Number_Node,info)
   
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Cauculate T-bare
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   !===================================================================================
    write(*,*)'Making Global Matrix'
   !===================================================================================
    call mkGSM_1(ND,Number_Element,Number_Node,eid,l_matrix,NODEX,XCOORD,YCOORD,g_matrix)
   !=================================================================================== 
    write(*,*)'Set boundary conditions'
   !===================================================================================
    call setBC(Number_Element,Number_Node,nd,Number_Node_for_fem,xcoord,ycoord,g_matrix,T_rhs,mnx,mxx,mny,mxy )
   !===================================================================================
    j = 1
    do i = 1, Number_Node
      if( g_matrix( i,i ) /= 0.0d0 ) then
         g_matrix2_id(j) = i
         j=j+1
      end if
    end do
        
    do j = 1, Number_Node_for_fem
       do i = 1, Number_Node_for_fem
          g_matrix2( i,j ) = g_matrix( g_matrix2_id( i ), g_matrix2_id( j ) )
       end do
       rhs2( j ) = T_rhs( g_matrix2_id( j ) )
    end do   
   !===================================================================================
    write(*,*)'Solving Linear Equation' 
   !===================================================================================
    call dgesv( Number_Node_for_fem,1,g_matrix2,Number_Node_for_fem,ipmkl,rhs2,Number_Node_for_fem,info)

    T_bare(:)=0.0d0

    do i=1,Number_Node_for_fem
     T_bare(g_matrix2_id(i))=rhs2(i)
    end do
   !======================================================== 
    write(*,*)'T-bare calculated'
   !======================================================== 
  
  !============================================================================================
   !対称化
    do i=1,Number_SElement
       SElement_coordinate(i,1)=int(xcoord(nodex_s(i,3)))
       SElement_coordinate(i,2)=int(ycoord(nodex_s(i,3)))
    end do   
   
    Reverse_SElement_coordinate(:,:)=0
 
    !set reverse_SElement_coordinate(座標から番号，縦横が1から始まるように値を調整)
    do i=1,Number_SElement
       i1=SElement_coordinate(i,1)+nelx/2
       i2=SElement_coordinate(i,2)+nely/2
       Reverse_SElement_coordinate(i1,i2)=i 
    end do 
 
    !create quarter of square element
    lo=1
    do i=nelx/2+1,nelx
      do j=nely/2+1,nely
         lo1=reverse_SElement_coordinate(i,j) ! 1st quadrant
         lo2=reverse_SElement_coordinate(nelx-i+1,j)!2nd quadrant
         lo3=reverse_SElement_coordinate(nelx-i+1,nely-j+1)!3rd quadrant
         lo4=reverse_SElement_coordinate(i,nely-j+1)!4th quadrant
       if(eid(lo1*2)==4)then
         SElement_Design_Area_quarter(lo,1)=lo1
         SElement_Design_Area_quarter(lo,2)=lo2
         SElement_Design_Area_quarter(lo,3)=lo3
         SElement_Design_Area_quarter(lo,4)=lo4
         lo=lo+1
       end if
      end do
    end do   

!============================================================================================
   !calculate psi_n
    do i=1,Number_Node
       npsi_n(i)=0.0d0
    end do

    do i=1,Number_Node
       npsi_n(i)=(abs(T_bare(i)-T_ref(i)))**2.0d0
    end do   

    psi_n=sum(npsi_n)

    flat_SElement_Design_Area_quarter=reshape(SElement_Design_Area_quarter,[Number_SEDesignArea])
  end if
!=============================================================================================

!==========================================================================================
!broadcast data
 call MPI_BCAST(T_ref,Number_Node,MPI_DOUBLE,0,MPI_COMM_WORLD,ierr) 
 call MPI_BCAST(g_matrix2_id,Number_Node_for_fem,MPI_INT,0,MPI_COMM_WORLD,ierr) 
 call MPI_BCAST(T_bare,Number_Node,MPI_DOUBLE,0,MPI_COMM_WORLD,ierr)
 call MPI_BCAST(flat_SElement_Design_Area_quarter,Number_SEDesignArea,MPI_INT,0,MPI_COMM_WORLD,ierr)
 call MPI_BCAST(psi_n,1,MPI_DOUBLE,0,MPI_COMM_WORLD,ierr) 
!==========================================================================================

!=============================================================================================  
  call MPI_Barrier(MPI_COMM_WORLD,ierr)
!=============================================================================================

!=============================================================================================
  if(process_id/=0) then
    SElement_Design_Area_quarter=reshape(flat_SElement_Design_Area_quarter,[Number_SEDesignArea/4,4])
  end if
!=============================================================================================

!============================================================================================
   !initialize
    Convergence_Ratio=1.0d0
    Convergence_Ratio_dowhile=1.0d0
    Optimization_step=0
    objective_function_average=0.0d0
    objective_function_pre=0.0d0
    Best_Objective_function=1.0d0
    lambda_va(:,:)=1.0d0
    rhs2(:)=0.0d0!clear rhs2
    Number_matrix_design_variable=Number_Design_Variable*Number_Sampling   

    allocate(flat_design_variable(Number_matrix_design_variable))

    !call omp_set_num_threads(1)
    !call mkl_set_num_threads(1)

!============================================================================================
call MPI_Barrier(MPI_COMM_WORLD,ierr)
!============================================================================================

!============================================================================================
  !CMA-ES      
100 if(process_id==0) then
      call cmaes(design_variable,Convergence_Ratio,&
                 'min', Type_CMA,Optimization_Step,Number_Design_Variable, lambda, 0.0d0, 1.0d0, fitness&
                ,average,Vec_D,mat_B,sigma)!, flg_min, type_cma) 
     
      write(*,*)'convergence ratio=', convergence_ratio
  
      if(Optimization_step>=1) then
        Convergence_Ratio_dowhile=abs(objective_function_average-objective_function_pre)
      end if
  
      objective_function_pre=objective_function_average
       
      write(*,*)'convergence error=', Convergence_Ratio_dowhile

      flat_design_variable=reshape(design_variable,[Number_matrix_design_variable])

    end if
!============================================================================================

!=====================================================================================================
!broadcast data
 call MPI_BCAST(objective_function_average,1,MPI_double,0,MPI_COMM_WORLD,ierr) 
 call MPI_BCAST(psi_n,1,MPI_double,0,MPI_COMM_WORLD,ierr)  
 call MPI_BCAST(flat_design_variable,Number_matrix_design_variable,MPI_double,0,MPI_COMM_WORLD,ierr) 
!=====================================================================================================

!============================================================================================
  if(process_id/=0) then
    design_variable=reshape(flat_design_variable,[Number_Design_Variable,Number_Sampling])
  end if
!=============================================================================================

!============================================================================================
call MPI_Barrier(MPI_COMM_WORLD,ierr)
!============================================================================================

!=============================================================================================
  !start fem
   !============================================================================================
    call mpi_proc_allocate(total_process,lambda,process_id,number_firsttask,tasks_per_process)
   !============================================================================================

    do sample=number_firsttask,number_firsttask+tasks_per_process!while(change>1.0d-2)
          do i=1,Number_Design_Variable
            density_filtered_vec(i)=design_variable(i,sample)
          end do
  
       !set lambda_va
       !do i=1,Number_Element
       !   lambda_va(i,sample)=1.0d0
       !end do
  
       do i=1,Number_Design_Variable
         do j=1,4
          lo5=SElement_Design_Area_quarter(i,j)
          lambda_va(lo5*2,sample)=1.0d0+(density_filtered_vec(i))
          lambda_va(lo5*2-1,sample)=1.0d0+(density_filtered_vec(i))
         end do 
       end do
         
     !===================================================================================
     write(*,*)'Making Global Matrix'
     !===================================================================================
     call mkGSM(ND,Number_Element,Number_Node,lambda,lambda_va,eid,l_matrix,NODEX,XCOORD,YCOORD,g_matrix,sample)
     !=================================================================================== 
     write(*,*)'Set boundary conditions'
     !===================================================================================
     call setBC(Number_Element,Number_Node,nd,Number_Node_for_fem,xcoord,ycoord,g_matrix,T_rhs,mnx,mxx,mny,mxy )
     !===================================================================================
     write(*,*)'change g_matrix for Solving Linear Equation' 
     !===================================================================================
     
     !j = 1
     !do i = 1, Number_Node
     !  if( g_matrix( i,i ) /= 0.0d0 ) then
     !     g_matrix2_id(j) = i
     !     j=j+1
     !  end if
     !end do
 
     do j = 1, Number_Node_for_fem
        do i = 1, Number_Node_for_fem
           g_matrix2( i,j ) = g_matrix( g_matrix2_id( i ), g_matrix2_id( j ) )
        end do
        rhs2( j ) = T_rhs( g_matrix2_id( j ) )
     end do   
  
     !===================================================================================
     write(*,*)'Solving Linear Equation' 
     !===================================================================================
     !call mkl_set_num_threads(1)
     call dgesv( Number_Node_for_fem,1,g_matrix2,Number_Node_for_fem,ipmkl,rhs2,Number_Node_for_fem,info)
      
      do i=1,Number_Node_for_fem
       T_rhs(g_matrix2_id(i))=rhs2(i)
      end do
  
    !============================================================================================
     !compute objective function
      do j=1,Number_Element
         npsi(j)=0.0d0
      end do
      
      do i=1,Number_Node
         npsi(i)=(abs(T_rhs(i)-T_ref(i)))**2.0d0
      end do
  
      psi=sum(npsi)/psi_n

      if(process_id/=0) then
         call MPI_Send(psi,1,MPI_INTEGER,0,0,MPI_COMM_WORLD,ierr)
         call MPI_Send(sample,1,MPI_INTEGER,0,0,MPI_COMM_WORLD,ierr)      
      end if

      if(process_id==0) then
         Objective_function(sample)=psi
         fitness(sample)=Objective_function(sample)
         do i=1,total_process-1
           call MPI_Recv(psi,1,MPI_DOUBLE,i,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
           call MPI_Recv(sample,1,MPI_INTEGER,i,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
           Objective_function(sample)=psi
           fitness(sample)=Objective_function(sample)
           write(*,*) 'Sample=',sample  
         end do
      end if
  
      !============================================================================================
      call MPI_Barrier(MPI_COMM_WORLD,ierr)
      !============================================================================================
    end do
 
!============================================================================================
    !find best result&output
    if(process_id==0) then 
      iw=102
       open(iw,file='temperature.txt',status='replace')
        do i=1,Number_Element
          if(eid(i)/=2) then
           do j=1,nd
             write(iw,*) xcoord(nodex(i,j)), ycoord(nodex(i,j)), T_rhs(nodex(i,j))
           end do
          end if
        end do
       close(iw)  
       
      iw=103
       open(iw,file='Convergence_ratio.dat',position='append')
         write(iw,*) Optimization_step,Convergence_ratio
       close(iw)
 
      objective_function_average=0.0d0
 
      do i=1,Number_Sampling
         objective_function_average=objective_function_average+objective_function(i)/100.0d0
      end do
 
      Best_Objective_function_Loc_Temp=minloc(objective_function,1)
 
      do i=1,Number_Element
        Density_plot(i)=lambda_va(i,Best_Objective_function_Loc_Temp)-1.0d0
      end do
      
      do e=1,Number_Element
         do i=1,3
            position_in_element(e,i,1)=xcoord(nodex(e,i))
            position_in_element(e,i,2)=ycoord(nodex(e,i))
        end do
     end do
 
     if( minval(Objective_function)<Best_Objective_function ) then 
       Best_Objective_function=minval(Objective_function)
       Best_Objective_function_Loc=minloc(Objective_function,1)
       !===================================================================================
        write(*,*)'Best result update'
       !===================================================================================
       iw=107
       open( iw,file='Best_Value.txt',status='replace' )
 
          write(iw,*) 'best result is in generation',Optimization_step,'sampling=',Best_Objective_function_Loc
          write(iw,*) 'best vale of Objective function is',Best_Objective_function
          write(iw,*) 'best design variables'
  
          do i=1,Number_Design_Variable
            write(iw,*) design_variable(i,Best_Objective_function_Loc)
          end do
       close(iw)  
 
       do i=1,Number_Element
          Best_lambda_va(i)=lambda_va(i,Best_Objective_function_Loc)    
       end do
       
       iw=108
       open(iw,file='Best_structure.ps',status='replace')
          write(iw,*)'%!PS-Adobe-3.0 EPSF-3.0'
          write(iw,*)'%%BoundingBox:',0,0,600,450
          write(iw,*)'/m {moveto} def'
          write(iw,*)'/l {moveto} def'
          write(iw,*)'1.0 setgray'
          write(iw,*)'gsave'
          write(iw,*)'0 0 moveto'
          write(iw,*)'600 0 lineto'
          write(iw,*)'600 450 lineto'
          write(iw,*)'0 450 lineto'
          write(iw,*)'closepath'
          write(iw,*)'fill'
          write(iw,*)'grestore'
          write(iw,*)'newpath'
          write(iw,*)'0.0 setgray'
           do e=1,Number_Element
             if(eid(e)==4) then
               write(iw,*) (600.0d0/dble(nelx))*position_in_element(e,1,1)+300.0d0,&
               (600.0d0/dble(nely))*position_in_element(e,1,2)+225.0d0, 'moveto'
               write(iw,*) (600.0d0/dble(nelx))*position_in_element(e,2,1)+300.0d0,&
               (600.0d0/dble(nely))*position_in_element(e,2,2)+225.0d0, 'lineto'
               write(iw,*) (600.0d0/dble(nelx))*position_in_element(e,3,1)+300.0d0,&
               (600.0d0/dble(nely))*position_in_element(e,3,2)+225.0d0, 'lineto'
               write(iw,*) 'closepath'
               write(iw,*) 1.0d0-Best_lambda_va(e),'setgray'
               write(iw,*) 'fill'
               write(iw,*) ''
             end if
           end do    
       close(iw)
     end if

    !===================================================================================
     write(*,*)'Result update'
    !===================================================================================
 
     write(filename,'(i5.5,"_Structure.ps")') Optimization_step
 
     write(*,*)'output.ps >>>>>' ,filename
     iw=104
     open(iw,file=filename,status='replace')
        write(iw,*)'%!PS-Adobe-3.0 EPSF-3.0'
        write(iw,*)'%%BoundingBox:',0,0,600,450
        write(iw,*)'/m {moveto} def'
        write(iw,*)'/l {moveto} def'
        write(iw,*)'1.0 setgray'
        write(iw,*)'gsave'
        write(iw,*)'0 0 moveto'
        write(iw,*)'600 0 lineto'
        write(iw,*)'600 450 lineto'
        write(iw,*)'0 450 lineto'
        write(iw,*)'closepath'
        write(iw,*)'fill'
        write(iw,*)'grestore'
        write(iw,*)'newpath'
        write(iw,*)'0.0 setgray'
         do e=1,Number_Element
           if(eid(e)==4) then
             write(iw,*) (600.0d0/dble(nelx))*position_in_element(e,1,1)+300.0d0,&
             (600.0d0/dble(nely))*position_in_element(e,1,2)+225.0d0, 'moveto'
             write(iw,*) (600.0d0/dble(nelx))*position_in_element(e,2,1)+300.0d0,&
             (600.0d0/dble(nely))*position_in_element(e,2,2)+225.0d0, 'lineto'
             write(iw,*) (600.0d0/dble(nelx))*position_in_element(e,3,1)+300.0d0,&
             (600.0d0/dble(nely))*position_in_element(e,3,2)+225.0d0, 'lineto'
             write(iw,*) 'closepath'
             write(iw,*) 1.0d0-(e),'setgray'
             write(iw,*) 'fill'
             write(iw,*) ''
           end if
         end do    
     close(iw)
 
     Optimization_step=Optimization_step+1
 
     write(*,*)'step',Optimization_step
    end if
!============================================================================================
call MPI_Barrier(MPI_COMM_WORLD,ierr)
!============================================================================================ 

!===================================================================================
   !収束判定
   if(process_id==0) then
      if(Convergence_Ratio_dowhile>=Convergence_Error)then
         goto 100
      end if     
   end if

!===================================================================================
 write(*,*)'Optimization Complete'
!===================================================================================

 deallocate(SElement_Design_Area,SElement_Design_Area_quarter)
 deallocate(g_matrix2_id)
 deallocate(g_matrix2,rhs2,ipmkl)
 deallocate(Design_Variable)
 deallocate(density_filtered_vec)
 deallocate(fitness)
 deallocate(Vec_D)
 deallocate(mat_B)
 deallocate(average)
 deallocate(flat_SElement_Design_Area_quarter)
 deallocate(flat_design_variable)

end program 

  subroutine Generate_Finite_Element(nelx,nely,mnx,mny,mxx,mxy,Number_Element,Number_SElement,Number_Node,&
  xcoord,ycoord,nodex,nodex_s,nesc,eid,r_designarea,r_ins)
  implicit none
  integer nelx,nely,Number_Element,Number_SElement,Number_Node,r_designarea,r_ins
  integer i,j,k,l,m,n,e1,e2
  integer :: nodex(Number_Element,3),nodex_s(Number_SElement,4),eid(Number_Element)
  real(8) mnx,mny,mxx,mxy
  real(8) dis
  real(8) :: xcoord(Number_Node),ycoord(Number_Node),nec(Number_Element,2),nesc(Number_SElement,2)

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
    end do
    
    do i = 1, Number_SElement
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
          nesc( n,1 ) = dble(i) + ( mnx - 0.5d0 )
          nesc( n,2 ) = dble(j) + ( mny - 0.5d0 )
       end do
    end do
    
    do i=1,nelx
      do j=1,nely
        e2=(i-1)*nely*2+j*2
        e1=e2-1
        nec(e2,1)=dble(i)+ (mnx-0.25d0)
        nec(e2,2)=dble(j)+ (mny-0.25d0)
        nec(e1,1)=dble(i)+ (mnx-0.75d0)
        nec(e1,2)=dble(j)+ (mny-0.75d0)
      end do
    end do
    
    do i = 1, Number_Element
       eid( i ) = 5
    end do
    
    do i = 1, Number_SElement
       dis = sqrt( nesc( i,1 )*nesc( i,1 ) + nesc( i,2 )*nesc( i,2 ) )
       if( dis <= dble(r_ins) ) then
          eid( i*2 ) = 2
          eid( i*2 - 1 ) = 2
       else if( dis <= dble(r_designarea) .and. dis > dble(r_ins) ) then
          eid( i*2 ) = 4
          eid( i*2 - 1 ) = 4
       end if
    end do
    return
  end subroutine Generate_Finite_Element

  subroutine mkGSM_0 (nd,Number_Element,Number_Node,nodex,xcoord,ycoord,g_matrix)
    implicit none
    integer i,j,k,l,ie
    integer nd,Number_Node,Number_Element
    integer :: nodex( Number_Element,nd )
    real(8) :: det
    real(8) :: x( nd ), y( nd ), A( nd ), B( nd ), l_matrix( nd,nd )
    real(8) :: xcoord( Number_Node ), ycoord( Number_Node ), g_matrix( Number_Node,Number_Node )

    g_matrix( :,: ) = 0.0d0
  
    do ie = 1, Number_Element
  
        do i = 1, nd
           x(i)=xcoord(nodex(ie,i))
           y(i)=ycoord(nodex(ie,i))
        end do
          
        B(1) = y(2)-y(3)
        B(2) = y(3)-y(1)
        B(3) = y(1)-y(2)
        A(1) = x(3)-x(2)
        A(2) = x(1)-x(3)
        A(3) = x(2)-x(1)
       
        det = A(3)*B(2) - B(3)*A(2)      
        do i = 1, nd
           B(i) = B(i) / det
           A(i) = A(i) / det
        end do
    
        do i = 1, nd
           do j = 1, nd
              l_matrix(i,j)=(B(i)*B(j)+A(i)*A(j))*det/2.0d0
           end do
        end do
                   
        do k=1,nd
           i=nodex(ie,k) 
           do l=1,nd
              j=nodex(ie,l)  
              g_matrix(i,j)=g_matrix(i,j)+l_matrix(k,l)
           end do
        end do      
     end do      
     return
  end subroutine mkGSM_0
  
  SUBROUTINE mkGSM_1 (ND,Number_Element,Number_Node,eid,l_matrix,NODEX,XCOORD,YCOORD,g_matrix)
    IMPLICIT none
    integer nd,ne,Number_Node,Number_Element
    integer :: nodex(Number_Element,nd),eid(Number_Element)
    real(8) det
    real(8) :: x(nd),y(nd),a(nd),b(nd),l_matrix(nd,nd),xcoord(Number_Node),ycoord(Number_Node),g_matrix(Number_Node,Number_Node)
    integer i,j,k,l,ie,n

    g_matrix (:,:) = 0.0D0

    do ie=1,Number_Element
      if(eid(ie)/=2) then
        
        do i=1,nd
          x(i)=xcoord(nodex(ie,i))
          y(i)=ycoord(nodex(ie,i))
        end do
        
        B(1)=y(2)-y(3)
        B(2)=y(3)-y(1)
        B(3)=y(1)-y(2)
        A(1)=x(3)-x(2)
        A(2)=x(1)-x(3)
        A(3)=x(2)-x(1)
        det=A(3)*B(2)-B(3)*A(2)
        
        DO I = 1 , ND
          B(I)=B(I)/det
          A(I)=A(I)/det
        END DO

        DO I=1,ND
          DO J=1,ND
            l_matrix(I,J)=((B(i)*B(j)+A(i)*A(j))*det)/2.0d0
          END DO
        END DO
        
        DO K=1,ND
          I=NODEX(IE,K)
          DO L=1,ND
            J=NODEX(IE,L)
            g_matrix(I,J)=g_matrix(I,J)+l_matrix(K,L)
          END DO
        END DO
      end if
    end do
    return
  end SUBROUTINE mkGSM_1

  SUBROUTINE mkGSM (ND,Number_Element,Number_Node,Number_Sampling,lambda_va,eid,l_matrix,NODEX,XCOORD,YCOORD,g_matrix,t)
     IMPLICIT none
     integer nd,ne,Number_Node,Number_Element
     integer t,Number_sampling!t means current stage 
     integer :: nodex(Number_Element,nd),eid(Number_Element)
     real(8) det
     real(8) :: x(nd),y(nd),a(nd),b(nd)
     real(8) :: l_matrix(nd,nd), xcoord(Number_Node), ycoord(Number_Node), g_matrix(Number_Node,Number_Node)
     real(8) :: lambda_va(Number_Element,Number_Sampling)
     integer i,j,k,l,ie,n
  
      g_matrix (:,:) = 0.0D0
     
     do ie=1,Number_Element
       n=eid(ie)
       if(n/=2) then
         
         do i=1,nd
           x(i)=xcoord(nodex(ie,i))
           y(i)=ycoord(nodex(ie,i))
         end do
         
         B(1)=y(2)-y(3)
         B(2)=y(3)-y(1)
         B(3)=y(1)-y(2)
         A(1)=x(3)-x(2)
         A(2)=x(1)-x(3)
         A(3)=x(2)-x(1)
         det=A(3)*B(2)-B(3)*A(2)
         
         DO I = 1 , ND
           B(I)=B(I)/det
           A(I)=A(I)/det
         END DO
         
         DO j=1,ND
           DO i=1,ND
             l_matrix(I,J)=((B(i)*B(j)+A(i)*A(j))*(lambda_va(ie,t))*det)/2.0d0
           END DO
         END DO
         
         DO K=1,ND
           I=NODEX(IE,K)
           DO L=1,ND
             J=NODEX(IE,L)
             g_matrix(I,J)=g_matrix(I,J)+l_matrix(K,L)
           END DO
         END DO
       end if
     end do
     return
   end SUBROUTINE mkGSM

   subroutine setBC(Number_Element,Number_Node,nd,Number_Node_for_fem,xcoord,ycoord,g_matrix,T_rhs,mnx,mxx,mny,mxy)
      implicit none
      integer :: i,j
      integer :: nd, Number_Element, Number_Node, Number_Node_for_fem
      real(8) :: mnx, mxx, mny, mxy 
      real(8) :: xcoord( Number_Node ), ycoord( Number_Node ), T_rhs( Number_Node )
      real(8) :: T_e( Number_Element )
      real(8) :: g_matrix(Number_Node,Number_Node)
      real(8) :: rhs2(Number_Node_for_fem)
      
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
              g_matrix( i,j ) = 0.0d0
           end do
              g_matrix( i,i ) = 1.0d0
              T_rhs( i ) = 0.0d0 
        else if( xcoord( i ) == mxx ) then        
           do j = 1, Number_Node
              g_matrix( i,j ) = 0.0d0
           end do
              g_matrix( i,i ) = 1.0d0
              T_rhs( i ) = 1.0d0
        end if
      end do
      return
   end subroutine setBC
     
   subroutine setBC_0(Number_Element,Number_Node,nd,xcoord,ycoord,g_matrix,T_rhs,mnx,mxx,mny,mxy)
      implicit none
      integer :: i,j
      integer :: nd, Number_Element, Number_Node
      real(8) :: mnx, mxx, mny, mxy 
      real(8) :: xcoord( Number_Node ), ycoord( Number_Node ), T_rhs( Number_Node )
      real(8) :: g_matrix(Number_Node,Number_Node)
      
      !Reset
      T_rhs(:) = 0.0d0
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
              g_matrix( i,j ) = 0.0d0
           end do
              g_matrix( i,i ) = 1.0d0
              T_rhs( i ) = 0.0d0 
        else if( xcoord( i ) == mxx ) then
           do j = 1, Number_Node
              g_matrix( i,j ) = 0.0d0
           end do
              g_matrix( i,i ) = 1.0d0
              T_rhs( i ) = 1.0d0
        end if
      end do
      return
   end subroutine setBC_0

   subroutine mpi_proc_allocate(total_process,Number_tasks,process_id,number_firsttask,tasks_per_process)
      implicit none
      integer :: total_process,Number_tasks,process_id,number_firsttask,tasks_per_process

      if(process_id==total_process+1) then
        number_firsttask=process_id*(Number_tasks/total_process)+1
        tasks_per_process=(Number_tasks/total_process)+mod(Number_tasks,total_process)
      else 
        number_firsttask=process_id*(Number_tasks/total_process)+1
        tasks_per_process=(Number_tasks/total_process)
      end if
      return
   end subroutine mpi_proc_allocate

   

          
               
               
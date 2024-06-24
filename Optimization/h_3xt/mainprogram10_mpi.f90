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
   integer,parameter :: r_ins=30
   real(8),parameter :: mnx=-(dble(nelx/2))
   real(8),parameter :: mxx=-mnx
   real(8),parameter :: mny=-(dble(nely/2))
   real(8),parameter :: mxy=-mny
!============================================================================================

!============================================================================================
   !fem関連の変数(要素、節点など)の定義
   integer,parameter :: nd=3
   integer,parameter :: Number_Node=(nelx+1)*(nely+1)
   integer Number_nonzero,Number_Node_for_fem!断熱領域の中にある節点を除いた値
   integer,parameter :: Width_Matrix_LHS=7 !g_matrixにおいて一行の中に値が代入された成分の数の最大値(有限要素不規則分布の場合はmaxelementshareonenode+1)
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
   real(8),parameter :: t_boundary_higher_side=1.0d0
   real(8),parameter :: t_boundary_lower_side=0.0d0
   real(8) psi_n,psi!目的関数計算用
   real(8) :: objective_function_average,objective_function_pre,Best_Objective_function
!============================================================================================

!============================================================================================
   !その他、プログラム処理用変数
   integer,parameter :: total_process=10!mpi関連
   integer,parameter :: size_flat_nodex=Number_Element*nd!mpi送受信用
   integer :: process_id, num_procs, namelen, rc, ierr, number_firsttask,number_lasttask,tasks_per_process!mpi関連
   integer :: process_current_sample
   integer num_threads!スレッド指定用
   integer info!配列計算用
   integer i,i1,i2,j,k,e,e1,e2,d1,d2,d3,d4,lo,lo1,lo2,lo3,lo4,lo5
   integer sample
   integer ir,iw!入出力関連
   integer convergence_flag ! 1 or 0 収束判定結果を示す
   character filename*128
!============================================================================================

!============================================================================================
   !fem関連配列
   integer :: nodex(nd,Number_Element),nodex_s(4,Number_SElement),eid(Number_Element)
   integer :: SElement_coordinate(Number_SElement,Number_SElement),Reverse_SElement_coordinate(nelx,nely)
   integer :: Number_Node_for_fem_count(Number_Node)
   integer,allocatable :: g_matrix2_id(:)
   integer,allocatable :: SElement_Design_Area(:),SElement_Design_Area_quarter(:,:)
   real(8) :: xcoord(Number_Node),ycoord(Number_Node),T_rhs(Number_Node)
   real(8) :: l_matrix(nd,nd,Number_Element),g_matrix(Number_Node,Number_Node)
   real(8) :: nec(Number_Element,2),nesc(Number_SElement,2)
!============================================================================================

!============================================================================================
   !最適化関連配列
   real(8) :: lambda_va(Number_Element,Number_Sampling),sample_lambda_va(Number_Element)
   real(8) :: Best_sample_lambda_va(Number_Element)!設計変数から等価変換した熱伝導率
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
   !integer :: current_sample(total_process-1)
   integer,allocatable :: flat_SElement_Design_Area_quarter(:)
   !real(8) :: current_psi(total_process-1)
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

    iw=101
    open(iw,file='FEM_Meshdata.txt',status='replace')
      write(iw,*) 1, "# Flag_Electrical_Insulation_BC"
      write(iw,*) Number_Node, "# The Number of Nodes"
      write(iw,*) Number_Element, "# The Number of Triangle Element"
      write(iw,*) 6, "# Maximum Number of Elements sharing one Nodes"
      write(iw,*) mnx, "Minimal Value of X"
      write(iw,*) mny, "Minimal Value of Y"
      write(iw,*) mxx, "Maximum Value of X"
      write(iw,*) mxy, "Maximum Value of Y"
      do j=1,Number_Node
        write(iw,*) j,xcoord(j),ycoord(j)
      end do
      do i=1,Number_Element
        write(iw,*) i,eid(i),nodex(1,i),nodex(2,i),nodex(3,i)
      end do
      !do l=1,nes
      !  write(iw,*) l,eid(l*2),nodex_s(l,1),nodex_s(l,2),nodex_s(l,3),nodex_s(l,4)
      !end do
    close(iw)
 
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
    nodex=reshape(flat_nodex,[nd,Number_Element])
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
         Number_Node_for_fem_count(nodex(j,i))=1
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

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Cauculate T-reference
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  if(process_id==0) then 
   !======================================================== 
    write(*,*)'Making local Matrix'
   !========================================================
    call mkLSM_0 (ND,Number_Element,Number_Node,eid,NODEX,XCOORD,YCOORD,l_matrix)
   !===================================================================================
    write(*,*)'Solving Linear Equation' 
   !===================================================================================
    call analyze_thermal_distribution(l_matrix, nodex, eid, &
           0, Number_Node, Number_Element, Width_Matrix_LHS, &
           xcoord, ycoord, t_boundary_higher_side, t_boundary_lower_side, &
           mxx,mxy,mnx,mny, &
           !===========================================================================
           T_ref)
   
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Cauculate T-bare
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   !===================================================================================
    write(*,*)'Making local Matrix'
   !===================================================================================
    call mkLSM_1 (ND,Number_Element,Number_Node,eid,NODEX,XCOORD,YCOORD,l_matrix)
   !===================================================================================
    write(*,*)'Solving Linear Equation' 
   !===================================================================================
    call analyze_thermal_distribution(l_matrix, nodex, eid, &
       1, Number_Node, Number_Element, Width_Matrix_LHS, &
       xcoord, ycoord, t_boundary_higher_side, t_boundary_lower_side, &
       mxx,mxy,mnx,mny, &
       !===========================================================================
       T_bare)

      iw=109
       open(iw,file='Tbare.txt',status='replace')
        do i=1,Number_Element
          if(eid(i)/=2) then
           do j=1,nd
             write(iw,*) xcoord(nodex(j,i)), ycoord(nodex(j,i)), T_bare(nodex(j,i))
           end do
          end if
        end do
       close(iw)  

   !======================================================== 
    write(*,*)'T-bare calculated'
   !======================================================== 
  
  !============================================================================================
   !対称化
    do i=1,Number_SElement
       SElement_coordinate(i,1)=int(xcoord(nodex_s(3,i)))
       SElement_coordinate(i,2)=int(ycoord(nodex_s(3,i)))
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

    do j=1,nd
       do i=1,Number_Element
          if(eid(i)/=2)then
             npsi_n(nodex(j,i))=(abs(T_bare(nodex(j,i))-T_ref(nodex(j,i))))**2.0d0
          end if
       end do
    end do
 

    psi_n=sum(npsi_n)

    flat_SElement_Design_Area_quarter=reshape(SElement_Design_Area_quarter,[Number_SEDesignArea])
  end if
!=============================================================================================

!==========================================================================================
!broadcast data
 call MPI_BCAST(T_ref,Number_Node,MPI_DOUBLE,0,MPI_COMM_WORLD,ierr) 
 !call MPI_BCAST(g_matrix2_id,Number_Node_for_fem,MPI_INT,0,MPI_COMM_WORLD,ierr) 
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
    objective_function(:)=1.0d0
    lambda_va(:,:)=1.0d0
    sample_lambda_va(:)=1.0d0
    convergence_flag=0
    Number_matrix_design_variable=Number_Design_Variable*Number_Sampling   

    allocate(flat_design_variable(Number_matrix_design_variable))

    !call omp_set_num_threads(1)
    !call mkl_set_num_threads(1)

!============================================================================================
100 call MPI_Barrier(MPI_COMM_WORLD,ierr)
!============================================================================================

    !CMA-ES      
 if(process_id==0) then
      call cmaes(design_variable,Convergence_Ratio,&
                 'min', Type_CMA,Optimization_Step,Number_Design_Variable, lambda, 0.0d0, 3.0d0, fitness&
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

!============================================================================================
call MPI_Barrier(MPI_COMM_WORLD,ierr)
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

!=============================================================================================
  !start fem
   !============================================================================================
    call mpi_proc_allocate(total_process,lambda,process_id,number_firsttask,tasks_per_process)
   !============================================================================================
    
    number_lasttask= number_firsttask + tasks_per_process-1

    write(*,*) process_id,'sample=', number_firsttask, '~', number_lasttask

    do sample=number_firsttask,number_lasttask!while(change>1.0d-2)
          do i=1,Number_Design_Variable
            density_filtered_vec(i)=design_variable(i,sample)
          end do
  
       !set sample_lambda_va
       !do i=1,Number_Element
       !   sample_lambda_va(i,sample)=1.0d0
       !end do
  
       do i=1,Number_Design_Variable
         do j=1,4
          lo5=SElement_Design_Area_quarter(i,j)
          sample_lambda_va(lo5*2)=1.0d0+(density_filtered_vec(i))
          sample_lambda_va(lo5*2-1)=1.0d0+(density_filtered_vec(i))
         end do 
       end do
          
     !===================================================================================
     write(*,*)'Making Local Matrix'
     !===================================================================================
     call mkLSM (ND,Number_Element,Number_Node,Number_Sampling,sample_lambda_va,eid,NODEX,XCOORD,YCOORD,sample,l_matrix)
     !===================================================================================
     write(*,*)'Solving Linear Equation' 
     !===================================================================================
     !call mkl_set_num_threads(1)
     call analyze_thermal_distribution(l_matrix, nodex, eid, &
       1, Number_Node, Number_Element, Width_Matrix_LHS, &
       xcoord, ycoord, t_boundary_higher_side, t_boundary_lower_side, &
       mxx,mxy,mnx,mny, &
       !===========================================================================
       T_rhs)

       write(*,*) 'process=', process_id,'Sample=',sample, 'complete'
  
    !============================================================================================
     !compute objective function
      npsi(:)=0.0d0
      
      do j=1,nd
        do i=1,Number_Element
           if(eid(i)/=2)then
              npsi(nodex(j,i))=(abs(T_rhs(nodex(j,i))-T_ref(nodex(j,i))))**2.0d0
           end if
        end do
      end do    
             
      psi=sum(npsi)/psi_n
         
!      if(process_id==0) then
!        write(*,*) 'Sample=',sample  
!      end if
      process_current_sample=sample

      if(process_id==0)then
        Objective_function(process_current_sample)=psi
        fitness(process_current_sample)=Objective_function(process_current_sample)
        do i=1,Number_Element
           lambda_va(i,process_current_sample)=sample_lambda_va(i)
        end do
      end if

      !============================================================================================
      call MPI_Barrier(MPI_COMM_WORLD,ierr)
      !============================================================================================ 

      !call MPI_Gather(psi,1,MPI_double,current_psi,1,MPI_double,0,MPI_COMM_WORLD,ierr)
      !call MPI_Gather(sample,1,MPI_INTEGER,current_sample,1,MPI_double,0,MPI_COMM_WORLD,ierr)
!
      !if(process_id==0)then
      !  sample=0
      !  do i=1,total_process-1
      !     sample=current_sample(i)
      !     Objective_Function(sample)
      !  end do

        if(process_id/=0) then
           call MPI_Send(psi,1,MPI_DOUBLE,0,0,MPI_COMM_WORLD,ierr)
           call MPI_Send(process_current_sample,1,MPI_INTEGER,0,0,MPI_COMM_WORLD,ierr)
           call MPI_Send(sample_lambda_va,Number_Element,MPI_DOUBLE,0,0,MPI_COMM_WORLD,ierr)
           write(*,*) 'Sample=',sample, 'status=data sent'
        end if

        if(process_id==0) then
          do i=1,total_process-1
           call MPI_Recv(psi,1,MPI_DOUBLE,i,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
           call MPI_Recv(process_current_sample,1,MPI_INTEGER,i,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
           call MPI_Recv(sample_lambda_va,Number_Element,MPI_DOUBLE,i,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
           Objective_function(process_current_sample)=psi
           fitness(process_current_sample)=Objective_function(process_current_sample)
           do j=1,Number_Element
              lambda_va(j,process_current_sample)=sample_lambda_va(j)
           end do
           write(*,*) 'Sample=',process_current_sample, 'status=data received'
          end do
        end if    

      !============================================================================================
      call MPI_Barrier(MPI_COMM_WORLD,ierr)
      !============================================================================================
    end do
 
!============================================================================================
call MPI_Barrier(MPI_COMM_WORLD,ierr)
!============================================================================================

!============================================================================================
    !find best result&output
    if(process_id==0) then 
      iw=102
       open(iw,file='temperature.txt',status='replace')
        do i=1,Number_Element
          if(eid(i)/=2) then
           do j=1,nd
             write(iw,*) xcoord(nodex(j,i)), ycoord(nodex(j,i)), T_rhs(nodex(j,i))
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
         objective_function_average=objective_function_average+objective_function(i)/dble(Number_Sampling)
      end do
 
      Best_Objective_function_Loc_Temp=minloc(objective_function,1)
 
      do i=1,Number_Element
        Density_plot(i)=lambda_va(i,Best_Objective_function_Loc_Temp)-1.0d0
      end do
      
      do e=1,Number_Element
         do i=1,3
            position_in_element(e,i,1)=xcoord(nodex(i,e))
            position_in_element(e,i,2)=ycoord(nodex(i,e))
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
          write(iw,*) 'best value of Objective function is',Best_Objective_function
          write(iw,*) 'best design variables'
  
          do i=1,Number_Design_Variable
            write(iw,*) design_variable(i,Best_Objective_function_Loc)
          end do
       close(iw)  
 
       do i=1,Number_Element
          Best_sample_lambda_va(i)=lambda_va(i,Best_Objective_function_Loc)    
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
               (450.0d0/dble(nely))*position_in_element(e,1,2)+225.0d0, 'moveto'
               write(iw,*) (600.0d0/dble(nelx))*position_in_element(e,2,1)+300.0d0,&
               (450.0d0/dble(nely))*position_in_element(e,2,2)+225.0d0, 'lineto'
               write(iw,*) (600.0d0/dble(nelx))*position_in_element(e,3,1)+300.0d0,&
               (450.0d0/dble(nely))*position_in_element(e,3,2)+225.0d0, 'lineto'
               write(iw,*) 'closepath'
               write(iw,*) 1.0d0-(Best_sample_lambda_va(e)-1.0d0)/3.0d0,'setgray'
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
             (450.0d0/dble(nely))*position_in_element(e,1,2)+225.0d0, 'moveto'
             write(iw,*) (600.0d0/dble(nelx))*position_in_element(e,2,1)+300.0d0,&
             (450.0d0/dble(nely))*position_in_element(e,2,2)+225.0d0, 'lineto'
             write(iw,*) (600.0d0/dble(nelx))*position_in_element(e,3,1)+300.0d0,&
             (450.0d0/dble(nely))*position_in_element(e,3,2)+225.0d0, 'lineto'
             write(iw,*) 'closepath'
             write(iw,*) 1.0d0-Density_plot(e)/3.0d0,'setgray'
             write(iw,*) 'fill'
             write(iw,*) ''
           end if
         end do    
     close(iw)
 
     Optimization_step=Optimization_step+1
 
     write(*,*)'step',Optimization_step
    end if

!===================================================================================
   !収束判定
   if(process_id==0) then
      if(Convergence_Ratio_dowhile<Convergence_Error)then
         !===================================================================================
          write(*,*)'Optimization Complete'
         !===================================================================================      
         convergence_flag=1
      end if     
   end if

   call MPI_BCAST(convergence_flag,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr) 

!============================================================================================
call MPI_Barrier(MPI_COMM_WORLD,ierr)
!============================================================================================

   if(convergence_flag==1)then 
      goto 101
   end if

   goto 100 

!============================================================================================
101 call MPI_Finalize(ierr)
!============================================================================================ 

 deallocate(SElement_Design_Area,SElement_Design_Area_quarter)
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
  integer :: nodex(3,Number_Element),nodex_s(4,Number_SElement),eid(Number_Element)
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
          nodex_s( 1,m ) = j + i*( nely + 1 )
          nodex_s( 2,m ) = j + i*( nely + 1 ) + ( nely + 1 )
          nodex_s( 3,m ) = j + i*( nely + 1 ) + ( nely + 2 )
          nodex_s( 4,m ) = j + i*( nely + 1 ) + 1
       end do
    end do
    
    do i = 1, Number_SElement
       nodex( 1,i*2 ) = nodex_s( 4,i )
       nodex( 2,i*2 ) = nodex_s( 2,i )
       nodex( 3,i*2 ) = nodex_s( 3,i )	
       nodex( 1,i*2-1 ) = nodex_s( 4,i )
       nodex( 2,i*2-1 ) = nodex_s( 1,i )
       nodex( 3,i*2-1 ) = nodex_s( 2,i )
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
    integer :: nodex( nd,Number_Element )
    real(8) :: det
    real(8) :: x( nd ), y( nd ), A( nd ), B( nd ), l_matrix( nd,nd )
    real(8) :: xcoord( Number_Node ), ycoord( Number_Node ), g_matrix( Number_Node,Number_Node )

    g_matrix( :,: ) = 0.0d0
  
    do ie = 1, Number_Element
  
        do i = 1, nd
           x(i)=xcoord(nodex(i,ie))
           y(i)=ycoord(nodex(i,ie))
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
           i=nodex(k,ie) 
           do l=1,nd
              j=nodex(l,ie)  
              g_matrix(i,j)=g_matrix(i,j)+l_matrix(k,l)
           end do
        end do      
     end do      
     return
  end subroutine mkGSM_0

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

   SUBROUTINE mkLSM (ND,Number_Element,Number_Node,Number_Sampling,sample_lambda_va,eid,NODEX,XCOORD,YCOORD,t,l_matrix)
   IMPLICIT none
   integer nd,Number_Node,Number_Element
   integer t,Number_sampling!t means current stage 
   integer :: nodex(nd,Number_Element),eid(Number_Element)
   real(8) det
   real(8) :: x(nd),y(nd),a(nd),b(nd)
   real(8) :: l_matrix(nd,nd,Number_Element), xcoord(Number_Node), ycoord(Number_Node)
   real(8) :: sample_lambda_va(Number_Element)
   integer i,j,ie
   
   do ie=1,Number_Element
     if(eid(ie)/=2) then
       
       do i=1,nd
         x(i)=xcoord(nodex(i,ie))
         y(i)=ycoord(nodex(i,ie))
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
           l_matrix(I,J,ie)=((B(i)*B(j)+A(i)*A(j))*(sample_lambda_va(ie))*det)/2.0d0
         END DO
       END DO
     end if
   end do
   return
   END SUBROUTINE mkLSM

   SUBROUTINE mkLSM_0 (ND,Number_Element,Number_Node,eid,NODEX,XCOORD,YCOORD,l_matrix)
   IMPLICIT none
   integer nd,Number_Node,Number_Element
   integer :: nodex(nd,Number_Element),eid(Number_Element)
   real(8) det
   real(8) :: x(nd),y(nd),a(nd),b(nd)
   real(8) :: l_matrix(nd,nd,Number_Element), xcoord(Number_Node), ycoord(Number_Node)
   integer i,j,ie

   l_matrix=0.0d0
   
   do ie=1,Number_Element
       
       do i=1,nd
         x(i)=xcoord(nodex(i,ie))
         y(i)=ycoord(nodex(i,ie))
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
           l_matrix(I,J,ie)=(B(i)*B(j)+A(i)*A(j)*det)/2.0d0
         END DO
       END DO
   end do
   return
   END SUBROUTINE mkLSM_0

   SUBROUTINE mkLSM_1 (ND,Number_Element,Number_Node,eid,NODEX,XCOORD,YCOORD,l_matrix)
   IMPLICIT none
   integer nd,Number_Node,Number_Element
   integer :: nodex(nd,Number_Element),eid(Number_Element)
   real(8) det
   real(8) :: x(nd),y(nd),a(nd),b(nd)
   real(8) :: l_matrix(nd,nd,Number_Element), xcoord(Number_Node), ycoord(Number_Node)
   integer i,j,ie

   l_matrix=0.0d0
   
   do ie=1,Number_Element
     if(eid(ie)/=2) then
       
       do i=1,nd
         x(i)=xcoord(nodex(i,ie))
         y(i)=ycoord(nodex(i,ie))
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
           l_matrix(I,J,ie)=(B(i)*B(j)+A(i)*A(j)*det)/2.0d0
         END DO
       END DO
     end if
   end do
   return
   END SUBROUTINE mkLSM_1

   

          
               
               
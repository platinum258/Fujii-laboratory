program thermalcloak
      !$ use omp_lib
      !   use parameters
      !   use Objective_function_subprogram
   use mpi

   implicit none
   integer :: my_id, num_procs, namelen, rc, ierr
   integer :: enum,countnum,num_threads,countstep,Number_Node,Number_Node_ins,node
   integer SE1
   integer ID_Element_Objective_Function_1
   !  integer npb,npec,epb
   integer INFO
   integer generation,situation
   integer Optimization_step, lambda!, mu
   integer :: flg_mi,g_opt,n_dv,popsize!,type_cma
   integer Number_Design_Variable,Number_matrix_design_variable,Number_SElement
   integer Number_EDesignArea, Number_SEDesignArea, Number_Element, Number_Eins, Number_SEins
   integer i,j,k,k1,k2,k3,l,m,m0,m2,m3,n,e,t,ir,iw,ie,iel,ier,iek,ik,e1,e2,n1,n2,n3,i1,i2,l1,l2,j1,j2,r1,r2,r3
   
   character(len=1), parameter :: Type_CMA = 'N'
   character filename*128
   
   real(8) MNX,MXX,MNY,MXY
   real(8) Convergence_Ratio_dowhile,det
   real(8) p,p1,p2,m1
   real(8) psi_n,psi,volume
   real(8) max_val!=1.0d0
   real(8) min_val!=0.0d0
   real(8) val_size,sigma
   

! real(8),parameter :: height=80.0d0
! real(8),parameter :: width=120.0d0
   !map size
   integer,parameter :: nelx=120
   integer,parameter :: nely=80
   integer,parameter :: design_area_lx=50
   integer,parameter :: design_area_ly=50

! integer,parameter :: Number_Node=nelx*nely*2
   integer,parameter :: nd=3
   integer,parameter :: penal=3,rmin=3
   !integer,parameter :: Number_Element=nelx*nely*2
   integer,parameter :: Number_sampling=100
   !integer,parameter :: Number_EDesignArea=4200
   !optimization step?
   !integer,parameter :: lambda=350 
   real(8),parameter :: Convergence_Error=1.0d-6 !収束しやすいよう1桁増した
   !real(8),parameter :: lambda_cu=391.0d0 無次元化するので不要
   real(8),parameter :: t_flat=1.0d0
   real(8),parameter :: vol_constraint=0.5d0

   real(8) :: Objective_function(Number_Sampling)
   real(8) :: sk(nd,nd)
   real(8) :: x(nd),y(nd),a(nd),b(nd)
   real(8) :: objective_function_average,objective_function_pre

   integer,allocatable :: nodex(:,:),nodex_s(:,:),eid(:)
   integer,allocatable :: SE1_count(:)
   integer,allocatable :: position_in_element(:,:,:)
   integer,allocatable :: Element_Design_Area(:)
   integer,allocatable :: SElement_Design_Area(:),Element_Design_Area_quarter(:,:)
   integer,allocatable :: rko_id(:)
   integer,allocatable :: SEcord(:,:),Reverse_SEcord(:,:)
  !  integer :: rkr_id(SE1)
  !  integer,allocatable :: edofvec(:),iH(:),jH(:),iK(:),jK(:)

   real(8),allocatable :: xcoord(:),ycoord(:)
   real(8),allocatable :: T_rhs(:)!,T_e(:)
   real(8),allocatable :: design_variable(:,:),dv_plot(:,:),rk(:,:)
   real(8),allocatable :: Vec_D(:),mat_B(:,:),average(:)
   real(8),allocatable :: convergence_ratio(:),fitness(:)
  !real(8) :: density_filtered_mat(Number_Element,lambda)
  !real(8) ipmkl(Number_Node)
   real(8),allocatable :: rko(:,:),rhso(:),ipmkl(:)
   real(8),allocatable :: density_filtered_vec(:),shape_function(:,:),xecord(:),yecord(:)
   real(8),allocatable :: density_vec(:),density_vec_s(:)!number of element in design area while square
   real(8),allocatable :: nec(:,:)!,SEcord(Number_SElement,2)
  !real(8) lambda_va(Number_Design_Variable)
   real(8),allocatable :: lambda_va(:)!設計変数から等価変換した熱伝導率
   real(8),allocatable :: T_ref(:),T_bare(:)
   real(8),allocatable :: epsi_n(:),epsi(:)
   real(8),allocatable :: T_e(:)
   real(8),allocatable :: E_str(:)

!start mpi
 call MPI_INIT(ierr)
 call MPI_COMM_RANK(MPI_COMM_WORLD,my_id,ierr)
 call MPI_COMM_SIZE(MPI_COMM_WORLD,num_procs,ierr)

   ir=1
   open(IR,file='MdataC.txt')
   read(IR,*)Number_Node
   allocate(xcoord(Number_Node),ycoord(Number_Node),T_rhs(Number_Node),rk(Number_Node,Number_Node))
   !read(IR,*)mepn
  !  read(IR,*)npb
  !  read(IR,*)npec
  !  read(IR,*)epb
   read(IR,*)Number_Element
   allocate(nodex(Number_Element,3),eid(Number_Element))
   read(IR,*)Number_SElement
   allocate(nodex_s(Number_SElement,4))
   read(IR,*)Number_EDesignArea
   read(IR,*)Number_SEDesignArea
   read(IR,*)Number_Eins
   read(IR,*)Number_SEins
   read(IR,*)ID_Element_Objective_Function_1
   read(IR,*)MNX
   read(IR,*)MXX
   read(IR,*)MNY
   read(IR,*)MXY
   read(IR,*) (node,XCOORD(node),YCOORD(node),j=1,Number_Node)
   read(ir,*) (iel,eid(iel),(nodex(iel,j),j=1,nd),i=1,Number_Element)
   read(ir,*) (iel,eid(iel*2),(nodex_s(iel,k),k=1,4),i=1,Number_SElement)
   close(ir)
   write(*,*)'input over'

   
   Number_Design_Variable=Number_SEDesignArea/4

    allocate(position_in_element(Number_Element,3,2))
    allocate(Element_Design_Area(Number_EDesignArea))
    allocate(SElement_Design_Area(Number_SEDesignArea),Element_Design_Area_quarter(Number_SEDesignArea/4,4))
    allocate(SEcord(Number_SElement,Number_SElement),Reverse_SEcord(nelx,nely))

    ! allocate(T_e(Number_EDesignArea))
    allocate(SE1_count(Number_Node))
    allocate(density_filtered_vec(Number_Design_Variable),shape_function(Number_Element,3),&
xecord(Number_Element),yecord(Number_Element))
    allocate(average(Number_Design_Variable))
    allocate(Vec_D(Number_Design_Variable),mat_B(Number_Design_Variable,Number_Design_Variable))
    allocate(density_vec(Number_EDesignArea),density_vec_s(Number_SEDesignArea))
    allocate(nec(Number_Element,2))!,SEcord(Number_SElement,2)
    allocate(lambda_va(Number_Element))!設計変数から等価変換した熱伝導率
    allocate(T_ref(Number_Element),T_bare(Number_Element))
    allocate(epsi_n(Number_Element),epsi(Number_Element))
    allocate(T_e(Number_Element))
    allocate(E_str(Number_Element))

   !count SE1
   do i=1,Number_Node
     SE1_count(i)=0
   end do

   do i=1,Number_Element
     if(eid(i)/=2)then
       do j=1,nd
        SE1_count(nodex(i,j))=1
       end do
     end if
   end do

   SE1=0
   do i=1,Number_Node
     if(SE1_count(i)==1)then 
       SE1=SE1+1
     end if
   end do

   write(*,*)"SE1=",SE1
   Number_Node_ins=Number_Node-SE1

   allocate(rko(SE1,SE1),rhso(SE1),ipmkl(SE1))
   allocate(rko_id(SE1))
   

  if(Number_Sampling>=4+int(3.0d0*log(dble(Number_Element)))) then
    lambda=Number_Sampling
  else
    lambda=4+int(3.0d0*log(dble(Number_Element)))
  end if

  allocate(design_variable(Number_Design_Variable,lambda))
  allocate(convergence_ratio(lambda))
  allocate(fitness(lambda))

  !set SElement cordinate(右上節点の座標に依存)
   do i=1,Number_SElement
      SEcord(i,1)=int(xcoord(nodex_s(i,3)))
      SEcord(i,2)=int(ycoord(nodex_s(i,3)))
   end do   
  
   do i=1,nelx
     do j=1,nely
       Reverse_SEcord(i,j)=0
     end do
   end do    

  !set reverse_SEcord(座標から番号，縦横が1から始まるように値を調整)
   do i=1,Number_SElement
      i1=SEcord(i,1)+60
      i2=SEcord(i,2)+40
      Reverse_SEcord(i1,i2)=i 
   end do 

   !create quarter of square element
   l2=1
   do i=nelx/2+1,nelx
     do j=nely/2+1,nely
        j1=reverse_SEcord(i,j) ! 1st quadrant
        r1=reverse_SEcord(nelx-i+1,j)!2nd quadrant
        r2=reverse_SEcord(nelx-i+1,nely-j+1)!3rd quadrant
        r3=reverse_SEcord(i,nely-j+1)!4th quadrant
      if(eid(j1*2)==4)then
        Element_Design_Area_quarter(l2,1)=j1
        Element_Design_Area_quarter(l2,2)=r1
        Element_Design_Area_quarter(l2,3)=r2
        Element_Design_Area_quarter(l2,4)=r3
        l2=l2+1
      end if
     end do
   end do   

   do i=1,Number_Element
      T_e(i)=0.0d0
   end do

   ir=2
   open(ir,file='Treferance.txt')
   do i=1,Number_Element
   read(ir,*) ier,T_ref(ier)
   end do
   close(ir)

   do i=1,Number_Element
   T_bare(i)=0.0d0
   end do

   ir=3
   open(ir,file='Tbare.txt')
   do j=1,Number_Element
     if(eid(j)/=2) then
      read(ir,*) T_bare(j)
     end if
   end do
   close(ir)
 
   !set element center
   e1=0
   e2=0
   do i=0,nelx-1
     do j=1,nely
       e2=i*160+j*2
       e1=e2-1
       nec(e2,1)=dble(i)-59.25d0
       nec(e2,2)=dble(j)-40.25d0
       nec(e1,1)=dble(i)-59.75d0
       nec(e1,2)=dble(j)-40.75d0
     end do
   end do

   ! set Element_Design_Area
  !  j2=1
  !  do i=1,Number_Element
    !  if(eid(i)==4) then
    !  Element_Design_Area(j2)=i
    !  j2=j2+1
    !  end if
  !  end do
 
      
   ! if(Number_Sampling>=4+int(3.0d0*log(dble(Number_Element)))) then
   !    lambda=Number_Sampling
   ! else 
   !    lambda=4+int(3.0d0*log(dble(Number_Element)))
   ! end if
   ! mu=lambda/2

   ! allocate(design_variable(Number_Element,lambda))
   ! allocate(convergence_ratio(lambda))
   ! allocate(fitness(lambda))
   ! allocate(average(Number_Element))
   ! allocate(Vec_D(Number_Element))
   ! allocate(mat_B(Number_Element,Number_Element))
   ! allocate(dv_plot(Number_Element,lambda))
   
   call omp_set_num_threads(4)
   call mkl_set_num_threads(4)

   ! Filter Preparation

   ! iH(:)=1
   ! jH(:)=1
   ! sH(:)=0.0d0
   
   ! z=0
   ! do i=1,nelx
   !   do j=1,nely
   !     e1=(i-1)*nely+j
   !      do k=max(i-(rmin-1),1),min(i+(rmin-1),nelx)
   !         do j2=max(j-(rmin-1),1), min(j+(rmin-1),nely)
   !            e2=(k-1)*nely+j2

   !            z=z+1
   !            iH(z)=e1
   !            jH(z)=e2
   !            sH(z)=max(0.0d0,rmin-sqrt((dble(i)-dble(k))**2+(dble(j)-dble(j2))**2))
   !         end do
   !      end do
   !   end do
   ! end do      

   ! H=0.0d0
   ! do i=1,nely*nelx*(2*(rmin-1)+1)**2
   !   H(iH(i),jH(i))=H(iH(i),jH(i))+sH(i)
   ! end do

   ! Hs=0.0d0
   ! do i=1,Number_Element
   !    do j=1, Number_Element
   !       Hs(i)=Hs(i)+H(i,j)
   !    end do
   ! end do

  ! initialization

  !  do j=1,nelx
  !    do i=1,nely
  !       design_variable(i,j)=vol_constraint
  !    end do
  !  end do

  !  do j=1,nelx
  !    do i=1,nely
  !      density_filtered_mat(i,j)=design_variable(i,j)
  !    end do
  !  end do
   
  !  change=1.0d0

   ! call function_inf(function_number,function_name,min_val,max_val)
   !calculate psi_n

!id0
if(my_id==0) then
    do i=1,Number_Element
       epsi_n(i)=0.0d0
    end do

    do i=1,Number_Element
       if(eid(i)/=2) then 
       epsi_n(i)=abs(T_bare(i)-T_ref(i))**2.0d0
       end if
    end do   

     psi_n=sum(epsi_n)

     val_size=max_val-min_val

     Convergence_Ratio_dowhile=1.0d0
     Optimization_step=0
     situation=0
     objective_function_average=0.0d0
     objective_function_pre=0.0d0
     Number_matrix_design_variable=Number_Design_Variable*100
   write(*,*) 'initialized'
end if 

call MPI_Barrier(MPI_COMM_WORLD,ierr)!他のスレッドが勝手に動かさぬようにする
  ! start CMA-ES(min0.0d0,max10.0d0)
       
  ! do while(Convergence_Ratio_dowhile>Convergence_Error)　!id0の最後に判断する
100 if(my_id==0) then  
    call cmaes(design_variable,Convergence_Ratio(Optimization_step),&
               'min', Type_CMA,Optimization_Step,Number_Design_Variable, lambda, 0.0d0, 1.0d0, fitness&
              ,average,Vec_D,mat_B,sigma)!, flg_min, type_cma) 

    write(*,*)'convergence ratio=', convergence_ratio(Optimization_step)

    if(Optimization_step>=1) then
      Convergence_Ratio_dowhile=abs(objective_function_average-objective_function_pre)
    end if

    objective_function_pre=objective_function_average
     
    write(*,*)'convergence error=', Convergence_Ratio_dowhile

    call MPI_BCAST(objective_function_average,1,MPI_double,0,MPI_COMM_WORLD,ierr) 
    call MPI_BCAST(psi_n,1,MPI_double,0,MPI_COMM_WORLD,ierr)    

end if 
!id0 から集団送信
 
!call MPI_Barrier(MPI_COMM_WORLD,ierr)!他のスレッドが勝手に動かさぬようにする

!id1~100
if(my_id/=0) then  
  !call MPI_BCAST(objective_function_average,1,MPI_double,0,MPI_COMM_WORLD,ierr) 
  !call MPI_BCAST(psi_n,1,MPI_double,0,MPI_COMM_WORLD,ierr)    
  !call MPI_BCAST(design_variable,Number_matrix_design_variable,MPI_double,0,MPI_COMM_WORLD,ierr)
  !  do t=1,lambda!while(change>1.0d-2)
        do i=1,Number_Design_Variable
        density_filtered_vec(i)=design_variable(i,my_id)
        end do

       !square to triangle
        ! do i=1,Number_SEDesignArea
          ! ik=i*2-1
          ! density_vec(i*2)=density_vec_s(i)
          ! density_vec(ik)=density_vec_s(i)
        ! end do              

     !set lambda_va
     do i=1,Number_Element
        lambda_va(i)=1.0d0
     end do

     do i=1,Number_Design_Variable
       do j=1,4
        iek=Element_Design_Area_quarter(i,j)
        lambda_va(iek*2)=1.0d0+(density_filtered_vec(i))
        lambda_va(iek*2-1)=1.0d0+(density_filtered_vec(i))
       end do 
     end do
    !write(*,*) my_id, 'start making GSM'
    call mkGSM(ND,Number_Element,Number_Node,lambda_va,eid,x,y,A,B,SK,NODEX,XCOORD,YCOORD,RK)
    !write(*,*) my_id, 'GSM complete'
    !Result matrix initialization
     do i=1,Number_Node
       T_rhs(i)=0.0d0
     end do
    !set Norman boundary
     do k=1,Number_Node
       if(ycoord(k)==mny .or. ycoord(k)==mxy) then
        T_rhs(k)=0.0d0
       end if
     end do
    !set Dirichlet boundary
     do k=1,Number_Node
       if(xcoord(k) .eq. mnx) then
        do j=1,Number_Node
          rk(k,j)=0.0d0
        end do
          rk(k,k)=1.0d0
          T_rhs(k)=0.0d0 

       else if(xcoord(k) .eq. mxx) then
       do l=1,Number_Node
          rk(k,l)=0.0d0
       end do
          rk(k,k)=1.0d0
          T_rhs(k)=1.0d0
       end if
     end do
    
    l1=1
    do i=1,Number_Node
      if(rk(i,i)/=0.0d0) then
       rko_id(l1)=i
       l1=l1+1
      end if
    end do

    do i=1,SE1
       do j=1,SE1
       rko(i,j)=rk(rko_id(i),rko_id(j))
       end do
      rhso(i)=T_rhs(rko_id(i))
    end do   
!write(*,*) my_id, 'start cauculate matrix'
call mkl_set_num_threads(num_threads)
call dgesv(SE1,1,rko,SE1,ipmkl,rhso,SE1,info)
!write(*,*) my_id, 'calculated'
    do i=1,SE1
     T_rhs(rko_id(i))=rhso(i)
    end do

    call etemp(Number_Element,Number_Node,nd,det,eid,xcoord,ycoord,nec,x,y,a,b,n1,n2,n3,nodex,shape_function,&
   T_rhs,T_e) 
    
    ! iw=106
    ! open(iw,file='etemperature.txt',status='replace')
    !  do i=1,Number_Element
      !  if(eid(i)/=2) then
          ! write(iw,*) nec(i,1), nec(i,2), T_e(i)
      !  end if
    !  end do
    ! close(iw)

   !compute objective function
    do j=1,Number_Element
       epsi(j)=0.0d0
    end do
    
    do i=1,Number_Element
       if(eid(i)==5) then
       epsi(i)=abs(T_e(i)-T_ref(i))**2.0d0
       end if
    end do

    psi=sum(epsi)/psi_n
    
    
    ! if(psi<=0.0d0) then
      ! situation=1
    ! end if

    volume=sum(density_filtered_vec)/dble(Number_Design_Variable)

    ! if(volume .gt. 0.5d0) then 
    !    fitness(t)=Objective_function(t)+exp(1.0d1*(volume/1.0d0)-0.5d0)-1.0d0
    ! else
       
    ! end if
    !目的関数、温度分布データなどをsend back

    write(*,*) 't=',my_id
end if  
  !  end do
call MPI_Barrier(MPI_COMM_WORLD,ierr)
if(my_id==Number_Sampling) then 
    iw=102!mpiによるエラーを回避するため一世代につき温度の出力を一回だけする
     open(iw,file='temperature.txt',status='replace')
      do i=1,Number_Element
        if(eid(i)/=2) then
         do j=1,nd
           write(iw,*) xcoord(nodex(i,j)), ycoord(nodex(i,j)), T_rhs(nodex(i,j))
         end do
        end if
      end do
     close(iw)
     iw=107
      open(iw,file='Design_Variable.txt',status='replace')
        write(iw,*) density_filtered_vec
      close(iw)      
end if     

call MPI_Barrier(MPI_COMM_WORLD,ierr)
call MPI_Gather(psi,1,MPI_double,Objective_Function,1,MPI_double,0,MPI_COMM_WORLD,ierr)
call MPI_Barrier(MPI_COMM_WORLD,ierr)

if(my_id==0) then !id0の処理再開
    
    do i=1,Number_Sampling
    fitness(t)=Objective_function(t)
    end do
     
     objective_function_average=0.0d0

     do i=1,Number_Sampling
        objective_function_average=objective_function_average+objective_function(i)/100.0d0
     end do

    iw=103!一世代の平均値を出力するようにした
     open(iw,file='Objective_Function.dat',position='append')
       write(iw,*) Optimization_step,objective_function_average
     close(iw)

     do i=1,Number_Element
       E_str(i)=lambda_va(i)-1.0d0
     end do
     
     do e=1,Number_Element
        do i=1,3
           position_in_element(e,i,1)=xcoord(nodex(e,i))
           position_in_element(e,i,2)=ycoord(nodex(e,i))
        end do
     end do

    write(filename,'(i5.5,"_Structure.ps")') Optimization_step

    write(*,*)'output.ps >>>>>' ,filename
    iw=104
    open(iw,file=filename,status='replace')
       write(iw,*)'%!PS-Adobe-3.0 EPSF-3.0'
       write(iw,*)'%%BoundingBox:',0,0,600,400
       write(iw,*)'/m {moveto} def'
       write(iw,*)'/l {moveto} def'
       write(iw,*)'1.0 setgray'
       write(iw,*)'gsave'
       write(iw,*)'0 0 moveto'
       write(iw,*)'600 0 lineto'
       write(iw,*)'600 400 lineto'
       write(iw,*)'0 400 lineto'
       write(iw,*)'closepath'
       write(iw,*)'fill'
       write(iw,*)'grestore'
       write(iw,*)'newpath'
       write(iw,*)'0.0 setgray'
        do e=1,Number_Element
        if(eid(e)==4) then
          write(iw,*) 5.0d0*position_in_element(e,1,1)+300.0d0,5.0d0*position_in_element(e,1,2)+200.0d0, 'moveto'
          write(iw,*) 5.0d0*position_in_element(e,2,1)+300.0d0,5.0d0*position_in_element(e,2,2)+200.0d0, 'lineto'
          write(iw,*) 5.0d0*position_in_element(e,3,1)+300.0d0,5.0d0*position_in_element(e,3,2)+200.0d0, 'lineto'
          write(iw,*) 'closepath'
          write(iw,*) 1.00-E_str(e),'setgray'
          write(iw,*) 'fill'
          write(iw,*) ''
        end if
        end do    
    close(iw)

    ! if(situation==1) then
      ! write(*,*) 'error'
    !  exit
    ! end if

    Optimization_step=Optimization_step+1

    write(*,*)'step',Optimization_step
end if

!収束判定
if(my_id==0) then
   if(Convergence_Ratio_dowhile>=Convergence_Error)then
      goto 100
   else 
      write(*,*) 'complet'
   end if     
end if

  ! write(*,*)'result output'

  ! iw=101
  ! open(iw,file='result.txt',status='new')
  ! write(iw,*) 'optimizationstep=',Optimization_step
  ! write(iw,*) 'psi=',psi
  ! do i=1,Number_Element
    ! if(eid(i)/=2) then
    ! write(iw,*) i,lambda_va(i),T_e(i)
    ! end if
  ! end do
  ! close(iw）

 deallocate(position_in_element)
 deallocate(Element_Design_Area)
 deallocate(SElement_Design_Area,Element_Design_Area_quarter)
 deallocate(rko_id)
 deallocate(SEcord,Reverse_SEcord)  
 deallocate(rko,rhso,ipmkl)
 deallocate(Design_Variable,SE1_count)
 deallocate(density_filtered_vec,shape_function,xecord,yecord)
 deallocate(convergence_ratio)
 deallocate(fitness)
 deallocate(Vec_D)
 deallocate(mat_B)
 deallocate(nodex)
 deallocate(nodex_s)
 deallocate(eid)
 deallocate(xcoord)
 deallocate(ycoord)
 deallocate(T_rhs)
 deallocate(rk)
 deallocate(average)
 deallocate(density_vec,density_vec_s)
 deallocate(nec)
 deallocate(lambda_va)
 deallocate(T_ref,T_bare)
 deallocate(epsi_n,epsi)
 deallocate(T_e)
 deallocate(E_str)
 
call MPI_FINALIZE(ierr)
!end mpi
end program 

  !temparature interpolation
  subroutine etemp(Number_Element,Number_Node,nd,det,eid,xcoord,ycoord,nec,x,y,a,b,n1,n2,n3,nodex,shape_function,&
  T_rhs,T_e)  
  implicit none
  integer i,n
  integer Number_Element,Number_Node,nd,n1,n2,n3
  integer :: nodex(Number_Element,nd),eid(Number_Element)
  real(8) det
  real(8) :: x(nd),y(nd),a(nd),b(nd), sk(nd,nd), xcoord(Number_Node), ycoord(Number_Node),nec(Number_Element,2)&
  ,T_rhs(Number_Node)
  real(8) :: shape_function(Number_Element,3),T_e(Number_Element)
  do n=1,Number_Element
     if(eid(n)/=2) then
     do i=1,nd
       x(i)=xcoord(nodex(n,i))
       y(i)=ycoord(nodex(n,i))
     end do
      !  B(1)=y(2)-y(3)
      !  B(2)=y(3)-y(1)
      !  B(3)=y(1)-y(2)
      !  A(1)=x(3)-x(2)
      !  A(2)=x(1)-x(3)
      !  A(3)=x(2)-x(1)
       det=(x(2)-x(1))*(y(3)-y(1))-(y(1)-y(2))*(x(1)-x(3))
       shape_function(n,1)=(x(2)*y(3)-x(3)*y(2))/det+nec(n,1)*(y(2)-y(3))/det+nec(n,2)*(x(3)-x(2))/det 
       shape_function(n,2)=(x(3)*y(1)-x(1)*y(3))/det+nec(n,1)*(y(3)-y(1))/det+nec(n,2)*(x(1)-x(3))/det 
       shape_function(n,3)=(x(1)*y(2)-x(2)*y(1))/det+nec(n,1)*(y(1)-y(2))/det+nec(n,2)*(x(2)-x(1))/det 
       n1=nodex(n,1)
       n2=nodex(n,2)
       n3=nodex(n,3)
       T_e(n)=shape_function(n,1)*T_rhs(n1)+shape_function(n,2)*T_rhs(n2)+shape_function(n,3)*T_rhs(n3)
     end if
  end do
  end  

  SUBROUTINE mkGSM (ND,Number_Element,Number_Node,lambda_va,eid,x,y,A,B,SK,NODEX,XCOORD,YCOORD,RK)
  IMPLICIT none
  integer nd,ne,Number_Node,Number_Element
  integer :: nodex(Number_Element,nd),eid(Number_Element)
  real(8) det
  real(8) :: x(nd),y(nd),a(nd),b(nd), sk(nd,nd), xcoord(Number_Node), ycoord(Number_Node), rk(Number_Node,Number_Node)
  real(8) :: lambda_va(Number_Element)
  integer i,j,k,l,ie,n

  DO I = 1 , Number_Node
  DO J = 1 , Number_Node
  RK ( I , J ) = 0.D0
  END DO
  END DO
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
  DO I=1,ND
  DO J=1,ND
  SK(I,J)=((B(i)*B(j)+A(i)*A(j))*(lambda_va(ie))*det)/2.0d0
  END DO
  END DO
  DO K=1,ND
  I=NODEX(IE,K)
  DO L=1,ND
  J=NODEX(IE,L)
  RK(I,J)=RK(I,J)+SK(K,L)
  END DO
  END DO
  end if
  end do
  return
  end


   

          
               
               

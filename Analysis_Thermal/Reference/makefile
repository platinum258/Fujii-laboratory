
#include /home/garuda_fujii/Library/MUMPS_ALL/include.h 

PROGS=a.out

FCFLAGS= $(FCFLAGS_O2)
#FCFLAGS= $(FCFLAGS_O2_dbg)
#FCFLAGS= $(FCFLAGS_Check)

MKLROOT= /opt/intel
MKLOPT=  -I$(MKLROOT)/include
#------------------------------
# Source File
# -----------------------------

SRCS= Main.f90

#------------------------------
# Compiler & Option 
#------------------------------

FC=ifort
FCFLAGS_O3_xHost= -warn all -O3 -xHOST $(MKLFLAGS)
FCFLAGS_O3_omp_vec_dbg= -warn all -CB -fpe0 -traceback -g -O3 -xHOST -unroll-aggressive -vec-report1 -openmp -openmp-report1 $(MKLFLAGS)
FCFLAGS_O3_omp_vec= -warn all -O3 -xHOST -vec-report1 -qopenmp $(MKLFLAGS)
FCFLAGS_warn_all= -warn all -O0 -CB -fpe0 -traceback -g $(MKLOPT) 
FCFLAGS_O3_omp= -warn all -g -O3 -xHOST -unroll-aggressive -vec-report1 -openmp -openmp-report1 $(MKLOPT)  
FCFLAGS_O3= -warn all -O3 -xHOST -unroll-aggressive -vec-report1 -g 
FCFLAGS_O2= -O2 $(MKLFLAGS) 
FCFLAGS_O2_dbg= -O2 -check all -warn all -std -gen_interfaces -fpe0 -ftrapuv -traceback -g $(MKLFLAGS) 
FCFLAGS_O0= -warn all -O0 $(MKLOPT) 
FCFLAGS_f77= -O3 -xHOST -ip -static-intel 
FCFLAGS_fast_omp= -warn all -fast -openmp -openmp-report1 $(MKLOPT)  
FCFLAGS_Check= -O0 -check uninit -check all -warn all -std -gen_interfaces -fpe0 -ftrapuv -traceback $(MKLFLAGS)

#FCFLAGS_O3= -O3 -warn all -xHOST -unroll-aggressive -ip -static-intel -vec-report2 -openmp -openmp-report2 $(MKLOPT) 
#FCFLAGS_O3= -warn all -O3 -xHOST -vec-report1 
#FCFLAGS_O3= -warn all -O3 -xHOST -unroll-aggressive -vec-report1 -g 
#FCFLAGS_O3= -warn all -O3 -xHOST -unroll-aggressive -vec-report2 -openmp -openmp-report2 $(MKLOPT) $(MKLOPT) 
#FCFLAGS_O3= -O3 -xHOST -warn all -unroll-aggressive -opt_report -vec-report2 
#FCFLAGS= -O2 #-xHOST -openmp -openmp-report2 
#FCFLAGS= -debug -openmp -openmp-report2 

# Debag
#FCFLAGS= -O3 
# Warning

# Investigate tuning performance
#FCFLAGS= -O3 -xHOST -openmp -openmp-report2 -pg 
#FCFLAGS= -O3 -xHOST -pg 
# gprof ./a.out gmon.out > gmon.log

#FCFLAGS= -O3 -openmp -openmp-report2 -p -debug

TMP=$(SRCS:.f90=.o)
#OBJS=$(SRCS:.f=.o)

#------------------------------
# Library 
#------------------------------

MKLLIB= -L$(MKLROOT)/lib/intel64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm

LIB= $(ZMUMPS) $(DMUMPS) $(SCOTCHLIB) $(METISLIB) $(MKLLIB)
#LIB= $(MKLLIB)

.SUFFIXES: .f90



OBJS= $(Main_Objs) 

Main_Objs= \
Refine_Finite_Elements.o \
Plot_Configuration.o \
Plot_Mesh_Configuration.o \
Plot_Result_Element_Triangle.o \
$(Plot_Result_Objs) \
Plot_Result_and_Configuration_Element_Triangle.o $(Plot_Result_and_Configuration_Element_Triangle_Objs) \
Compute_Position_and_Element_PV.o $(Compute_Position_and_Element_PV_Objs) \
Compute_Poynting_Vector.o \
Plot_HF_T_OC.o \
Plot_Contour_T_OC.o \
Plot_Gradient_T_OC.o \
Plot_Gradient_Arrow_T_OC.o \
Detect_Element_inwhich_Point_Exists.o \
Plot_Single_Arrow.o \
Plot_Stream_Line.o Plot_Stream_Line_Symmetry.o \
Plot_DTFD_OC.o \
$(Plot_HF_T_OC_Objs) \
Plot_OC_PV.o \
Counter_Number_Node_Reference.o \
Operate_Data_Reference_Field.o $(Operate_Data_Reference_Field_Objs) \
Compute_Objective_Function.o \
Compute_Objective_Function_Flux.o \
Analyze_Steady_State_Heat_Conduction.o \
$(Analyze_Steady_State_Heat_Conduction_OBJS) \
Compute_Mode_Volume.o \
Compute_Q_Factor.o $(Compute_Q_Factor_OBJS) 

Plot_Result_Objs= \
Plot_Element_InHomogeneous.o \
Plot_Element_InHomogeneous_3_Node.o \
Compute_Position_Point_Plot_Result.o \
Compute_Position_Point_Plot_Scalor.o \
Plot_Gradation_Element.o \
Plot_Gradation_Element_Monocromatic.o \
Create_PlotData_Element_Homogeneous.o \
Plot_Color_Var_Upper_Side.o \
Plot_Color_Var_Right_Side.o \
Set_RGB.o

Plot_Result_and_Configuration_Element_Triangle_Objs= \
Plot_Element_InHomogeneous.o Plot_Element_InHomogeneous_3_Node.o Detect_Boundary_Data.o Set_RGB.o $(Plot_Result_Objs) $(Basic_Objs)

Compute_Position_and_Element_PV_Objs= \
Judge_Cross_Product.o $(Basic_Objs)

Plot_HF_T_OC_Objs= \
Plot_Arrow_Postscript.o \
Set_RGB.o \
$(Plot_Result_and_Configuration_Element_Triangle_Objs) $(Basic_Objs)

Compute_Line_Integral_Poynting_Vector_Objs= \
Compute_Unit_Normal_Vector_All.o $(Compute_Unit_Normal_Vector_All_Objs) \
Compute_Basis_Function_Triangle_Element.o \
$(Basic_Objs)

Compute_Unit_Normal_Vector_All_Objs= \
Compute_Unit_Normal_Vector_Single.o \
$(Basic_Objs)

          Operate_Data_Reference_Field_Objs= \
          Transform_Integer_2_Character.o \
          Transform_Real_2_Character.o $(Transform_Real_2_Character_Objs) \
          $(Basic_Objs)

               Transform_Real_2_Character_Objs= \
               Obtain_Number_Power_Decimal_Part.o
               $(Basic_Objs)

Analyze_Steady_State_Heat_Conduction_OBJS= \
Format_Local_Global_DMUMPS.o Format_Global_aa_DMUMPS.o \
Compute_Unit_Normal_Vector.o Compute_Incident_Wave.o \
Implement_Dirichlet_BC_real_CSRF.o Check_DMUMPS_Format.o \
Detect_Node_and_Element_Triangle_on_Thermal_Insulator_BC.o \
Detect_Node_and_Element_Triangle_for_Boundary_Integral.o \
Change_Direction_Vector_for_Unit_Normal_Vector.o \
MA57.o MA57_Component.o \
$(Basic_Objs)

Compute_Q_Factor_OBJS= \
Compute_Line_Integral_Poynting_Vector.o $(Compute_Line_Integral_Poynting_Vector_Objs) \
Detect_Line_on_Waveguide.o

Basic_Objs= match.o Output_Error.o Parameters.o

#------------------------------
# Compile 
#------------------------------

all:$(PROGS) 

$(PROGS): Main.o $(OBJS)
	$(FC) $(FCFLAGS) -o $@ $^ $(MKLLIB) $(LIB) 

Main.o: Main.f90 $(Main_Objs) 
	$(FC) $(FCFLAGS) -c $<

Refine_Finite_Elements.o: Refine_Finite_Elements.f90 
	$(FC) $(FCFLAGS) -c $<

Plot_Configuration.o : Plot_Configuration.f90 $(Basic_Objs)
	$(FC) $(FCFLAGS) -c $<

Plot_Mesh_Configuration.o : Plot_Mesh_Configuration.f90 $(Basic_Objs)
	$(FC) $(FCFLAGS) -c $<

Plot_Result_Element_Triangle.o : Plot_Result_Element_Triangle.f90 Plot_Element_InHomogeneous.o Plot_Element_InHomogeneous_3_Node.o $(Basic_Objs)
	$(FC) $(FCFLAGS) -c $<

Plot_Result_and_Configuration_Element_Triangle.o : Plot_Result_and_Configuration_Element_Triangle.f90 $(Plot_Result_and_Configuration_Element_Triangle_Objs)
	$(FC) $(FCFLAGS) -c $<

Judge_Cross_Product.o : Judge_Cross_Product.f90
	$(FC) $(FCFLAGS) -c $<

Compute_Position_and_Element_PV.o : Compute_Position_and_Element_PV.f90 $(Compute_Position_and_Element_PV_Objs)
	$(FC) $(FCFLAGS) -c $<

Compute_Poynting_Vector.o : Compute_Poynting_Vector.f90 Compute_Basis_Function_Triangle_Element.o
	$(FC) $(FCFLAGS) -c $<

Plot_HF_T_OC.o : Plot_HF_T_OC.f90 $(Plot_HF_T_OC_Objs)
	$(FC) $(FCFLAGS) -c $<

Plot_Contour_T_OC.o : Plot_Contour_T_OC.f90 $(Plot_HF_T_OC_Objs) 
	$(FC) $(FCFLAGS) -c $<

Plot_Gradient_T_OC.o : Plot_Gradient_T_OC.f90 Plot_Stream_Line.o Plot_Stream_Line_Symmetry.o $(Basic_Objs)
	$(FC) $(FCFLAGS) -c $<

Plot_Gradient_Arrow_T_OC.o: Plot_Gradient_Arrow_T_OC.f90 Plot_Arrow_Postscript.o Detect_Element_inwhich_Point_Exists.o Plot_Single_Arrow.o $(Basic_Objs)
	$(FC) $(FCFLAGS) -c $<

Detect_Element_inwhich_Point_Exists.o : Detect_Element_inwhich_Point_Exists.f90 Judge_Cross_Product.o $(Basic_Objs)
	$(FC) $(FCFLAGS) -c $<

Plot_Single_Arrow.o : Plot_Single_Arrow.f90 $(Basic_Objs)
	$(FC) $(FCFLAGS) -c $<

Plot_Stream_Line.o : Plot_Stream_Line.f90 $(Basic_Objs)
	$(FC) $(FCFLAGS_Check) -c $<
#	$(FC) $(FCFLAGS) -c $<

Plot_Stream_Line_Symmetry.o : Plot_Stream_Line_Symmetry.f90 $(Basic_Objs)
	$(FC) $(FCFLAGS) -c $<

Plot_DTFD_OC.o : Plot_DTFD_OC.f90 $(Plot_HF_T_OC_Objs)
	$(FC) $(FCFLAGS) -c $<

Plot_OC_PV.o : Plot_OC_PV.f90 $(Plot_HF_T_OC_Objs)
	$(FC) $(FCFLAGS) -c $<

Plot_Arrow_Postscript.o : Plot_Arrow_Postscript.f90 $(Basic_Objs) 
	$(FC) $(FCFLAGS) -c $<

Set_RGB.o : Set_RGB.f90 $(Basic_Objs)
	$(FC) $(FCFLAGS) -c $<

Detect_Boundary_Data.o : Detect_Boundary_Data.f90 $(Basic_Objs)
	$(FC) $(FCFLAGS) -c $<

Plot_Element_InHomogeneous.o : Plot_Element_InHomogeneous.f90 $(Basic_Objs)
	$(FC) $(FCFLAGS) -c $<

Plot_Element_InHomogeneous_3_Node.o : Plot_Element_InHomogeneous_3_Node.f90 $(Basic_Objs)
	$(FC) $(FCFLAGS) -c $<

Compute_Position_Point_Plot_Result.o : Compute_Position_Point_Plot_Result.f90 $(Basic_Objs)
	$(FC) $(FCFLAGS) -c $<

Compute_Position_Point_Plot_Scalor.o : Compute_Position_Point_Plot_Scalor.f90 $(Basic_Objs)
	$(FC) $(FCFLAGS) -c $<

Plot_Gradation_Element.o : Plot_Gradation_Element.f90 $(Basic_Objs)
	$(FC) $(FCFLAGS) -c $<

Plot_Gradation_Element_Monocromatic.o : Plot_Gradation_Element_Monocromatic.f90 $(Basic_Objs)
	$(FC) $(FCFLAGS) -c $<

Create_PlotData_Element_Homogeneous.o : Create_PlotData_Element_Homogeneous.f90 $(Basic_Objs)
	$(FC) $(FCFLAGS) -c $<

Plot_Color_Var_Upper_Side.o : Plot_Color_Var_Upper_Side.f90 $(Basic_Objs)
	$(FC) $(FCFLAGS) -c $<

Plot_Color_Var_Right_Side.o : Plot_Color_Var_Right_Side.f90 $(Basic_Objs)
	$(FC) $(FCFLAGS) -c $<

Counter_Number_Node_Reference.o : Counter_Number_Node_Reference.f90 $(Basic_Objs)
	$(FC) $(FCFLAGS) -c $<

Operate_Data_Reference_Field.o : Operate_Data_Reference_Field.f90 $(Operate_Data_Reference_Field_Objs)
	$(FC) $(FCFLAGS) -c $<

Transform_Integer_2_Character.o : Transform_Integer_2_Character.f90 $(Basic_Objs)
	$(FC) $(FCFLAGS) -c $<

Transform_Real_2_Character.o : Transform_Real_2_Character.f90 $(Transform_Real_2_Character_Objs)
	$(FC) $(FCFLAGS) -c $<

Obtain_Number_Power_Decimal_Part.o : Obtain_Number_Power_Decimal_Part.f90 $(Basic_Objs)
	$(FC) $(FCFLAGS) -c $<

Compute_Objective_Function.o : Compute_Objective_Function.f90 $(Basic_Objs)
	$(FC) $(FCFLAGS) -c $<

Compute_Objective_Function_Flux.o : Compute_Objective_Function_Flux.f90 $(Basic_Objs)
	$(FC) $(FCFLAGS) -c $<

Analyze_Steady_State_Heat_Conduction.o : Analyze_Steady_State_Heat_Conduction.f90 $(Analyze_Steady_State_Heat_Conduction_OBJS) $(ZMUMPSINCLUDE)
	$(FC) $(FCFLAGS) -c $<

Detect_Node_and_Element_Triangle_on_Thermal_Insulator_BC.o : Detect_Node_and_Element_Triangle_on_Thermal_Insulator_BC.f90 $(Basic_Objs)
	$(FC) $(FCFLAGS) -c $<

Detect_Node_and_Element_Triangle_for_Boundary_Integral.o : Detect_Node_and_Element_Triangle_for_Boundary_Integral.f90 $(Basic_Objs)
	$(FC) $(FCFLAGS) -c $<

Change_Direction_Vector_for_Unit_Normal_Vector.o : Change_Direction_Vector_for_Unit_Normal_Vector.f90 $(Basic_Objs)
	$(FC) $(FCFLAGS) -c $<

Compute_Mode_Volume.o : Compute_Mode_Volume.f90 
	$(FC) $(FCFLAGS) -c $<

Compute_Q_Factor.o : Compute_Q_Factor.f90 $(Compute_Q_Factor_OBJS)
	$(FC) $(FCFLAGS) -c $<

Compute_Line_Integral_Poynting_Vector.o : Compute_Line_Integral_Poynting_Vector.f90 
	$(FC) $(FCFLAGS) -c $<

Compute_Unit_Normal_Vector_All.o : Compute_Unit_Normal_Vector_All.f90 $(Compute_Unit_Normal_Vector_All_Objs)
	$(FC) $(FCFLAGS) -c $<

Compute_Unit_Normal_Vector_Single.o : Compute_Unit_Normal_Vector_Single.f90 $(Basic_Objs)
	$(FC) $(FCFLAGS) -c $<

Compute_Basis_Function_Triangle_Element.o : Compute_Basis_Function_Triangle_Element.f90 $(Basic_Objs)
	$(FC) $(FCFLAGS) -c $<

Detect_Line_on_Waveguide.o : Detect_Line_on_Waveguide.f90
	$(FC) $(FCFLAGS) -c $<

Format_Global_aa_DMUMPS.o : Format_Global_aa_DMUMPS.f90 $(Basic_Objs)
	$(FC) $(FCFLAGS) -c $<

Format_Local_Global_DMUMPS.o : Format_Local_Global_DMUMPS.f90 $(Basic_Objs)
	$(FC) $(FCFLAGS) -c $<

Compute_Incident_Wave.o : Compute_Incident_Wave.f90 $(Basic_Objs)
	$(FC) $(FCFLAGS) -c $<

Compute_Unit_Normal_Vector.o : Compute_Unit_Normal_Vector.f90 $(Basic_Objs)
	$(FC) $(FCFLAGS) -c $<

Implement_Dirichlet_BC_real_CSRF.o : Implement_Dirichlet_BC_real_CSRF.f90 $(Basic_Objs) 
	$(FC) $(FCFLAGS) -c $<

Check_DMUMPS_Format.o : Check_DMUMPS_Format.f90 $(Basic_Objs) 
	$(FC) $(FCFLAGS) -c $<

MA57.o : MA57.f90 MA57_Component.o $(Basic_Objs)
	$(FC) $(FCFLAGS) -c $<

MA57_Component.o : MA57_Component.f $(Basic_Objs)
	$(FC) $(FCFLAGS) -c $<

match.o : match.f90
	$(FC) -c $<

Output_Error.o : Output_Error.f90 
	$(FC) $(FCFLAGS) -c $<

Parameters.o: Parameters.f90
	$(FC) -c $<

#------------------------------
# Dependency 
#------------------------------

$(MUMPSHOME)/include/zmumps_struc.h: $(MUMPSHOME)/include/zmumps_root.h
$(MUMPSHOME)/include/dmumps_struc.h: $(MUMPSHOME)/include/dmumps_root.h

#------------------------------
# Clean 
#------------------------------

cleanF:
	rm -f fort.2? fort.3? fort.4? fort.5? fort.6? fort.7? fort.??? fort.???? gmon.out prof.txt  

cleanO:
	rm -f *.mod *.o a.out

clean:
	rm -f fort.?? fort.??? fort.???? *.dat *.txt *.fld *.inp *.el2 *ls2 filename_log *.core *.fort? ??2d00?? fort.100000 fort.* tmp* *mod.f90 Ea2d* Es2d* Et2d* Se2d* *.lsf *.cgp *.ps FEM_Data 

copy:
	cp ../NmakeInput/fort.1? .
	cp ../NmakeInput/fort.2? .
	cp ../NmakeInput/fort.3? .



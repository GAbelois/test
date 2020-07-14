# This is the main makefile.
# Use with GNU make.
# Relies on makedepf90 to do dependency lists:
# You can obtain this free program from:
#   http://www.helsinki.fi/~eedelman/makedepf90.html
# New $ARCH entries may be required, depending on your compiler.
# Also requires fpx3 fortran preprocessor available at
#   http://wwwuser.gwdg.de/~jbehren/fpx3.html

### user configurable options #####

FC          = mpiifort
#OPTFLAGS    = -O3 -ip -ftz -warn all -shared-intel -mcmodel=large
OPTFLAGS    = -O3 -ip -ftz -warn nointerfaces -shared-intel -mcmodel=large
#OPTFLAGS = -O0 -check bounds -g -debug all -traceback
WARNFLAGS   = -w
#### End User cenfigurable options ###

FFLAGS =  $(OPTFLAGS)
FLIBS =
EXECS =  test

INCLPATH = -I/usr/include
#LIBPATH = -L/usr/lib64
LIBPATH = -L/glade/u/home/chaoy/opt
LIBS = $(LIBPATH) -lfftw3f -lfftw3 -lm
#LIBPATH = -L/opt/fftw/2.1.5-intel/lib
#LIBS = $(LIBPATH) -ldrfftw -ldfftw -lm

# Directory for the .o files
OPATH = obj
# Directory for the .mod files, if your compiler generates them
# May want to just make this 'obj' as well
MPATH = mod

LDFLAGS = -threads
MODDIR = -I$(OPATH) -module $(MPATH)
#FFLAGS += $(MODDIR) $(INCLPATH)

FOBJ=000_types.o 001_param.o 002_fft.o 003_sim_param.o 004_derivatives.o 005_sgs.o 006_test_filter.o 007_topbc.o \
008_bottombc.o 009_intermediate.o 010_scalars.o 011_scalars2.o 012_io.o 013_intp.o 014_stokes_drift.o  015_canopy_param.o 016_ibm.o \
101_main.o 102_01_ic.o 102_initial.o 103_flosol.o 104_wallstress.o 105_dns_stress.o 106_01_std_dynamic.o \
106_02_scaledep_dynamic.o 106_03_lagrange_Sdep.o 106_04_lagrange_Ssim.o 106_05_interpolag_Ssim.o \
106_06_interpolag_Sdep.o 106_sgs_stag.o 107_divstress.o 108_convec.o 108_01_padd.o 109_press_stag_array.o 110_project.o \
111_scalar_dynamic.o 111_01_box_filter.o 111_02_dealias.o 201_canopy.o 202_parallel.o 203_fftn_mpi.o 204_intp_field.o 301_output.o 901_debug.o

all: $(EXECS)

$(EXECS): $(FOBJ)
	$(FC) $(FFLAGS) $(LDFLAGS) $(FOBJ) -o $(EXECS) $(LIBS)
000_types.o : ../code/000_types.f90
	$(FC) -c ../code/000_types.f90 $(FFLAGS)
001_param.o : ../code/001_param.f90 000_types.o
	$(FC) -c ../code/001_param.f90 $(FFLAGS)
002_fft.o : ../code/002_fft.f90 000_types.o 001_param.o
	$(FC) -c ../code/002_fft.f90 $(FFLAGS)
003_sim_param.o : ../code/003_sim_param.f90 000_types.o
	$(FC) -c ../code/003_sim_param.f90 $(FFLAGS)
004_derivatives.o : ../code/004_derivatives.f90 000_types.o 001_param.o 002_fft.o 
	$(FC) -c ../code/004_derivatives.f90 $(FFLAGS)
005_sgs.o : ../code/005_sgs.f90 000_types.o 001_param.o
	$(FC) -c ../code/005_sgs.f90 $(FFLAGS)
006_test_filter.o : ../code/006_test_filter.f90 000_types.o 001_param.o 002_fft.o 
	$(FC) -c ../code/006_test_filter.f90 $(FFLAGS)
007_topbc.o : ../code/007_topbc.f90 000_types.o 001_param.o
	$(FC) -c ../code/007_topbc.f90 $(FFLAGS)
008_bottombc.o : ../code/008_bottombc.f90 000_types.o 001_param.o 003_sim_param.o
	$(FC) -c ../code/008_bottombc.f90 $(FFLAGS)
009_intermediate.o : ../code/009_intermediate.f90 000_types.o 001_param.o
	$(FC) -c ../code/009_intermediate.f90 $(FFLAGS)
010_scalars.o : ../code/010_scalars.f90 000_types.o 001_param.o 003_sim_param.o 004_derivatives.o 005_sgs.o 006_test_filter.o 007_topbc.o 008_bottombc.o 014_stokes_drift.o 015_canopy_param.o 
	$(FC) -c ../code/010_scalars.f90 $(FFLAGS)
011_scalars2.o : ../code/011_scalars2.f90 000_types.o 001_param.o 009_intermediate.o 010_scalars.o 003_sim_param.o 005_sgs.o 008_bottombc.o
	$(FC) -c ../code/011_scalars2.f90 $(FFLAGS)
012_io.o : ../code/012_io.f90 000_types.o 001_param.o 002_fft.o 003_sim_param.o 005_sgs.o 008_bottombc.o 009_intermediate.o 010_scalars.o 011_scalars2.o
	$(FC) -c ../code/012_io.f90 $(FFLAGS)
013_intp.o:../code/013_intp.f90 000_types.o
	$(FC) -c ../code/013_intp.f90 $(FFLAGS) 
014_stokes_drift.o : ../code/014_stokes_drift.f90 000_types.o 001_param.o 
	$(FC) -c ../code/014_stokes_drift.f90 $(FFLAGS)
015_canopy_param.o : ../code/015_canopy_param.f90 000_types.o
	$(FC) -c ../code/015_canopy_param.f90 $(FFLAGS)
016_ibm.o : ../code/016_ibm.f90 000_types.o 001_param.o
	$(FC) -c ../code/016_ibm.f90 $(FFLAGS)
101_main.o : ../code/101_main.f90 000_types.o 001_param.o 003_sim_param.o 012_io.o
	$(FC) -c ../code/101_main.f90 $(FFLAGS)
102_initial.o : ../code/102_initial.f90 001_param.o 002_fft.o 006_test_filter.o 009_intermediate.o 011_scalars2.o 012_io.o 016_ibm.o
	$(FC) -c ../code/102_initial.f90 $(FFLAGS)
102_01_ic.o : ../code/102_01_ic.f90 000_types.o 001_param.o 003_sim_param.o 008_bottombc.o 012_io.o 015_canopy_param.o 
	$(FC) -c ../code/102_01_ic.f90 $(FFLAGS)
103_flosol.o : ../code/103_flosol.f90 001_param.o 003_sim_param.o 004_derivatives.o 007_topbc.o 009_intermediate.o 010_scalars.o 011_scalars2.o 012_io.o  
	$(FC) -c ../code/103_flosol.f90 $(FFLAGS)
104_wallstress.o : ../code/104_wallstress.f90 000_types.o 001_param.o 003_sim_param.o 006_test_filter.o 008_bottombc.o 014_stokes_drift.o 
	$(FC) -c ../code/104_wallstress.f90 $(FFLAGS)
105_dns_stress.o : ../code/105_dns_stress.f90 000_types.o 001_param.o 003_sim_param.o 
	$(FC) -c ../code/105_dns_stress.f90 $(FFLAGS)
106_sgs_stag.o : ../code/106_sgs_stag.f90 000_types.o 001_param.o 003_sim_param.o 005_sgs.o 006_test_filter.o
	$(FC) -c ../code/106_sgs_stag.f90 $(FFLAGS)
106_01_std_dynamic.o : ../code/106_01_std_dynamic.f90 000_types.o 001_param.o 003_sim_param.o 006_test_filter.o
	$(FC) -c ../code/106_01_std_dynamic.f90 $(FFLAGS)
106_02_scaledep_dynamic.o : ../code/106_02_scaledep_dynamic.f90 000_types.o 001_param.o 003_sim_param.o 005_sgs.o 006_test_filter.o
	$(FC) -c ../code/106_02_scaledep_dynamic.f90 $(FFLAGS)
106_03_lagrange_Sdep.o : ../code/106_03_lagrange_Sdep.f90 000_types.o 001_param.o 003_sim_param.o 005_sgs.o 006_test_filter.o
	$(FC) -c ../code/106_03_lagrange_Sdep.f90 $(FFLAGS)
106_04_lagrange_Ssim.o : ../code/106_04_lagrange_Ssim.f90 000_types.o 001_param.o 003_sim_param.o 005_sgs.o 006_test_filter.o
	$(FC) -c ../code/106_04_lagrange_Ssim.f90 $(FFLAGS)
106_05_interpolag_Ssim.o : ../code/106_05_interpolag_Ssim.f90 000_types.o 001_param.o 003_sim_param.o 005_sgs.o 006_test_filter.o
	$(FC) -c ../code/106_05_interpolag_Ssim.f90 $(FFLAGS)
106_06_interpolag_Sdep.o : ../code/106_06_interpolag_Sdep.f90 000_types.o 001_param.o 003_sim_param.o 005_sgs.o 006_test_filter.o
	$(FC) -c ../code/106_06_interpolag_Sdep.f90 $(FFLAGS)
107_divstress.o : ../code/107_divstress.f90 000_types.o 001_param.o 004_derivatives.o 
	$(FC) -c ../code/107_divstress.f90 $(FFLAGS)
108_convec.o : ../code/108_convec.f90 000_types.o 001_param.o 002_fft.o 003_sim_param.o 009_intermediate.o 014_stokes_drift.o
	$(FC) -c ../code/108_convec.f90 $(FFLAGS)
108_01_padd.o : ../code/108_01_padd.f90 000_types.o 001_param.o 
	$(FC) -c ../code/108_01_padd.f90 $(FFLAGS)
109_press_stag_array.o : ../code/109_press_stag_array.f90 000_types.o 001_param.o 002_fft.o 003_sim_param.o 008_bottombc.o 009_intermediate.o
	$(FC) -c ../code/109_press_stag_array.f90 $(FFLAGS)
110_project.o : ../code/110_project.f90 000_types.o 001_param.o 003_sim_param.o 012_io.o 015_canopy_param.o 
	$(FC) -c ../code/110_project.f90 $(FFLAGS)
111_scalar_dynamic.o : ../code/111_scalar_dynamic.f90 000_types.o 001_param.o 005_sgs.o 010_scalars.o        
	$(FC) -c ../code/111_scalar_dynamic.f90 $(FFLAGS)
111_01_box_filter.o : ../code/111_01_box_filter.f90 000_types.o 001_param.o  
	$(FC) -c ../code/111_01_box_filter.f90 $(FFLAGS)
111_02_dealias.o : ../code/111_02_dealias.f90 000_types.o 001_param.o 002_fft.o
	$(FC) -c ../code/111_02_dealias.f90 $(FFLAGS)
201_canopy.o : ../code/201_canopy.f90 000_types.o 001_param.o 003_sim_param.o 015_canopy_param.o 
	$(FC) -c ../code/201_canopy.f90 $(FFLAGS)
202_parallel.o : ../code/202_parallel.f90 000_types.o 001_param.o 003_sim_param.o 012_io.o 
	$(FC) -c ../code/202_parallel.f90 $(FFLAGS)
203_fftn_mpi.o : ../code/203_fftn_mpi.f90 000_types.o 001_param.o 
	$(FC) -c ../code/203_fftn_mpi.f90 $(FFLAGS)
204_intp_field.o : ../code/204_intp_field.f90 000_types.o 001_param.o 003_sim_param.o 005_sgs.o 008_bottombc.o 010_scalars.o 012_io.o 013_intp.o 
	$(FC) -c ../code/204_intp_field.f90 $(FFLAGS)
301_output.o : ../code/301_output.f90 000_types.o 001_param.o 003_sim_param.o 005_sgs.o 010_scalars.o 014_stokes_drift.o 
	$(FC) -c ../code/301_output.f90 $(FFLAGS)
901_debug.o : ../code/901_debug.f90 000_types.o 001_param.o 003_sim_param.o 005_sgs.o 014_stokes_drift.o
	$(FC) -c ../code/901_debug.f90 $(FFLAGS)


clean:
	rm -f *.o *.mod *.f90 *.log ./$(EXECS) PI* core $(EXECS)

#compile:
#	$(FC) $(FFLAGS) -c ../code/*.f90 

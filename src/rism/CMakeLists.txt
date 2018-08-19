#AMBER buildfile for RISM (Reference Interaction Site Model)
project(RISM Fortran)

set(RISM1D_SOURCES rism1d.F90)

set(RISM_LIB_SOURCES rism1d_c.F90 solvmdl_c.F90 rism1d_potential_c.F90
	rism1d_closure_c.F90 rism1d_kh_c.F90 rism1d_hnc_c.F90 rism1d_py_c.F90 
	rism1d_mv0_c.F90 rism1d_psen_c.F90 rism1d_nub_c.F90 rism1d_polyt_c.F90 rism1d_devb_c.F90 
	quaternion.F90 array_util.F90 rism_util.F90 rism_report_c.F90 
	rism3d_grid_c.F90 rism3d_closure_c.F90 
	rism3d_kh_c.F90 rism3d_hnc_c.F90 rism3d_psen_c.F90 
	rism3d_c.F90 rism3d_potential_c.F90 rism3d_restart.F90 
	rism3d_xyzv.F90 rism3d_opendx.F90 rism3d_ccp4.F90 
	rism3d_solvent_c.F90 rism3d_solute_c.F90 
	pubfft.F90 rism3d_fft.F90 rism_parm.F90 rism_io.F90 
	mdiis_orig_c.F90 mdiis_blas_c.F90 mdiis_blas2_c.F90 mdiis_c.F90 
	fce_c.F90 erfcfun.F90 safemem.F90 blend.F90 rism_timer_c.F90 
	getopts_c.F90 bspline.F90  fftw3.F90)
	
set(RISM_MPI_SOURCES quaternion.F90 array_util.F90 rism_util.F90 rism_report_c.F90 
	rism3d_grid_c.F90 rism3d_closure_c.F90 
	rism3d_kh_c.F90 rism3d_hnc_c.F90 rism3d_psen_c.F90 
	rism3d_c.F90 rism3d_potential_c.F90 rism3d_restart.F90 
	rism3d_xyzv.F90 rism3d_opendx.F90 rism3d_ccp4.F90 
	rism3d_solvent_c.F90 rism3d_solute_c.F90 
	pubfft.F90 rism3d_fft.F90 rism_parm.F90 rism_io.F90 
	mdiis_orig_c.F90 mdiis_blas_c.F90 mdiis_blas2_c.F90 mdiis_c.F90 
	fce_c.F90 erfcfun.F90 safemem.F90 blend.F90 rism_timer_c.F90 
	bspline.F90 fftw3.F90)
	
set(RISM_THERMO_SOURCES rism3d.thermo.F90 getopts_c.F90 array_util.F90)

set(RISM_ORAVE_SOURCES rism3d.orave.F90 getopts_c.F90 array_util.F90)

set(VOLSLICE_SOURCES volslice.F90)

set(SINGLE_POINT_SOURCES rism3d.snglpnt.c)

set(RISM_INTERFACE_SOURCES amber_rism_interface.F90)

apply_optimization_declarations(${RISM_LIB_SOURCES} ${RISM_MPI_SOURCES} ${RISM_THERMO_SOURCES} ${RISM_INTERFACE_SOURCES})	


#---------------------------------------------------------------------------------------------------------------------------------------------------------------------

#rism 3D library

#note: librism must be static for the time being due to undefined references to timer_start and timer_stop
add_library(rism STATIC ${RISM_LIB_SOURCES})
config_module_dirs(rism ${RISMLIB_MOD_DIR} ${AMBER_COMMON_MOD_DIR})

target_link_libraries(rism fftw netlib amber_common xblas)
install_libraries(rism)
make_pic_if_needed(rism)

#rism *sander* interface
add_library(sander_rism_interface OBJECT amber_rism_interface.F90)
set_property(TARGET sander_rism_interface PROPERTY COMPILE_DEFINITIONS SANDER RISM_CRDINTERP)
config_module_dirs(sander_rism_interface ${RISMLIB_SANDER_INTERFACE_DIR} ${RISMLIB_MOD_DIR} ${AMBER_COMMON_MOD_DIR}) #this goes into the sander module dir
add_dependencies(sander_rism_interface rism amber_common) 
make_pic_if_needed(sander_rism_interface)

#rism *sff* interface
add_library(sff_rism_interface OBJECT amber_rism_interface.F90)
set_property(TARGET sff_rism_interface PROPERTY COMPILE_DEFINITIONS NAB)
config_module_dirs(sff_rism_interface ${RISMLIB_SFF_INTERFACE_DIR} ${RISMLIB_MOD_DIR} ${AMBER_COMMON_MOD_DIR})
add_dependencies(sff_rism_interface rism amber_common) 
make_pic_if_needed(sff_rism_interface)

# rism3d.snglpnt
add_executable(rism3d.snglpnt ${SINGLE_POINT_SOURCES})
set_property(TARGET rism3d.snglpnt PROPERTY COMPILE_OPTIONS ${OPT_CFLAGS})
target_link_libraries(rism3d.snglpnt sff rism fftw netlib amber_common xblas)
install(TARGETS rism3d.snglpnt DESTINATION ${BINDIR})

#rism 1D executable
add_executable(rism1d ${RISM1D_SOURCES})
config_module_dirs(rism1d ${RISM1D_MOD_DIR} ${RISMLIB_MOD_DIR} ${AMBER_COMMON_MOD_DIR})

#rism thermo executable
add_executable(rism3d.thermo ${RISM_THERMO_SOURCES})
config_module_dirs(rism3d.thermo ${RISMTHERMO_MOD_DIR} ${RISMLIB_MOD_DIR} ${AMBER_COMMON_MOD_DIR})

#rism orientational averaging executable
add_executable(rism3d.orave ${RISM_ORAVE_SOURCES})
config_module_dirs(rism3d.orave ${RISMORAVE_MOD_DIR} ${RISMLIB_MOD_DIR} ${AMBER_COMMON_MOD_DIR})

#volslice executable
add_executable(volslice ${VOLSLICE_SOURCES})
config_module_dirs(volslice ${VOLSLICE_MOD_DIR} ${RISMLIB_MOD_DIR} ${AMBER_COMMON_MOD_DIR})

#handle link libraries
targets_link_libraries(rism1d rism3d.thermo rism3d.orave volslice LIBRARIES rism)

install(TARGETS rism1d rism3d.thermo rism3d.orave volslice DESTINATION ${BINDIR})

# MPI version
# --------------------------------------------------------------------
if(MPI)
	make_mpi_version(rism rism_mpi LANGUAGES Fortran SWAP_SOURCES ${RISM_LIB_SOURCES} TO ${RISM_MPI_SOURCES} INSTALL)
	config_module_dirs(rism_mpi ${RISMLIB_MPI_MOD_DIR} ${AMBER_COMMON_MOD_DIR})
	target_link_libraries(rism_mpi fftw_mpi)
	
	make_mpi_version(sander_rism_interface sander_rism_interface_mpi LANGUAGES Fortran)
	config_module_dirs(sander_rism_interface_mpi ${RISMLIB_SANDER_INTERFACE_MPI_DIR} ${RISMLIB_MPI_MOD_DIR} ${AMBER_COMMON_MOD_DIR}) #this goes into the sander module dir
	add_dependencies(sander_rism_interface rism_mpi)
	
	make_mpi_version(sff_rism_interface sff_rism_interface_mpi LANGUAGES Fortran)
	config_module_dirs(sff_rism_interface_mpi ${RISMLIB_SFF_INTERFACE_MPI_DIR} ${RISMLIB_MPI_MOD_DIR} ${AMBER_COMMON_MOD_DIR})
	add_dependencies(sff_rism_interface_mpi rism_mpi)
	
# rism3d.snglpnt.MPI  (assumes make_mpi_version works about the same 
#                     for both C and NAB)
	make_mpi_version(rism3d.snglpnt rism3d.snglpnt.MPI LANGUAGES C INSTALL)
	set_property(TARGET rism3d.snglpnt.MPI PROPERTY COMPILE_OPTIONS ${OPT_CFLAGS})

endif()

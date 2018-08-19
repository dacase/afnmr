project(etc C Fortran)

#list of all programs built by this file
#used for setting optimization flags
set(C_PROGRAMS elsize molsurf)
set(F_PROGRAMS resp tinker_to_amber new_crd_to_dyn new_to_old_crd fix_new_inpcrd_vel)

#perl scripts to install
set(PERL_PROGS process_mdout.perl process_minout.perl mdout2pymbar.pl)

#shell scripts to install
set(SHELL_PROGS sgldinfo.sh sgldwt.sh)
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------

# feel like using source lists here would make things LESS clear, not more.
add_executable(elsize elsize.c)
target_link_libraries(elsize m)

add_nab_executable(molsurf NAB_SOURCES molsurf.nab)
target_include_directories(molsurf PUBLIC ${NAB_HEADER_DIR})

add_executable(resp resp.F lapack.F)

add_executable(tinker_to_amber tinker_to_amber.F90 cspline.F)

add_executable(new_crd_to_dyn new_crd_to_dyn.F90)
target_link_libraries(new_crd_to_dyn amber_common)

add_executable(new_to_old_crd new_to_old_crd.F90)
target_link_libraries(new_to_old_crd amber_common)

add_executable(fix_new_inpcrd_vel fix_new_inpcrd_vel.F90)
target_link_libraries(fix_new_inpcrd_vel amber_common)

if(MPI)
	# used by the tests to count how many MPI processes are in use
	add_executable(numprocs numprocs.F90)
	target_link_libraries(numprocs mpi_fortran)
endif()

#optimization flags
set_property(TARGET ${C_PROGRAMS} PROPERTY COMPILE_OPTIONS ${NO_OPT_CFLAGS})
set_property(TARGET ${F_PROGRAMS} PROPERTY COMPILE_OPTIONS ${NO_OPT_FFLAGS})

#install everything
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------
install(TARGETS ${C_PROGRAMS} ${F_PROGRAMS} DESTINATION ${BINDIR})
install(PROGRAMS ${SHELL_PROGS} DESTINATION ${BINDIR})

if(BUILD_PERL)
	install(PROGRAMS ${PERL_PROGS} DESTINATION ${BINDIR})
endif()


if(MPI) 
	#install numprocs
	#---------------------------------------------------------------------------------------------------------------------------------------------------------------------

	# I'm not too wild about installing this executable to the source directory, but there's more or less no other way to get it to work
	# other than changing every single test script to point to a different file path.  
	
	# not allowed to install to absolute paths when making a Windows installer, so we do this to quash the errors.
	install(CODE "set(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION FALSE)")
	install(CODE "set(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION FALSE)")
	
	install(TARGETS numprocs DESTINATION ${CMAKE_SOURCE_DIR}/test )
	install(TARGETS numprocs DESTINATION ${CMAKE_SOURCE_DIR}/AmberTools/test)
	
	# the "Tests" component doesn't exist if test installation is disabled, so in that case, use serial instead.
	if(INSTALL_TESTS)
		set(NUMPROCS_COMPONENT Tests)
	else()
		set(NUMPROCS_COMPONENT Serial)
	endif()
	
	install(TARGETS numprocs DESTINATION ${CMAKE_INSTALL_POSTFIX}test COMPONENT ${NUMPROCS_COMPONENT})
	install(TARGETS numprocs DESTINATION ${CMAKE_INSTALL_POSTFIX}AmberTools/test COMPONENT ${NUMPROCS_COMPONENT})

endif()

#install python stuff
install_python_library()
     

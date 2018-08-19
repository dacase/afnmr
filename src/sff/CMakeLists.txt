project(sff)

#main SFF library
set(SFF_LIB_C_SOURCES conjgrad.c lmodC.c memutil.c nblist.c newton.c
	nmode.c prm.c rand2.c sasad.c sff.c xminC.c AmberNetcdf.c )

#smaller library of some sff functions with Fortran name mangling.
#used in some other amber programs, like sander and sqm.
set(SFF_FORTRAN_C_SOURCES lmodC.c rand2.c xminC.c) #sff_fortran C language sources

#This file is part of arpack.
#But, it is called by sff, and it calls the hessvec function in sff.
#So, it creates a circular dependency, wich is No Good with shared libraries.
#So, we add it to the sff and sff_fortran libraries
set(DSARPACK_SOURCES ../arpack/dsarpack.f)

if(NOT ARPACK_HAS_ARSECOND)
	# If arpack is missing the arsecond function, build it inside sff so that sff.c can use it
	list(APPEND SFF_LIB_C_SOURCES ../arpack/wallclock.c)
endif()

#------------------------------------------------------------------------------------------

#deal with the lexed code
set(LEXED_FILES_DIR "${CMAKE_BINARY_DIR}/AmberTools/src/sff")
set(LEX_INTERMEDIATES_DIR "${LEXED_FILES_DIR}/intermediates")
file(MAKE_DIRECTORY ${LEX_INTERMEDIATES_DIR})

set(INTERMEDIATE_LEXED_FILE "${LEX_INTERMEDIATES_DIR}/yy.mm_options.c")

set(FINAL_LEXED_FILE "${LEXED_FILES_DIR}/lex.mm_options.c")

flex_target(mm_options_lex mm_options.l ${INTERMEDIATE_LEXED_FILE} COMPILE_FLAGS -Pmmo)

buildtime_file_replace(${INTERMEDIATE_LEXED_FILE} ${FINAL_LEXED_FILE} TO_REPLACE stdout REPLACEMENT nabout)
add_custom_target(sff_lex ALL DEPENDS ${FINAL_LEXED_FILE}) #give CMake a reason to do the lexing

#------------------------------------------------------------------------------------------
check_all_symbols(float.h)

# use internal c9x-complex if needed
if(c9x-complex_INTERNAL)
	include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../c9x-complex)
endif()

#------------------------------------------------------------------------------------------

set_property(SOURCE ${SFF_LIB_C_SOURCES} ${SFF_FORTRAN_C_SOURCES} PROPERTY COMPILE_FLAGS "${OPT_CFLAGS_SPC}")
set_property(SOURCE ${SFF_FORTRAN_CXX_SOURCES} ${SFF_LIB_CXX_SOURCES} PROPERTY COMPILE_FLAGS "${OPT_CXXFLAGS_SPC}")

add_definitions(-Dflex)
include_directories(${LEXED_FILES_DIR}) #the lexed file is #include-ed into another c file

if(USE_FFT)
	add_definitions(-DRISMSFF)
	set(RISM_INTERFACE $<TARGET_OBJECTS:sff_rism_interface>)
	set(RISM_INTERFACE_MPI $<TARGET_OBJECTS:sff_rism_interface_mpi>)
else()
	set(RISM_INTERFACE "")
	set(RISM_INTERFACE_MPI "")
endif()
	
#------------------------------------------------------------------------------------------

# dsarpack needs to be in a seperate target so that we can avoif including SFF's debug.h, which conflicts with Arpack's debug.h

add_library(dsarpack_obj OBJECT ${DSARPACK_SOURCES})
target_compile_options(dsarpack_obj PRIVATE ${OPT_FFLAGS})
target_include_directories(dsarpack_obj PRIVATE ${CMAKE_SOURCE_DIR}/AmberTools/src/arpack)
make_pic_if_needed(dsarpack_obj)
 
add_library(sff STATIC ${SFF_LIB_C_SOURCES} $<TARGET_OBJECTS:dsarpack_obj> ${RISM_INTERFACE})
target_link_libraries(sff netcdf libpbsa netlib)
add_dependencies(sff sff_lex)
make_pic_if_needed(sff)
target_include_directories(sff INTERFACE .)



#these files have alternate code (enabled with the SQM flag) for SANDER and SQM.
add_library(sff_fortran ${SFF_FORTRAN_C_SOURCES} $<TARGET_OBJECTS:dsarpack_obj>)
set_property(TARGET sff_fortran PROPERTY COMPILE_DEFINITIONS SQM)
set_property(SOURCE mme_initC.c embed.c PROPERTY COMPILE_DEFINITIONS flex)
target_link_libraries(sff_fortran netlib)

add_executable(hcp_getpdb hcp_getpdb.c)

install(TARGETS hcp_getpdb DESTINATION ${BINDIR})
install_libraries(sff sff_fortran)

# install header
install(FILES sff.h DESTINATION ${INCDIR})

# parallelization
# --------------------------------------------------------------------

if(MPI)
	# NOTE: sff_mpi has a Fortran LANGUAGE because sff_rism_interface is written in Fortran
	make_mpi_version(sff sff_mpi LANGUAGES C Fortran INSTALL SWAP_SOURCES ${RISM_INTERFACE} TO ${RISM_INTERFACE_MPI})
	make_mpi_version(sff_fortran sff_fortran_mpi LANGUAGES C INSTALL)
	
	if(USE_FFT)
		target_link_libraries(sff_mpi rism_mpi)
	endif()
endif()

# wait until here to link to serial rism, so sff_mpi can link to parallel rism
if(USE_FFT)
	target_link_libraries(sff rism)
endif()

if(OPENMP)
	make_openmp_version(sff sff_omp LANGUAGES C INSTALL)
	target_compile_definitions(sff_omp PRIVATE OPENMP)
endif()

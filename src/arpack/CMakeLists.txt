#AMBER CMake buildfile for arpack
project(arpack C Fortran)

set(FORTRAN_SOURCE ivout.f dvout.f dgetv0.f dlaqrb.f dnaitr.f dnapps.f 
	dnaup2.f dnaupd.f dnconv.f dneigh.f dmout.f 
	dneupd.f dngets.f dsaitr.f dsapps.f dsaup2.f 
	dsaupd.f dsconv.f dseigt.f dsesrt.f dseupd.f 
	dsgets.f dsortc.f dsortr.f dstatn.f dstats.f dstqrb.f)
	
#note: dsarpack.f has been moved to sff_fortran to avoid a circular dependency between arpack and sff_fortran

set(C_SOURCE wallclock.c)

#language-specific compile flags
set_property(SOURCE ${FORTRAN_SOURCE} PROPERTY COMPILE_OPTIONS ${OPT_FFLAGS})
set_property(SOURCE ${C_SOURCE} PROPERTY COMPILE_OPTIONS ${OPT_CFLAGS})

add_library(arpack ${FORTRAN_SOURCE} ${C_SOURCE})

if(mkl_ENABLED)
	target_link_libraries(arpack ${MKL_FORTRAN_LIBRARIES})
else()
	target_link_libraries(arpack lapack blas)
endif()

install_libraries(arpack)
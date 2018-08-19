#AMBER CMake buildfile for leap
project(leap C)

add_subdirectory(src)

if(BUILD_GUI)
	add_subdirectory(pixmap)
endif()

if(NOT BUILD_HOST_TOOLS)
	if(TARGET_WINDOWS)
		# use batch file
		install(PROGRAMS tleap.bat DESTINATION ${BINDIR})
	else()
		install(PROGRAMS tleap DESTINATION ${BINDIR})
	endif()
endif()

if(BUILD_GUI)
	install(PROGRAMS xleap DESTINATION ${BINDIR})
endif()

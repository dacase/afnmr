
if(BUILD_GUI)
	include_directories(${X11_INCLUDE_DIR})
	
	add_subdirectory(Wc)
	add_subdirectory(Xmu)
	add_subdirectory(Xpm)
	add_subdirectory(Xraw)
endif()


add_subdirectory(leap)

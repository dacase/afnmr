#Compared to the original Makefile, I'm kind of reformatting stuff here, but I'm trying to keep the same categories
#ah well, at least it's well-commented


#generated sourcefiles
set(YACC_PARSER_FILE ${CMAKE_CURRENT_BINARY_DIR}/y.tab.c)
set(HELPTEXT_FILE ${CMAKE_CURRENT_BINARY_DIR}/helptext.c)

set(COMMON_SOURCE basics.c sysdepend.c stringExtra.c varArray.c getline.c avl.c)

set(GENERAL_SOURCE pdb_format.c pdb_read.c pdb_sprntf.c pdb_sscanf.c pdb_write.c
	vector.c zMatrix.c sort.c bag.c hash.c dictionary.c database.c
	nVector.c ring.c matrix.c fortran.c displayer.c)

set(OBJECT_HIERARCHY_SOURCE assoc.c atom.c byteArray.c 
	collection.c container.c internal.c list.c loop.c molecule.c 
	oDouble.c oInteger.c oString.c objekt.c parmSet.c residue.c 
	unit.c unitio.c tripos.c graphUtil.c select.c)

#  xaLeap and leap specific objects and tools
set(LEAP_SOURCE amber.c build.c elements.c library.c chirality.c
	minimizer.c model.c parmLib.c pdbFile.c mol2File.c 
	variables.c help.c octree.c mol3File.c tools.c
	commands.c mathop.c block.c restraint.c hybrid.c
	${YACC_PARSER_FILE} ${HELPTEXT_FILE})

set(TERMINAL_INTERFACE_SOURCE tLeap.c)

#  X-Windows LEaP interface specifics
set(XLEAP_SOURCE xTank.c xAction.c x3d.c xBasics.c)
 
#  X-Windows ATHENA widget set specifics
set(ATHENA_SOURCE xaLeapc.c xaUnitEditor.c
	xaTable.c xaAtomTable.c  XrawRegistr.c xaCommand.c xaTools.c 
	xaAtomParmTable.c xaBondParmTable.c xaAngleParmTable.c xaParmEditor.c
	xaTorsionParmTable.c xaImproperParmTable.c xaHBondParmTable.c)

set(MAIN_RESOURCE_FILE ${CMAKE_CURRENT_SOURCE_DIR}/xaLeap.rm4)
set(PROCESSED_RESOURCE_FILE ${CMAKE_CURRENT_BINARY_DIR}/XaLeap_wcl)

#Other files included into the main resource file.
#Only used to update the resource file when these files change.
set(RESOURCE_FILE_SOURCE xaLeap.rm4 xaDialog.rm4 xaTorsionParmTable.rm4 
	xaUnitEditor.rm4 xaCommand.rm4 xaAtomTable.rm4 xaUtilMessageFilter.rm4
	xaAtomParmTable.rm4 xaBondParmTable.rm4 xaAngleParmTable.rm4 
	xaImproperParmTable.rm4 xaHBondParmTable.rm4 xaParmEditor.rm4)

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Find a boolean type

if(BUILD_GUI)
	# X11/Xmd.h contains a BOOL definition that we have to check for
	set(CMAKE_REQUIRED_INCLUDES ${X11_INCLUDE_DIR})
endif()

check_all_includes(C sys/param.h sys/dir.h dirent.h string.h strings.h pwd.h stdbool.h sys/time.h Windows.h X11/Xmd.h)

check_all_symbols(ctype.h toupper tolower)

# look for BOOL type

set(BOOL_IN_X11_XMD FALSE)

if(HAVE_WINDOWS_H)
	set(CMAKE_EXTRA_INCLUDE_FILES Windows.h)
elseif(BUILD_GUI)
	set(CMAKE_EXTRA_INCLUDE_FILES X11/Xmd.h)
	set(BOOL_IN_X11_XMD TRUE)
else()
	set(CMAKE_EXTRA_INCLUDE_FILES "")
endif()
check_all_types(BOOL)


# look for _Bool type
if(HAVE_STDBOOL_H)
	# _Bool might be in stdbool.h, so include it if it exists
	set(CMAKE_EXTRA_INCLUDE_FILES stdbool.h)
else()
	set(CMAKE_EXTRA_INCLUDE_FILES "")
endif()
check_all_types(_Bool)
set(CMAKE_EXTRA_INCLUDE_FILES "")

# --------------------------------------------------------------------

configure_file(config.h.in ${CMAKE_CURRENT_BINARY_DIR}/config.h)

# tell basics.h to use config.h
add_definitions(-DCMAKE)

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------

#no optimization flags at all are given in the original makefile.  I'll assume it should be optimized.
add_compile_options(${OPT_CFLAGS})

include_directories(. ${CMAKE_CURRENT_BINARY_DIR})

#if we are building host tools, only build utilMakeHelp
if(BUILD_HOST_TOOLS)
	add_executable(utilMakeHelp utilMakeHelp.c)
	if(libm_ENABLED)
		target_link_libraries(utilMakeHelp m)
	endif()
	install(TARGETS utilMakeHelp DESTINATION ${BINDIR})
else() #build the leap programs
	include_directories($<TARGET_PROPERTY:netcdf,INTERFACE_INCLUDE_DIRECTORIES>)
	
	#---------------------------------------------------------------------------------------------------------------------------------------------------------------------

	byacc_target(leap_yacc parser.y ${CMAKE_CURRENT_BINARY_DIR})
	
	#---------------------------------------------------------------------------------------------------------------------------------------------------------------------
	add_library(leap_common_obj OBJECT ${COMMON_SOURCE})
	add_library(leap_general_obj OBJECT ${GENERAL_SOURCE})
	add_library(leap_hierarchy_obj OBJECT ${OBJECT_HIERARCHY_SOURCE})
	add_library(leap_leap_obj OBJECT ${LEAP_SOURCE})
	add_library(leap_terminal_obj OBJECT ${TERMINAL_INTERFACE_SOURCE})
	
	if(BUILD_GUI)
		add_library(leap_xleap_obj OBJECT ${XLEAP_SOURCE})
		add_library(leap_athena_obj OBJECT ${ATHENA_SOURCE})
	endif()
	
	#utility executables
	#---------------------------------------------------------------------------------------------------------------------------------------------------------------------
	
	# every leap program uses the same one or two libraries
	if(libm_ENABLED)
		set(LEAP_LIBRARIES netcdf libm)
	else()
		set(LEAP_LIBRARIES netcdf)
	endif()
	
	if(NOT USE_HOST_TOOLS)
		add_executable(utilMakeHelp utilMakeHelp.c)
		target_link_libraries(utilMakeHelp ${LEAP_LIBRARIES})
	endif()
	
	# These utility executables are turned off since they aren't used for anything.
	
	#add_executable(utilConvertWaterBox $<TARGET_OBJECTS:leap_common_obj> $<TARGET_OBJECTS:leap_general_obj> $<TARGET_OBJECTS:leap_hierarchy_obj> $<TARGET_OBJECTS:leap_leap_obj> utilConvertWaterBox.c)
	
	#add_executable(utilDescribeTopology $<TARGET_OBJECTS:leap_common_obj> $<TARGET_OBJECTS:leap_general_obj> $<TARGET_OBJECTS:leap_hierarchy_obj> $<TARGET_OBJECTS:leap_leap_obj> utilDescribeTopology.c)
	
	#add_executable(utilCrd2Off $<TARGET_OBJECTS:leap_common_obj> $<TARGET_OBJECTS:leap_general_obj> $<TARGET_OBJECTS:leap_hierarchy_obj> $<TARGET_OBJECTS:leap_leap_obj> utilCrd2Off.c)
	
	#add_executable(mol2_to_off $<TARGET_OBJECTS:leap_common_obj> $<TARGET_OBJECTS:leap_general_obj> $<TARGET_OBJECTS:leap_hierarchy_obj> $<TARGET_OBJECTS:leap_leap_obj> mol2_to_off.c)
	
	#add_executable(utilLib2Pdb $<TARGET_OBJECTS:leap_common_obj> $<TARGET_OBJECTS:leap_general_obj> $<TARGET_OBJECTS:leap_hierarchy_obj> $<TARGET_OBJECTS:leap_leap_obj> utilLib2Pdb.c)
	
	#targets_link_libraries(utilConvertWaterBox utilDescribeTopology utilCrd2Off mol2_to_off utilLib2Pdb LIBRARIES ${LEAP_LIBRARIES})
	
	#---------------------------------------------------------------------------------------------------------------------------------------------------------------------
	set(HELPTEXT_SOURCEFILE ${CMAKE_CURRENT_SOURCE_DIR}/helptext.text)
	add_custom_command(OUTPUT ${HELPTEXT_FILE}
		COMMAND ${RUNNABLE_utilMakeHelp} ${HELPTEXT_SOURCEFILE} ${HELPTEXT_FILE}
		DEPENDS ${HELPTEXT_SOURCEFILE})
	
	add_executable(teLeap $<TARGET_OBJECTS:leap_common_obj> $<TARGET_OBJECTS:leap_hierarchy_obj> $<TARGET_OBJECTS:leap_general_obj>
		$<TARGET_OBJECTS:leap_leap_obj> $<TARGET_OBJECTS:leap_terminal_obj>)

	target_link_libraries(teLeap ${LEAP_LIBRARIES})
	
	install(TARGETS teLeap DESTINATION ${BINDIR})
	
	if(BUILD_GUI)	
		add_custom_command(
			OUTPUT ${PROCESSED_RESOURCE_FILE}
			COMMAND ${M4}
			ARGS "-I${CMAKE_CURRENT_SOURCE_DIR}" ${MAIN_RESOURCE_FILE} > ${PROCESSED_RESOURCE_FILE} VERBATIM
			DEPENDS ${RESOURCE_FILE_SOURCE}
			VERBATIM
			COMMENT "Creating XLeap resource file")
			
		add_custom_target(xleap_resources ALL DEPENDS ${PROCESSED_RESOURCE_FILE})
		
		add_executable(xaLeap $<TARGET_OBJECTS:leap_common_obj> $<TARGET_OBJECTS:leap_hierarchy_obj> $<TARGET_OBJECTS:leap_general_obj>
		$<TARGET_OBJECTS:leap_leap_obj> $<TARGET_OBJECTS:leap_xleap_obj> $<TARGET_OBJECTS:leap_athena_obj>)
		target_link_libraries(xaLeap Xaw wcleap Xpm Xmu netcdf ${X11_Xt_LIB} ${X11_Xext_LIB} ${X11_SM_LIB}  ${X11_X11_LIB} ${X11_Xau_LIB} ${LEAP_LIBRARIES})
		
		install(TARGETS xaLeap DESTINATION ${BINDIR})
		install(FILES ${PROCESSED_RESOURCE_FILE} DESTINATION ${DATADIR})
	endif()
endif()

#AMBER CMake buildfile for cifparse
project(cifparse)

set(LEXED_FILES_DIR ${CMAKE_CURRENT_BINARY_DIR})

# Flex and bison
flex_target(CifLex cifparse.l "${LEXED_FILES_DIR}/lex.cif.c" COMPILE_FLAGS -Pcifp)
byacc_target(CifYacc cifparse.y ${LEXED_FILES_DIR} GENERATE_HEADER COMPILE_FLAGS -pcifp FILE_PREFIX cifp)
add_flex_bison_dependency(CifLex CifYacc)

#build the library

if(${CMAKE_C_COMPILER_ID} STREQUAL GNU)
	add_compile_options(-Wno-format)
endif()

add_compile_options(${NO_OPT_CFLAGS})

add_library(cifparse cifparse.c "${LEXED_FILES_DIR}/lex.cif.c" "${LEXED_FILES_DIR}/cifp.tab.c" "${LEXED_FILES_DIR}/cifp.tab.h")
target_include_directories(cifparse PUBLIC . ${LEXED_FILES_DIR})
install_libraries(cifparse)
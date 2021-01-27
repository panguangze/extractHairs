cmake_minimum_required(VERSION 3.4)

project(ExtractHAIRs)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_BUILD_TYPE DEBUG)
set(CMAKE_CXX_FLAGS_DEBUG "$ENV{CXX_FLAGS} -O0 -Wall -ggdb")

set(CMAKE_INCLUDE_CURRENT_DIR ON)

aux_source_directory(. EXTRACT_SRC)
aux_source_directory(./htslib HTS_SRC)
aux_source_directory(./samtools SAM_SRC)

set(LIBSAM libbam.a)
set(LIBHTS libhts.a)

add_custom_command(OUTPUT ${LIBHTS}
                   COMMAND make lib-static
                   WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/htslib
                   COMMENT "===================Making static htslib=======")
                   
add_custom_command(OUTPUT ${LIBSAM}
                   COMMAND make lib
                   WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/samtools
                   COMMENT "===================Making static libbam=======")

add_executable(ExtractHAIRs ${EXTRACT_SRC})

target_link_libraries(ExtractHAIRs /usr/local/lib/libhts.so z) 
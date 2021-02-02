cmake_minimum_required(VERSION 3.4)

set(CMAKE_CXX_COMPILE g++)
set(CMAKE_C_COMPILE gcc)

project(ExtractHAIRs)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_BUILD_TYPE RELEASE)
set(CMAKE_CXX_FLAGS_DEBUG "$ENV{CXX_FLAGS} -O0 -Wall -ggdb")
set(CMAKE_CXX_FLAGS_RELEASE "$ENV{CXXFLAGS} -O3 ")

set(CMAKE_INCLUDE_CURRENT_DIR ON)

aux_source_directory(. EXTRACT_SRC)
aux_source_directory(./htslib HTS_SRC)
aux_source_directory(./samtools SAM_SRC)

set(LIBSAM libbam.a)
set(LIBHTS libhts.a)


add_executable(ExtractHAIRs ${EXTRACT_SRC})

target_link_libraries(ExtractHAIRs /usr/local/lib/libhts.so z) 
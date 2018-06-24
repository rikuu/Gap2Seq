project(Gap2Seq)

cmake_minimum_required(VERSION 3.1)

################################################################################
# Define cmake modules directory
################################################################################
set (GATB_CORE_HOME  ${PROJECT_SOURCE_DIR}/thirdparty/gatb-core/gatb-core)
set (CMAKE_MODULE_PATH ${GATB_CORE_HOME}/cmake)

################################################################################
# THIRD PARTIES
################################################################################

# We don't want to install some GATB-CORE artifacts
#SET (GATB_CORE_EXCLUDE_TOOLS     1)
#SET (GATB_CORE_EXCLUDE_TESTS     1)
#SET (GATB_CORE_EXCLUDE_EXAMPLES  1)

# GATB CORE
include (GatbCore)

################################################################################
# Check for lambda expressions
################################################################################
INCLUDE(CheckCXXSourceCompiles)

IF (CMAKE_COMPILER_IS_GNUXX)
  SET(CMAKE_REQUIRED_FLAGS "-std=c++0x")
ELSEIF ("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
  SET(CMAKE_REQUIRED_FLAGS "-std=c++11 -stdlib=libc++")
ENDIF ()

CHECK_CXX_SOURCE_COMPILES(
"int main(int argc, char *args[])
{
int ret=0;
return ([&ret]() -> int { return ret; })();
}
" CXX_HAS_LAMBDA)

IF (NOT CXX_HAS_LAMBDA)
  MESSAGE (FATAL_ERROR "Lambda expressions not available. Use a newer C++ compiler (e.g. GCC Version 4.5 or greater)")
ENDIF (NOT CXX_HAS_LAMBDA)

################################################################################
# TOOL
################################################################################

# we get compilation definitions from the gatb-core part
add_definitions (${gatb-core-flags})

# htslib
set (HTSLIB_HOME ${PROJECT_SOURCE_DIR}/thirdparty/htslib/)
link_directories (${HTSLIB_HOME})

include_directories (SYSTEM src  ${gatb-core-includes} ${HTSLIB_HOME})

# we define the files to be compiled
file (GLOB  GAP2SEQ_SOURCES  src/Gap2Seq.cpp src/main.cpp)
file (GLOB  GAPCUTTER_SOURCES  src/GapCutter.cpp)
file (GLOB  GAPMERGER_SOURCES  src/GapMerger.cpp)
file (GLOB  READFILTER_SOURCES  src/ReadFilter.cpp)

# we define the artifact to be built: the project binary
add_executable ("Gap2Seq-core"  ${GAP2SEQ_SOURCES})
add_executable ("GapCutter"  ${GAPCUTTER_SOURCES})
add_executable ("GapMerger"  ${GAPMERGER_SOURCES})
add_executable ("ReadFilter"  ${READFILTER_SOURCES})

# we define which libraries to be linked with project binary
target_link_libraries ("Gap2Seq-core"  ${gatb-core-libraries})
target_link_libraries ("GapCutter"  ${gatb-core-libraries})
target_link_libraries ("GapMerger"  ${gatb-core-libraries})
target_link_libraries ("ReadFilter"  ${gatb-core-libraries} hts)

# we add a compiler flag
set (CXX_WARNINGS "-Wall -Wextra")
set (GCC_COMPILE_MAP_FLAGS "-std=c++0x -g")
set (CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} ${GCC_COMPILE_MAP_FLAGS} ${CXX_WARNINGS}")

# we add the main executable script
configure_file ("${CMAKE_SOURCE_DIR}/src/Gap2Seq.py" "${CMAKE_CURRENT_BINARY_DIR}/Gap2Seq" COPYONLY)

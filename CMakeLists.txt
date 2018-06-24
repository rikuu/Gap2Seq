project(Gap2Seq)

cmake_minimum_required(VERSION 3.1.0)

################################################################################
# Define the version number.
################################################################################
SET (GAP2SEQ_VERSION_MAJOR 3)
SET (GAP2SEQ_VERSION_MINOR 1)
SET (GAP2SEQ_VERSION_PATCH 0)

SET (GAP2SEQ_VERSION ${GAP2SEQ_VERSION_MAJOR}.${GAP2SEQ_VERSION_MINOR}.${GAP2SEQ_VERSION_PATCH})

################################################################################
# THIRD PARTIES
################################################################################

# We don't want to install some GATB-CORE artifacts
SET (GATB_CORE_EXCLUDE_TOOLS     1)
SET (GATB_CORE_EXCLUDE_TESTS     1)
SET (GATB_CORE_INCLUDE_EXAMPLES  0)

SET (HTSLIB_HOME ${PROJECT_SOURCE_DIR}/thirdparty/htslib/)
SET (GATB_CORE_HOME  ${PROJECT_SOURCE_DIR}/thirdparty/gatb-core/gatb-core)
SET (CMAKE_MODULE_PATH ${GATB_CORE_HOME}/cmake)

# GATB CORE
include (GatbCore)

################################################################################
# TOOL
################################################################################

add_definitions (${gatb-core-flags})

include_directories (SYSTEM src ${gatb-core-includes} ${HTSLIB_HOME})

link_directories (${gatb-core-extra-libraries-path})
link_directories (${HTSLIB_HOME})

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

set (CXX_WARNINGS "-Wall -Wextra")
#set (GCC_COMPILE_MAP_FLAGS "-std=c++0x -g")
set (CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} ${GCC_COMPILE_MAP_FLAGS} ${CXX_WARNINGS}")

# we add the main executable script
configure_file ("${CMAKE_SOURCE_DIR}/src/Gap2Seq.py" "${CMAKE_CURRENT_BINARY_DIR}/Gap2Seq" COPYONLY)

project(Gap2Seq)

cmake_minimum_required(VERSION 3.1.0)

################################################################################
# THIRD PARTIES
################################################################################

set(GATB_CORE_HOME ${PROJECT_SOURCE_DIR}/thirdparty/gatb-core/gatb-core)

# GATB-core includes parts of boost (critically not graph) which likely
# conflicts with whatever version we can find
find_package(Boost REQUIRED COMPONENTS graph)
include_directories(SYSTEM ${Boost_INCLUDE_DIRS})
file(REMOVE_RECURSE ${GATB_CORE_HOME}/thirdparty/boost)

set(HTSLIB_HOME ${PROJECT_SOURCE_DIR}/thirdparty/htslib)
include_directories(BEFORE SYSTEM ${HTSLIB_HOME})
link_directories(${HTSLIB_HOME})

# GATB CORE
set(GATB_CORE_EXCLUDE_TOOLS 1)
set(GATB_CORE_EXCLUDE_TESTS 1)
set(GATB_CORE_INCLUDE_EXAMPLES 0)

set(CMAKE_MODULE_PATH ${GATB_CORE_HOME}/cmake)

include(GatbCore)
add_definitions(${gatb-core-flags})
include_directories(${gatb-core-includes})
link_directories(${gatb-core-extra-libraries-path})

################################################################################
# TOOL
################################################################################

add_executable("Gap2Seq-core" src/Gap2Seq.cpp src/main.cpp)
add_executable("GapCutter" src/GapCutter.cpp)
add_executable("GapMerger" src/GapMerger.cpp)
add_executable("ReadFilter" src/ReadFilter.cpp)

target_link_libraries("Gap2Seq-core" ${gatb-core-libraries})
target_link_libraries("GapCutter" ${gatb-core-libraries})
target_link_libraries("GapMerger" ${gatb-core-libraries})
target_link_libraries("ReadFilter" ${gatb-core-libraries} hts)

set(CXX_WARNINGS "-Wall") # -Wextra
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CXX_WARNINGS}")

# we add the main executable script
configure_file("${CMAKE_SOURCE_DIR}/src/Gap2Seq.py" "${CMAKE_CURRENT_BINARY_DIR}/Gap2Seq" COPYONLY)

cmake_minimum_required(VERSION 3.12...3.27)
project(livpropa CXX)

# Python and SWIG
find_package(Python COMPONENTS Interpreter Development)
find_package(SWIG REQUIRED)

# general options
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
set(CRPropa_SWIG_FILE "crpropa.i")
set(SWIG_MODE_FLAG "")

# use SWIG's builtin
option(ENABLE_SWIG_BUILTIN "Use SWIG builtin option" ON) # Segmentation errors on OSX if false
if(ENABLE_SWIG_BUILTIN)
  set(CRPropa_SWIG_FILE "crpropa-builtin.i")
  set(SWIG_MODE_FLAG "-builtin")
endif(ENABLE_SWIG_BUILTIN)


# get information on CRPropa
execute_process(COMMAND ${Python_EXECUTABLE} "${CMAKE_CURRENT_SOURCE_DIR}/python/findCRPropa.py" swig_interface OUTPUT_VARIABLE crpropaSwigPath)
execute_process(COMMAND ${Python_EXECUTABLE} "${CMAKE_CURRENT_SOURCE_DIR}/python/findCRPropa.py" install_prefix OUTPUT_VARIABLE CRPropa_INSTALL_PREFIX)

find_path(CRPropa_SWIG_INTERFACE_FILES ${CRPropa_SWIG_FILE} HINTS share/crpropa share/crpropa/swig_interface python/ ${crpropaSwigPath})
find_path(CRPropa_INCLUDE_PATH CRPropa.h HINTS crpropa/ include/ include/crpropa ${CRPropa_INSTALL_PREFIX}/include/)
find_library(CRPropa_LIBRARY NAMES crpropa libcrpropa HINTS crpropa/ lib/crpropa crpropa/lib ${CRPropa_INSTALL_PREFIX}/lib/)

message(STATUS "CRPropa install prefix: ${CRPropa_INSTALL_PREFIX}")
message(STATUS "CRPropa SWIG interface files: ${CRPropa_SWIG_INTERFACE_FILES}")
message(STATUS "CRPropa include path : ${CRPropa_INCLUDE_PATH}")
message(STATUS "CRPropa library: ${CRPropa_LIBRARY}")


# Python won't be found for lower versions of cmake due to changes in cmake's FindPython file:
#   https://cmake.org/cmake/help/latest/module/FindPython.html#module:FindPython
include_directories(${Python_INCLUDE_DIRS})

# Python
if(Python_VERSION VERSION_GREATER 3.9)
  execute_process(COMMAND ${Python_EXECUTABLE} -c "import sys; import sysconfig; sys.stdout.write(sysconfig.get_path(\'stdlib\'))" OUTPUT_VARIABLE Python_SITE_PACKAGES OUTPUT_STRIP_TRAILING_WHITESPACE)
else()
  execute_process(COMMAND ${Python_EXECUTABLE} -c "import sys; from distutils import sysconfig; sys.stdout.write(sysconfig.get_python_lib(1,0,prefix='${CMAKE_INSTALL_PREFIX}'))" OUTPUT_VARIABLE Python_SITE_PACKAGES OUTPUT_STRIP_TRAILING_WHITESPACE)
endif()

message(STATUS "Python version: ${Python_VERSION}")
message(STATUS "Python headers: ${Python_INCLUDE_DIRS}")
message(STATUS "Python library: ${Python_LIBRARIES}, ${Python_LIBRARY}")
message(STATUS "Python site packages: ${Python_SITE_PACKAGES}")


# kiss library (provided by CRPropa)
include_directories(include ${CRPropa_INSTALL_PREFIX}/libs/kiss)

# Eigen library
include_directories(include external/eigen3)
set(Eigen_INCLUDE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/external/eigen3)


# SWIG plugin
file(GLOB_RECURSE CRPROPA_SWIG_INPUTS python/*.i)
set_source_files_properties(${CMAKE_CURRENT_BINARY_DIR}/livpropa_wrap.cxx PROPERTIES GENERATED true)
add_custom_command(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/livpropa_wrap.cxx
  COMMAND swig ${SWIG_MODE_FLAG} -c++ -python ${PY3} -I${CRPropa_SWIG_INTERFACE_FILES} -I${CRPropa_INCLUDE_PATH} -I${CMAKE_SOURCE_DIR}/include -I${Eigen_INCLUDE_PATH} -I${Python_INCLUDE_DIRS} ${CRPROPA_SWIG_DEFINES} -dirprot -o ${CMAKE_CURRENT_BINARY_DIR}/livpropa_wrap.cxx -outdir ${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/python/livpropa.i
)
include_directories(include ${CMAKE_CURRENT_BINARY_DIR}/include ${CRPropa_INCLUDE_PATH} ${Python_INCLUDE_DIRS} ${CMAKE_CURRENT_SOURCE_DIR}/../include)


# build interface and source in one library
add_library(livpropa MODULE 
	${CMAKE_CURRENT_BINARY_DIR}/livpropa_wrap.cxx
	${CMAKE_CURRENT_SOURCE_DIR}/src/Data.cc
	${CMAKE_CURRENT_SOURCE_DIR}/src/Kinematics.cc
	${CMAKE_CURRENT_SOURCE_DIR}/src/PairProduction.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/InverseComptonScattering.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/PhotonDecay.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src/VacuumCherenkov.cc
  )
target_link_libraries(livpropa ${CRPropa_LIBRARY} ${Python_LIBRARY} ${Python_LIBRARIES})
set_target_properties(livpropa PROPERTIES PREFIX "")
set_target_properties(livpropa PROPERTIES COMPILE_FLAGS "-w")
set_target_properties(livpropa PROPERTIES OUTPUT_NAME "_livpropa")


# ----------------------------------------------------------------------------
# Tests
# ----------------------------------------------------------------------------

option(ENABLE_TESTING "Build tests and enable test target" ON)
if(ENABLE_TESTING)
  set(GTEST_INCLUDE_PATH ${CRPropa_INSTALL_PREFIX}/../libs/gtest/include)
  include_directories(${GTEST_INCLUDE_PATH})

  if(APPLE)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DGTEST_USE_OWN_TR1_TUPLE=1")
  endif(APPLE)

  link_directories(${CRPropa_INSTALL_PREFIX}/libs/gtest)
  message(STATUS "Google testing suite (from CRPropa): ${CRPropa_INSTALL_PREFIX}/libs/gtest/include")
  enable_testing()
endif(ENABLE_TESTING)


# ----------------------------------------------------------------------------
# Install
# ----------------------------------------------------------------------------

add_definitions(-DLIVpropa_INSTALL_PREFIX="${CMAKE_INSTALL_PREFIX}")

install(TARGETS livpropa DESTINATION lib)

install(DIRECTORY external/eigen3/Eigen DESTINATION include)
install(DIRECTORY include/ DESTINATION include FILES_MATCHING PATTERN "*.h")
install(DIRECTORY ${CMAKE_BINARY_DIR}/include/ DESTINATION include FILES_MATCHING PATTERN "*.h")
install(DIRECTORY ${CMAKE_BINARY_DIR}/data/ DESTINATION share/livpropa/ PATTERN ".git" EXCLUDE)

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/data DESTINATION ${CMAKE_INSTALL_PREFIX}/.)


cmake_minimum_required(VERSION 3.12...3.27)

project(livpropa CXX)

# Python and SWIG
find_package(Python COMPONENTS Interpreter Development)
find_package(SWIG REQUIRED)

# general options
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
set(crpropaSwigFile "crpropa.i")
set(swigModeFlag "")

# use SWIG's builtin
option(ENABLE_SWIG_BUILTIN "Use SWIG builtin option" ON) # Segmentation errors on OSX if false
if(ENABLE_SWIG_BUILTIN)
  set(crpropaSwigFile "crpropa-builtin.i")
  set(swigModeFlag "-builtin")
endif(ENABLE_SWIG_BUILTIN)


# get information on CRPropa
execute_process(COMMAND ${PYTHON_EXECUTABLE} "${CMAKE_CURRENT_SOURCE_DIR}/python/findCRPropa.py" swig_interface OUTPUT_VARIABLE crpropaSwigPath)
execute_process(COMMAND ${PYTHON_EXECUTABLE} "${CMAKE_CURRENT_SOURCE_DIR}/python/findCRPropa.py" install_prefix OUTPUT_VARIABLE crpropaInstallPrefix)

find_path(CRPROPA_SWIG_INTERFACE_FILES ${crpropaSwigFile} HINTS share/crpropa share/crpropa/swig_interface python/ ${crpropaSwigPath})
find_path(CRPROPA_INCLUDE_PATH CRPropa.h HINTS crpropa/ include/ include/crpropa ${crpropaInstallPrefix}/include/)
find_library(CRPROPA_LIBRARY NAMES crpropa libcrpropa HINTS crpropa/ lib/crpropa crpropa/lib ${crpropaInstallPrefix}/lib/)

message(STATUS "CRPropa install prefix: ${crpropaInstallPrefix}")
message(STATUS "CRPropa SWIG interface files: ${CRPROPA_SWIG_INTERFACE_FILES}")
message(STATUS "CRPropa include path : ${CRPROPA_INCLUDE_PATH}")
message(STATUS "CRPropa library: ${CRPROPA_LIBRARY}")



if(CMAKE_VERSION VERSION_LESS "3.12")
  # Use built-in Python finder; otherwise use CMake's
  include(python/Python.cmake)
  set(PYTHON_VERSION ${PYTHON_VERSION_STRING})
else()
  set(PYTHON_INCLUDE_PATH ${Python_INCLUDE_DIRS})
  set(PYTHON_VERSION ${Python_VERSION})
  set(PYTHON_LIBRARY ${Python_LIBRARY})
endif()
include_directories(${PYTHON_INCLUDE_PATH})

# Python
if(PYTHON_VERSION VERSION_GREATER 3.9)
  execute_process(COMMAND ${PYTHON_EXECUTABLE} -c "import sys; import sysconfig; sys.stdout.write(sysconfig.get_path(\'stdlib\'))" OUTPUT_VARIABLE PYTHON_SITE_PACKAGES OUTPUT_STRIP_TRAILING_WHITESPACE)
else()
  execute_process(COMMAND ${PYTHON_EXECUTABLE} -c "import sys; from distutils import sysconfig; sys.stdout.write(sysconfig.get_python_lib(1,0,prefix='${CMAKE_INSTALL_PREFIX}'))" OUTPUT_VARIABLE PYTHON_SITE_PACKAGES OUTPUT_STRIP_TRAILING_WHITESPACE)
endif()

message(STATUS "Python version: ${PYTHON_VERSION}")
message(STATUS "Python headers: ${PYTHON_INCLUDE_PATH}")
message(STATUS "Python library: ${PYTHON_LIBRARY}")
message(STATUS "Python site packages: ${PYTHON_SITE_PACKAGES}")


# kiss library (provided by CRPropa)
include_directories(include ${crpropaInstallPrefix}/libs/kiss)

# Eigen library
include_directories(include external/eigen3)

# SWIG plugin
file(GLOB_RECURSE CRPROPA_SWIG_INPUTS python/*.i)
set_source_files_properties(${CMAKE_CURRENT_BINARY_DIR}/livpropa_wrap.cxx PROPERTIES GENERATED true)
add_custom_command(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/livpropa_wrap.cxx
  COMMAND swig ${swigModeFlag} -c++ -python ${PY3} -I${CRPROPA_SWIG_INTERFACE_FILES} -I${CRPROPA_INCLUDE_PATH} -I${CMAKE_SOURCE_DIR}/include -I${Python_INCLUDE_DIRS} ${CRPROPA_SWIG_DEFINES} -dirprot -o ${CMAKE_CURRENT_BINARY_DIR}/livpropa_wrap.cxx -outdir ${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/python/livpropa.i
)
include_directories(include ${CMAKE_CURRENT_BINARY_DIR}/include ${CRPROPA_INCLUDE_PATH} ${Python_INCLUDE_DIRS} ${CMAKE_CURRENT_SOURCE_DIR}/../include)


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
target_link_libraries(livpropa ${CRPROPA_LIBRARY} ${PYTHON_LIBRARY})
set_target_properties(livpropa PROPERTIES PREFIX "")
set_target_properties(livpropa PROPERTIES COMPILE_FLAGS "-w")
set_target_properties(livpropa PROPERTIES OUTPUT_NAME "_livpropa")





# ----------------------------------------------------------------------------
# Tests
# ----------------------------------------------------------------------------
option(ENABLE_TESTING "Build tests and enable test target" ON)
if(ENABLE_TESTING)
  set(GTEST_INCLUDE_PATH ${CRPROPA_INSTALL_PREFIX}/../libs/gtest/include)
  include_directories(${GTEST_INCLUDE_PATH})
  link_directories(${CRPROPA_INSTALL_PREFIX}/libs/gtest)
  message(STATUS "Google testing suite (from CRPropa): ${CRPROPA_INSTALL_PREFIX}/libs/gtest/include")
  enable_testing()
endif(ENABLE_TESTING)


# ----------------------------------------------------------------------------
# Install
# ----------------------------------------------------------------------------
add_definitions(-DLIVPROPA_INSTALL_PREFIX="${CMAKE_INSTALL_PREFIX}")
install(TARGETS livpropa DESTINATION lib)
install(DIRECTORY include/ DESTINATION include FILES_MATCHING PATTERN "*.h")
install(DIRECTORY ${CMAKE_BINARY_DIR}/include/ DESTINATION include FILES_MATCHING PATTERN "*.h")

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/data DESTINATION ${CMAKE_INSTALL_PREFIX}/.)


install(DIRECTORY ${CMAKE_BINARY_DIR}/data/ DESTINATION share/livpropa/ PATTERN ".git" EXCLUDE)
# install(DIRECTORY libs/kiss/include/ DESTINATION include)



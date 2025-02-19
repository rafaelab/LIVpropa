cmake_minimum_required(VERSION 3.14)
project(livpropa CXX)

message(STATUS "Building LIVropa...")


# Set default build-type to release to enable performance improvements
if(NOT CMAKE_BUILD_TYPE)
	set(CMAKE_BUILD_TYPE Release)
endif(NOT CMAKE_BUILD_TYPE)
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")


# compilation flags
set(CMAKE_CXX_FLAGS "-std=c++14")

# Enable extra warnings on debug builds
SET(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS} -Wextra")

# # workaround for LLVM > 15
# if(CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang" OR CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
# 	if(CMAKE_CXX_COMPILER_VERSION VERSION_GREATER_EQUAL 15.0)
# 		set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D_LIBCPP_ENABLE_CXX17_REMOVED_UNARY_BINARY_FUNCTION")
# 	endif(CMAKE_CXX_COMPILER_VERSION VERSION_GREATER_EQUAL 15.0)
# endif(CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang" OR CMAKE_CXX_COMPILER_ID STREQUAL "Clang")

# display compiler information
message(STATUS "C++ compiler: ${CMAKE_CXX_COMPILER} (version ${CMAKE_CXX_COMPILER_VERSION})")
message(STATUS "C++ compilation flags: ${CMAKE_CXX_FLAGS}")


# find CRPropa
include(cmake/FindCRPropa.cmake)


# Eigen library
include_directories(include external/eigen3)
set(Eigen_INCLUDE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/external/eigen3)

# assemble includes (using own Eigen version)
include_directories(include 
	${CRPropa_kiss_INCLUDE_DIR}
	${CRPropa_HepID_INCLUDE_DIR}
	${CRPropa_INCLUDE_DIR}  
	# external/eigen3
	${CMAKE_CURRENT_BINARY_DIR}/include 
	${CMAKE_CURRENT_SOURCE_DIR}/../include
	)


# source files for compilation
set(SOURCE_FILES
	${CMAKE_CURRENT_SOURCE_DIR}/src/Data.cc
	${CMAKE_CURRENT_SOURCE_DIR}/src/Histogram.cc
	${CMAKE_CURRENT_SOURCE_DIR}/src/Sampler.cc
	${CMAKE_CURRENT_SOURCE_DIR}/src/Weighter.cc
	${CMAKE_CURRENT_SOURCE_DIR}/src/Kinematics.cc
	${CMAKE_CURRENT_SOURCE_DIR}/src/PairProduction.cc
	${CMAKE_CURRENT_SOURCE_DIR}/src/InverseComptonScattering.cc
	${CMAKE_CURRENT_SOURCE_DIR}/src/PhotonDecay.cc
	${CMAKE_CURRENT_SOURCE_DIR}/src/VacuumCherenkov.cc
	)


# decide on whether to use Python or not
option(ENABLE_PYTHON "Create python library via SWIG" On)
message(STATUS "Build Python bindings via SWIG: ${ENABLE_PYTHON}")

# compile with Python support
if(ENABLE_PYTHON)
	# find Python
	find_package(Python COMPONENTS Interpreter Development NumPy)

	if(Python_FOUND)
		message(STATUS "Python version: ${Python_VERSION}")
		message(STATUS "Python executable: ${Python_EXECUTABLE}")
		message(STATUS "Python headers: ${Python_INCLUDE_DIRS}")
		message(STATUS "Python library: ${Python_LIBRARY}")
		message(STATUS "Python site-packges: ${Python_SITE_PACKAGES}")
	else(Python_FOUND)
		set(Python_HELPER_EXECUTABLE "Python_EXECUTABLE (path to python executable)")
		set(Python_HELPER_INCLUDE_DIRS "Python_INCLUDE_DIRS (directory where Python.h is located)")
		set(Python_HELPER_LIBRARY "Python_LIBRARY (path to file libpython.so or libpython.dylib [in OSX]; sometimes contains version before extension)")
		message(STATUS 	"Python: NOT FOUND. Try setting the relevant flags manually: \n\t${Python_HELPER_EXECUTABLE} \n\t${Python_HELPER_LIBRARY} \n\t${Python_HELPER_INCLUDE_DIRS}")
	endif(Python_FOUND)


	if(Python_NumPy_FOUND)
		set(CMAKE_SWIG_FLAGS -DWITHNUMPY ${CRP})
		list(APPEND CRPROPA_SWIG_DEFINES -DWITHNUMPY)
		include_directories(${Python_NumPy_INCLUDE_DIRS})
		message(STATUS "NumPy headers: ${Python_NumPy_INCLUDE_DIRS} (version ${Python_NumPy_VERSION})")
	elseif(Python_NumPy_FOUND)
		message(STATUS "NumPy: NOT FOUND. The code might work just fine without NumPy, but some features in Python might not be straightfoward to use.")
	endif(Python_NumPy_FOUND)

	# find SWIG
	find_package(SWIG REQUIRED)

	# only support versions of SWIG higher than 4.0
	if(SWIG_VERSION VERSION_LESS 4.0)
		message(STATUS "SWIG version too old (${SWIG_VERSION} < 4.0). Code will not be compile with Python support.")
		set(SWIG_FOUND False)
	endif()

	if(SWIG_FOUND)
		message(STATUS "SWIG installation directory: ${SWIG_DIR}")
		message(STATUS "SWIG executable: ${SWIG_EXECUTABLE}")
	else(SWIG_FOUND)
		message(STATUS "SWIG: NOT FOUND (cannot build with Python bindings)")
	endif(SWIG_FOUND)

	# # add Python headers
	# include_directories(include ${Python_INCLUDE_DIRS})


	# SWIG plugin
	set_source_files_properties(${CMAKE_CURRENT_BINARY_DIR}/livpropa_wrap.cxx PROPERTIES GENERATED True)
	add_custom_command(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/livpropa_wrap.cxx
		COMMAND swig ${SWIG_MODE_FLAG} -c++ -python ${PY3} 
		-I${CRPropa_SWIG_INTERFACE_FILE}
		-I${CMAKE_SOURCE_DIR}/include 
		-I${Python_INCLUDE_DIRS}
		-I{Eigen_INCLUDE_PATH}
		-I${CRPropa_SWIG_PATH}  
		-I${CRPropa_INCLUDE_DIR} ${CRPropa_SWIG_DEFINES} 
		-dirprot -o ${CMAKE_CURRENT_BINARY_DIR}/livpropa_wrap.cxx 
		-outdir ${CMAKE_CURRENT_BINARY_DIR} 
		${CMAKE_CURRENT_SOURCE_DIR}/python/livpropa.i
	)
	include_directories(include 
		${CMAKE_CURRENT_BINARY_DIR}/include 
		${CRPropa_INCLUDE_PATH} 
		${Python_INCLUDE_DIRS} 
		${CMAKE_CURRENT_SOURCE_DIR}/../include
		)

	
	# build interface and source in one library
	add_library(livpropa MODULE ${CMAKE_CURRENT_BINARY_DIR}/livpropa_wrap.cxx ${SOURCE_FILES})

	target_link_libraries(livpropa ${CRPropa_LIBRARY} ${Python_LIBRARIES})

# compile code to work in C++ only, without Python support
elseif()

	# build interface and source in one library
	add_library(livpropa MODULE ${SOURCE_FILES})
	target_link_libraries(livpropa ${CRPropa_LIBRARY})
endif()

set_target_properties(livpropa PROPERTIES PREFIX "")
set_target_properties(livpropa PROPERTIES COMPILE_FLAGS "-w")
set_target_properties(livpropa PROPERTIES OUTPUT_NAME "_livpropa")


# option(ENABLE_TESTING "Build tests and enable test target" On)
# if(ENABLE_TESTING)
# 	set(GTEST_INCLUDE_PATH ${CRPropa_INSTALL_PREFIX}/../libs/gtest/include)
# 	include_directories(${GTEST_INCLUDE_PATH})
# 	link_directories(${CRPropa_INSTALL_PREFIX}/libs/gtest)
# 	message(STATUS "Google testing suite (from CRPropa): ${CRPropa_INSTALL_PREFIX}/libs/gtest/include")
# 	enable_testing()
# endif(ENABLE_TESTING)


# Install
add_definitions(-DLIVpropa_INSTALL_PREFIX="${CMAKE_INSTALL_PREFIX}")
add_definitions(-DLIVpropa_DATA_PATH="${CMAKE_INSTALL_PREFIX}/data")
install(TARGETS livpropa DESTINATION lib)
install(DIRECTORY external/eigen3/Eigen DESTINATION include)
install(DIRECTORY include/ DESTINATION include FILES_MATCHING PATTERN "*.h")
install(DIRECTORY ${CMAKE_BINARY_DIR}/include/ DESTINATION include FILES_MATCHING PATTERN "*.h")
install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/python/livpropa.i DESTINATION swig_interface)
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/data DESTINATION ${CMAKE_INSTALL_PREFIX}/.)


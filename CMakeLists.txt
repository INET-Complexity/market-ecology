CMAKE_MINIMUM_REQUIRED(VERSION 3.12)

SET(CMAKE_MODULE_PATH              "${CMAKE_SOURCE_DIR}/cmake")

PROJECT(market_ecology)

SET(CMAKE_CXX_STANDARD              17)
SET(CMAKE_CXX_STANDARD_REQUIRED     ON)
SET(CMAKE_VERBOSE_MAKEFILE          ON)
SET(CMAKE_COLOR_MAKEFILE            ON)

SET(CONFIGURATION_SHARED ON CACHE BOOL "Build ESL shared library")
ADD_SUBDIRECTORY(${CMAKE_CURRENT_SOURCE_DIR}/ESL)   # make CMake aware of the shared library
SET(CONFIGURATION_SHARED ON CACHE BOOL "Build ESL shared library")

IF(WIN32)
    ADD_DEFINITIONS(-D_NO_W32_PSEUDO_MODIFIERS )
ENDIF()

IF(WIN32)
    set(Boost_USE_STATIC_LIBS   ON)
    SET(CMAKE_EXE_LINKER_FLAGS " -static")
ENDIF()

OPTION(WITH_MPI "Enable MPI using Boost.MPI" OFF)
IF(WITH_MPI)
    MESSAGE("Enabling MPI")
    FIND_PACKAGE(MPI REQUIRED)
    ADD_DEFINITIONS(-DWITH_MPI)
ENDIF()

OPTION(WITH_PYTHON "Enables export of Python bindings using Boost.Python" OFF)
IF(WITH_PYTHON)
    ADD_DEFINITIONS(-DWITH_PYTHON)
    ADD_DEFINITIONS(-DBOOST_PYTHON_STATIC_LIB)
ENDIF()

OPTION(WITH_QL  "Enables QuantLib" OFF)
IF(WITH_QL)
    MESSAGE("Enabling QuantLib")
    ADD_DEFINITIONS(-DWITH_QL)
ENDIF()

#set(CMAKE_POSITION_INDEPENDENT_CODE ON)

IF(CMAKE_BUILD_TYPE MATCHES RELEASE OR ESL_BUILD_RELEASE MATCHES RELEASE)
    INCLUDE(CheckCXXCompilerFlag)
    CHECK_CXX_COMPILER_FLAG("-march=native" COMPILER_SUPPORTS_MARCH_NATIVE)
    if(COMPILER_SUPPORTS_MARCH_NATIVE)
        SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")
    endif()
ENDIF(CMAKE_BUILD_TYPE MATCHES RELEASE OR ESL_BUILD_RELEASE MATCHES RELEASE)

SET(CMAKE_C_FLAGS_DEBUG "-g -Og")#   --coverage -fprofile-arcs -ftest-coverage") #  -fanalyzer
SET(CMAKE_CXX_FLAGS_DEBUG "-g -Og ")#    --coverage -fprofile-arcs -ftest-coverage")

SET(CMAKE_C_FLAGS_RELEASE "-O3  ")
SET(CMAKE_CXX_FLAGS_RELEASE "-O3  ")

ADD_DEFINITIONS(-DADEPT_RECORDING_PAUSABLE)
ADD_DEFINITIONS(-DADEPT_STORAGE_THREAD_SAFE) #ADD_DEFINITIONS(-DADEPT_STACK_THREAD_UNSAFE)
#ADD_DEFINITIONS(-DADEPT_TRACK_NON_FINITE_GRADIENTS)
#ADD_DEFINITIONS(-DADEPT_NO_AUTOMATIC_DIFFERENTIATION)

################################################################################

IF(WITH_PYTHON)
    FIND_PACKAGE(PythonLibs 3)  # for building bindings
    FIND_PACKAGE(PythonInterp 3)# for running tests
ENDIF()

################################################################################


IF(WITH_PYTHON)
    LIST(APPEND Boost_LIBRARIES python)
ELSE()
ENDIF()

IF(WITH_MPI)
    LIST(APPEND Boost_LIBRARIES mpi)
    IF(WITH_PYTHON)
        LIST(APPEND Boost_LIBRARIES mpi_python38)
    ENDIF()
ENDIF()

################################################################################
#SET(Boost_USE_STATIC_LIBS ON)
LIST(APPEND Boost_LIBRARIES
        program_options
        date_time
        serialization
        unit_test_framework
        thread
        #log
        #log_setup
        )
#ADD_DEFINITIONS(-DBOOST_LOG_DYN_LINK)
#ADD_DEFINITIONS(-DBOOST_ALL_DYN_LINK)
FIND_PACKAGE(Boost 1.65 COMPONENTS ${Boost_LIBRARIES} REQUIRED)
################################################################################
IF(WIN32)
    LIST(APPEND ALL_INCLUDE_DIRECTORIES ${PROJECT_SOURCE_DIR} SYSTEM  )
ELSE()
    LIST(APPEND ALL_INCLUDE_DIRECTORIES ${PROJECT_SOURCE_DIR})
ENDIF()

IF(WITH_MPI)
    list(APPEND ALL_INCLUDE_DIRECTORIES ${MPI_INCLUDE_PATH})
ENDIF()

IF(WITH_PYTHON)
    list(APPEND ALL_INCLUDE_DIRECTORIES ${PYTHON_INCLUDE_DIRS})
ENDIF()



INCLUDE_DIRECTORIES(${ALL_INCLUDE_DIRECTORIES} ${CMAKE_CURRENT_SOURCE_DIR}/ESL)

################################################################################


LIST(APPEND ALL_LINKED_LIBRARIES pthread esl)

IF(WITH_PYTHON)
    LIST(APPEND ALL_LINKED_LIBRARIES ${PYTHON_LIBRARIES})
ENDIF()

IF(WITH_MPI)
    LIST(APPEND ALL_LINKED_LIBRARIES ${MPI_CXX_LIBRARIES} )
ENDIF()

IF(WITH_QL)
    list(APPEND ALL_LINKED_LIBRARIES QuantLib ${QuantLib_LIBRARIES})
ENDIF()


IF(WIN32)
    LINK_LIBRARIES(${ALL_LINKED_LIBRARIES} -static -static-libgcc -static-libstdc++)
ELSE()
    LINK_LIBRARIES(${ALL_LINKED_LIBRARIES} dl )
ENDIF()


################################################################################

file(GLOB_RECURSE ESL_SOURCES "ESL/esl/**.cpp" "ESL/esl/**.hpp")

SET(SOURCE_FILES main.cpp
        traded_company.cpp traded_company.hpp
        fund.cpp fund.hpp

        strategies/noise_trader.cpp strategies/noise_trader.hpp

        strategies/fundamental_value/dividend_discount.hpp
        strategies/fundamental_value/dividend_discount.cpp

        strategies/technical/momentum.hpp
        strategies/technical/momentum.cpp

        strategies/kelly_bettor.hpp
        strategies/kelly_bettor.cpp

        strategies/constant_demand.hpp
        strategies/constant_demand.cpp

        strategies/fundamental_value/mean_reverting_noise.cpp
        strategies/fundamental_value/mean_reverting_noise.hpp

        experiment/e_3_grid.cpp
        experiment/e_3_grid.hpp

        experiment/e_4_trajectories.cpp
        experiment/e_4_trajectories.hpp

        experiment/e_5.cpp
        experiment/e_5.hpp

        experiment/e_6.cpp
        experiment/e_6.hpp
        )

#list(TRANSFORM ESL_SOURCES PREPEND "ESL/")

ADD_EXECUTABLE(market_ecology ${SOURCE_FILES} ${ESL_SOURCES})

################################################################################


FILE(GLOB TEST_SRCS RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} test*/test_*.cpp)

FOREACH(test_src ${TEST_SRCS})
    MESSAGE("TEST " ${test_src})
    GET_FILENAME_COMPONENT(test_name ${test_src} NAME_WE)

    ADD_EXECUTABLE(${test_name} ${test_src})
    TARGET_LINK_LIBRARIES(${test_name} ${Boost_LIBRARIES} ${ESL_TARGET_NAME})
    SET_TARGET_PROPERTIES(${test_name} PROPERTIES
            RUNTIME_OUTPUT_DIRECTORY  ${CMAKE_BINARY_DIR}/test/)

    # test_mpi_* and test_mpi_python_* need to be run using a MPI job runner
    IF(${test_name} MATCHES "^(test_mpi_)")
        IF(WITH_MPI)
            MESSAGE("\t MPI TEST")
            ADD_TEST(NAME ${test_name}
                    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/test/
                    COMMAND mpirun -np 2 ${CMAKE_BINARY_DIR}/test/${test_name})
        ENDIF()
    ELSE()
        ADD_TEST(NAME ${test_name}
                WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/test/
                COMMAND ${CMAKE_BINARY_DIR}/test/${test_name})
    ENDIF()
ENDFOREACH(test_src)
